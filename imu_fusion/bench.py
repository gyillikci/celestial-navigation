''' Real-time benchmark: batch re-solve vs incremental streaming.

    Shows the core real-time argument:
      * batch Levenberg-Marquardt re-optimises the WHOLE trajectory each fix, so
        its per-fix latency grows with the number of shots so far;
      * the fixed-lag smoother marginalises old keyframes, so its per-update
        latency is BOUNDED (flat) regardless of trip length;
    while the current-position accuracy of the two matches.

    Timings are wall-clock on the host CPU (single-threaded); an iPhone 17 Pro
    (A19 Pro) is comparable or faster single-thread. All randomness is seeded.

    Run:  python -m imu_fusion.bench
    (c) 2026.  MIT License (see LICENSE file).
'''

import time
import random
from dataclasses import replace

import numpy as np

import imu_fusion.celestial_factor_graph as fg
from .scenario import build_scenario
from .realtime import StreamingEstimator

FULL_SC = dict(horizon_mode="fused", use_azimuth=True, heading_source="optical",
               use_parallactic=True, sun_spots=True)
FULL_SV = dict(use_imu=True, use_azimuth=True, use_parallactic=True)


def _timed_batch_solve(sc):
    t0 = time.perf_counter()
    fg.solve(sc, **FULL_SV)
    return (time.perf_counter() - t0) * 1000.0     # ms


def latency_vs_length(regime, shot_counts, seeds=3):
    ''' Per-fix latency (ms) for batch (full re-solve of n shots) vs streaming
        (the n-th incremental update), as n grows. '''
    batch_ms = {n: [] for n in shot_counts}
    stream_ms = {n: [] for n in shot_counts}
    nmax = max(shot_counts)
    want = set(shot_counts)
    for s in range(seeds):
        base = build_scenario(regime, random.Random(700 + s), n_shots=nmax,
                              **FULL_SC)
        # batch: re-solve the first n keyframes from scratch
        for n in shot_counts:
            sc_n = replace(base, keyframes=base.keyframes[:n])
            batch_ms[n].append(_timed_batch_solve(sc_n))
        # streaming: one pass, time each update; record at the wanted n
        est = StreamingEstimator(base, **{k: FULL_SV[k] for k in
                                          ("use_imu", "use_azimuth",
                                           "use_parallactic")})
        kfs = base.keyframes
        for k, kf in enumerate(kfs):
            t0 = time.perf_counter()
            est.add_keyframe(kf)
            dt = (time.perf_counter() - t0) * 1000.0
            if (k + 1) in want:
                stream_ms[k + 1].append(dt)
            if FULL_SV["use_imu"] and k < len(kfs) - 1:
                est.integrate_imu(kf.imu_to_next)
    agg = lambda d: {n: (float(np.mean(v)), float(np.std(v))) for n, v in d.items()}
    return agg(batch_ms), agg(stream_ms)


def jacobian_speedup(regime, n_shots=14, seeds=4):
    ''' Solve time with analytic vs reduced finite-difference Jacobians. '''
    out = {}
    for name, analytic in (("analytic", True), ("reduced-FD", False)):
        fg.USE_ANALYTIC_JAC = analytic
        ts = []
        for s in range(seeds):
            sc = build_scenario(regime, random.Random(800 + s), n_shots=n_shots,
                                **FULL_SC)
            t0 = time.perf_counter()
            fg.solve(sc, **FULL_SV)
            ts.append((time.perf_counter() - t0) * 1000.0)
        out[name] = float(np.mean(ts))
    fg.USE_ANALYTIC_JAC = True
    return out


def accuracy_parity(regime, n_shots=14, seeds=6):
    ''' Current-position (final keyframe) error: batch vs streaming. '''
    from .realtime import solve_streaming
    bf, sf = [], []
    for s in range(seeds):
        sc = build_scenario(regime, random.Random(900 + s), n_shots=n_shots,
                            **FULL_SC)
        bf.append(fg.solve(sc, **FULL_SV)["final_err_km"])
        sf.append(solve_streaming(sc)["final_err_km"])
    return float(np.mean(bf)), float(np.mean(sf))


def run(shot_counts=(2, 6, 10, 16, 22, 30)):
    ''' Assemble the benchmark data structure for the study/report. '''
    data = {"shot_counts": list(shot_counts), "latency": {}, "jacobian": {},
            "parity": {}}
    for r in ("land", "sea", "air"):
        b, s = latency_vs_length(r, shot_counts)
        data["latency"][r] = {"batch": b, "stream": s}
        data["jacobian"][r] = jacobian_speedup(r)
        data["parity"][r] = accuracy_parity(r)
    return data


if __name__ == "__main__":
    d = run()
    for r in ("land", "sea", "air"):
        lat = d["latency"][r]
        nmax = d["shot_counts"][-1]
        print(f"[{r}]")
        print(f"  batch  @ {nmax} shots : {lat['batch'][nmax][0]:.1f} ms/fix")
        print(f"  stream @ {nmax} shots : {lat['stream'][nmax][0]:.1f} ms/update")
        print(f"  jacobian analytic {d['jacobian'][r]['analytic']:.0f} ms vs "
              f"reduced-FD {d['jacobian'][r]['reduced-FD']:.0f} ms")
        bf, sf = d["parity"][r]
        print(f"  final-error parity: batch {bf:.2f} km, stream {sf:.2f} km")
