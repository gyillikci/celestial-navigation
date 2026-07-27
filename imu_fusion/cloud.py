''' Cloud occlusion of the tracked Sun / Moon.

    Cloud hits two things at once: the SIGHT for the obscured body (no line of
    position; a cloud edge would also corrupt the centroid/limb, so those frames
    are rejected) and the VISUAL ANCHOR (gyro calibration stops, attitude coasts
    -- see visual_anchor.coast_attitude_arcmin).

    Cloud does not flicker independently frame to frame; it drifts across the sky
    in patches, so occlusion comes in RUNS.  This module models each body's
    clear/obscured state as a two-state Markov chain over the shot times,
    parameterised by a clear-sky fraction and a mean cloud-passage duration.

    All randomness flows through an injected random.Random for reproducibility.

    (c) 2026.  MIT License (see LICENSE file).
'''

from dataclasses import dataclass


@dataclass(frozen=True)
class CloudSpec:
    ''' Sky-condition model. '''
    clear_fraction: float = 1.0        # long-run fraction of time a body is clear
    mean_passage_s: float = 25.0       # mean duration of a cloud passage / clear run
    correlated_bodies: float = 0.5     # 0 = independent skies, 1 = Sun&Moon share


def _markov_clear(times, clear_fraction, mean_passage_s, rng):
    ''' A clear/obscured boolean per time via a two-state Markov chain whose
        stationary clear-probability is `clear_fraction` and whose mean dwell in
        each state is ~`mean_passage_s`. '''
    if clear_fraction >= 1.0:
        return [True] * len(times)
    if clear_fraction <= 0.0:
        return [False] * len(times)
    out = []
    clear = rng.random() < clear_fraction
    prev_t = times[0] if times else 0.0
    # Per-second switch rates giving the target dwell time and clear-fraction.
    for i, t in enumerate(times):
        dt = max(0.0, t - prev_t) if i > 0 else 0.0
        # Probability of a state switch over dt (exponential dwell).
        if clear:
            # clear -> obscured; dwell in clear ~ mean_passage_s * f/(1-f) skew
            rate = (1.0 - clear_fraction) / mean_passage_s
        else:
            rate = clear_fraction / mean_passage_s
        p_switch = 1.0 - pow(2.718281828, -rate * dt) if dt > 0 else 0.0
        if rng.random() < p_switch:
            clear = not clear
        out.append(clear)
        prev_t = t
    return out


def body_clear_flags(times, bodies, cloud: CloudSpec, rng):
    ''' Return {body: [clear?  per shot]} for the given shot times. '''
    if cloud is None:
        return {b: [True] * len(times) for b in bodies}
    flags = {}
    base = _markov_clear(times, cloud.clear_fraction, cloud.mean_passage_s, rng)
    for b in bodies:
        if cloud.correlated_bodies >= 1.0:
            flags[b] = list(base)
        else:
            own = _markov_clear(times, cloud.clear_fraction,
                                cloud.mean_passage_s, rng)
            # Blend the shared sky and the body's own sky.
            flags[b] = [base[i] if rng.random() < cloud.correlated_bodies
                        else own[i] for i in range(len(times))]
    return flags


def availability(flags) -> float:
    ''' Fraction of (body, shot) pairs that are clear. '''
    total = sum(len(v) for v in flags.values())
    clear = sum(sum(1 for c in v if c) for v in flags.values())
    return clear / total if total else 0.0


if __name__ == "__main__":
    import random
    ts = [i * 8.0 for i in range(30)]
    for cf in (1.0, 0.7, 0.4, 0.1):
        fl = body_clear_flags(ts, ("Sun", "Moon"),
                              CloudSpec(clear_fraction=cf), random.Random(1))
        print(f"clear_fraction={cf}:  availability={availability(fl):.2f}  "
              f"Sun runs={''.join('C' if c else '.' for c in fl['Sun'])}")
