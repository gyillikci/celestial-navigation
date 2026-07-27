''' Least-rotation capture trigger ("smart shutter") for the sighting app.

    A hand-held phone is never perfectly still.  Its angular rate |omega| and
    residual linear acceleration wander with hand tremor, boat swell or aircraft
    turbulence.  Because the synthetic horizon (see `iphone_model`) degrades
    sharply with motion, the app should fire the camera automatically at the
    *calmest* instants -- local minima of a stillness metric -- instead of on a
    naive periodic/manual shutter.  This yields cleaner IMU-photo pairs.

    This module (a) simulates a realistic per-regime disturbance trace and
    (b) selects shutter instants two ways -- gated (least-rotation) and ungated
    (periodic) -- so the study can quantify the benefit.

    (c) 2026.  MIT License (see LICENSE file).
'''

from dataclasses import dataclass
from math import sin, pi

from .iphone_model import KinematicState


@dataclass(frozen=True)
class DisturbanceProfile:
    ''' Per-regime hand/platform disturbance model.

        The trace is a sum of a slow oscillation (swell / sway) and fast tremor
        (hand / vibration), for both angular rate [rad/s] and residual linear
        acceleration [m/s^2].
    '''
    name: str
    omega_base: float
    omega_swell_amp: float
    omega_tremor: float
    accel_base: float
    accel_swell_amp: float
    accel_tremor: float
    swell_period_s: float


# Representative profiles.  "Sea" has strong slow swell; "air" has smaller slow
# motion but persistent higher-frequency vibration; "land" is a braced hand.
PROFILES = {
    "land": DisturbanceProfile("land", 0.020, 0.010, 0.060,
                               0.030, 0.020, 0.120, 3.0),
    "sea":  DisturbanceProfile("sea", 0.120, 0.350, 0.120,
                               0.250, 1.300, 0.300, 6.0),
    "air":  DisturbanceProfile("air", 0.060, 0.060, 0.140,
                               0.250, 0.200, 0.600, 4.0),
}


def simulate_trace(profile: DisturbanceProfile, duration_s: float,
                   rate_hz: float, rng) -> tuple[list, list, list]:
    ''' Simulate |omega|(t) and |a_lin|(t) for a regime.

        Returns (times, omega, accel) as lists.  `rng` is a random.Random so the
        trace is reproducible (no wall-clock randomness).
    '''
    n = int(duration_s * rate_hz)
    times, omega, accel = [], [], []
    # AR(1) tremor states for temporal correlation.
    w_tremor = 0.0
    a_tremor = 0.0
    alpha = 0.85
    for i in range(n):
        t = i / rate_hz
        swell = sin(2 * pi * t / profile.swell_period_s)
        # |swell| gives the rotational magnitude peaking twice per period.
        w_tremor = alpha * w_tremor + (1 - alpha) * rng.gauss(0, profile.omega_tremor)
        a_tremor = alpha * a_tremor + (1 - alpha) * rng.gauss(0, profile.accel_tremor)
        w = abs(profile.omega_base + profile.omega_swell_amp * abs(swell) + w_tremor)
        a = abs(profile.accel_base + profile.accel_swell_amp * abs(swell) + a_tremor)
        times.append(t)
        omega.append(w)
        accel.append(a)
    return times, omega, accel


def stillness(omega: float, accel: float, g: float = 9.80665) -> float:
    ''' Combined stillness cost (lower = calmer).  Weighs angular rate and the
        tilt-inducing linear acceleration on comparable footing (a/g is the tilt
        angle in rad; omega is rad/s over a ~short exposure).
    '''
    return omega + accel / g


def find_shutter_instants(times: list, omega: list, accel: list,
                          n_shots: int, min_gap_s: float,
                          gated: bool) -> list[tuple[float, KinematicState]]:
    ''' Select `n_shots` shutter instants.

        gated=True  : least-rotation -- the calmest local minima of the
                      stillness metric, respecting a minimum spacing.
        gated=False : naive periodic shutter -- evenly spaced in time,
                      regardless of motion (what a manual/timer capture does).

        Returns a list of (time_s, KinematicState) at the chosen instants.
    '''
    rate_hz = 1.0 / (times[1] - times[0])
    min_gap = int(min_gap_s * rate_hz)

    if not gated:
        # Evenly spaced indices across the trace.
        step = max(1, len(times) // (n_shots + 1))
        idxs = [min((k + 1) * step, len(times) - 1) for k in range(n_shots)]
    else:
        cost = [stillness(w, a) for w, a in zip(omega, accel)]
        # Local minima of the stillness cost.
        minima = [i for i in range(1, len(cost) - 1)
                  if cost[i] <= cost[i - 1] and cost[i] < cost[i + 1]]
        minima.sort(key=lambda i: cost[i])
        chosen: list[int] = []
        for i in minima:
            if all(abs(i - j) >= min_gap for j in chosen):
                chosen.append(i)
            if len(chosen) == n_shots:
                break
        idxs = sorted(chosen)

    return [(times[i], KinematicState(ang_rate=omega[i], lin_accel=accel[i]))
            for i in idxs]


def summarise_regime(regime: str, rng, n_shots: int = 12,
                     duration_s: float = 120.0) -> str:
    ''' Compare mean stillness of gated vs ungated captures for a regime. '''
    prof = PROFILES[regime]
    t, w, a = simulate_trace(prof, duration_s, 50.0, rng)
    g = find_shutter_instants(t, w, a, n_shots, 3.0, gated=True)
    u = find_shutter_instants(t, w, a, n_shots, 3.0, gated=False)
    mean = lambda pts: sum(stillness(s.ang_rate, s.lin_accel) for _, s in pts) / len(pts)
    return (f"{regime:5}  gated stillness={mean(g):.3f}   "
            f"ungated stillness={mean(u):.3f}   "
            f"improvement x{mean(u) / mean(g):.1f}")


if __name__ == "__main__":
    import random
    for r in ("land", "sea", "air"):
        print(summarise_regime(r, random.Random(42)))
