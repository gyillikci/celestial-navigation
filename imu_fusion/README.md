<!---
    © 2026.  MIT License (see LICENSE file).
-->

# iPhone 17 Pro daytime Sun + Moon sighting — IMU synthetic horizon + factor graph

A reproducible **modelling study**: can a hand-held iPhone 17 Pro fix its
position by photographing the **Sun and the Moon together in daylight**, using
the phone's **inertial sensors as the horizon** and a **GTSAM factor graph** to
beat the error down? Studied for **sea, land and air**.

The premise is real — the daytime Moon and the Sun are both easy phone subjects
(see the reference photos: a half/gibbous Moon in a blue sky, and the Sun as a
clean disk). What a phone lacks that a sextant has is a **horizon**: at 10 000 m
or on a hazy sea there is none. So we take the horizon from **gravity**: the
IMU's attitude gives the local vertical, turning the phone into a *digital
theodolite*. Each photo yields the measured **altitude** of the Sun and of the
Moon; two altitudes are two lines of position; many shots, linked by IMU
dead-reckoning, are fused into one trajectory with a covariance.

## The idea in one paragraph

The phone is a digital theodolite. Its **IMU** defines the gravity-horizontal
plane (a synthetic *artificial horizon*, so there is **no sea-horizon dip** to
correct), and its **camera** gives the direction to a body (the pixel of the
disk centre). Together → a measured altitude. A single hand-held sight is weak:
the synthetic horizon is only good to ~6′ braced, tens of arc-minutes while
moving (a sextant reaches ~1–2′). The win is **statistical** — fuse many *calm*
shots. An in-app **least-rotation shutter** fires the camera at the quietest
instants (local minima of |ω|) to keep each IMU-photo pair clean.

## What's here

| File | Role |
|---|---|
| `astro.py` | Body geographic position + predicted altitude/azimuth. Reuses `starfix`'s real Sun/Moon ephemeris. Local ENU ↔ lat/lon. |
| `iphone_model.py` | Representative iPhone-class IMU + tele-camera noise; how motion corrupts the gravity horizon. |
| `realtime.py` | Streaming estimator (`gtsam_unstable.IncrementalFixedLagSmoother`) — bounded per-shot latency for on-device use; `bench.py` benchmarks it vs batch. |
| `ultrawide_horizon.py` | Optical horizon from the ultrawide lens (fired with the tele): an acceleration-immune tilt reference, fused with the IMU. Reuses `starfix.get_dip_of_horizon`. |
| `optical_attitude.py` | Orientation from the tele-resolved disk: Moon bright-limb PA, illuminated fraction, Sun P-angle, and the parallactic angle *q* → magnetometer-free heading + position line. The *differential* Sun−Moon orientation (`differential_orientation_sigma_deg`) is the genuinely horizon-free line (shared platform roll cancels). |
| `visual_anchor.py` | Tracked sunspot / disk feature as a star-tracker anchor: bounds gyro drift and gives an acceleration-immune attitude/vertical (moving or stationary) — rescues the fix when the optical horizon is unavailable. Also the cloud-outage coast models. |
| `cloud.py` | Temporally-correlated cloud occlusion (Markov passages): drops obscured sights and coasts the anchor — graceful degradation and a coast-time budget. |
| `capture_trigger.py` | Per-regime hand/platform disturbance model and the least-rotation "smart shutter" (gated vs. periodic). |
| `scenario.py` | Sea/land/air ground truth from the ephemeris + noisy Sun+Moon measurements + IMU stream. |
| `validate_ephemeris.py` | **Ground-truth check**: recomputes the Sun/Moon GHA/Dec from an *independent* engine (astropy/ERFA, or Skyfield/JPL if a kernel is present) and reports the residual vs the almanac. Also loads a **Stellarium** CSV export as a third witness (`stellarium_reference.md`). |
| `mission_plan.py` | **Mission planner**: for a route A→B, forecasts the Sun/Moon sky and answers the operational question *"is the Moon available, and when?"* — the both-bodies-up windows, a SUN-ONLY fallback flag, and warnings (near-new/full Moon, near-zenith Sun). CLI: `python -m imu_fusion.mission_plan --from Istanbul --to Ankara --next`. |
| `celestial_factor_graph.py` | GTSAM graph: celestial altitude `CustomFactor` per sight, `ImuFactor` between shots, priors, LM solve, `Marginals` covariance. |
| `baseline.py` | Incumbent `starfix` single-epoch two-LOP fix, for comparison. |
| `run_study.py` | Runs the whole study → `results/` figures, `RESULTS.md`, `dashboard.html`. |
| `../test/test_imu_fusion.py` | Fast regression tests. |

## Running it

```bash
pip install -r imu_fusion/requirements.txt     # gtsam wheel + numpy<2 + matplotlib + pandas
python -m imu_fusion.run_study                  # writes RESULTS.md, dashboard.html, results/*.png
python -m unittest test.test_imu_fusion         # tests
# quick component demos:
python -m imu_fusion.iphone_model
python -m imu_fusion.capture_trigger
python -m imu_fusion.scenario
python -m imu_fusion.validate_ephemeris    # independent ground-truth check (needs astropy or skyfield)
```

## The unified fusion

One factor graph fuses **every observable at once** (this is `FULL_SC`/`FULL_SV`
in `run_study.py`):

| Observable | Sensor | Role |
|---|---|---|
| Altitude of Sun & Moon | tele + horizon reference | position lines |
| Horizon reference | IMU gravity ⊕ ultrawide optical horizon | sets altitude covariance |
| Azimuth of Sun & Moon | tele-disk heading (magnetometer-free) | position lines |
| Parallactic angle *q* (per body) | tele-disk orientation vs. vertical | position line (needs the horizon) |
| Differential Δq (Sun−Moon) | difference of the two disks' orientations | **horizon-free** line (shared roll cancels) |
| IMU preintegration + bias | IMU | links & smooths the trajectory |
| Least-rotation gating | gyro | picks clean shutter instants |
| Coarse position prior | last known / DR | disambiguation only |

Per keyframe: a pose `X(i)` (position + attitude), velocity `V(i)`, shared bias
`B`; six celestial factors per shot (alt/az/q × Sun/Moon); IMU factors between
shots. A **leave-one-out ablation** (`exp_ablation`, `fig_ablation.png`) shows
what each observable is worth inside the combined graph. Deliberately out:
lunar-distance timing (a phone already has a clock), and refraction / Moon
parallax (handled by `starfix.Sight` in a field build). Per-shot attitude is
marginalised into the celestial factor covariances.

## Headline results

Canonical epoch: **2026-03-24 12:00 UTC near Greenwich**, first-quarter Moon,
Sun and Moon ~94° apart in azimuth (a well-conditioned fix). RMS horizontal
position error, mean over 8 seeds — see [`RESULTS.md`](RESULTS.md) for the full
tables and figures.

| Regime | starfix single-fix | Factor graph, no IMU | **Factor graph + IMU** |
|---|---|---|---|
| Land (stationary) | ~17 km | ~17 km | **~4 km** |
| Sea (vessel + swell) | ~240 km | ~160 km | **~25 km** |
| Air (aircraft) | ~75 km | ~75 km | **~8 km** |

- Fusing many Sun+Moon shots with IMU cuts error **4–10×** vs a single fix.
- The **least-rotation shutter** cuts per-shot horizon noise 3–6× → ≈2× lower
  fix error **on land and in the air**; **at sea it is a wash** (even the
  calmest swell instant is too tilted for the IMU gravity horizon alone).
- **The ultrawide camera's optical horizon** (captured simultaneously with the
  tele body shot) is immune to acceleration and **rescues the moving
  platforms**: sea **~27 km → ~2 km**, air **~9 km → ~2 km**, bringing them to
  land-class accuracy. On land there is no true sea horizon, so it falls back to
  the IMU.

  | Regime | IMU horizon | ultrawide horizon | **fused** |
  |---|---|---|---|
  | Land | ~4 km | ~4 km (n/a) | **~4 km** |
  | Sea | ~27 km | ~2 km | **~2 km** |
  | Air | ~9 km | ~2 km | **~2 km** |

- **The tele lens is an instrument, not a pointer.** Resolving the disk — the
  Moon's bright limb, the Sun's sunspot P-angle — gives an absolute celestial
  orientation: a **magnetometer-free heading** (~0.35° vs ~1° for the phone
  compass) that makes the Sun/Moon **azimuth lines of position usable**, plus a
  horizon-free **parallactic-angle** position line. On a weak (IMU) horizon this
  cuts the sea fix ~28 → ~12 km and the air fix ~12 → ~5 km. The Moon's limb is
  the workhorse; the Sun's P-angle needs a solar filter and visible spots.
- Accuracy tracks **azimuth separation** of the two bodies (best near 90°).

## Real-time on-device

`realtime.py` runs the same fusion **incrementally** with a fixed-lag smoother
so per-shot latency is **bounded** regardless of voyage length. Benchmark
(`python -m imu_fusion.bench`): batch re-solve reaches **~330 ms/fix at 30
shots and grows**, while streaming stays **flat at ~5–6 ms/update** — same
current-position accuracy (within ~0.02 km). Fixes are low-rate (a gated shot
every few seconds); IMU dead-reckoning gives the continuous position between
them. The main saving is doing IMU preintegration **once online** instead of
re-preintegrating every leg per solve; analytic Jacobians + a reduced two-DOF
finite difference (only east/north move a celestial factor) remove the rest.
The image processing (disk centroid, horizon-line fit, limb/spot detect) runs
on the Neural Engine/GPU and is out of scope here.

## Modelling assumptions & honest limitations

- IMU/camera specs are **representative order-of-magnitude** values, not Apple
  datasheet figures.
- **Spherical Earth**; observer positions geocentric. Because truth *and*
  estimate share this convention, the reported error is a pure estimation error
  (unaffected by the geocentric-vs-geodetic figure offset). The Sun/Moon
  **ephemeris is the real one** from `starfix` (nautical-almanac CSVs), and it is
  **independently cross-checked**: `validate_ephemeris.py` recomputes the GHA/Dec
  from a separate engine (astropy/ERFA, or Skyfield/JPL, or a Stellarium export)
  and finds agreement to a few arc-seconds — well under the arc-minute
  measurement noise. That check caught a real sign bug in the almanac reader near
  Dec ≈ 0; see the *Ground-truth validation* section of `RESULTS.md`.
- Altitudes are **geometric** (refraction/parallax/semidiameter treated as
  pre-corrected); `starfix.Sight` already models those and would be layered in
  for a field build.
- Straight, constant-velocity legs with a known level attitude; manoeuvres
  (turns, climb) are left as an extension.
- Batch Levenberg–Marquardt is used; `gtsam.ISAM2` would give the same graph a
  real-time incremental form on-device.
- This is a **simulation of data gathering + fusion**, not iOS app code. The
  natural on-device seam already exists in the repo:
  `az_calc.get_latlon_for_solar_obs(azimuth, altitude, timestamp)` and the
  `Sight(object_name, set_time, measured_alt, …)` constructor.
  [`ios_app_architecture.md`](ios_app_architecture.md) maps every validated module
  onto a concrete iPhone 17 Pro app (frameworks, capture/solver data flow, the
  on-device calibration protocol for the tilt floor, Xcode layout).
