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
| `ultrawide_horizon.py` | Optical horizon from the ultrawide lens (fired with the tele): an acceleration-immune tilt reference, fused with the IMU. Reuses `starfix.get_dip_of_horizon`. |
| `capture_trigger.py` | Per-regime hand/platform disturbance model and the least-rotation "smart shutter" (gated vs. periodic). |
| `scenario.py` | Sea/land/air ground truth from the ephemeris + noisy Sun+Moon measurements + IMU stream. |
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
```

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

- Accuracy tracks **azimuth separation** of the two bodies (best near 90°).

## Modelling assumptions & honest limitations

- IMU/camera specs are **representative order-of-magnitude** values, not Apple
  datasheet figures.
- **Spherical Earth**; observer positions geocentric. Because truth *and*
  estimate share this convention, the reported error is a pure estimation error
  (unaffected by the geocentric-vs-geodetic figure offset). The Sun/Moon
  **ephemeris is the real one** from `starfix` (nautical-almanac CSVs).
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
