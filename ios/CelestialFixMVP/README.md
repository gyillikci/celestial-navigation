# CelestialFixMVP — a buildable iPhone slice of the daytime celestial fix

The smallest **end-to-end, dependency-free** slice of the larger iPhone 17 Pro
celestial-navigation app (spec: [`../../imu_fusion/ios_app_architecture.md`](../../imu_fusion/ios_app_architecture.md)):

> Point the **telephoto** camera at the Moon, tap **Capture**, and get a
> position fix with a real 1σ covariance — from **one** photo plus the phone's
> gravity vector, no GTSAM, no ultrawide, no Sun, no network.

It exists to prove the capture → metrology → sight → fix chain runs on a real
device, and to be the scaffold the full app grows into. Every heavy module
(differential Sun−Moon Δq, the fixed-lag GTSAM smoother, the ultrawide horizon,
the Metal rotational-NCC roll) is deliberately **out of scope here** — see the
map below for where each one lands.

---

## ⚠️ Honest caveats — read before trusting a number

- **I have not compiled, run, signed, or tested this.** It was written on Linux
  with no macOS/Xcode toolchain. Treat it as a careful first draft: the math
  mirrors the validated Python in `imu_fusion/`, but Swift compile errors,
  CoreMotion axis conventions, and camera behaviour must be shaken out **on the
  device**. The `Tests/` target (⌘U) is the first thing to run.
- **Refraction + topocentric parallax ARE applied** (`CelestialMath.swift`,
  ported from `imu_fusion/corrections.py`): Bennett refraction and exact parallax
  from the Moon's distance — the ~0.5–1° parallax that used to dominate the
  budget is now removed. These transforms are validated against IAU ERFA and
  astropy to **sub-arcminute** in `validate_ephemeris.py`, and the Swift port is
  checked against the exported `GoldenVectors.swift` in the test target.
- **Low-precision ephemeris — the remaining accuracy gap.** `MoonEphemeris.swift`
  is still a *truncated* Meeus ch.47 lunar theory — good to a **few arc-minutes**
  (≈ a few km on the ground), enough to demonstrate the pipeline, **not** the
  arc-second ephemeris the production app needs. The next step swaps it for a
  compact table baked from the validated `imu_fusion/astro.py` / `starfix`
  pipeline (cross-checked by `validate_ephemeris.py`). Until then this MVP is a
  plumbing demo, not a navigation instrument.
- **One sight + a prior, not a redundant fix.** A single altitude line is
  under-determined; the GPS/last-known **prior** pins the along-line direction.
  So the reported covariance is honestly bounded by the prior. Two bodies (Sun +
  Moon) or two Moon sights minutes apart are what remove that dependence — that's
  the full app, not this slice.
- **Night / high-contrast Moon assumed** for the limb seed (bright disk on dark
  sky). The daytime low-contrast path uses the gradient sky-limb RANSAC seed in
  `disk_metrology.py`; only the full-circle gradient refine is ported here.
- **The 0.1° tilt floor is an assumption, not a measurement.** `altitudeSigmaDeg`
  hard-codes it. Measure the real synthetic-horizon σ on the device with the
  calibration protocol in the architecture spec and feed it back.

---

## Build & run

Requires macOS + Xcode 16 (iOS 17 SDK) and a physical iPhone (the camera and
CoreMotion don't exist in the Simulator). The project is generated with
[XcodeGen](https://github.com/yonaskolb/XcodeGen) so the repo carries no
`.xcodeproj` blob.

```bash
brew install xcodegen           # once
cd ios/CelestialFixMVP
xcodegen generate               # writes CelestialFixMVP.xcodeproj
open CelestialFixMVP.xcodeproj
```

Then in Xcode:

1. Select the **CelestialFixMVP** scheme.
2. Signing & Capabilities → set your **Team** (automatic signing).
   The bundle id is `com.example.celestialfix.CelestialFixMVP`; change the prefix
   in `project.yml` if that collides.
3. Pick your connected iPhone as the run destination and **⌘R**.
4. Grant the **camera**, **motion**, and **location** prompts on first launch
   (usage strings live in `project.yml` → generated `Info.plist`).
5. **⌘U** runs `Tests/PipelineTests.swift` — device-independent math checks
   (ephemeris sanity, sub-point altitude = 90°, the fix satisfies the sight and
   moves off the prior, sub-pixel limb recovers a synthetic full disk, and the
   full-circle centre beats the lit-blob centroid on a half phase). Run these
   first; they need no camera.

No Swift package dependencies — only Apple frameworks (SwiftUI, AVFoundation,
CoreMotion, CoreLocation, MapKit, simd).

## How one capture becomes a fix

```
tele frame (BGRA) ─► toGray ─┐
                             ├─► DiskMetrology.subpixelLimb ─► full-circle centre + radius
CMDeviceMotion.gravity ──────┘                                        │
                                                                      ▼
MoonEphemeris.position(t) ─► dec, GHA, distance ─► angularRadius ─► arcsec/pixel (plate scale)
                                                                      │
        gravity + Moon-centre pixel offset + plate scale ─► CelestialMath.altitudeSightDeg
                                                                      │
                          PositionFix.solve (2-DOF Gauss-Newton: one altitude LOP + prior)
                                                                      ▼
                                          lat, lon, ENU covariance, 1σ major axis
```

The plate scale comes from the Moon's **known angular size** ÷ the fitted radius
in pixels — so no camera-intrinsics calibration is needed. The altitude is the
IMU boresight elevation plus the Moon centre's vertical pixel offset from the
principal point.

## Simulation module → MVP file map

| Validated Python (`imu_fusion/`)                     | MVP file                        | Status in MVP |
|------------------------------------------------------|---------------------------------|---------------|
| `disk_metrology.py` (sub-pixel limb, full circle)    | `Sources/DiskMetrology.swift`   | ported (CPU; gradient refine only, threshold seed) |
| `astro.py` predicted altitude/azimuth                | `Sources/CelestialMath.swift`   | ported |
| `corrections.py` refraction + parallax + SD          | `Sources/CelestialMath.swift`   | ported, golden-checked |
| `astro.py` / `starfix` ephemeris                     | `Sources/MoonEphemeris.swift`   | **truncated** stand-in (few-arcmin, geocentric) |
| `celestial_factor_graph.py` intercept-method LOP     | `Sources/PositionFix.swift`     | 2-DOF hand-rolled GN (GTSAM replaces at scale) |
| `iphone_model.py` synthetic horizon (gravity)        | `CelestialMath.altitudeSightDeg`, `CaptureController.gravity` | ported |
| capture / preview                                    | `Sources/CaptureController.swift`, `ContentView.swift` | tele single-frame |
| — glue / view model —                                | `Sources/FixViewModel.swift`    | new |
| `optical_attitude.py` differential Sun−Moon Δq       | —                               | **out of scope** (full app) |
| `ultrawide_horizon.py` optical horizon               | —                               | **out of scope** |
| `lunar_orientation.py` Metal rotational-NCC roll     | —                               | **out of scope** |
| `realtime.py` IncrementalFixedLagSmoother (GTSAM)    | —                               | **out of scope** |

## What "done" looks like for this MVP

A green ⌘U test run, then on-device: aim at the Moon at night, tap Capture, and
see a pin drop within a few km of where you're standing (bounded by the prior
and the truncated ephemeris). That confirms the chain end-to-end; everything
after that is the architecture spec's build order.
