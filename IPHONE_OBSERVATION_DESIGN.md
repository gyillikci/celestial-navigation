<!---
    Design document for an iPhone-based celestial observation front-end
    feeding the celestial-navigation (starfix.py) reduction engine.
-->

# iPhone Celestial Observation — Design & Accuracy Study

A design study for an **iPhone 17 Pro** app that observes the **Sun and Moon**
through the built-in telephoto camera plus a **Reeflex 300–600 mm** clip-on
lens, extracts precise astronomical observables (altitude, orientation /
P-angle, bright-limb angle, phase), and hands them to the existing
[`starfix.py`](./starfix.py) sight-reduction engine to produce a position fix
with an error ellipse.

> **Scope.** This document covers the observation front-end and its error
> budget. The *reduction* half (lines of position, intersections, refraction,
> parallax, WGS-84 handling, Monte-Carlo error simulation) already exists in
> this repository and is reused, not redesigned.

## Table of Contents

1. [The one thing that determines accuracy](#the-limit)
2. [What a "sighting" is](#sighting)
3. [Optical chain and plate scale](#optics)
4. [Error budget](#budget)
5. [How to increase accuracy on iPhone 17 Pro](#improve)
6. [Artificial horizon — the key unlock](#artificial-horizon)
7. [Sun vs Moon specifics](#bodies)
8. [App architecture and data flow](#architecture)
9. [Reuse map into starfix.py](#reuse)
10. [Things easily missed](#missed)
11. [Suggested build order](#roadmap)

---

## 1. The one thing that determines accuracy <a name="the-limit"></a>

**Every optical measurement in this system — limb fitting, sunspot P-angle,
Moon crater/texture matching, phase — is arcsecond-class and is _not_ the
bottleneck. The bottleneck is the vertical reference: how well the phone knows
"straight down" at the instant of the shot.**

The reason is the fundamental conversion factor of celestial navigation:

```
1 arcminute of altitude error  =  1 nautical mile (1.85 km) of position error
```

A handheld iPhone derives its gravity vector from the accelerometer. In motion
that vector is good to perhaps ±6–18 arcmin, i.e. **±11–33 km**. The sunspot
method measures orientation to **±1 arcsec** — roughly 400× better — but
orientation mainly constrains *longitude via time*, while *altitude* constrains
the position circle, and altitude is IMU-limited.

> **Consequence for effort allocation:** do not keep refining limb detection
> past arcsecond precision. Spend the effort on the vertical
> (see [§6](#artificial-horizon) and [§5](#improve)). Until the vertical is
> solved, better optics buy nothing.

This is consistent with the toolkit's own Monte-Carlo result
([`starfixdata_stat_1_mc.py`](./starfixdata_stat_1_mc.py)):
*2 arcmin altitude + 2 s time → ~5 km position sigma.*

---

## 2. What a "sighting" is <a name="sighting"></a>

One sighting is a single timestamped observation that reduces to one line of
position. Each record carries:

| Field | Source | Achievable precision |
|---|---|---|
| UTC timestamp | GPS/NTP time + calibrated capture latency | ±10–50 ms (±0.2–0.8″ of hour angle) |
| Body | user / ephemeris | — |
| Disk-center pixel | **full-limb ellipse fit** (not edge sighting) | ±0.1–0.3″ |
| Plate scale (″/px) | **self-calibrated from known semidiameter** each frame | ~0.1 % |
| Camera attitude (quaternion) + gravity vector | CoreMotion `deviceMotion` | ±6–18 arcmin handheld — *the limiter* |
| Roll vs celestial north | sunspot match (Sun) / crater-texture match (Moon) | ±1″ |
| Bright-limb PA (χ), phase | ellipse + terminator fit vs ephemeris | arcsec-class |
| Temp, pressure | sensor / manual | for Bennet refraction |
| Assumed position | last fix / GPS-denied seed | — |

The corrected altitude `Ho` and a coarse azimuth then feed
`Sight` → `SightCollection` → intersection + Monte-Carlo ellipse.

**Advantage over a sextant:** because the whole disk is fit, the geometric
center is obtained directly — **no semidiameter correction, no index error** —
and the known angular diameter becomes a ruler that **self-calibrates the plate
scale in every frame**.

---

## 3. Optical chain and plate scale <a name="optics"></a>

Working numbers for a ~600 mm-equivalent configuration on the phone sensor:

- **Field of view** ≈ 3.4° × 2.3° at 600 mm-equiv (≈ 6.8° at 300 mm-equiv).
  The Sun/Moon disk is ~0.5° (1800″), so it comfortably fits with margin at
  both ends of the zoom.
- **Plate scale** ≈ 1.5″/px with the full disk spanning ~1200 px on an
  ~8000-px-wide sensor. Comparable to the Dwarf 2 reference (2.8″/px) used in
  the sunspot accuracy study.
- **Diffraction limit** of a ~75 mm aperture (600 mm at ~f/8 catadioptric):
  1.22·λ/D ≈ **1.85″** — well matched to the plate scale.
- **Practical resolution** is seeing-limited in daytime to ~2–5 arcsec, so the
  system is slightly oversampled — good for sub-pixel centroiding.

The Reeflex is a low-cost clip-on: expect distortion, decentering, vignetting,
and chromatic aberration. Mitigations: keep the disk near frame center,
self-calibrate scale from the semidiameter, and fit distortion once against a
star field.

---

## 4. Error budget <a name="budget"></a>

Ordered by impact:

1. **Vertical / IMU tilt — DOMINANT.**
   Handheld ±6–18 arcmin → **±11–33 km**.
   Braced + burst-averaged ±1–3 arcmin → **±2–6 km** (matches the 5 km MC result).
   Star-calibrated ±5–15 arcsec → **±0.15–0.5 km**.
2. **Refraction.** ~1 arcmin at 45° altitude, but grows and gets uncertain by
   several arcmin below ~15°. **Stay above ~20–30°.** Model with Bennet
   (already in `get_refraction`) using real temp/pressure.
3. **Timing / capture latency.** GPS time is ~ms, but the camera pipeline
   (rolling shutter + buffering) can add 10–100 ms if uncharacterized →
   up to ±1.5″. Calibrate once.
4. **Plate-scale & distortion of the Reeflex.** Arcsec-class if the disk drifts
   off-center; controlled by self-calibration + one-time distortion fit.
5. **Atmospheric dispersion** (blue above red by arcsec at low altitude) biases
   the limb. Use one color channel (green) or model it. Another reason to stay high.
6. **Parallax.** Moon up to ~1°(!), Sun ~8.8″. Already handled by
   `get_vertical_parallax` — just wire it in.
7. **Limb/feature centroiding** — sub-arcsec, effectively negligible. This is
   the part it is easy to *over*-engineer.
8. **Gyro drift** — only matters for sequential inter-body pointing
   (~0.01–0.1°/slew).
9. **Magnetometer heading** — ±1–2°, useless for precision. Azimuth should come
   from the sky solution, not the compass.

**Realistic outcomes**

| Setup | Expected position error |
|---|---|
| Naive handheld | ±10–20 km |
| Braced + averaged + limb-fit + good refraction | ±3–6 km |
| Star-calibrated IMU on a tripod | ±0.3–1 km |
| + multi-sight geometry (removes single-sight ambiguity) | ±0.5–2 km, better constrained |

---

## 5. How to increase accuracy on iPhone 17 Pro <a name="improve"></a>

### A. Beat the vertical (biggest wins — do these first)

1. **Artificial horizon** — see [§6](#artificial-horizon). The cleanest path
   from km-class to ~100 m-class.
2. **Night star-field calibration** of the accelerometer bias, carried into
   daytime. Image stars of known altitude, solve boresight + gravity offset,
   store per-orientation bias. (~±5″ path.)
3. **Stillness-gated burst averaging.** Brace on tripod/beanbag; average the
   accelerometer over 2–5 s only when |a|≈g and |gyro|≈0. Cuts noise by √N.
4. **Two-position (face-flip) sights.** Measure, rotate 180° about the
   boresight, measure again; averaging cancels mounting/bias offset the way
   reversing a real instrument does.

### B. Exploit the optics (already most of the way there)

5. **Full-limb ellipse fit** for disk center + self-calibrated plate scale.
6. **Sunspot differential-rotation / template match** for celestial-north roll
   (verified ±1″). Extend the same idea to the **Moon** via crater matching
   against a libration-correct rendered Moon (reuse the Stellarium P/B₀/L₀
   texture-mapping code; add lunar libration + phase + bright-limb PA χ).
7. **Fit the bright limb only** on the Moon — the terminator is fuzzy and the
   limb is mountainous; the illuminated limb is the clean arc, and χ is a
   strong ephemeris-comparable orientation observable.
8. **Lucky imaging** in continuous mode: stack and select the sharpest frames
   to fight daytime seeing.

### C. Geometry & procedure

9. **Stay above 20–30° altitude** (refraction + dispersion).
10. **Multi-sight**: several Sun sights across the day for good azimuth crossing
    angles, or Sun + daytime Moon simultaneously. `SightCollection` already
    intersects these; the MC script gives the ellipse.
11. **Characterize capture-to-timestamp latency once** (photograph a
    GPS-disciplined clock, or measure the AVFoundation presentation-timestamp
    offset) and log per-frame UTC.

---

## 6. Artificial horizon — the key unlock <a name="artificial-horizon"></a>

Image the Sun/Moon **and its reflection** in a tray of still water or oil. The
angle between the direct and reflected images equals **2× the true altitude,
measured purely optically in the image plane — it bypasses the IMU entirely**.

```
altitude = (angular separation between direct and reflected disk) / 2
```

This is the classic aviation bubble-sextant / mercury-horizon technique. On this
hardware it is the single most effective route from kilometre-class to
~100 m-class accuracy, because it removes the dominant error term
([§1](#the-limit)) rather than merely shrinking it.

With a telephoto the direct and reflected images are captured either in one
wide framing or sequentially with gyro transfer between the two pointings. Both
disks are then found by the same limb-fit used everywhere else, so the extra
machinery is small.

---

## 7. Sun vs Moon specifics <a name="bodies"></a>

**Sun**
- Requires a **certified objective solar filter** in front of the Reeflex
  (see [§10](#missed)). White-light shows sunspots; H-alpha shows more
  features but changes the limb.
- Orientation from **sunspot differential rotation** (±1″, verified).
- Clean full disk → best limb fit and plate-scale self-calibration.

**Moon**
- Orientation from **crater/texture matching** against a rendered Moon that
  **includes libration**, phase, and bright-limb PA χ.
- Fit the **bright limb only**; the terminator side is not a clean circle.
- **Parallax up to ~1°** must be corrected (`get_vertical_parallax`).
- Usable in daylight alongside the Sun for crossing geometry.

Both bodies reduce through the same pipeline; only the orientation observable
and the limb-fit masking differ.

---

## 8. App architecture and data flow <a name="architecture"></a>

```
 ┌──────────────── iPhone 17 Pro (Swift) ────────────────┐
 │  AVFoundation  ── ProRAW/RAW frames + CMSampleBuffer   │
 │                    presentation timestamps             │
 │  CoreMotion    ── gravity vector + attitude @100 Hz    │
 │  CLLocation    ── GPS time (truth ref / GPS-denied off)│
 │        │                                               │
 │        ▼                                               │
 │  Metal/OpenCV  ── limb ellipse fit, feature match,     │
 │                    synthetic Sun/Moon render           │
 │        │                                               │
 │        ▼    reduced observables (alt, az, time,        │
 │             P-angle, χ, phase, temp/pressure)          │
 └────────┼──────────────────────────────────────────────┘
          ▼
   starfix.py reduction  (Sight → SightCollection →
   get_intersections → Monte-Carlo error ellipse)
          ▼
   position fix + error ellipse + logged RAW/IMU for reprocessing
```

- **Capture:** AVFoundation ProRAW/RAW with `CMSampleBuffer` presentation
  timestamps; CoreMotion `deviceMotion` at 100 Hz for gravity + attitude
  quaternion; `CLLocation` only for GPS *time* (and as ground truth to
  validate, or disabled to simulate a GPS-denied fix).
- **Vision:** Metal/OpenCV for the ellipse fit and feature matching; render the
  synthetic Sun/Moon from the existing Stellarium orientation math.
- **Logging:** persist RAW frames + IMU + timestamps so every fix can be
  reprocessed and each fix ships with a Monte-Carlo error ellipse.

---

## 9. Reuse map into starfix.py <a name="reuse"></a>

| Need | Existing symbol |
|---|---|
| Line of position (small circle) | `Circle`, `get_circle_for_angle` |
| Multi-sight fix | `SightCollection`, `SightPair`, `get_intersections` |
| Refraction | `get_refraction` (Bennet) |
| Dip (only if a real horizon is used) | `get_dip_of_horizon` |
| Moon/Sun parallax | `get_vertical_parallax`, `get_geocentric_alt` |
| Datum handling | `LatLonGeodetic` / `LatLonGeocentric` (WGS-84) |
| Ephemeris | `Almanac`, `get_mr_item`, machine-readable almanac in `sample_data/` |
| Error ellipse | `starfixdata_stat_1_mc.py` (Monte-Carlo) |

> **Portability note.** The Android app embeds this Python via
> buildozer/kivy; that path does **not** work on iOS. Plan to either port the
> core reduction (`Sight`, `Circle`, `get_intersections`, refraction/parallax)
> to Swift, or run it as a paired service. This is the one place the toolkit
> does not lift-and-shift.

---

## 10. Things easily missed <a name="missed"></a>

- **Solar safety, non-negotiable:** a certified objective filter
  (Baader/ND5 or H-alpha) **in front of the Reeflex**, not a software trick —
  otherwise the sensor (and an eye through the finder) is at risk.
- **The IMU is the ceiling.** Every arcsec won optically is wasted until the
  vertical is solved. Prioritize the artificial horizon and star-cal.
- **Thermal:** long continuous capture heats the phone → read-noise up, possible
  throttling, and the clip-on lens flexes/defocuses as it warms. Refocus and
  watch for drift.
- **Manual focus is critical** and the Reeflex is manual — use live limb-edge
  sharpness (max gradient) as the focus metric.
- **Datum consistency:** use `LatLonGeodetic` vs geocentric + WGS-84 — do not mix.
- **Libration** must be in the Moon render or crater matching biases; the Moon's
  limb topography is arcsec-level, so weight the fit accordingly.
- **Extra targets:** Venus and Jupiter are reachable in twilight/daytime with
  this focal length and add crossing geometry.

---

## 11. Suggested build order <a name="roadmap"></a>

1. **Sensor-stack Monte-Carlo first.** Extend `starfixdata_stat_1_mc.py` to
   model *this* stack (IMU tilt + timing + limb + P-angle) and print the
   predicted error ellipse. This quantifies, before any app code, how much the
   artificial horizon and star-cal actually buy.
2. **Capture + logging app** (RAW + IMU + UTC), no reduction yet.
3. **Limb ellipse fit** → disk center + self-calibrated plate scale.
4. **Orientation observables** (sunspot match; Moon crater/χ match).
5. **Artificial horizon** altitude path.
6. **Wire to `starfix.py`** (ported or paired) for the fix + ellipse.
7. **Night star-field IMU calibration** and multi-sight geometry.

The step-1 Monte-Carlo is the highest-value next action: it turns the accuracy
claims in this document into concrete, hardware-specific numbers.
