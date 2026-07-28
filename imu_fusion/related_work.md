<!--- © 2026.  MIT License (see LICENSE file). -->

# Related work: terrain-aided navigation, factor graphs and DEMs

Prior art for `terrain_resection.py` / `terrain_factors.py` — what exists, what
is safe to reuse, and what is deliberately not reused.

**Licences below were checked by fetching the licence file from each repository
on 2026-07-28.** A repository with *no licence file* grants **no rights** — it is
reading material only, however permissive it looks. Re-check before depending on
anything here, and note the difference between *reading an algorithm* and
*linking code into a shipped app*.

---

## The gap

**No public repository does what this one does**: fuse a photographed **horizon
silhouette** matched against SRTM into a **factor graph** alongside celestial
sights and IMU.

The adjacent fields each have one half of it:

| Field | Observable | Estimator | Has our combination? |
|---|---|---|---|
| Terrain-referenced nav (TERCOM/TRN) | elevation strip **under** the vehicle | correlation + KF/ESKF | DEM yes, factor graph no |
| Bathymetric SLAM | sonar **submaps** of the seabed | **pose graph / factor graph** | graph yes, horizon no |
| Visual geo-localisation | **skyline** from a photo | CNN retrieval / matching | horizon yes, graph no |
| DEM-anchored visual SLAM | VO + DEM surface | **factor graph anchoring** | closest in spirit |

The important distinction: TERCOM matches a **1-D profile flown over**;
we match a **2-D silhouette seen from a fixed point**. Different geometry,
different failure modes — the TERCOM literature's accuracy figures do not
transfer.

---

## Usable — permissive licences, verified

| Repository | What it is | Licence |
|---|---|---|
| [ignaciotb/bathymetric_slam](https://github.com/ignaciotb/bathymetric_slam) | Bathymetric graph SLAM for AUVs: submap registration → pose-graph optimisation | **BSD-3-Clause** |
| [borglab/gtsam](https://github.com/borglab/gtsam) | The factor-graph solver this project already uses | **BSD** |
| [lewisgibson/tercom-missile-guidance](https://github.com/lewisgibson/tercom-missile-guidance) | C++20 terrain-referenced-navigation library + SITL simulator | **MIT** |
| [jblindsay/whitebox-tools](https://github.com/jblindsay/whitebox-tools) | Rust GIS toolbox incl. viewshed / line-of-sight (see `hw_acceleration.md`) | **MIT** |

### The one worth studying closely

**`bathymetric_slam` is the closest mature analogue to our problem.** Underwater
terrain-aided navigation solved, years ago, the thing we are doing: sense local
terrain, match it to a map, turn the match into a graph constraint. Its
architecture is the transferable part:

> It does **not** try to localise absolutely against the map at every step.
> It accumulates local **submaps**, registers map-to-map, and feeds each
> registration into the pose graph as a **relative constraint**.

Applied here that argues for treating a **stitched panorama as a submap** and
the skyline match as a **constraint with a covariance**, rather than as a
per-frame absolute fix — which is exactly the direction the
`terrain_resection` → `terrain_factors` split already takes (search to
*identify*, graph to *estimate*), and an argument for pushing further that way.

---

## Reference only — GPL or unlicensed

Read the algorithms; do not link the code.

| Repository | Licence | Note |
|---|---|---|
| [rpng/open_vins](https://github.com/rpng/open_vins) | **GPL-3.0** | Excellent factor-graph / MSCKF engineering; GPL blocks a closed app |
| [OSGeo/grass](https://github.com/OSGeo/grass) `r.viewshed` | **GPL-2.0** | The reference IO-efficient viewshed |
| [mzahana/tercom_nav](https://github.com/mzahana/tercom_nav) | **no licence file** | TERCOM + ESKF for fixed-wing UAVs |
| [alti3/python-tercom](https://github.com/alti3/python-tercom) | **no licence file** | Readable Python TERCOM demo |
| [YFS90/GNSS-Denied-UAV-Geolocalization](https://github.com/YFS90/GNSS-Denied-UAV-Geolocalization) | **no licence file** | Terrain-weighted constraint optimisation |
| [ayushmankumar7/TERCOM-python](https://github.com/ayushmankumar7/TERCOM-python) | **no licence file** | TERCOM + DSMAC via OpenCV |
| [smarc-project/UWExploration](https://github.com/smarc-project/UWExploration) | **no licence file** | Underwater exploration stack |
| [tombh/total-viewsheds](https://github.com/tombh/total-viewsheds) | **no licence file** | Cache-efficient total viewshed, CPU SIMD + GPU |

Four of five TERCOM repositories found have no licence at all.

---

## Papers

**Factor graphs for navigation**
- Taylor & Gross, *Factor Graphs for Navigation Applications: A Tutorial*,
  NAVIGATION 71(3), 2024 — the canonical treatment; factor graphs as a
  generalisation of the Kalman filter, with the sparse-linear-algebra view.
  <https://navi.ion.org/content/71/3/navi.653>

**DEM as a graph constraint**
- *Visual SLAM with DEM Anchoring for Lunar Surface Navigation*,
  arXiv 2603.17229 — augments a pose graph with **anchoring factors derived
  from a reference DEM** to bound VO drift. Structurally the same idea as
  `terrain_factors.py`; worth reading for how the anchor is weighted against
  odometry.

**Bathymetric factor-graph SLAM**
- Bichucher et al., *Bathymetric factor graph SLAM with sparse point cloud
  alignment*, OCEANS 2015. <https://ieeexplore.ieee.org/document/7404433/>
- Torroba et al., *Active Bathymetric SLAM for autonomous underwater
  exploration*, Ocean Engineering, 2023.
- VanMiddlesworth et al., *Mapping 3D Underwater Environments with Smoothed
  Submaps*, FSR 2013.
  <https://www.cs.cmu.edu/~kaess/pub/VanMiddlesworth13fsr.pdf>

**Terrain-referenced navigation**
- *An advanced ESKF-based terrain contour matching (TERCOM) method for tracking
  an aerial vehicle using a low-cost digital elevation map*, 2025.
- *A Robust Terrain Aided Navigation Using the Rao-Blackwellized Particle Filter
  Trained by LSTM Networks*, Sensors, 2018.

**Skyline / horizon geo-localisation** (see also `hw_acceleration.md`)
- Baatz et al., *Large Scale Visual Geo-Localization of Images in Mountainous
  Terrain*, ECCV 2012.
- Liu et al., *CMLocate: A cross-modal automatic visual geo-localization
  framework for a natural environment without GNSS information*, IET Image
  Processing, 2023 — DEM-rendered **skyline database** + semantic segmentation
  + CNN matching; reports **49 m** average positioning error over a 203 km² site.
- *A New Method of Improving the Azimuth in Mountainous Terrain by Skyline
  Matching*, PFG, 2020 — consistent with our finding that skyline matching
  recovers **azimuth** far better than position.

---

## What this project takes, and what it does not

**Taken**
- GTSAM as the estimator (BSD), same as the celestial stack.
- The bathymetric-SLAM idea of *match → constraint with covariance* rather than
  *match → absolute fix*, reflected in the `terrain_resection` (identify) →
  `terrain_factors` (estimate) split.
- The skyline-database idea from the geo-localisation literature, recorded as a
  future optimisation in `hw_acceleration.md` §2.4.

**Not taken**
- **TERCOM implementations** — the licences are unusable *and* the observable is
  wrong for us (profile-under-vehicle, not horizon-from-a-point).
- **Published accuracy figures** — the 49 m of CMLocate comes with a rendered
  skyline database, learned segmentation and a specific test site. Our own
  measured numbers are in `terrain_resection.py`'s docstring and are what this
  project claims.
- **GPL code**, on the assumption the iOS app may ship closed. If it ships
  open-source that constraint relaxes and `r.viewshed`/OpenVINS become options.
