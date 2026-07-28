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
| [TouqeerAhmad/skyline_detection](https://github.com/TouqeerAhmad/skyline_detection) | **BSD-like but NON-COMMERCIAL, source-form only** | Mountainous skyline extraction. Looks permissive, clause 3 forbids commercial use — see the horizon section below |

Four of five TERCOM repositories found have no licence at all.

---

## Horizon- and skyline-based navigation specifically

The closest field to what this project does, searched separately. **The
conclusion is that there is no usable open-source horizon-based *localisation*
system** — the technique exists in papers, in shipped products, and in patents,
but not in code you can build on.

### The one on-target repository, and its licence trap

[TouqeerAhmad/skyline_detection](https://github.com/TouqeerAhmad/skyline_detection)
— *Resource Efficient Mountainous Skyline Extraction using Shallow Learning*,
IJCNN 2021. Structure-tensor-selected linear filters plus dynamic programming
(shortest path through a multistage graph) to find the sky/mountain boundary.
Explicitly targeted at **resource-constrained platforms: mobile phones,
planetary rovers, UAVs** — i.e. exactly our deployment target, and exactly the
skyline-extraction stage that is the weakest link in our real-photo results.

**Its licence looks permissive and is not.** It opens
`Copyright (c) 2021 Touqeer Ahmad. All rights reserved.` and then reads like
BSD, but:

- **clause 3**: redistribution is *"permitted only for non-commercial research
  collaboration and demonstration purposes"* — **not** open source in the usual
  sense;
- it grants **source-form redistribution only** (there is no binary clause),
  which is a further problem for a compiled iOS app.

So: **read the paper, reimplement the algorithm, do not ship the code.** The
method is simple enough to reimplement from the paper in Accelerate or Metal,
which is the recommended path.

### What exists only as papers, products or patents

| Work | What it is | Available as |
|---|---|---|
| **PeakLens** | FCN skyline extraction aligned against a DEM panorama from GPS + compass — architecturally the same as ours, shipped on Android | product, no source found |
| **CMLocate** (Liu et al., IET IP 2023) | DEM-rendered skyline database + segmentation + CNN matching; **49 m** mean error over 203 km² | paper |
| Baatz et al., ECCV 2012 | Large-scale visual geo-localisation in mountainous terrain | paper |
| PFG 2020 | Improving **azimuth** by skyline matching — matches our finding that skylines fix heading far better than position | paper |
| US 9165217 / 9292766 | Ground-level photo geolocation using digital elevation | **patents** |
| US 8311285 | Localising in urban environments from omni-directional skyline images | **patent** |

The patent coverage is worth noting before this becomes a product; it is outside
what this study can assess.

### Consequence for this project

The two stages need different sourcing:

- **Skyline extraction from the photo** — reimplement the IJCNN shallow-learning
  + dynamic-programming approach, or use Vision / Core ML segmentation. This is
  the stage most likely to improve the real-photo numbers, since the panorama
  registration and extraction (not the DEM) were diagnosed as the limit in
  `terrain_resection.py`'s docstring.
- **Skyline → position** — `terrain_resection.py` + `terrain_factors.py`. There
  is nothing off-the-shelf to adopt.

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
