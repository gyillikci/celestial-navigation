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

---

## Video pans and the vertical: what I should have read first

Added after attempting a resection from a panning iPhone video (Bodrum, 559
frames, 93° of yaw).  Two of the three things I derived from scratch are
textbook, and the third is a *known* degeneracy I walked straight into.  Recorded
here so the next attempt starts from the literature.

**Caveat on sourcing.** Every full text below returned HTTP 403 through this
session's egress proxy, so these are summaries from search results, not from
reading the primary papers.  Verify before relying on a number.

### 1. Focal length from a pan is Hartley 1994, not a new idea

- Hartley, *Self-Calibration from Multiple Views with a Rotating Camera*, ECCV
  1994. <https://users.cecs.anu.edu.au/~hartley/Papers/calibration/eccv94/calib.pdf>
  Calibration from point matches alone across ≥3 images taken from one point with
  different orientations; no knowledge of the orientations required.

What I did — measuring the `1 + (x/f)²` variation in across-frame image motion —
is a special case of this, and it worked (f to 0.4% of the 4× nominal).  The part
I should have known in advance: **rotation about a single axis is a recognised
"critical motion"**, degenerate for full calibration, though the focal length
specifically can still be recovered.  Kahl et al. and Zisserman et al. refine
when single-axis rotation is and is not critical.

That is precisely the outcome I measured empirically and then explained as if it
were a discovery: the pan gave me the scale and took away the vertical.

### 2. The skyline literature does NOT globally align a curve

- Baatz, Saurer et al., *Large Scale Visual Geo-Localization of Images in
  Mountainous Terrain*, ECCV 2012.
  <https://link.springer.com/chapter/10.1007/978-3-642-33709-3_37>
- Saurer, Baatz et al., *Image Based Geo-localization in the Alps*, IJCV.

Method: a **bag of curvelets** — *local* skyline shape descriptors aggregated
across the whole skyline and matched against a large database of panoramic
skylines rendered offline from DEMs, with each descriptor's viewing direction
stored for on-the-fly geometric verification in an inverted-file search.  No
prior camera position or field of view needed.  Reported **88% of 200+ images
localised within 1 km**.

Two corrections to how this project has been working:

- **Local descriptors, not a global fit.**  My solver aligns one long curve with
  two global nuisance offsets.  That is exactly the estimator that collapses when
  the panorama carries slow drift — which is what happened at Bodrum.  A bag of
  local curvelets is invariant to drift that a global alignment cannot survive.
- **1 km is the field's success criterion**, not 100 m.  The Tahoe result
  (≈1 km) is therefore *at* the state of the art rather than disappointing, and I
  should stop treating a kilometre as a failure.

### 3. Roll, tilt and FOV are SAMPLED in the literature, not solved

The profile-matching work sweeps camera roll explicitly (e.g. −6° to +6°) and
samples tilt and field of view, reporting that roll has a significant effect on
match quality and that sampling it recovers images the baseline cannot localise.

This project fits those as free offsets instead, which is what let them absorb
the signal — the same nuisance-absorption failure documented throughout
`RESULTS.md`, arrived at independently three times.

### 4. The vertical comes from the horizon — this is a solved sub-problem

- Grelsson & Felsberg, *Highly Accurate Attitude Estimation via Horizon
  Detection*.  Canny plus probabilistic Hough voting for an approximate attitude
  and horizon line, then refinement by **registering the extracted edge pixels
  against the geometrical horizon from SRTM3 at an approximate position**.
- Grelsson et al., *GPS-level accurate camera localization with HorizonNet*,
  Journal of Field Robotics, 2020.
  <https://onlinelibrary.wiley.com/doi/abs/10.1002/rob.21929>
- Dumble & Gibbens, *Horizon Profile Detection for Attitude Determination* —
  extracts the actual horizon *profile shape* for visual attitude determination.
- Carrio et al., *Attitude estimation using horizon detection in thermal images*,
  IJARS 2018.

This is the direct answer to the Bodrum failure.  The pan cannot supply its own
vertical; the *horizon* supplies it, and the sea horizon is plainly visible in
that very clip.  The published order is **attitude from the horizon first, then
position** — not position and attitude jointly from a drifting panorama.

### 5. Skyline extraction: continuity, not per-column independence

- *Comparison of Semantic Segmentation Approaches for Horizon/Sky Line
  Detection*, arXiv:1805.08105. <https://arxiv.org/abs/1805.08105>
- Survey: horizon detection for maritime video surveillance.
  <https://ouci.dntb.gov.ua/en/works/7WGXR6O7/>

The classical baselines enforce **continuity along the skyline** (shortest-path /
dynamic programming across columns) rather than deciding each column
independently.  `skyline_extract` decides per column and then repairs with
`drop_straight_runs` / `drop_outliers`; the mislocks onto rooflines and the
waterline in the Bodrum frames are the predictable failure of that design, and a
DP formulation would be the principled fix.

### What to do differently next time

1. Take attitude from the horizon **before** attempting position (Grelsson).
2. Match **local curvelets**, not a globally aligned curve (Baatz).
3. **Sample** roll/tilt rather than fitting them as free offsets.
4. Extract the skyline with a **continuity-constrained** method, not per column.
5. Judge success against the field's **1 km** benchmark.
