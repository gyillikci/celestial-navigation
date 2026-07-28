<!--- © 2026.  MIT License (see LICENSE file). -->

# Hardware-accelerating the terrain resection

How to make `terrain_resection.py` fast enough to run on an iPhone 17 Pro.

Everything marked **[measured]** was timed in this repository on **one core of a
2.1 GHz Xeon** (single-threaded Python + NumPy). Everything marked **[estimate]**
is an extrapolation and has **not** been run on device.

---

## 1. The workload, measured

One *candidate position* costs three things:

| Stage | Cost **[measured]** | Shape |
|---|---|---|
| `render_skyline`, full 360°, fine (0.05°, 50 m steps, 45 km) | **2074 ms** | 7200 rays × 897 samples = 6.46 M lookups |
| `render_skyline`, full 360°, search (0.10°, 80 m, 40 km) | **401 ms** | 3600 × 498 = 1.79 M lookups |
| `render_skyline`, 90° window (0.10°, 80 m, 40 km) | **78 ms** | 900 × 498 = 0.45 M lookups |
| `skyline_peaks` | 6.3 ms | 1-D scan |
| `best_match`, blind (55 scales × 720 headings × 2) | **1294 ms** | 79 200 hypotheses |
| `best_match`, pruned (±2° heading, known focal) | **0.5 ms** | 32 hypotheses |
| `solve_landmark_fix` (GTSAM, 7 factors + covariance) | **0.3 ms** | — |

That is **≈4.5 M DEM lookups/second** single-threaded. The ray-march is ~95 % of
a pruned candidate's cost, so **that is the only thing worth accelerating.**

Real runs over the Bitez case (51 km², 159 candidates scored): **287 s**
**[measured]**. The pruned 2.6 km² run: **17.9 s** for 64 candidates
**[measured]**.

---

## 2. Do the algorithmic work first — it beats the hardware

Ranked by measured or expected payoff. The first one is already implemented.

### 2.1 Prune the search space — **~2600× [measured]**

`resect_with_priors` already does this. A phone knows its dead-reckoned
position, its magnetometer heading and its own lens:

| Configuration | Hypotheses | `best_match` |
|---|---|---|
| Blind | 79 200 | 1294 ms |
| + magnetometer ±20° | 8 800 | 150 ms |
| + known focal length | 160 | 3.1 ms |
| + magnetometer ±2° | **32** | **0.5 ms** |

It also shrinks the *render*: with a heading you only sweep the arc the camera
saw (90° instead of 360° → **78 ms instead of 401 ms** **[measured]**).

**No GPU work competes with this.** Do it before anything else.

### 2.2 Hierarchical ray-marching (maximum mipmap) — **[estimate] 5–20×**

The inner loop marches a ray in fixed 50–80 m steps and takes the maximum
elevation angle. That is exactly height-field ray casting, and the standard
acceleration is a **maximum mipmap**: a pyramid where each texel of level *n*
holds the *maximum* elevation of its 2×2 children. A ray can then skip whole
blocks whose maximum cannot possibly raise the running horizon angle.

- Tevs, Ihrke & Seidel, *Maximum Mipmaps for Fast, Accurate, and Scalable
  Dynamic Height Field Rendering* (MPI) — the canonical reference.
- Build cost is trivial (one pass per tile, cached with the tile) and it helps
  **CPU and GPU equally**.

This is the single biggest algorithmic win left, and it is independent of Metal.

### 2.3 Early termination and monotone bounds — **[estimate] 1.5–3×**

The horizon angle is monotone non-decreasing along a ray only in the maximum;
once the remaining terrain's *maximum possible* angle (from the mipmap, or a
global max-elevation bound) falls below the running maximum, the ray can stop.
Near-field terrain frequently dominates, so many rays terminate early.

### 2.4 Precomputed horizon maps — **[estimate] large, but only if reused**

If the same area is queried repeatedly (a voyage, a mission plan), precompute a
horizon profile per grid node offline and interpolate at query time. This turns
resection into a lookup. Cost: storage. This is what the published
skyline-matching systems do — they build a *skyline database* rather than
ray-marching per candidate.

---

## 3. Hardware paths on iOS

### 3.1 Metal compute — the primary path

The ray-march is close to an ideal GPU workload: thousands of **independent**
rays, pure texture sampling, one `max` reduction per ray, no divergence beyond
early-exit.

Mapping:

- **DEM tile → `MTLTexture`**, `.r16Sint` or `.r32Float`, with a
  **max-mipmap chain** (§2.2). Hardware bilinear sampling replaces the manual
  interpolation in `DemTiles.elevation`.
- **One threadgroup per candidate position**, **one thread per azimuth ray**.
  A 90° window at 0.1° is 900 rays — a good threadgroup-per-candidate size.
- **Reduce along the ray in-thread** (a running `max`), so no cross-thread
  reduction is needed; write one float per ray.
- The whole candidate grid becomes one dispatch: 64 candidates × 900 rays =
  57 600 threads, trivial for the GPU.

Apple references: [Metal sample code](https://developer.apple.com/metal/sample-code/),
[Metal Performance Shaders](https://developer.apple.com/documentation/metalperformanceshaders),
[Metal Feature Set Tables](https://developer.apple.com/metal/Metal-Feature-Set-Tables.pdf).

**[estimate]** A modern mobile GPU sustains billions of texture samples/s
against the ~4.5 M/s measured here. Even allowing an order of magnitude for
dispatch overhead and cache behaviour, a pruned grid should drop from seconds to
**well under a second**. This is an extrapolation, not a benchmark.

### 3.2 Accelerate / SIMD — the cheap CPU win

If Metal is more integration than you want initially, the same loop vectorises
on the CPU with **`simd`** types or **vDSP**, plus `DispatchQueue.concurrentPerform`
across candidates. **[estimate]** 4–8× on a 6-core A-series, no GPU code.
This is the pragmatic first port of `render_skyline` to Swift.

### 3.3 What the Neural Engine is *not* for

The ANE accelerates fixed-function neural graphs, not irregular ray marching.
The one place it fits is **skyline extraction from the photograph** — a
semantic-segmentation model (sky vs terrain) via Core ML / Vision, which is how
the published systems get robust skylines in haze. That is a different stage
from the resection and should not be conflated with it.

---

## 4. Reference implementations and libraries

This problem is the GIS **viewshed / line-of-sight** problem. Licenses below
were checked against the repositories; **verify again before depending on any of
them**, and note the difference between *reading an algorithm* and *linking a
library into a shipped app*.

| Project | What it gives you | Licence | Usable in a closed app? |
|---|---|---|---|
| [WhiteboxTools](https://github.com/jblindsay/whitebox-tools) | Rust GIS toolbox with viewshed / LOS; clean, readable reference | **MIT** (verified) | **Yes** |
| [GRASS GIS](https://github.com/OSGeo/grass) `r.viewshed` | The reference IO-efficient viewshed implementation | **GPL-2.0** (verified) | **No** — read it, don't link it |
| [total-viewsheds / CacheTVS](https://github.com/tombh/total-viewsheds) | Cache-efficient total-viewshed, CPU SIMD + GPU kernels | **No LICENSE file found** | **No** — treat as unlicensed |
| Apple **Metal** / **MPS** | The GPU path itself | Apple SDK | Yes |
| Apple **Accelerate** / **simd** | Vectorised CPU fallback | Apple SDK | Yes |

Papers worth reading before writing the kernel:

- **R2 viewshed on CUDA** — Osterman et al., *An IO-efficient parallel
  implementation of an R2 viewshed algorithm for large terrain maps on a CUDA
  GPU*, *IJGIS* (2014). Reported speedups in the tens; the IO-efficiency work
  is the transferable part.
- **XDraw viewshed on GPU** — Zhao et al., *Parallel Computing* (2015).
- **Maximum mipmaps** — Tevs, Ihrke & Seidel (§2.2).

**Caveat on published speedups.** Figures like "28–925×" in the viewshed
literature are GPU-vs-*single-threaded-CPU* on *total-viewshed* workloads
(visibility from every cell), which is a much larger and more regular problem
than our few-hundred-candidate sparse search. Do not budget against them.

---

## 5. Recommended order of work

1. **Port `render_skyline` to Swift with `simd` + `concurrentPerform`.** Cheapest
   real speedup, no GPU code, and it establishes correctness against the Python
   reference and the golden vectors.
2. **Add the max-mipmap** (§2.2). Algorithmic, helps both CPU and GPU, and is
   independent of the port.
3. **Move the ray-march to a Metal compute kernel** (§3.1) once the CPU version
   is correct and there is a regression test to compare against.
4. **Only then** consider a precomputed skyline database (§2.4), which trades
   storage for compute and changes the data-shipping story.

Keep `terrain_resection.py` as the executable specification: any accelerated
implementation must reproduce its output on the synthetic fixtures in
`test_imu_fusion.TestSyntheticSkyline` and, when tiles are present,
`TestRealTerrainResection`.

---

## 6. What not to bother with

- **Accelerating `best_match`.** Already 0.5 ms when pruned **[measured]**.
- **Accelerating the factor-graph solve.** 0.3 ms **[measured]**.
- **A finer candidate grid instead of a finer render.** Measured on the Bitez
  case: a 200 m grid with a coarsened render scored *worse* (431 m) than a 445 m
  grid with a finer render (297 m). **Render resolution buys more than grid
  resolution** — spend the compute there.
