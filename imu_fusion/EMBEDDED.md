# Running the skyline fix on a hardware-accelerated embedded device

Every number here was measured on this repository's own kernels. Where a claim
is an estimate rather than a measurement it says so. The scripts are in the
session scratchpad (`embed_num.py`, `embed_mip.py`, `embed_mip2.py`); the
numbers they produce are reproduced below because the scripts are not committed.

## 1. Where the time actually goes

The forward model is `terrain_resection.render_skyline`. Stripped of its guards
it is four lines:

```python
A, D = meshgrid(radians(azs), ds)            # azimuth x range
H     = dem.elevation(lat + D*cos(A)*k1,     # <-- bilinear gather
                      lon + D*sin(A)*k2)
alpha = degrees((H - cam)/1000/D - D/(2*R_eff))
profile = alpha.max(axis=1)                  # <-- reduction along the ray
```

Measured on 403 200 samples from a real tile:

| stage | time | share |
|---|---|---|
| `dem.elevation` bilinear gather | 88 ms | **97 %** |
| angle arithmetic + max-reduction | 3 ms | 3 % |

Scattered gather throughput is **1.4–3.8 M samples/s** on this CPU. One photo at
the sweep's settings costs

| pass | grid | samples/candidate |
|---|---|---|
| coarse | 1200 az × 111 range | 0.13 M |
| fine | 4500 az × 298 range | 1.34 M |

which for 121 coarse + 14 fine candidates is **35.0 M DEM samples per photo** —
and at the measured gather rate that is the entire ~10 s/photo the sweep spends.
Nothing else in the pipeline is worth optimising until this is fixed. Extraction
is 118 ms and scoring is milliseconds.

**The kernel is a bilinear texture fetch feeding a max-reduction.** That is
literally what a GPU texture unit plus a wavefront reduction does in hardware,
which is why this problem is a good fit for an accelerator — but see §2, because
the obvious accelerator is the wrong one.

**Correction: "97 % gather" is a fact about numpy, not about the algorithm.**
The table above times `dem.elevation` against the angle arithmetic. Breaking
`dem.elevation` open shows that most of it is not the gather at all:

| inside `dem.elevation` (0.113 M samples, 16.2 ms) | time | share |
|---|---|---|
| float64 index arithmetic + per-tile boolean masking | 9.7 ms | **60 %** |
| 4-tap fetch + bilinear blend | 3.1–3.3 ms | **19–20 %** |
| remainder (allocation, dispatch) | ~3.2 ms | ~20 % |

So the *memory traffic* is under a fifth of the call, and the 60 % is Python-layer
overhead that does not exist in a CUDA kernel at all. Two consequences. First,
this CPU profile systematically misleads about a GPU port: the thing to size is
the memory traffic (§6), not the wall clock measured here. Second, any effort
spent making the *CPU* path faster should go at the index arithmetic first, not
the fetch — which is the opposite of what the 97 % figure suggests.

## 2. Precision: fp32 is required, fp16 is fatal

Rendering the same site in three precisions and comparing against float64:

| precision | max error | p99 error |
|---|---|---|
| fp32 | **0.0038′** | 0.0032′ |
| fp16 | **29.2′** | 20.4′ |

Residuals that decide a fix are 10–80′. An fp16 render therefore injects an
error the size of the signal it is trying to measure, and the position solution
is destroyed rather than degraded.

The cause is dynamic range, not rounding. The term `(H - cam)/d` spans roughly
four decades between the near field (`d = 0.3 km`, terrain 4 km up) and the far
field, and fp16 carries about three decimal digits.

The practical consequence is sharp: **NPU and DSP accelerators whose fast path is
fp16 or int8 cannot run this kernel as written.** A GPU running fp32 texture math
(Metal, Vulkan, CUDA) can. If an fp16 unit is the only accelerator available the
kernel must first be reformulated to bound its dynamic range — marching in a
normalised height/range coordinate, or splitting near and far field with separate
scalings — and that reformulation must be re-validated against the table above,
not assumed.

Related, and in the same direction: SRTM's int16 metre quantisation is worth
**3.44′ at 1 km** but only **0.076′ at 45 km**. The near field is where both
precision problems live, because both are the `1/d` factor.

## 3. What does NOT work: decimating the DEM

At 0.08° azimuth step adjacent rays are 1.4 m apart at 1 km and 63 m apart at
45 km, so the far field looks massively oversampled against SRTM's 30 m posting.
Reading a decimated DEM beyond 8 km should therefore be free. It is not:

| site | subsample 2× | subsample 4× | max-pool 2× | max-pool 4× |
|---|---|---|---|---|
| 46.55 N | 0.00′ | 0.00′ | 0.00′ | 0.00′ |
| 46.02 N | 29.4′ | 49.9′ | 32.9′ | 59.3′ |
| 47.05 N | 26.6′ | 18.8′ | 28.0′ | 46.9′ |

Tens of arcminutes — again the size of the signal. Max-pooling does not rescue
it: it removes the *loss* of crests at the cost of *raising* them, trading one
bias for another.

The reason is structural. A skyline is not an average of terrain, it is an
extremum over a thin set — one ridge crest per azimuth. Decimation is a
low-pass operation and the quantity being computed is not band-limited. The
first site scores 0.00′ only because its skyline is formed entirely inside
0.3 km, so the far-field LOD never engages.

**Keep the DEM at full 1-arcsec posting at every range.** This is a hard
memory-footprint constraint, not a tunable.

## 4. What DOES work: conservative early-out (exact, 68–99 % saved)

The profile is `max_d alpha(d)`. If the tallest terrain remaining beyond range
`d` cannot beat the best angle already found, the rest of that ray is dead. With
`suf[i] = max_{j>=i} H[j]`, an upper bound on everything from `d_i` outward is

```
opt[i] = (suf[i] - cam)/d_i      - d_i/(2 R_eff)     when suf[i] > cam
       = (suf[i] - cam)/d_max    - d_i/(2 R_eff)     otherwise
```

taking the nearest range for a positive height rise and the farthest for a
negative one, so it bounds from above in both cases. Measured:

| site | range march skipped | mean stop | skyline formed at (p50/p90/max) | profile error |
|---|---|---|---|---|
| 46.55 N | **99 %** | 0.6 km | 0.3 / 0.3 / 0.3 km | **0.0000′** |
| 46.02 N | **79 %** | 9.9 km | 5.2 / 11.8 / 45.0 km | **0.0000′** |
| 47.05 N | **68 %** | 14.8 km | 6.4 / 19.2 / 42.4 km | **0.0000′** |

Exact, not approximate: the truncated render reproduces the full render bit for
bit while touching a third to a hundredth of the samples. It works because the
skyline is usually formed in the near field (p50 = 0.3–6.4 km) even though it
occasionally reaches to 45 km — so a fixed `d_max` cut would be wrong, and an
adaptive bound is not.

That table is what a *perfect* bound achieves, computed from the finished ray.
A real implementation has to get the bound from somewhere cheaper, so this is
now implemented and measured: `render_skyline_early` builds a **max-pyramid**
over each DEM tile (`DemTiles.max_tile`, ~113 × 113 floats against 25.9 MB for
the tile) and marches in range blocks, dropping each ray as its bound falls
below its own best-so-far.

| case | full | early | speedup | max error |
|---|---|---|---|---|
| alpine high 46.55 N | 0.379 s | 0.049 s | **7.7×** | 0.0000′ |
| alpine 46.02 N | 0.407 s | 0.115 s | **3.6×** | 0.0000′ |
| valley 47.05 N | 0.421 s | 0.161 s | **2.6×** | 0.0000′ |
| low camera 46.8 N | 0.391 s | 0.052 s | **7.5×** | 0.0000′ |
| lake, water clamped 47.0 N | 0.326 s | 0.027 s | **11.9×** | 0.0000′ |

Bit-identical output, 2.6–11.9× less work. The pyramid bound is looser than the
exact ray suffix — a square block is not a ray — which is why the realised
speedup is below the 68–99 % ceiling above.

**Counting the probes, which a first version of this section did not.** The
early-out is not free: it reads a bound table before it can skip anything, and
those reads are DEM lookups too. Per photo (121 coarse + 14 fine candidates):

| | bilinear gathers (26 MB tile) | pyramid probes (49 KB) |
|---|---|---|
| plain march | 35.1 M | 0 |
| early-out, probe stride 1 | 10.5 M | **31.2 M** |
| early-out, probe stride capped by geometry | 10.7 M | **12.2 M** |

At stride 1 the probes cost three times the march they save, and by raw lookup
count the early-out is a net **loss** (0.84×). It still won on wall clock only
because a probe is a single nearest-neighbour read of a 49 KB array while a
gather is four taps into 26 MB — measured 7.18 vs 3.89 M lookups/s, i.e. only
1.85× cheaper, not free.

The waste was that the probe table was built at the march's own azimuth step. At
45 km adjacent rays are 63 m apart and a pyramid cell is 0.7–1.0 km, so the bound
was being sampled 10–15× more finely than it can possibly resolve. Probing every
`stride`-th azimuth fixes it — but `stride` **cannot be a free parameter**. It is
capped so the lateral gap between neighbouring probe rays at `d_max` stays under
half a pyramid cell; past that, terrain can sit between two probe rays in a cell
neither touches. Removing the cap costs 144.6′. Note the cap runs opposite to
intuition: a *coarse* azimuth step forces a *smaller* stride (2 for the coarse
pass, 7 for the fine one), because the gap grows with the step.

For an accelerator the distinction matters more than on CPU. The probes touch
49 KB — L1/L2-resident, effectively zero DRAM traffic — so the number that sizes
memory bandwidth is the **10.7 M bilinear gathers**, not the 22.9 M total.

Note this is the same max-pooling that failed in §3. As a **bound** it is sound;
as a **substitute for the data** it is not. The two uses are not in tension:
§3 replaced the DEM with the pool and lost crests, §4 keeps the DEM and uses the
pool only to prove that a stretch of ray cannot matter.

**The subtlety that cost two wrong answers.** A first implementation was
*almost* exact — it failed on 2 of 7 sites by 1.2′ and 11.1′. The cause is
geometric: a ray can clip the CORNER of a pyramid cell between two probe points,
so an undilated bound never sees the peak inside it. The fix is to dilate the
pyramid 3 × 3, so any cell the ray touches is covered by a cell a probe lands in.
Anyone porting this to an accelerator's native mip chain has to reproduce that
dilation — and must build the chain with **max**, never the default average,
since a bound that can sit below the terrain is not a bound.

Two lessons for the test suite came out of it. Smooth analytic terrain
(`SyntheticDem`) passes with the dilation removed, so it cannot defend this
code; the regression test writes a small rough `.hgt` and reads it through
`DemTiles`, where removing the dilation breaks 6 of 8 configurations. And of the
two defensive fixes made at the time, only the dilation was load-bearing — the
probe-index off-by-one was correct reasoning about a gap that the actual
probe/block alignment never opened. It is kept because other settings would open
it, but it fixed nothing observed.

## 5. Memory

SRTM 1-arcsec is 3601 × 3601 int16 = **25.9 MB per 1° tile**; this repo holds 54
tiles = 1.4 GB. But the working set for one fix is bounded: a 25 km² search box
plus a 45 km march radius spans about 95 × 95 km, so **4 tiles ≈ 104 MB** covers
any single photo, and §3 says it cannot be compressed by decimation.

Options, in order of preference:

1. **Stream by tile with an LRU.** Access is a radial fan from a known centre,
   so the tile set is predictable before the march starts — prefetch is easy and
   there are no surprises mid-kernel.
2. **Store as int16, convert in the sampler.** Already the on-disk format; do
   not expand to float in memory, that is a free 2× loss.
3. **Bake a skyline database offshore** (the Baatz/Saurer route). This is the
   only option that removes the DEM from the device entirely, and it converts
   layer 4 from compute into lookup. It costs offline precompute over the whole
   operating area and is the right answer for a fixed corridor, the wrong one
   for go-anywhere.

## 6. Concurrency

The search is embarrassingly parallel over `(candidate, azimuth)` — 121 × 1200 =
145 000 independent rays in the coarse pass, each ending in a private
max-reduction with no cross-ray communication. There is no reduction tree, no
atomics, no synchronisation until the per-candidate residual at the very end.

Occupancy is therefore not the problem; **memory divergence is**.

**Order the march so a warp holds ADJACENT AZIMUTHS AT ONE RANGE**, not adjacent
ranges on one ray. An earlier revision of this document said the opposite; it was
reasoned, not measured, and the measurement reverses it. Counting the distinct
memory granules a 32-sample warp touches over all four bilinear taps, on the real
fine-pass index arrays:

| ordering | 32 B sectors / warp (360° mean) | ≈ bytes/sample |
|---|---|---|
| ray-major (adjacent ranges on one ray) | 70.0 | ~70 B |
| **range-major (adjacent azimuths, one range)** | **29.2** | **~29 B** |
| address-sorted, unreachable ideal | 7.2 | ~7 B |
| uniformly shuffled | 72.0 | ~72 B |

Range-major wins at all 19 azimuths sampled — 4.72× at due N/S, 1.54× at due
E/W, **2.40× averaged over a full fan**. The damning row is the last one:
ray-major at az 0° and 45° scores 72.00 against a random shuffle's 71.97, so the
ordering previously recommended here is indistinguishable from random access.

The mechanism is the one the old text had inverted. Azimuth neighbours
**converge** toward the camera — their spacing is `d · az_step`, 1.4 m at 1 km
and 63 m at 45 km — while range neighbours stay a fixed 150 m apart at every
range. So range-major locality is best exactly where the early-out (§4) keeps
most of the samples: the near field. Per range band, due north, 32 B sectors per
warp: 0.1–1 km 13.9 ray / **2.8** range; 3–8 km 72.0 / **6.2**; 8–20 km 72.0 /
**11.5**.

Two caveats on this. The claim above is an exact transaction count on the real
index arrays, **not** a GPU timing — there is no GPU in this environment. And on
*this CPU* the change is worth only 2–23%, because (see §1) the fetch is under a
fifth of `dem.elevation`; the traffic saving is real but numpy overhead hides it.
`render_skyline` currently builds its grid with `meshgrid(indexing='ij')`, i.e.
ray-major, the worse order — harmless here, wrong on a GPU.

The second correction is smaller and geometric: at 46.5 °N an SRTM 1-arcsec post
is 30.9 m N–S but only **21.3 m E–W**, because the posting is 1 arcsec in
longitude too. Anything reasoning about "how many posts apart" has to carry the
`cos(latitude)`.

## 7. Latency budget, and the honest part

Per photo, on this CPU: extraction 118 ms, render ~10 s, scoring ~4 ms. Removing
the render is the whole problem, and §4 alone takes 68–99 % of it before any
hardware is involved.

### Sizing a hardware-accelerated target (worked example: Jetson Orin NX)

Combining §1, §4 and §6, one photo with a 25 km² box costs **10.7 M bilinear
gathers** (after the early-out) at **~29 B/sample** if the march is ordered
range-major — about **310 MB** of DRAM traffic. Against Orin NX's 102.4 GB/s
LPDDR5 that is **~3 ms at peak, ~5 ms at 60 % achieved**. The other two ceilings
are far away: ~15 FLOP/sample × 10.7 M ÷ ~1.9 TFLOPS fp32 ≈ 0.08 ms of
arithmetic, and the 12.2 M pyramid probes touch 49 KB so they cost issue slots,
not bandwidth.

Order the march ray-major instead and the traffic is ~70 B/sample — 750 MB,
7–12 ms. That is the whole practical consequence of §6.

**The conclusion is an inversion, and it is the point of this section.** On this
x86 box the render is 10 s and extraction is 118 ms — 99 % render. Put the render
on the GPU and it becomes single-digit milliseconds while extraction, still on
the CPU, stays in the 150–400 ms range. **Extraction becomes ~95 % of the
latency.** Every profile taken on a workstation says "optimise the search"; the
device says "optimise the extractor". Whoever ports this should port the DP
extractor to C/NEON first and leave the renderer alone.

Two things could still move this and are NOT established here: whether CUDA's
hardware bilinear filtering (fixed-point filter weights) is precise enough given
that fp16 *arithmetic* already costs 29′ (§2) — if not, the blend must be done
manually in fp32, which raises the ALU term but not the bandwidth term, so the
millisecond conclusion likely survives; and the power mode, since clocks and
memory bandwidth both drop at 10 W. The bandwidth figures above are peak-derived
with a single efficiency assumption, not measured on the part.

What an accelerator cannot fix is that **the fix rate is not limited by compute**
(see `RESULTS.md`): a 65 m grid produced no better answer than a 750 m grid, and
the sweep's separation values say most scenes never clear the fix gate at all. A
faster renderer buys more candidates per second; it does not buy discrimination
between them. The device should spend its budget on the §4 early-out to make
each fix cheap, and on rejecting scenes early — the `extraction_quality` gate is
O(W) and costs nothing — rather than on searching harder.
