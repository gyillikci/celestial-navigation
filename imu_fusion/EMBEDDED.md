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
| alpine high 46.55 N | 0.422 s | 0.118 s | **3.6×** | 0.0000′ |
| alpine 46.02 N | 0.404 s | 0.199 s | **2.0×** | 0.0000′ |
| valley 47.05 N | 0.432 s | 0.286 s | **1.5×** | 0.0000′ |
| low camera 46.8 N | 0.707 s | 0.120 s | **5.9×** | 0.0000′ |
| lake, water clamped 47.0 N | 0.319 s | 0.092 s | **3.5×** | 0.0000′ |

Bit-identical output, 1.5–5.9× less work. The pyramid bound is looser than the
exact ray suffix — a square block is not a ray — which is why the realised
speedup is below the 68–99 % ceiling above.

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

Occupancy is therefore not the problem; **memory divergence is**. Neighbouring
rays at long range land 63 m apart and walk in different directions, so the
gather scatters across the tile. Ordering the march so that a warp/wavefront
holds *adjacent ranges on one ray* rather than *one range across many rays*
keeps each fetch group inside a cache line. This is a layout decision, and on
the measured 97 %/3 % split it is the decision that sets the achieved rate.

## 7. Latency budget, and the honest part

Per photo, on this CPU: extraction 118 ms, render ~10 s, scoring ~4 ms. Removing
the render is the whole problem, and §4 alone takes 68–99 % of it before any
hardware is involved.

What an accelerator cannot fix is that **the fix rate is not limited by compute**
(see `RESULTS.md`): a 65 m grid produced no better answer than a 750 m grid, and
the sweep's separation values say most scenes never clear the fix gate at all. A
faster renderer buys more candidates per second; it does not buy discrimination
between them. The device should spend its budget on the §4 early-out to make
each fix cheap, and on rejecting scenes early — the `extraction_quality` gate is
O(W) and costs nothing — rather than on searching harder.
