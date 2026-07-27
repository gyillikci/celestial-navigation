# Results — iPhone 17 Pro Sun+Moon daytime sighting with IMU + factor graph

*Epoch 2026-03-24 12:00 UTC, observer near Greenwich (51.5°N, 0°). Sun+Moon ~94° apart in azimuth. 8 seeds per cell; values are mean ± std of RMS position error in km.*

## Key findings

1. **Fusing many Sun+Moon shots with IMU cuts error 4–10×** versus a single-epoch two-body fix — from tens/hundreds of km down to 4 km (land), 8 km (air) and 25 km (sea).

2. **A single hand-held phone sight is a weak instrument.** The synthetic horizon's tilt error is ~6′ braced on land but grows to tens of arc-minutes at sea/air (vs ~1–2′ for a marine sextant), so the whole value is in fusing many shots, not in any one shot.

3. **The least-rotation shutter clearly helps on land and in the air** (≈2× lower error) by cutting the per-shot horizon noise 3–6×. **At sea it is a wash** — even the calmest swell instant is too tilted, so a handheld phone needs a gimbal/mount there; the swell floor, not the shutter, dominates.

4. **Geometry matters:** the fix is well-conditioned only when the Sun and Moon are well separated in azimuth (~90° here, a first-quarter Moon); near-parallel lines of position degrade it.

## 1. Factor graph vs. the incumbent single-fix

| Regime | starfix single-fix (per-epoch RMS) | Factor graph, no IMU | **Factor graph + IMU** |
|---|---|---|---|
| Land (stationary) | 16.7 ± 1.5 | 16.7 ± 1.5 | **4.1 ± 2.1** |
| Sea (vessel + swell) | 241.0 ± 35.2 | 157.1 ± 102.4 | **25.4 ± 14.6** |
| Air (aircraft) | 75.5 ± 19.7 | 75.1 ± 19.3 | **7.7 ± 4.4** |

![main](results/fig_main.png)

## 2. Least-rotation capture trigger

Modelled synthetic-horizon altitude noise at the chosen shutter instants (arc-minutes):

| Regime | gated (least-rotation) | ungated (periodic) |
|---|---|---|
| Land (stationary) | 6.4′ | 18.1′ |
| Sea (vessel + swell) | 95.9′ | 402.6′ |
| Air (aircraft) | 22.2′ | 138.3′ |

Resulting position error (factor graph + IMU):

| Regime | gated | ungated |
|---|---|---|
| Land (stationary) | 5.7 ± 2.4 | 7.7 ± 4.8 |
| Sea (vessel + swell) | 26.2 ± 14.5 | 23.1 ± 12.5 |
| Air (aircraft) | 8.6 ± 3.1 | 16.4 ± 9.6 |

![gating](results/fig_gating.png)

## 3. Error vs. number of fused shots

![convergence](results/fig_convergence.png)

## 4. The trigger in action

![trigger](results/fig_trigger.png)

## 5. Fix with covariance (error ellipse)

![ellipse](results/fig_ellipse.png)

## Model assumptions

- Representative phone-class MEMS IMU and periscope tele camera (see `iphone_model.py`); numbers are order-of-magnitude, not datasheet values.
- Spherical Earth; observer positions geocentric; Sun/Moon ephemeris reused verbatim from `starfix` (real nautical almanac).
- Altitudes are geometric (refraction/parallax/semidiameter treated as pre-corrected). A ~30 km dead-reckoning offset seeds the coarse prior, so reported accuracy comes from the sky, not the prior.
- The synthetic horizon comes from the phone's gravity vector (no sea-horizon dip); its dominant error is IMU tilt, which the least-rotation shutter minimises.
