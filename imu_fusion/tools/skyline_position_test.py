''' Can a photographed skyline fix your position?  A controlled test.

    Three Theodolite frames from one spot near Istanbul, each carrying its own
    GPS.  The GPS is therefore ground truth, so the question can be asked
    properly: inject a known error into the position prior, search around the
    WRONG prior, and see how far the recovered position lands from the truth.

    THE MODEL.  Position (lat, lon) plus two nuisance parameters shared by all
    frames -- one compass bias and one pitch bias, because both are properties of
    the device rather than of the shot.  Per-frame azimuth, roll and pitch come
    from the app.

    WHY THE NUISANCE PARAMETERS MATTER SO MUCH.  They absorb exactly the two
    first-order signals a position shift produces:

      * moving ACROSS the line of sight swings every feature's bearing by nearly
        the same amount -- which is what a compass bias looks like;
      * moving ALONG it raises or lowers the whole horizon together -- which is
        what a pitch bias looks like.

    What survives is only the DIFFERENCE between near and far features, and for
    this scene (ridges at 49 and 72 km) that is 0.371 deg/km laterally but just
    0.0073 deg/km radially.  So the geometry is a line of position, not a fix:
    the error ellipse is ~50x longer along the sight line than across it.

    THE ANSWER, on this data: no.  See `RESULTS.md`.  The residual floor is
    ~0.14 deg and it is SYSTEMATIC -- smooth in azimuth, correlated across the
    839 samples -- so it does not average down, and a false minimum 15 km away
    scores better than the truth.  Fixing the pitch bias at its calibrated value
    does not rescue it, which is the useful part: the binding constraint is
    skyline-match accuracy, not the free parameters.

        python -m imu_fusion.tools.skyline_position_test --half 6 --step 0.3

    (c) 2026.  MIT License (see LICENSE file).
'''

import argparse

import numpy as np

from ..terrain_resection import DemTiles, render_skyline
from .theodolite_skyline import Frame, extract_skyline

EYE_M = 1.6                     # camera above the ground the DEM reports
KM_PER_DEG_LAT = 111.32


def build_observations(frames):
    ''' Each frame -> (azimuth offsets, elevations) in the app's own frame. '''
    obs = []
    for fr in frames:
        y = extract_skyline(fr.path)
        x = np.arange(fr.w, dtype=float)
        ok = np.isfinite(y)
        x, y = x[ok][::5], y[ok][::5]
        tilt = np.tan(np.radians(fr.roll_deg)) * (x - fr.xc)
        obs.append(dict(name=fr.time or fr.path, az0=fr.az_deg,
                        daz=(x - fr.xc) / fr.f,
                        el=-((y - tilt) - fr.yc) / fr.f + fr.pitch_deg))
    return obs


def residual(obs, dem, lat, lon, az_slack=2.0, az_step=0.05, el_bias=None,
             az_lo=122.0, az_hi=150.0, d_max_km=95.0):
    ''' rms elevation residual (deg) of every frame at one candidate site.

        `el_bias=None` solves the pitch bias; giving it fixes the horizon's
        absolute height and is far more constraining -- in principle.
    '''
    g = float(dem.elevation(np.array([lat]), np.array([lon]))[0])
    azs, els = render_skyline(dem, lat, lon, g + EYE_M, az_start=az_lo,
                              az_end=az_hi, az_step=0.05, d_max_km=d_max_km,
                              d_step_km=0.25)
    best = None
    for dz in np.arange(-az_slack, az_slack + 1e-9, az_step):
        r = np.concatenate([o['el'] - np.interp(o['az0'] + dz + o['daz'], azs, els)
                            for o in obs])
        b = r.mean() if el_bias is None else -el_bias
        s = float(np.sqrt(((r - b) ** 2).mean()))
        if best is None or s < best[0]:
            best = (s, float(dz), float(r.mean()))
    return best


def surface(obs, dem, lat, lon, half_km, step_km, **kw):
    ''' Residual over a square grid centred on (lat, lon). '''
    offs = np.arange(-half_km, half_km + 1e-9, step_km)
    km_lon = KM_PER_DEG_LAT * np.cos(np.radians(lat))
    Z = np.full((len(offs), len(offs)), np.nan)
    for i, dn in enumerate(offs):
        for j, de in enumerate(offs):
            Z[i, j] = residual(obs, dem, lat + dn / KM_PER_DEG_LAT,
                               lon + de / km_lon, **kw)[0]
    return Z, offs


def prior_error_sweep(Z, offs, levels=(0.0, 0.5, 1.0, 2.0, 5.0, 10.0)):
    ''' Offset the prior by each level in eight directions, search a box around
        the WRONG prior, and report how far the winner is from the truth.

        The box is +/-1.5x the injected error, so the truth is always inside it:
        this measures whether the skyline can FIND the truth, not whether the
        search was lucky enough to be pointed at it.
    '''
    dirs = [(1, 0), (.707, .707), (0, 1), (-.707, .707),
            (-1, 0), (-.707, -.707), (0, -1), (.707, -.707)]
    near = lambda v: int(np.argmin(np.abs(offs - v)))
    out = []
    for E in levels:
        half = max(1.5 * E, 1.0)
        errs = []
        for dn, de in dirs:
            cn, ce = E * dn, E * de
            i0, i1 = near(cn - half), near(cn + half)
            j0, j1 = near(ce - half), near(ce + half)
            sub = Z[i0:i1 + 1, j0:j1 + 1]
            if not np.isfinite(sub).any():
                errs.append(np.nan)
                continue
            k = np.unravel_index(np.nanargmin(sub), sub.shape)
            errs.append(float(np.hypot(offs[i0 + k[0]], offs[j0 + k[1]])))
            if E == 0.0:
                break
        if E == 0.0:
            errs = errs * len(dirs)
        out.append(dict(injected_km=E, box_km=half, errors_km=np.array(errs),
                        mean_km=float(np.nanmean(errs)),
                        worst_km=float(np.nanmax(errs))))
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--half", type=float, default=6.0, help="grid half-width, km")
    p.add_argument("--step", type=float, default=0.3, help="grid spacing, km")
    p.add_argument("--dem", default="imu_fusion/dem")
    p.add_argument("--uploads", required=True,
                   help="directory holding the TH0000{20,21,22}.jpg frames")
    a = p.parse_args(argv)
    import os
    u = a.uploads.rstrip("/") + "/"
    frames = [Frame(u + "TH000020.jpg", 40.904778, 29.209650, 315, 2471, -1.4, -0.5, time="13:07:32"),
              Frame(u + "TH000021.jpg", 40.904444, 29.209346, 310, 2400, -0.2, -0.5, time="13:07:53"),
              Frame(u + "TH000022.jpg", 40.904450, 29.209346, 310, 2329, -0.7, -0.2, time="13:07:58")]
    for f in frames:
        if not os.path.exists(f.path):
            raise SystemExit(f"missing {f.path}")
    dem = DemTiles(a.dem)
    obs = build_observations(frames)
    lat = float(np.mean([f.lat for f in frames]))
    lon = float(np.mean([f.lon for f in frames]))
    s, dz, off = residual(obs, dem, lat, lon)
    print(f"{sum(len(o['daz']) for o in obs)} samples; at the true position "
          f"rms {s:.4f} deg, compass bias {dz:+.2f} deg, pitch bias {-off:+.2f} deg")
    Z, offs = surface(obs, dem, lat, lon, a.half, a.step)
    k = np.unravel_index(np.nanargmin(Z), Z.shape)
    print(f"best match in a +/-{a.half:.0f} km box: N{offs[k[0]]:+.1f} E{offs[k[1]]:+.1f} km "
          f"({np.hypot(offs[k[0]], offs[k[1]]):.1f} km from truth), rms {Z[k]:.4f}")
    print(f"\n{'injected':>9}{'box':>7}{'mean':>8}{'worst':>8}")
    for r in prior_error_sweep(Z, offs):
        print(f"{r['injected_km']:8.1f}k{r['box_km']:7.1f}{r['mean_km']:8.1f}{r['worst_km']:8.1f}")
    print("\n(km from the TRUE position to the recovered minimum)")


if __name__ == "__main__":
    main()
