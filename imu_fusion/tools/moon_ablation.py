''' What each ingredient of the lunar render is actually worth.

    The disk match went from a schematic ten-blob cartoon to a
    Stellarium-textured, librated, terminator-carrying render.  That is four
    changes at once, and "the match improved" is not a finding unless it says
    WHICH change bought what.  This measures them one at a time, on a real
    photograph, two ways:

      * **pattern match** -- whole-disk normalised cross-correlation, with the
        in-plane rotation re-optimised for every condition so orientation is
        never the excuse for a bad score;
      * **geometric match** -- the tie-point residual, with disk centre, radius
        and rotation always free, so the fit gets every chance to absorb a wrong
        libration before it is charged for one.

    The second is the honest test.  Libration is 0.96 correlated with the disk
    centre (see `lunar_match`), so a large part of a libration error simply
    translates the pattern and re-centring hides it; what survives is the
    genuine, irreducible cost.

    Run it:

        python -m imu_fusion.tools.moon_ablation PHOTO.JPG \\
            --time 2026-07-28T20:45:00 --lat 41.0082 --lon 28.9784

    (c) 2026.  MIT License (see LICENSE file).
'''

import argparse

import numpy as np

from ..disk_metrology import subpixel_limb
from ..lunar_geometry import topocentric_libration, geocentric_libration
from ..lunar_texture import render, subsolar_point, parallactic_angle, find_texture
from ..lunar_orientation import (render_moon, standardize_disk, _rotate,
                                 _masked_ncc)
from ..lunar_match import tie_points, fit_geometry, project

N = 260


def pattern_ablation(gray, time_iso, obs_lat, obs_lon):
    ''' Whole-disk NCC for each modelling step, rotation re-optimised each time. '''
    fit = subpixel_limb(gray)
    tgt, ctr, r_out = standardize_disk(gray, (fit["cx"], fit["cy"]), fit["R"], out=N)
    cx, cy = ctr
    yy, xx = np.mgrid[0:N, 0:N]
    mask = ((xx - cx) ** 2 + (yy - cy) ** 2) < (0.92 * r_out) ** 2

    lib = topocentric_libration(time_iso, obs_lat, obs_lon)
    geo = geocentric_libration(time_iso)
    sun = subsolar_point(time_iso)
    LIB, GEO = (lib["lon"], lib["lat"]), (geo["lon"], geo["lat"])
    SUN = (sun["lon"], sun["lat"])
    pq = lib["pole_pa"] - parallactic_angle(time_iso, obs_lat, obs_lon)
    rf = r_out / (N / 2.0)

    def best(ref):
        ang = np.arange(-180.0, 180.0, 1.0)
        n = [_masked_ncc(tgt, _rotate(ref, a, cx, cy), mask) for a in ang]
        k = int(np.argmax(n))
        return float(n[k]), -float(ang[k])

    def texture(L, S):
        return render(size=N, radius_frac=rf, libration=L, subsolar=S,
                      rotation_deg=0.0)[0]

    def schematic(L):
        return render_moon(size=N, radius_frac=rf, libration=L,
                           pole_pa_deg=0.0, roll_deg=0.0)[0]

    conditions = [
        ("schematic maria, no libration (the code as it stood)", schematic((0.0, 0.0))),
        ("schematic maria + libration", schematic(LIB)),
        ("Stellarium texture, no libration, no terminator", texture((0.0, 0.0), None)),
        ("Stellarium texture + geocentric libration", texture(GEO, None)),
        ("Stellarium texture + topocentric libration", texture(LIB, None)),
        ("Stellarium texture + topocentric libration + terminator", texture(LIB, SUN)),
    ]
    out = []
    for name, ref in conditions:
        ncc, rot = best(ref)
        held = float(_masked_ncc(tgt, _rotate(ref, -pq, cx, cy), mask))
        out.append(dict(name=name, ncc=ncc, rot=rot, ncc_at_predicted_rot=held))
    return out


def geometric_ablation(gray, time_iso, obs_lat, obs_lon, params=None):
    ''' Tie-point rms under each libration model, everything else free. '''
    lib = topocentric_libration(time_iso, obs_lat, obs_lon)
    geo = geocentric_libration(time_iso)
    sun = subsolar_point(time_iso)
    SUN = (sun["lon"], sun["lat"])

    if params is None:
        from ..lunar_match import solve_photo
        params = solve_photo(gray, time_iso, obs_lat, obs_lon)["fit"]["params"]
    p = np.array(params, float)

    # One fixed set of observations, measured at the converged geometry, so every
    # condition below is scored against the SAME pixels.
    ties = tie_points(gray, p[0], p[1], p[2], (p[4], p[5]), SUN, p[3])
    sel = np.array([[t["lon"], t["lat"]] for t in ties])
    obs = np.array([[t["x"], t["y"]] for t in ties])

    def resid(q):
        out = np.empty(2 * len(sel))
        for k, (lo, la) in enumerate(sel):
            x, y, _ = project(lo, la, (q[4], q[5]), q[3], q[0], q[1], q[2])
            out[2 * k], out[2 * k + 1] = x - obs[k, 0], y - obs[k, 1]
        return out

    def solve(free_libration, libration):
        q = p.copy()
        if not free_libration:
            q[4], q[5] = libration
        idx = [0, 1, 2, 3] + ([4, 5] if free_libration else [])
        for _ in range(60):
            r0 = resid(q)
            J = np.empty((len(r0), len(idx)))
            for k, i in enumerate(idx):
                z = q.copy(); z[i] += 1e-3
                J[:, k] = (resid(z) - r0) / 1e-3
            d = np.linalg.lstsq(J, -r0, rcond=None)[0]
            for k, i in enumerate(idx):
                q[i] += d[k]
            if np.abs(d).max() < 1e-9:
                break
        return q, float(np.sqrt((resid(q) ** 2).mean()))

    cases = [("libration forced to (0,0) — the old default", False, (0.0, 0.0)),
             ("libration fixed at the geocentric value", False,
              (geo["lon"], geo["lat"])),
             ("libration fixed at the topocentric value", False,
              (lib["lon"], lib["lat"])),
             ("libration solved from the pixels", True, None)]
    out = []
    for name, free, val in cases:
        q, rms = solve(free, val)
        out.append(dict(name=name, rms_px=rms, lon=q[4], lat=q[5]))

    # how far features move if libration is ignored, AFTER re-centring
    q0, _ = solve(False, (0.0, 0.0))
    q1, _ = solve(True, None)
    d = [np.hypot(*(np.subtract(project(lo, la, (q0[4], q0[5]), q0[3], q0[0], q0[1], q0[2])[:2],
                                project(lo, la, (q1[4], q1[5]), q1[3], q1[0], q1[1], q1[2])[:2])))
         for lo, la in sel]
    return out, dict(median_px=float(np.median(d)), max_px=float(max(d)),
                     n_ties=len(ties))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("photo")
    p.add_argument("--time", required=True)
    p.add_argument("--lat", type=float, required=True)
    p.add_argument("--lon", type=float, required=True)
    a = p.parse_args(argv)
    if find_texture() is None:
        raise SystemExit("Stellarium's lunar albedo map is not installed")
    from PIL import Image
    gray = np.asarray(Image.open(a.photo).convert("L"), float)

    print(f"{'pattern match (whole-disk NCC, rotation re-optimised)':<58}"
          f"{'NCC':>8}{'rot':>8}{'at P-q':>9}")
    for r in pattern_ablation(gray, a.time, a.lat, a.lon):
        print(f"  {r['name']:<56}{r['ncc']:>8.4f}{r['rot']:>+8.2f}"
              f"{r['ncc_at_predicted_rot']:>9.4f}")

    rows, spread = geometric_ablation(gray, a.time, a.lat, a.lon)
    print(f"\n{'geometric match (centre, radius, rotation always free)':<58}"
          f"{'rms px':>8}")
    for r in rows:
        print(f"  {r['name']:<56}{r['rms_px']:>8.3f}"
              f"   lib ({r['lon']:+.2f}, {r['lat']:+.2f})")
    print(f"\nignoring libration displaces features by median "
          f"{spread['median_px']:.1f} px, worst {spread['max_px']:.1f} px "
          f"({spread['n_ties']} tie-points; 1 deg = 5.6 px at this scale)")


if __name__ == "__main__":
    main()
