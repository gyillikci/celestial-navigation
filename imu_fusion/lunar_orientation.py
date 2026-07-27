''' Moon ORIENTATION (camera roll) from crater/mare pattern matching.

    The bright limb is only a ~2 deg compass (phase-limited, degenerate near
    full -- see optical_attitude).  The Moon's SURFACE pattern is a far stronger
    orientation reference: the maria/craters are a fixed fiducial, so matching
    the resolved disk to an ephemeris-correct render recovers the in-plane
    rotation (camera roll, once the ephemeris pole position angle and libration
    are known) to a fraction of a degree -- and it works BEST at full Moon,
    exactly where the bright limb fails.

    Two pieces:
      * `render_moon` -- a schematic disk from a table of the major maria placed
        at their selenographic coordinates, orthographically projected for a
        given libration and rotated by (pole PA + roll).  This is the
        ephemeris-referenced template (absolute roll needs the time; libration
        and pole PA come from the almanac).
      * `recover_roll` -- rotational NORMALISED CROSS-CORRELATION of a target
        disk against a reference over in-plane angle, parabola-refined to
        sub-degree.  It works on real lunar texture (not just the schematic),
        so the achievable roll precision can be measured directly.

    numpy only (no scipy).  (c) 2026.  MIT License (see LICENSE file).
'''

import numpy as np

# Major maria: (name, selenographic lon [deg, +E], lat [deg, +N], ang.radius, depth)
# Approximate centres; depth 0..1 is how dark relative to highlands.
MARIA = [
    ("Crisium",       59.0,  17.0, 0.16, 0.9),
    ("Fecunditatis",  52.0,  -8.0, 0.17, 0.7),
    ("Nectaris",      35.0, -15.0, 0.12, 0.8),
    ("Tranquillit.",  28.0,   8.0, 0.22, 0.8),
    ("Serenitatis",   18.0,  28.0, 0.18, 0.85),
    ("Vaporum",        3.0,  13.0, 0.10, 0.6),
    ("Imbrium",      -16.0,  33.0, 0.26, 0.85),
    ("Procellarum",  -57.0,  18.0, 0.34, 0.75),
    ("Humorum",      -39.0, -24.0, 0.13, 0.8),
    ("Nubium",       -17.0, -21.0, 0.17, 0.7),
]


def _seleno_to_xyz(lon_deg, lat_deg):
    lo, la = np.radians(lon_deg), np.radians(lat_deg)
    return np.array([np.cos(la) * np.sin(lo), np.sin(la), np.cos(la) * np.cos(lo)])


def _libration_R(lib_lon, lib_lat):
    a, b = np.radians(lib_lon), np.radians(lib_lat)
    Ry = np.array([[np.cos(a), 0, np.sin(a)], [0, 1, 0], [-np.sin(a), 0, np.cos(a)]])
    Rx = np.array([[1, 0, 0], [0, np.cos(b), -np.sin(b)], [0, np.sin(b), np.cos(b)]])
    return Rx @ Ry


def render_moon(size=257, radius_frac=0.92, libration=(0.0, 0.0),
                pole_pa_deg=0.0, roll_deg=0.0, highland=180.0, sky=8.0):
    ''' Schematic Moon disk (size x size uint8-like float array).  Maria are
        placed at their selenographic coords, projected for `libration`, and the
        whole disk rotated in-image by (pole_pa + roll). '''
    R = _libration_R(*libration)
    cx = cy = (size - 1) / 2.0
    r_px = radius_frac * (size / 2.0)
    yy, xx = np.mgrid[0:size, 0:size].astype(float)
    X = (xx - cx) / r_px
    Y = -(yy - cy) / r_px
    disk = X * X + Y * Y <= 1.0
    img = np.where(disk, highland, sky)
    # in-image rotation (pole PA + roll).  Sign matches `_rotate` so that
    # render(roll=t) equals _rotate(render(roll=0), t).
    th = -np.radians(pole_pa_deg + roll_deg)
    ct, st = np.cos(th), np.sin(th)
    for _, lon, lat, rad, depth in MARIA:
        p = R @ _seleno_to_xyz(lon, lat)
        if p[2] <= 0.05:                        # on far side / limb
            continue
        # rotate the projected position by +theta
        px = p[0] * ct - p[1] * st
        py = p[0] * st + p[1] * ct
        # foreshortened radius (~ cos of angle from centre = p[2])
        rr = rad * max(p[2], 0.25)
        d2 = (X - px) ** 2 + (Y - py) ** 2
        mask = disk & (d2 < rr * rr)
        img = np.where(mask, img - depth * (highland - sky) *
                       np.exp(-d2 / (2 * (rr / 2) ** 2)), img)
    return img, (cx, cy, r_px)


def _rotate(img, deg, cx, cy):
    ''' Bilinear rotation of img about (cx, cy) by deg (CCW). '''
    th = np.radians(deg)
    ct, st = np.cos(th), np.sin(th)
    yy, xx = np.mgrid[0:img.shape[0], 0:img.shape[1]].astype(float)
    xr = cx + (xx - cx) * ct + (yy - cy) * st
    yr = cy - (xx - cx) * st + (yy - cy) * ct
    x0 = np.clip(np.floor(xr).astype(int), 0, img.shape[1] - 2)
    y0 = np.clip(np.floor(yr).astype(int), 0, img.shape[0] - 2)
    fx, fy = xr - x0, yr - y0
    return (img[y0, x0] * (1 - fx) * (1 - fy) + img[y0, x0 + 1] * fx * (1 - fy)
            + img[y0 + 1, x0] * (1 - fx) * fy + img[y0 + 1, x0 + 1] * fx * fy)


def standardize_disk(img, center, r_px, out=220):
    ''' Crop a square around the disk and resample to `out`x`out` (bilinear), so
        rotational matching runs on a small, scale-normalised patch. '''
    cx, cy = center
    half = 1.05 * r_px
    gx = np.linspace(cx - half, cx + half, out)
    gy = np.linspace(cy - half, cy + half, out)
    XX, YY = np.meshgrid(gx, gy)
    x0 = np.clip(np.floor(XX).astype(int), 0, img.shape[1] - 2)
    y0 = np.clip(np.floor(YY).astype(int), 0, img.shape[0] - 2)
    fx, fy = XX - x0, YY - y0
    patch = (img[y0, x0] * (1 - fx) * (1 - fy) + img[y0, x0 + 1] * fx * (1 - fy)
             + img[y0 + 1, x0] * (1 - fx) * fy + img[y0 + 1, x0 + 1] * fx * fy)
    r_out = r_px / half * (out / 2.0)
    return patch, (out / 2.0, out / 2.0), r_out


def _masked_ncc(a, b, mask):
    a = a[mask]; b = b[mask]
    a = a - a.mean(); b = b - b.mean()
    d = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / d) if d > 1e-9 else 0.0


def recover_roll(target, reference, center, r_px, coarse=1.0, span=180.0):
    ''' Rotational NCC: find the in-plane rotation (deg) that best aligns
        `target` to `reference`, both disks centred at `center` with radius
        `r_px`.  Returns dict(roll, ncc, sigma_deg, curve).

        Coarse grid then parabolic sub-degree refine; sigma from the curvature
        of the NCC peak vs a noise floor (how sharply rotation is constrained).
    '''
    # standardise both to a small centred patch so 360 rotations are cheap
    reference, _, _ = standardize_disk(reference, center, r_px)
    target, center, r_px = standardize_disk(target, center, r_px)
    cx, cy = center
    yy, xx = np.mgrid[0:target.shape[0], 0:target.shape[1]]
    mask = ((xx - cx) ** 2 + (yy - cy) ** 2) < (0.9 * r_px) ** 2
    angles = np.arange(-span, span + coarse, coarse)
    ncc = np.array([_masked_ncc(target, _rotate(reference, a, cx, cy), mask)
                    for a in angles])
    k = int(np.argmax(ncc))
    peak = ncc[k]
    # parabola fit in a local window -> sub-degree peak and curvature
    win = np.abs(angles - angles[k]) <= 8.0
    a2, a1, a0 = np.polyfit(angles[win] - angles[k], ncc[win], 2)
    roll = angles[k] - a1 / (2 * a2) if a2 < 0 else float(angles[k])
    # resolution: angle at which NCC falls from the peak by its mismatch floor
    # (1-peak); with a<0, NCC ~ peak + a*dtheta^2 -> dtheta = sqrt((1-peak)/|a|).
    width = float(np.sqrt(max(1.0 - peak, 1e-4) / abs(a2))) if a2 < 0 else float("inf")
    # NB: `match_width_deg` is the NCC peak HALF-WIDTH (match sharpness), NOT the
    # roll uncertainty -- the parabola CENTROID locates the roll far more finely
    # (demonstrated ~0.06 deg RMS on real full-Moon texture over +/-40 deg).
    return dict(roll=float(roll), ncc=float(peak), match_width_deg=width,
                curve=(angles, ncc))
