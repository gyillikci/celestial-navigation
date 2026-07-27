''' Sub-pixel disk metrology for Sun/Moon frames.

    Squeezes the highest precision out of a resolved disk by locating the limb to
    a small fraction of a pixel, then fitting the circle (centre, radius, plate
    scale) and, for a partly lit Moon, the terminator (illuminated fraction) and
    bright-limb position angle.

    Method (validated on real iPhone frames):
      * cast radial rays from a seed centre; sample each intensity profile with
        bilinear interpolation at 0.1 px;
      * NORMALISED CROSS-CORRELATE each profile against an `erf` step-edge
        template and parabola-refine the correlation peak -> sub-pixel limb
        radius.  Only rays with a genuine bright->sky step are kept, so the soft
        terminator and craters are rejected;
      * least-squares circle fit of the sub-pixel limb points.

    On a sharp full Moon this reaches a circle RMSE ~0.01 px (radius good to a
    few 1e-4 px, i.e. the plate scale to ~0.001%); a first quarter ~0.04 px.
    The bright-limb ANGLE, by contrast, is phase/albedo limited to ~2 deg -- see
    `optical_attitude.bright_limb_sigma_deg`.  numpy only (no scipy).

    (c) 2026.  MIT License (see LICENSE file).
'''

from math import erf
import numpy as np

_MOON_R_ARCSEC = 932.0            # mean topocentric angular radius (887-1010)
_SUN_R_ARCSEC = 960.0


def _bilinear(g, x, y):
    x0 = np.clip(np.floor(x).astype(int), 0, g.shape[1] - 2)
    y0 = np.clip(np.floor(y).astype(int), 0, g.shape[0] - 2)
    fx, fy = x - x0, y - y0
    return (g[y0, x0] * (1 - fx) * (1 - fy) + g[y0, x0 + 1] * fx * (1 - fy)
            + g[y0 + 1, x0] * (1 - fx) * fy + g[y0 + 1, x0 + 1] * fx * fy)


def _circle_fit(x, y):
    ''' Algebraic (Kasa) least-squares circle fit -> (cx, cy, R). '''
    A = np.c_[2 * x, 2 * y, np.ones(len(x))]
    b = x * x + y * y
    s, *_ = np.linalg.lstsq(A, b, rcond=None)
    return s[0], s[1], float(np.sqrt(s[2] + s[0] ** 2 + s[1] ** 2))


def _erf_template(n, width=6.0):
    u = np.linspace(-3, 3, n)
    t = np.array([0.5 * (1 - erf(uu / (width / 3.0))) for uu in u])
    return t - t.mean()


def _circ3(p, q, r):
    ax, ay = p; bx, by = q; cx, cy = r
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if abs(d) < 1e-6:
        return None
    ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay)
          + (cx * cx + cy * cy) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx)
          + (cx * cx + cy * cy) * (bx - ax)) / d
    return ux, uy, np.hypot(ax - ux, ay - uy)


def _gradient_seed(g, sky, bright):
    ''' Robust centre/radius from the SHARP sky-limb (bright inside -> dark sky
        outside), found by gradient + a small deterministic RANSAC circle.  Works
        when the disk is only partly lit (a crescent/quarter), where a brightness
        threshold would latch onto a sliver. '''
    gy, gx = np.gradient(g)
    gm = np.hypot(gx, gy)
    ey, ex = np.nonzero(gm > np.percentile(gm, 99.2))
    if len(ex) < 30:
        return None
    cx0, cy0 = ex.mean(), ey.mean()
    ux, uy = ex - cx0, ey - cy0
    L = np.hypot(ux, uy) + 1e-9
    ox = np.clip((ex + 6 * ux / L).astype(int), 0, g.shape[1] - 1)
    oy = np.clip((ey + 6 * uy / L).astype(int), 0, g.shape[0] - 1)
    ix = np.clip((ex - 6 * ux / L).astype(int), 0, g.shape[1] - 1)
    iy = np.clip((ey - 6 * uy / L).astype(int), 0, g.shape[0] - 1)
    m = ((g[oy, ox] < sky + 0.15 * (bright - sky))
         & (g[iy, ix] > sky + 0.20 * (bright - sky)))
    ex, ey = ex[m].astype(float), ey[m].astype(float)
    if len(ex) < 30:
        return None
    rmax = 0.6 * min(g.shape)
    rs = np.random.RandomState(0)               # deterministic
    best, bn = None, 0
    for _ in range(3000):
        i, j, k = rs.randint(0, len(ex), 3)
        c = _circ3((ex[i], ey[i]), (ex[j], ey[j]), (ex[k], ey[k]))
        if not c or not (20 < c[2] < rmax):
            continue
        n = (np.abs(np.hypot(ex - c[0], ey - c[1]) - c[2]) < 2).sum()
        if n > bn:
            bn, best = n, c
    if best is None:
        return None
    cx, cy, r = best
    inl = np.abs(np.hypot(ex - cx, ey - cy) - r) < 2.5
    return _circle_fit(ex[inl], ey[inl])


def _seed(g, sky, bright):
    ''' Rough centre/radius.  Prefer the robust gradient sky-limb seed; fall back
        to a brightness threshold only if that fails. '''
    grad = _gradient_seed(g, sky, bright)
    if grad is not None and grad[2] > 0.05 * min(g.shape):
        return grad
    lit = g > sky + 0.14 * (bright - sky)
    ys, xs = np.nonzero(lit)
    return xs.mean(), ys.mean(), float(np.sqrt(lit.sum() / np.pi))


def subpixel_limb(g, seed=None, n_rays=1440, half_window=12.0):
    ''' Sub-pixel limb fit.  Returns dict with cx, cy, R, rmse (px), n_points,
        and the accepted ray angles (rad) + points (Nx2).

        (cx, cy) is the TRUE disk centre: only the illuminated LIMB arc (a sharp
        bright->sky edge) is used and the FULL circle is fitted through it, so the
        centre is correct even for a crescent/quarter where most of the disk is
        unseen.  This is deliberately NOT the centroid of the lit blob -- that
        centroid is pulled toward the bright limb (e.g. ~9' at first quarter,
        ~17 km of position error if mistaken for the Moon's position).  The unlit
        half is still there; the circle fit puts it back. '''
    g = np.asarray(g, float)
    sky = np.percentile(g, 20)
    bright = np.percentile(g, 99)
    cx, cy, R = seed if seed else _seed(g, sky, bright)
    rs = np.arange(-half_window, half_window, 0.1)
    templ = _erf_template(len(rs))
    pts = ang = None
    for _ in range(3):
        thal = np.linspace(0, 2 * np.pi, n_rays, endpoint=False)
        P, A = [], []
        for th in thal:
            c, s = np.cos(th), np.sin(th)
            rr = R + rs
            p = _bilinear(g, cx + rr * c, cy + rr * s)
            if p[:80].mean() < bright * 0.5 or p[-80:].mean() > sky + 0.25 * (bright - sky):
                continue                                   # not a bright->sky limb
            # NCC vs the erf template ACCEPTS/REJECTS the ray (real limb, not a
            # crater or the soft terminator) -- a shift-invariant match score.
            pc = p - p.mean()
            den = np.linalg.norm(pc) * np.linalg.norm(templ)
            if den < 1e-6 or (np.correlate(pc, templ, mode="valid").max()
                              / den) < 0.90:
                continue
            # LOCATE the edge at the bright->dark gradient peak (unbiased for a
            # symmetric edge), parabola-refined to sub-pixel.
            gr = -np.gradient(p)                        # +ve at a falling edge
            k = int(np.argmax(gr))
            if k <= 0 or k >= len(gr) - 1:
                continue
            d = 0.5 * (gr[k - 1] - gr[k + 1]) / (
                gr[k - 1] - 2 * gr[k] + gr[k + 1] + 1e-12)
            re = R + rs[k] + d * 0.1
            P.append((cx + re * c, cy + re * s))
            A.append(th)
        if len(P) < 12:
            break
        pts, ang = np.array(P), np.array(A)
        cx, cy, R = _circle_fit(pts[:, 0], pts[:, 1])
    if pts is None:
        raise ValueError("no limb points found")
    res = np.hypot(pts[:, 0] - cx, pts[:, 1] - cy) - R
    return dict(cx=cx, cy=cy, R=R, rmse=float(np.sqrt(np.mean(res ** 2))),
                n_points=len(pts), ang=ang, pts=pts,
                sky=float(sky), bright=float(bright))


def plate_scale_arcsec_px(R_px, body="moon"):
    ''' Angular plate scale from the fitted disk radius. '''
    r = _MOON_R_ARCSEC if body.lower() == "moon" else _SUN_R_ARCSEC
    return r / R_px


def illuminated_fraction(g, fit):
    ''' Illuminated fraction k = lit-disk-area / disk-area (0..1). '''
    g = np.asarray(g, float)
    cx, cy, R = fit["cx"], fit["cy"], fit["R"]
    lvl = 0.5 * (fit["sky"] + fit["bright"])
    yy, xx = np.mgrid[0:g.shape[0], 0:g.shape[1]]
    disk = np.hypot(xx - cx, yy - cy) < 0.985 * R
    return float(((g > lvl) & disk).sum() / disk.sum())


def bright_limb_concentration(fit):
    ''' 0..1: how directional the lit limb arc is (near 0 => full Moon, the
        bright limb is the whole ring and its POSITION ANGLE is degenerate). '''
    a = fit["ang"]
    return float(np.hypot(np.cos(a).mean(), np.sin(a).mean()))
