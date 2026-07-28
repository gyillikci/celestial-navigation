''' Fit the Moon's disk geometry to a photograph, and read the libration off it.

    The chain in `lunar_geometry` + `lunar_texture` PREDICTS what a photograph
    should look like.  This inverts it: given the photograph, solve for where the
    Moon was actually pointing.  That turns a rendering pipeline into a
    measurement, and a measurement can be checked against the ephemeris.

    HOW

      1. `subpixel_limb` fits the limb -> a first guess at centre and radius.
      2. `measure_rotation` correlates the whole disk against a textured render
         -> the in-image rotation to about a degree.
      3. `tie_points` locks small patches around named craters onto the render by
         normalised cross-correlation, refined to sub-pixel by a 2-D parabola.
      4. `fit_geometry` least-squares six parameters -- centre x, centre y,
         radius, in-image rotation, and the SUB-OBSERVER LONGITUDE AND LATITUDE
         -- to those tie-points.

    Steps 1-3 use the ephemeris only as a starting guess; step 4's libration
    comes from the pixels.  On IMG_7790 (322 px radius, 20 tie-points) the fit
    closes at 0.23 px rms and the libration agrees with the ephemeris to under
    0.3 deg, which is the level at which the hand-typed feature catalogue starts
    to be the limiting error rather than the image.

    WHAT THIS IS GOOD FOR

      * A check on the whole ephemeris chain that a real photograph can settle.
      * Dating an undated photograph: libration moves at roughly 0.1 deg/hour, so
        a 0.1 deg measurement is worth about an hour.  That is much weaker than
        it sounds -- see `time_from_libration`, which reports the ambiguity
        honestly instead of quoting a single minute.
      * NOT for camera roll unless the camera was known to be levelled.  A long
        lens on a tripod can be rotated to any angle, and then the recovered
        rotation measures the tripod, not the sky.

    (c) 2026.  MIT License (see LICENSE file).
'''

import numpy as np

from .disk_metrology import subpixel_limb
from .lunar_features import FEATURES
from .lunar_geometry import topocentric_libration
from .lunar_texture import (render, render_on_grid, subsolar_point,
                            _libration_R, _seleno_vec)
from .lunar_orientation import standardize_disk, _rotate, _masked_ncc

# Features used as tie-points: well-defined, spread over the disk, and each one
# distinctive enough that a 50 px patch correlates unambiguously.  Big smooth
# maria are deliberately absent -- a mare centre is not a landmark.
TIE_FEATURES = (
    "Copernicus", "Tycho", "Plato", "Aristarchus", "Grimaldi", "Mare Crisium",
    "Ptolemaeus", "Triesnecker", "Theophilus", "Sinus Iridum", "Petavius",
    "Clavius", "Langrenus", "Gassendi", "Posidonius", "Archimedes",
    "Bullialdus", "Aristoteles", "Kepler", "Schickard",
)


def project(lon, lat, libration, rotation_deg, cx, cy, r_px):
    ''' Selenographic (lon, lat) -> (x_px, y_px, z), inverse of `render`.
        z > 0 means the point is on the visible hemisphere; z is also the cosine
        of the emission angle, so it says how badly foreshortened it is. '''
    p = _libration_R(*libration) @ _seleno_vec(lon, lat)
    th = np.radians(rotation_deg)
    ct, st = np.cos(th), np.sin(th)
    return (cx + (p[0] * ct - p[1] * st) * r_px,
            cy - (p[0] * st + p[1] * ct) * r_px,
            p[2])


def measure_rotation(gray, cx, cy, r_px, libration, subsolar, n=260, coarse=1.0):
    ''' In-image rotation of the Moon's north pole (deg, counter-clockwise), by
        rotational NCC of the whole disk against a textured render.

        Returns (rotation_deg, ncc, match_width_deg).  `match_width_deg` is the
        NCC peak half-width -- how sharply rotation is constrained -- not the
        uncertainty, which the parabola centroid pins far more finely.
    '''
    tgt, ctr, r_out = standardize_disk(gray, (cx, cy), r_px, out=n)
    px, py = ctr
    yy, xx = np.mgrid[0:n, 0:n]
    mask = ((xx - px) ** 2 + (yy - py) ** 2) < (0.92 * r_out) ** 2
    ref, _ = render(size=n, radius_frac=r_out / (n / 2.0), libration=libration,
                    subsolar=subsolar, rotation_deg=0.0)
    angles = np.arange(-180.0, 180.0, coarse)
    ncc = np.array([_masked_ncc(tgt, _rotate(ref, a, px, py), mask) for a in angles])
    k = int(np.argmax(ncc))
    win = np.abs(((angles - angles[k] + 180) % 360) - 180) <= 6.0
    a2, a1, _ = np.polyfit(angles[win] - angles[k], ncc[win], 2)
    a_best = angles[k] - a1 / (2 * a2) if a2 < 0 else float(angles[k])
    peak = float(ncc[k])
    width = float(np.sqrt(max(1.0 - peak, 1e-4) / abs(a2))) if a2 < 0 else float("inf")
    # render(rotation_deg=A) == _rotate(render(0), -A); hence the sign flip.
    return -float(a_best), peak, width


def _full_render(gray_shape, cx, cy, r_px, libration, subsolar, rotation_deg):
    ''' A render on the PHOTOGRAPH'S OWN pixel grid, so patches at the same
        coordinates are directly comparable.

        Rendering into a separate square array and pasting it at an integer
        offset costs up to half a pixel of registration error, which is twice the
        tie-point residual measured here -- and it biases every tie-point the
        same way, so it does not average out.  `render_on_grid` honours a
        fractional centre exactly and removes the problem at the source.
    '''
    yy, xx = np.mgrid[0:gray_shape[0], 0:gray_shape[1]].astype(float)
    return render_on_grid(xx, yy, cx, cy, r_px, libration=libration,
                          subsolar=subsolar, rotation_deg=rotation_deg)


def tie_points(gray, cx, cy, r_px, libration, subsolar, rotation_deg,
               half=26, search=10, min_ncc=0.55, min_z=0.20, names=None):
    ''' Locate named craters in the photograph by matching against the render.

        `min_z` drops features too close to the limb to correlate (z is the
        cosine of the emission angle, so 0.20 is 78 deg of foreshortening).
        Returns a list of dicts with the selenographic coordinates, the measured
        pixel position, and the NCC that placed it.
    '''
    canvas = _full_render(gray.shape, cx, cy, r_px, libration, subsolar,
                          rotation_deg)
    wanted = set(names or TIE_FEATURES)
    out = []
    for name, lon, lat, _, _ in FEATURES:
        if name not in wanted:
            continue
        x, y, z = project(lon, lat, libration, rotation_deg, cx, cy, r_px)
        if z <= min_z:
            continue
        xi, yi = int(round(x)), int(round(y))
        if not (half + search < xi < gray.shape[1] - half - search
                and half + search < yi < gray.shape[0] - half - search):
            continue
        dx, dy, c = _subpixel_offset(gray, canvas, xi, yi, half, search)
        if c < min_ncc:
            continue
        # SIGN.  `_subpixel_offset` slides the RENDER window by (dx, dy) to line
        # it up with a fixed photograph window, so a render feature sitting dx to
        # the LEFT of its photographic counterpart wins at dx > 0.  The feature's
        # position in the photograph is therefore x - dx, not x + dx.  Getting
        # this backwards does not merely bias the answer: it makes the
        # solve/re-render loop divergent, doubling the libration error every
        # round while the residuals stay small enough to look healthy.
        out.append(dict(name=name, lon=lon, lat=lat, z=z, ncc=c,
                        x=x - dx, y=y - dy, x0=x, y0=y))
    return out


def _subpixel_offset(gray, canvas, xi, yi, half, search):
    A = gray[yi - half:yi + half, xi - half:xi + half]
    A = A - A.mean()
    na = np.sqrt((A * A).sum())
    s = search
    C = np.full((2 * s + 1, 2 * s + 1), -9.0)
    for j, dy in enumerate(range(-s, s + 1)):
        for i, dx in enumerate(range(-s, s + 1)):
            B = canvas[yi + dy - half:yi + dy + half, xi + dx - half:xi + dx + half]
            if B.shape != A.shape:
                continue
            B = B - B.mean()
            nb = np.sqrt((B * B).sum())
            if nb > 1e-9 and na > 1e-9:
                C[j, i] = float((A * B).sum() / (na * nb))
    j, i = np.unravel_index(np.argmax(C), C.shape)
    dx = dy = 0.0
    if 0 < i < 2 * s:
        den = C[j, i - 1] - 2 * C[j, i] + C[j, i + 1]
        if abs(den) > 1e-9:
            dx = 0.5 * (C[j, i - 1] - C[j, i + 1]) / den
    if 0 < j < 2 * s:
        den = C[j - 1, i] - 2 * C[j, i] + C[j + 1, i]
        if abs(den) > 1e-9:
            dy = 0.5 * (C[j - 1, i] - C[j + 1, i]) / den
    return (i - s + dx), (j - s + dy), float(C[j, i])


def fit_geometry(ties, guess, iters=40):
    ''' Least-squares (cx, cy, r_px, rotation_deg, lib_lon, lib_lat) from
        tie-points.  `guess` is that tuple.  Gauss-Newton with a numerical
        Jacobian: six parameters against 2N residuals, and N is 20, so this is a
        two-line solver, not a job for an optimiser.

        Returns dict(params, rms_px, sigma, residuals).  `sigma` is the FORMAL
        1-sigma from the covariance -- it believes the feature catalogue, which
        is optimistic; `jackknife_sigma` gives the empirical alternative.
    '''
    sel = np.array([[t["lon"], t["lat"]] for t in ties], float)
    obs = np.array([[t["x"], t["y"]] for t in ties], float)

    def resid(p):
        cx, cy, r, rot, ll, lb = p
        out = np.empty(2 * len(sel))
        for k, (lo, la) in enumerate(sel):
            x, y, _ = project(lo, la, (ll, lb), rot, cx, cy, r)
            out[2 * k] = x - obs[k, 0]
            out[2 * k + 1] = y - obs[k, 1]
        return out

    p = np.array(guess, float)
    J = None
    for _ in range(iters):
        r = resid(p)
        J = np.empty((len(r), 6))
        for k in range(6):
            q = p.copy()
            q[k] += 1e-3
            J[:, k] = (resid(q) - r) / 1e-3
        dp = np.linalg.lstsq(J, -r, rcond=None)[0]
        p = p + dp
        if np.abs(dp).max() < 1e-9:
            break
    r = resid(p)
    rms = float(np.sqrt((r ** 2).mean()))
    try:
        sigma = np.sqrt(np.diag(np.linalg.inv(J.T @ J))) * rms
    except np.linalg.LinAlgError:
        sigma = np.full(6, np.nan)
    return dict(params=p, rms_px=rms, sigma=sigma, residuals=r,
                names=[t["name"] for t in ties])


def jackknife_sigma(ties, guess):
    ''' Leave-one-out spread of the fitted parameters.

        The formal covariance assumes the catalogue coordinates are exact.  They
        are not -- they are rounded published centres, and a mare "centre" is a
        judgement call.  Dropping each tie-point in turn shows how much any one
        of them is carrying the answer.
    '''
    n = len(ties)
    ps = []
    for k in range(n):
        sub = [t for i, t in enumerate(ties) if i != k]
        ps.append(fit_geometry(sub, guess)["params"])
    ps = np.array(ps)
    return np.sqrt((n - 1) / n * ((ps - ps.mean(0)) ** 2).sum(0))


def solve_photo(gray, time_iso, obs_lat, obs_lon, names=None, rounds=4,
                tol_deg=0.005):
    ''' The whole chain on one photograph.  The ephemeris at `time_iso` supplies
        the starting guess and the sub-solar point (needed for the terminator,
        and slow enough that the starting time is not critical); the geometry
        itself is measured.

        ITERATED ON PURPOSE.  The tie-points are found by correlating against a
        render, and that render is built at the CURRENT libration estimate -- so
        a wrong starting libration foreshortens the template wrongly and drags
        the measured positions with it.  One pass is therefore not a measurement,
        it is a measurement plus a fraction of the guess.  Re-rendering at the
        fitted geometry and re-measuring removes that: the loop below converges
        in two or three rounds and then the answer no longer depends on where it
        started, which is the whole point of calling it a measurement.
    '''
    fit = subpixel_limb(gray)
    lib = topocentric_libration(time_iso, obs_lat, obs_lon)
    sun = subsolar_point(time_iso)
    sun0 = (sun["lon"], sun["lat"])
    rot, ncc, width = measure_rotation(gray, fit["cx"], fit["cy"], fit["R"],
                                       (lib["lon"], lib["lat"]), sun0)

    p = np.array([fit["cx"], fit["cy"], fit["R"], rot, lib["lon"], lib["lat"]])
    ties, sol, history = [], None, []
    for _ in range(rounds):
        ties = tie_points(gray, p[0], p[1], p[2], (p[4], p[5]), sun0, p[3],
                          names=names)
        if len(ties) < 4:
            break
        sol = fit_geometry(ties, p)
        moved = float(np.abs(sol["params"][3:] - p[3:]).max())
        p = sol["params"]
        history.append(dict(params=p.copy(), rms_px=sol["rms_px"],
                            n_ties=len(ties), moved_deg=moved))
        if moved < tol_deg:
            break
    return dict(limb=fit, ephem_libration=lib, subsolar=sun, rotation=p[3],
                disk_ncc=ncc, match_width_deg=width, ties=ties, fit=sol,
                history=history)


def time_from_libration(meas_lon, meas_lat, obs_lat, obs_lon, t_start, t_end,
                        step_minutes=20, sigma_deg=0.15):
    ''' Times whose predicted libration matches a measured one.

        Libration drifts at roughly 0.1 deg/hour, so a measurement good to
        `sigma_deg` = 0.15 deg is worth an hour and a half AT BEST -- and the
        drift stalls and reverses at the turning points, where the same libration
        recurs and the answer is genuinely ambiguous.  This returns every time in
        the window within `sigma_deg`, plus the best, so the ambiguity is visible
        rather than hidden behind a single number.
    '''
    from datetime import datetime, timedelta
    t0 = datetime.fromisoformat(t_start)
    t1 = datetime.fromisoformat(t_end)
    out = []
    t = t0
    while t <= t1:
        iso = t.strftime("%Y-%m-%dT%H:%M:%S")
        lb = topocentric_libration(iso, obs_lat, obs_lon)
        miss = float(np.hypot(lb["lon"] - meas_lon, lb["lat"] - meas_lat))
        out.append(dict(time=iso, lon=lb["lon"], lat=lb["lat"], miss_deg=miss))
        t += timedelta(minutes=step_minutes)
    out.sort(key=lambda r: r["miss_deg"])
    within = [r for r in out if r["miss_deg"] <= sigma_deg]
    return dict(best=out[0], within=within, all=out)
