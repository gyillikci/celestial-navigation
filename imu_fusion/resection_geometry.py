''' Will this view fix my position?  Answer it BEFORE running the search.

    Two real datasets in this study were searched exhaustively before the answer
    became obvious, and in both cases the geometry could have said so in
    milliseconds:

      * **Istanbul, 8x zoom over the Sea of Marmara.**  A 20 km search returned a
        false minimum 15 km from the truth WITH A LOWER RESIDUAL than the truth.
        Not a bug -- the information was not there.
      * **Denver, a stock photograph of the Front Range.**  The skyline matched to
        1.11 arcmin, 7.6x better, and the best skyline-only cell was still
        decisively wrong.  What fixed it was one angle to a NEAR object.

    THE MECHANISM, which this module makes computable.  Any real solve carries
    two nuisance parameters -- a compass bias and a pitch bias -- because neither
    a phone magnetometer nor a phone horizon is trustworthy at the arcminute
    level.  Those two absorb exactly the two FIRST-ORDER signals a position shift
    produces:

      * moving ACROSS the line of sight swings every feature's bearing by nearly
        the same amount, which is what a compass bias looks like;
      * moving ALONG it raises the whole horizon together, which is what a pitch
        bias looks like.

    What survives is only the DIFFERENCE between near and far features.  So the
    single number that decides whether a view is usable is the SPREAD of range
    among the things you can identify -- not the sharpness of the picture, not
    the number of pixels, and not the quality of the DEM.

    A worked contrast, both computed by `sensitivity` below:

        Istanbul, ridges at 49 and 72 km   lateral 0.371 deg/km   radial 0.0073
        Denver, tower at 15, peak at 92 km lateral 3.2   deg/km   radial 0.042

    An order of magnitude, from including something close.

    HONEST LIMIT.  These are first-order sensitivities about a nominal position:
    they predict the shape and scale of the error ellipse, not the probability
    that a global search lands in the right basin.  They also assume the residual
    is independent between samples, which on real photographs it is NOT -- see
    `effective_samples`, which exists precisely because that assumption inflated
    a real uncertainty by more than an order of magnitude.

    (c) 2026.  MIT License (see LICENSE file).
'''

import math

import numpy as np

from .landfall import K_REFRACTION, effective_radius_km

EARTH_R_KM = 6371.0


# --------------------------------------------------------------------------- #
# Bearings and ranges
# --------------------------------------------------------------------------- #

def bearing_deg(lat0, lon0, lat1, lon1):
    ''' Initial great-circle bearing from 0 to 1, degrees from north. '''
    p1, p2 = math.radians(lat0), math.radians(lat1)
    dl = math.radians(lon1 - lon0)
    return math.degrees(math.atan2(
        math.sin(dl) * math.cos(p2),
        math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl))) % 360.0


def distance_km(lat0, lon0, lat1, lon1):
    p1, p2 = math.radians(lat0), math.radians(lat1)
    dl = math.radians(lon1 - lon0)
    h = (math.sin((p2 - p1) / 2) ** 2
         + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2)
    return 2 * EARTH_R_KM * math.asin(math.sqrt(min(1.0, h)))


def elevation_angle_deg(distance_km_, height_m, eye_m=0.0, k=K_REFRACTION):
    ''' Apparent elevation of a point, with curvature and refraction. '''
    Re = effective_radius_km(k)
    return (math.degrees((height_m - eye_m) / (distance_km_ * 1000.0))
            - math.degrees(distance_km_ / (2 * Re)))


# --------------------------------------------------------------------------- #
# The thing that actually decides it
# --------------------------------------------------------------------------- #

def sensitivity(features, eye_m, k=K_REFRACTION):
    ''' How much the observables move per kilometre of observer motion.

        `features` : iterable of (distance_km, height_m).  Bearings are not
        needed -- what matters is the range spread.

        `eye_m` is DELIBERATELY REQUIRED.  Defaulting it to zero silently
        computes the elevation of every feature as seen from sea level, which
        for a 1770 m tower viewed from a 1693 m hilltop overstates the radial
        sensitivity by a factor of twenty -- and the answer still looks
        plausible, which is worse than a crash.

        Returns dict with, for lateral and radial motion:
          raw       the per-feature sensitivity, deg/km, of the closest feature
          absorbed  the SPREAD across features, deg/km -- what survives once a
                    common compass/pitch bias has eaten the mean
          ratio     raw/absorbed, i.e. how much leverage the biases destroy

        LATERAL moves bearing: 1 km across the line of sight at range d swings
        the bearing by degrees(1/d).  RADIAL moves elevation: closing 1 km on a
        peak raises it by d(elev)/d(range).

        The absorbed figures are the honest ones.  Quoting the raw lateral
        sensitivity is how a single-ridge view comes to look 50x better than it
        is.
    '''
    feats = [(float(d), float(h)) for d, h in features]
    if len(feats) < 1:
        raise ValueError("need at least one feature")
    lat_s = np.array([math.degrees(1.0 / d) for d, _ in feats])
    rad_s = np.array([abs(elevation_angle_deg(d - 0.5, h, eye_m, k)
                          - elevation_angle_deg(d + 0.5, h, eye_m, k))
                      for d, h in feats])
    lat_abs = float(lat_s.max() - lat_s.min()) if len(feats) > 1 else 0.0
    rad_abs = float(rad_s.max() - rad_s.min()) if len(feats) > 1 else 0.0
    return dict(
        ranges_km=[d for d, _ in feats],
        lateral_raw=float(lat_s.max()), lateral_absorbed=lat_abs,
        radial_raw=float(rad_s.max()), radial_absorbed=rad_abs,
        lateral_ratio=(float(lat_s.max()) / lat_abs) if lat_abs > 0 else float("inf"),
        radial_ratio=(float(rad_s.max()) / rad_abs) if rad_abs > 0 else float("inf"))


def focal_absorbed_radial(features, eye_m, k=K_REFRACTION):
    ''' Radial sensitivity that survives an UNKNOWN FOCAL LENGTH, deg/km.

        A third nuisance parameter, and the most destructive of the three.  If
        the pixels-per-degree is not known -- a cropped image, a screenshot, a
        stock photograph, anything without EXIF -- then every measured angle
        carries an unknown common scale.  Moving the camera toward the terrain
        scales all the apparent angles up together, which is *precisely* what
        reducing the focal length does, so the two are degenerate to first order
        and only the departure from proportionality survives.

        Measured on the Lake Tahoe wallpaper: with the focal length free, moving
        the camera EIGHT KILOMETRES changed the skyline residual by 0.1 arcmin
        while the fitted scale simply rescaled to compensate.  Position was not
        merely imprecise, it was absent.

        This is a FIRST-ORDER statement and must not be read as "unsolvable".
        Bounding the scale physically -- see `focal_bounds_from_relief` -- leaves
        enough second-order signal to localise the same scene to 0.3 km on a 19
        degree field, where the same search with the scale free railed and
        localised nothing.  Quantify the degeneracy with this function; do not
        let it talk you into an unrestricted search.

        Computed as the residual of regressing each feature's d(elevation)/
        d(range) against its elevation -- the component a common scale cannot
        explain.
    '''
    feats = [(float(d), float(h)) for d, h in features]
    if len(feats) < 3:
        return 0.0
    e = np.array([elevation_angle_deg(d, h, eye_m, k) for d, h in feats])
    r = np.array([elevation_angle_deg(d - 0.5, h, eye_m, k)
                  - elevation_angle_deg(d + 0.5, h, eye_m, k) for d, h in feats])
    denom = float(e @ e)
    if denom <= 0:
        return 0.0
    lam = float(r @ e) / denom          # best common scale
    return float(np.std(r - lam * e))


def focal_bounds_from_relief(relief_px, terrain_reliefs_deg, margin=0.25):
    ''' Bracket the unknown pixels-per-degree from the skyline's own relief.

        `focal_absorbed_radial` says an unknown scale destroys radial position.
        It does not say the scale is unknowable.  A search that lets the focal
        length run free will rail it against whichever bound flatters the wrong
        answer -- on the Tahoe wallpaper, f free over 60-520 px/deg pinned at 510
        and returned scattered positions at 2.2-2.7 arcmin, no minimum at all.
        Bounding it PHYSICALLY, on the other hand, cost nothing and let the same
        search converge to 400 m.

        The bound is free because the extraction already measured it.  The
        vertical span of the extracted crest, in pixels, divided by the angular
        span the terrain can plausibly subtend, in degrees, IS the scale.  The
        angular span comes from the DEM: render the horizon at a handful of
        candidate positions and take the spread of (max elevation - min).

        `relief_px`         peak-to-trough span of the extracted skyline, pixels
        `terrain_reliefs_deg`  angular reliefs rendered across candidates, deg
        `margin`            fractional widening, to cover positions not sampled

        Returns (f_lo, f_hi) in pixels per degree.  This is a BRACKET, not an
        estimate: it is only as good as the assumption that the extracted crest
        and the rendered arc cover the same features, so mask foreground before
        measuring the relief -- a near rock point will dominate the span and
        drive the bracket low, which biases the position outward.

        The two sides are not equally useful.  The UPPER bound is what does the
        work, because railing to a long focal length is the failure mode, and it
        is set by the SMALLEST relief in `terrain_reliefs_deg` -- a robust number,
        since it only asserts that no candidate position sees less than that.  The
        lower bound is set by the largest relief and is usually far too loose to
        bind.  Over the whole Tahoe candidate set the bracket came out 5-205
        px/deg -- a 40-fold span, lower bound 30x below the answer -- and the
        search still converged to within 0.3 km of the hand-tuned result.  A loose
        bracket is not a weak one, so widen rather than tune: the value comes from
        excluding the unphysical long-focal rail, not from being tight.
    '''
    relief_px = float(relief_px)
    r = np.asarray([float(v) for v in terrain_reliefs_deg], dtype=float)
    r = r[np.isfinite(r) & (r > 0)]
    if relief_px <= 0 or r.size == 0:
        raise ValueError('need a positive pixel relief and at least one '
                         'positive terrain relief')
    lo_deg = float(r.min()) * (1.0 - margin)
    hi_deg = float(r.max()) * (1.0 + margin)
    if lo_deg <= 0:
        raise ValueError('margin too large: lower relief bound is non-positive')
    # bigger angular relief for the same pixels means a SHORTER focal length,
    # so the degree bounds invert on the way through.
    return relief_px / hi_deg, relief_px / lo_deg


def waterline_range(depression_deg, eye_m, k=K_REFRACTION, d_max_km=60.0,
                    n=6000):
    ''' Range to a far shoreline from how far its waterline sits below level.

        The far edge of a lake or bay is a point at the observer's OWN water
        level, so its depression is pure geometry -- no DEM, no landmark
        identification, no scale.  That makes it the one observable in a terrain
        photograph that can break the focal-length degeneracy the way Denver's
        near tower did, which is why it is worth having.

        It comes with a trap.  The depression is the sum of an eye-height term
        that FALLS with range (h/d) and a curvature term that GROWS with it
        (d/2R_eff), so it is not monotonic: it reaches a minimum magnitude at
        sqrt(2 h R_eff) and increases on both sides.  Two consequences, both
        returned:

        * at that range the derivative is zero and the measurement carries NO
          range information at all -- a blind spot, 4.8 km for a standing
          observer at 1.6 m, 12 km from a 10 m deck;
        * every depression except the extremum matches TWO ranges, one inside
          the blind range and one outside, and nothing in the angle itself says
          which.

        Returns dict(ranges_km, blind_range_km, blind_depression_deg,
        sensitivity_deg_per_km) with one entry in `ranges_km` per solution, in
        increasing order.  An empty `ranges_km` means the depression is smaller
        than anything achievable at this eye height -- usually a sign the
        horizon reference, not the range, is wrong.
    '''
    eye_m = float(eye_m)
    if eye_m <= 0:
        raise ValueError('eye height must be above the water surface')
    d = np.linspace(d_max_km / n, d_max_km, n)
    a = np.array([elevation_angle_deg(x, 0.0, eye_m, k) for x in d])
    i = int(np.argmax(a))               # least depressed
    want = float(depression_deg)
    out = []
    for lo, hi in ((0, i + 1), (i, n)):
        seg_d, seg_a = d[lo:hi], a[lo:hi]
        if seg_d.size < 2:
            continue
        j = int(np.argmin(np.abs(seg_a - want)))
        # only accept a genuine crossing, not the nearest point on a branch
        # that never reaches the requested angle
        if seg_a.min() - 1e-9 <= want <= seg_a.max() + 1e-9:
            out.append(float(seg_d[j]))
    out = sorted(set(round(v, 6) for v in out))
    grad = [float((elevation_angle_deg(v + 0.05, 0.0, eye_m, k)
                   - elevation_angle_deg(v - 0.05, 0.0, eye_m, k)) / 0.1)
            for v in out]
    return dict(ranges_km=out, blind_range_km=float(d[i]),
                blind_depression_deg=float(a[i]), sensitivity_deg_per_km=grad)


def position_dilution(features, sigma_deg, eye_m, k=K_REFRACTION,
                      biases_free=True, focal_free=False):
    ''' Expected 1-sigma position error, km, from a given angular accuracy.

        `sigma_deg` is the accuracy of the MATCH, not of a pixel -- on real
        photographs it is dominated by systematics (extraction bias, DEM error,
        a drifting pitch bias), so use the residual you actually achieve, not
        the pixel scale.

        With `biases_free` (the realistic case) the absorbed sensitivities are
        used; set it False to see what a perfectly calibrated compass and
        horizon would buy.

        `focal_free=True` additionally admits that the pixels-per-degree is
        unknown -- a crop, a screenshot, a stock photograph.  Do not skip this
        when it applies: an unknown scale is degenerate with radial position and
        it is far more damaging than either bias.

        Returns semi-axes in km and their ratio -- an ellipse much longer than
        it is wide means a LINE of position, and one line is not a fix.
    '''
    s = sensitivity(features, eye_m, k=k)
    lat = s["lateral_absorbed"] if biases_free else s["lateral_raw"]
    rad = s["radial_absorbed"] if biases_free else s["radial_raw"]
    if focal_free:
        # An unknown pixels-per-degree scales every angle together, which is what
        # a radial move does -- so only the non-proportional part is left.
        rad = min(rad, focal_absorbed_radial(features, eye_m, k=k))
        s = dict(s, radial_focal_free=rad)
    across = sigma_deg / lat if lat > 0 else float("inf")
    along = sigma_deg / rad if rad > 0 else float("inf")
    return dict(across_km=across, along_km=along,
                elongation=(along / across) if across > 0 else float("inf"),
                sensitivity=s)


def verdict(features, sigma_deg, eye_m, target_km=1.0, **kw):
    ''' A blunt go/no-go, with the numbers that justify it. '''
    d = position_dilution(features, sigma_deg, eye_m, **kw)
    ok = d["across_km"] <= target_km and d["along_km"] <= target_km
    line = d["across_km"] <= target_km < d["along_km"]
    if ok:
        msg = (f"fix to about {max(d['across_km'], d['along_km']):.2f} km "
               f"({d['across_km']:.2f} across x {d['along_km']:.2f} along)")
    elif line:
        msg = (f"LINE OF POSITION only: {d['across_km']:.2f} km across, "
               f"{d['along_km']:.0f} km along. Needs a second bearing.")
    else:
        msg = (f"not usable: {d['across_km']:.1f} x {d['along_km']:.0f} km. "
               f"Range spread is {min(d['sensitivity']['ranges_km']):.0f}"
               f"-{max(d['sensitivity']['ranges_km']):.0f} km; add a nearer "
               f"landmark.")
    return dict(usable=ok, line_of_position=line, message=msg, **d)


# --------------------------------------------------------------------------- #
# One angle between two known points -> a circle of position
# --------------------------------------------------------------------------- #

def circle_of_position(p_near, p_far, angle_deg, centre, span_km=40.0,
                       step_km=0.25):
    ''' Locus of observers seeing `angle_deg` between two known points.

        This is the classic horizontal-angle station pointer, and it is what
        actually solved the Denver photograph: the angle between a downtown tower
        at 15 km and Longs Peak at 92 km swings about 3 degrees per kilometre of
        lateral movement, so one measured angle pins the observer to a strip a
        couple of hundred metres wide.

        `angle_deg` is signed: bearing(far) - bearing(near) at the observer,
        wrapped to +/-180.

        SOLVED NUMERICALLY, NOT BY THE INSCRIBED-ANGLE CIRCLE.  The textbook
        construction -- radius |AB| / (2 sin theta) -- is a PLANAR theorem, and
        bearings are spherical.  Across a 78 km chord the two disagree by about
        0.3 degrees, which at 3 deg/km is a 100 metre bias, and worse, the
        constructed arc runs through one of the two landmarks where the bearing
        degenerates and the angle swings through 180.  Evaluating the true
        spherical angle on a grid and contouring it has neither problem and costs
        milliseconds.

        Returns (lats, lons) along the locus within `span_km` of `centre`.
        One locus is one LINE of position -- cross it with a second for a fix,
        as in `landfall.two_landmark_circle_fix`.
    '''
    want = ((float(angle_deg) + 180.0) % 360.0) - 180.0
    la_n, lo_n = float(p_near[0]), float(p_near[1])
    la_f, lo_f = float(p_far[0]), float(p_far[1])
    clat, clon = float(centre[0]), float(centre[1])
    ky = 110.574
    kx = 111.320 * math.cos(math.radians(clat))
    n = max(8, int(2 * span_km / step_km) + 1)
    dn = np.linspace(-span_km, span_km, n)
    lats = clat + dn / ky
    lons = clon + dn / kx

    def sep(la, lo):
        return ((bearing_deg(la, lo, la_f, lo_f)
                 - bearing_deg(la, lo, la_n, lo_n) + 180.0) % 360.0) - 180.0

    G = np.empty((n, n))
    for i, la in enumerate(lats):
        for j, lo in enumerate(lons):
            G[i, j] = sep(la, lo)
    F = G - want

    def refine(a_la, a_lo, b_la, b_lo, iters=40):
        # BISECT, do not interpolate.  The angle changes by ~3 deg per km here,
        # so a single linear step across a 250 m cell still leaves ~0.6 deg of
        # error -- half the signal.  Bisection costs 40 bearing evaluations and
        # lands on the locus exactly.
        fa = sep(a_la, a_lo) - want
        for _ in range(iters):
            m_la, m_lo = 0.5 * (a_la + b_la), 0.5 * (a_lo + b_lo)
            fm = sep(m_la, m_lo) - want
            if fa * fm <= 0:
                b_la, b_lo = m_la, m_lo
            else:
                a_la, a_lo, fa = m_la, m_lo, fm
        return 0.5 * (a_la + b_la), 0.5 * (a_lo + b_lo)

    out_la, out_lo = [], []
    # crossings along each row, then each column: a connected locus without
    # needing a contouring library
    for i in range(n):
        srow = F[i]
        sign = np.sign(srow)
        idx = np.where((sign[:-1] * sign[1:] < 0) & (np.abs(srow[:-1]) < 90))[0]
        for j in idx:
            la_, lo_ = refine(lats[i], lons[j], lats[i], lons[j + 1])
            out_la.append(la_); out_lo.append(lo_)
    for j in range(n):
        scol = F[:, j]
        sign = np.sign(scol)
        idx = np.where((sign[:-1] * sign[1:] < 0) & (np.abs(scol[:-1]) < 90))[0]
        for i in idx:
            la_, lo_ = refine(lats[i], lons[j], lats[i + 1], lons[j])
            out_la.append(la_); out_lo.append(lo_)
    la_arr, lo_arr = np.array(out_la), np.array(out_lo)
    if not len(la_arr):
        return la_arr, lo_arr
    # FINAL VERIFICATION.  A crossing found next to one of the two landmarks is
    # spurious: the bearing to it swings through 180 there, so the sign of the
    # residual flips for a reason that has nothing to do with the locus.  Rather
    # than special-case the geometry, re-evaluate every surviving point and keep
    # only those that genuinely reproduce the angle.
    keep = np.array([abs((((sep(a, b) - want) + 180.0) % 360.0) - 180.0) < 1e-3
                     for a, b in zip(la_arr, lo_arr)])
    la_arr, lo_arr = la_arr[keep], lo_arr[keep]
    if not len(la_arr):
        return la_arr, lo_arr
    d = np.hypot((la_arr - clat) * ky, (lo_arr - clon) * kx)
    m = d <= span_km
    la_arr, lo_arr = la_arr[m], lo_arr[m]
    order = np.argsort(la_arr)
    return la_arr[order], lo_arr[order]


def angle_from_pixels(x_near, x_far, px_per_deg):
    ''' Signed angle between two features from their image columns. '''
    return (float(x_far) - float(x_near)) / float(px_per_deg)


# --------------------------------------------------------------------------- #
# Correlated residuals: why more pixels stop helping
# --------------------------------------------------------------------------- #

def effective_samples(x, residual, max_lag_frac=0.25):
    ''' Number of INDEPENDENT samples in a residual series.

        On a real photograph the skyline residual is not white.  It is smooth in
        azimuth -- extraction bias, haze-dependent edge placement, DEM error, a
        pitch bias that drifts between frames -- and correlated error is exactly
        what a position shift looks like, so it never averages down.

        Measured on the Istanbul frames: 839 samples, 0.14 deg residual.  If that
        had been white noise the lateral position would have come out to ~15 m.
        It came out to nothing.  The effective count was of order the number of
        FRAMES, not of pixels.

        Estimates the integrated autocorrelation time on a uniform grid and
        returns dict(n, n_eff, correlation_length, sigma_inflation).  Multiply a
        naive sigma by `sigma_inflation` to get an honest one.
    '''
    x = np.asarray(x, float)
    r = np.asarray(residual, float)
    ok = np.isfinite(x) & np.isfinite(r)
    x, r = x[ok], r[ok]
    n = len(r)
    if n < 16:
        return dict(n=n, n_eff=float(n), correlation_length=0.0,
                    sigma_inflation=1.0)
    order = np.argsort(x)
    x, r = x[order], r[order]
    grid = np.linspace(x[0], x[-1], n)
    r = np.interp(grid, x, r)
    r = r - r.mean()
    var = float((r * r).mean())
    if var <= 0:
        return dict(n=n, n_eff=float(n), correlation_length=0.0,
                    sigma_inflation=1.0)
    max_lag = max(2, int(n * max_lag_frac))
    tau = 1.0
    for lag in range(1, max_lag):
        c = float((r[:-lag] * r[lag:]).mean()) / var
        if c <= 0.0:
            break
        tau += 2.0 * c * (1.0 - lag / n)
    tau = max(tau, 1.0)
    n_eff = max(1.0, n / tau)
    return dict(n=n, n_eff=float(n_eff), correlation_length=float(tau),
                sigma_inflation=float(math.sqrt(n / n_eff)))


# --------------------------------------------------------------------------- #
# Straight from a DEM
# --------------------------------------------------------------------------- #

def scene_report(dem, lat, lon, eye_m, az_start, az_end, sigma_deg=0.05,
                 az_step=0.05, d_max_km=150.0, d_step_km=0.3,
                 k=K_REFRACTION, min_prominence=0.04, extra_landmarks=()):
    ''' Look at what a viewpoint can actually see, and rule on it.

        Ray-marches the sector, records the RANGE of whatever sets the horizon at
        each azimuth, and feeds the spread to `verdict`.  Call this before
        planning a shoot: it costs a second and it is the difference between a
        panorama that fixes your position and one that cannot.

        `extra_landmarks` is (lat, lon, height_m) for things a DEM does not
        contain -- towers, masts, chimneys, a lighthouse.  This matters more than
        it sounds.  On the Denver photograph the terrain alone scored 0.49 km
        across by 15 km along; adding one downtown tower at 15 km against peaks
        at 92 km improved the lateral leverage by nearly a hundredfold, and that
        tower is what actually solved the picture.  If your near object is
        man-made, the DEM cannot see it and you must say so here.
    '''
    from .terrain_resection import render_skyline, skyline_peaks
    from .tools.fetch_peak_names import summit_position
    Re = effective_radius_km(k)
    azs, els = render_skyline(dem, lat, lon, eye_m, az_start=az_start,
                              az_end=az_end, az_step=az_step,
                              d_max_km=d_max_km, d_step_km=d_step_km)
    peaks = skyline_peaks(azs, els, window=int(1.0 / az_step),
                          min_prominence=min_prominence)
    d = np.arange(1.0, d_max_km, d_step_km)
    feats = []
    for a, _ in peaks:
        la = np.empty(len(d)); lo = np.empty(len(d))
        for j, dd in enumerate(d):
            la[j], lo[j] = summit_position(lat, lon, a, dd)
        h = dem.elevation(la, lo)
        e = np.degrees((h - eye_m) / (d * 1000.0)) - np.degrees(d / (2 * Re))
        i = int(np.argmax(e))
        feats.append(dict(azimuth=float(a), distance_km=float(d[i]),
                          height_m=float(h[i]), elevation_deg=float(e[i])))
    for la_, lo_, h_ in extra_landmarks:
        dd = distance_km(lat, lon, la_, lo_)
        feats.append(dict(azimuth=bearing_deg(lat, lon, la_, lo_),
                          distance_km=dd, height_m=float(h_),
                          elevation_deg=elevation_angle_deg(dd, h_, eye_m, k),
                          man_made=True))
    if not feats:
        return dict(features=[], relief_deg=float(np.ptp(els)),
                    usable=False, line_of_position=False,
                    message="no prominent skyline features in this sector")
    v = verdict([(f["distance_km"], f["height_m"]) for f in feats],
                sigma_deg, eye_m, k=k)
    v.update(features=feats, relief_deg=float(np.ptp(els)),
             sector_deg=float(az_end - az_start))
    return v
