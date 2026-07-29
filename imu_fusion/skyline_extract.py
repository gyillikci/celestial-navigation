''' Pull the terrain/sky boundary out of a photograph, robustly.

    Three real photographs in this study needed three different answers, and the
    differences were not cosmetic -- each one silently destroyed the previous
    method:

      * **Hazy ridges over the Sea of Marmara.**  The ridge is ~30 DN darker than
        the sky in RED but barely 8 DN in blue, because the haze washing it out
        is itself blue.  A luminance detector sees almost nothing; the red
        channel sees it easily.
      * **The Front Range above Denver.**  Here a single channel fails in BOTH
        directions at once -- the snowfields are BRIGHTER than the sky (R 250 vs
        175) and the hazed foothills are BLUER (B-R 123 vs 68).  Any threshold
        picks one and misses the other.  What works is modelling the sky and
        flagging departure in any direction.
      * **An AR app's overlay on top of either.**  Translucent panels cancel out
        if the sky reference is per-column, but their sharp horizontal top edges
        are indistinguishable from a ridge to a "first row that darkens"
        detector -- same sign, comparable contrast.

    So this module offers both detectors behind one interface, plus the
    corrections that took three debugging rounds to find.  `extract` defaults to
    the sky-model method, which is the more general of the two.

    THE THREE TRAPS, all of which fail SILENTLY:

      1. `np.convolve(..., mode='same')` zero-pads.  The first smoothed samples
         of every column therefore dive toward zero and look exactly like a
         skyline edge at index 0 -- which any "not at the very top" guard then
         rejects, discarding the entire column.  It presents as terrain that
         simply cannot be detected.  `smooth_columns` edge-replicates instead.
      2. Masking every column that carries overlay furniture threw away two
         thirds of one frame, leaving two disjoint 1.5-degree windows -- and over
         1.5 degrees a distant ridge is a straight line, so the azimuth fit had
         nothing to bite on and ran to the edge of its search range.  Reject what
         is opaque, found in the data; keep the rest.
      3. A panel edge returns the IDENTICAL sub-pixel row for every one of the
         hundreds of columns it spans.  Terrain 50 km away never does that.
         Left in, such runs are worse than missing data: long, smooth and
         confident enough that least squares weights them heavily.

    (c) 2026.  MIT License (see LICENSE file).
'''

import numpy as np


def load_rgb(path):
    from PIL import Image
    return np.asarray(Image.open(path).convert("RGB"), dtype=float)


def smooth_columns(band, n):
    ''' Smooth each column of `band` by an n-tap boxcar, replicating the edges.

        Not a nicety.  Zero-padding here is the difference between a working
        detector and one that returns nothing at all -- see trap 1 above.
    '''
    if n <= 1:
        return band.astype(float)
    h = n // 2
    pad = np.vstack([np.repeat(band[:1], h, 0), band, np.repeat(band[-1:], h, 0)])
    ker = np.ones(n) / n
    out = np.empty((band.shape[0], band.shape[1]), float)
    for x in range(band.shape[1]):
        out[:, x] = np.convolve(pad[:, x], ker, mode="valid")
    return out


def _subpixel_cross(col, i, level):
    ''' Linear crossing of `level` between samples i-1 and i. '''
    a, b = col[i - 1], col[i]
    return 0.0 if b == a else (level - a) / (b - a)


def sky_model_edge(rgb, y_fit=(20, 110), y_search=(80, 400), thresh=13.0,
                   smooth=5):
    ''' Skyline as the first row whose COLOUR departs from the extrapolated sky.

        The sky is fitted per column as a linear function of row over `y_fit` --
        which captures the vertical gradient that a flat "median of the top"
        reference does not -- and extrapolated downward.  The boundary is where
        the RGB distance from that prediction first exceeds `thresh`.

        Departure in any direction counts, so snow brighter than the sky and haze
        bluer than it are both caught by the same test.  That is the property no
        single-channel threshold has.
    '''
    H, W, _ = rgb.shape
    y0, y1 = y_fit
    ys = np.arange(y0, min(y1, H), dtype=float)
    A = np.vstack([np.ones_like(ys), ys]).T
    pinv = np.linalg.pinv(A)
    ya, yb = y_search[0], min(y_search[1], H)
    yy = np.arange(ya, yb, dtype=float)
    B = np.vstack([np.ones_like(yy), yy]).T
    out = np.full(W, np.nan)
    for x in range(W):
        coef = pinv @ rgb[y0:y0 + len(ys), x, :]
        d = np.sqrt(((rgb[ya:yb, x, :] - B @ coef) ** 2).sum(axis=1))
        if smooth > 1:
            h = smooth // 2
            d = np.convolve(np.r_[np.repeat(d[:1], h), d, np.repeat(d[-1:], h)],
                            np.ones(smooth) / smooth, mode="valid")
        hit = np.where(d > thresh)[0]
        if len(hit) and hit[0] > 0:
            i = hit[0]
            out[x] = ya + (i - 1) + _subpixel_cross(d, i, thresh)
    return out


def channel_drop_edge(rgb, y_band=(330, 660), drop=13.0, smooth=9,
                      channel=0, sky_min=105.0, sky_rows=22):
    ''' Skyline as the first row where one CHANNEL falls `drop` below the sky.

        Built for hazed ridges, where red carries nearly four times the contrast
        of luminance.  The sky reference is taken PER COLUMN from the top of the
        band, which is what lets this run straight through a translucent overlay:
        a uniform grey wash lifts the reference and the ridge together and
        cancels.  `sky_min` rejects columns where an opaque dark object (a pole,
        a mast) has taken the sky reference down with it.
    '''
    C = rgb[:, :, channel]
    ya, yb = y_band[0], min(y_band[1], C.shape[0])
    band = smooth_columns(C[ya:yb, :], smooth)
    n = C.shape[1]
    out = np.full(n, np.nan)
    for x in range(n):
        col = band[:, x]
        sky = np.median(col[:sky_rows])
        if sky < sky_min:
            continue
        hit = np.where(col < sky - drop)[0]
        if len(hit) and hit[0] > 3:
            i = hit[0]
            out[x] = ya + (i - 1) + _subpixel_cross(col, i, sky - drop)
    return out


def drop_straight_runs(y, tol=0.6, min_run=45):
    ''' Discard detections lying on a perfectly horizontal line.

        An overlay panel's top edge and a ridge look the same to an edge
        detector.  What separates them is not contrast but SHAPE: the panel
        returns the same sub-pixel row for every column it spans, and at a few
        hundred pixels per degree half a pixel is a tenth of an arcminute of
        relief -- which real terrain tens of kilometres away never holds.
    '''
    out = np.asarray(y, float).copy()
    good = np.isfinite(out)
    if good.sum() < min_run:
        return out
    idx = np.where(good)[0]
    order = np.argsort(out[idx])
    sv, idx = out[idx][order], idx[order]
    i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] - sv[i] <= tol:
            j += 1
        if j - i + 1 >= min_run:
            out[idx[i:j + 1]] = np.nan
        i = j + 1
    return out


def drop_outliers(y, jump=6.0, win=45, rounds=2):
    ''' Remove detections that jump away from their neighbours (poles, masts,
        tower tops).  Two passes: the first clears the worst spikes so the
        running median the second pass uses is trustworthy. '''
    out = np.asarray(y, float).copy()
    n = len(out)
    for _ in range(rounds):
        med = np.full(n, np.nan)
        for x in range(n):
            w = out[max(0, x - win):min(n, x + win + 1)]
            w = w[np.isfinite(w)]
            if len(w) >= 8:
                med[x] = np.median(w)
        bad = np.isfinite(out) & np.isfinite(med) & (np.abs(out - med) > jump)
        out[bad] = np.nan
    return out


def extract(path_or_rgb, method="sky_model", straight=True, outliers=True,
            **kw):
    ''' Terrain/sky boundary per column, sub-pixel, NaN where not found.

        method="sky_model"    the general one; use it unless you know better
        method="channel_drop" hazed ridge under a translucent overlay

        Returns a float array of length = image width.
    '''
    rgb = load_rgb(path_or_rgb) if isinstance(path_or_rgb, str) else np.asarray(
        path_or_rgb, float)
    if method == "sky_model":
        y = sky_model_edge(rgb, **kw)
    elif method == "channel_drop":
        y = channel_drop_edge(rgb, **kw)
    else:
        raise ValueError(f"unknown method {method!r}")
    if outliers:
        y = drop_outliers(y)
    if straight:
        y = drop_straight_runs(y)
    return y


def to_angles(y, az_centre, px_per_deg, el_centre=0.0, roll_deg=0.0,
              width=None, xc=None, yc=None):
    ''' Pixel rows -> (azimuth, elevation) in degrees.

        Small-angle/gnomonic-lite: exact enough for the few-degree fields this
        study works in, and it keeps the inverse trivially available.  `roll_deg`
        is removed first, because over a narrow field a degree of camera roll
        tilts the frame by as much as the entire terrain signal.
    '''
    y = np.asarray(y, float)
    n = width if width is not None else len(y)
    x = np.arange(len(y), dtype=float)
    if xc is None:
        xc = (n - 1) / 2.0
    if yc is None:
        yc = np.nanmedian(y)
    tilt = np.tan(np.radians(roll_deg)) * (x - xc)
    az = az_centre + (x - xc) / px_per_deg
    el = -((y - tilt) - yc) / px_per_deg + el_centre
    return az, el
