''' Measure a device's pitch bias against the DEM horizon, from AR-app frames.

    The study's whole premise is that an iPhone's synthetic horizon carries an
    attitude error, and that an outside reference is needed to bound it.  This
    measures that error on a real device, in the field, with a real
    instrument-app screenshot -- and the outside reference is the terrain.

    THE INPUT is a frame from an AR theodolite app (here Theodolite, on an
    iPhone, 8x zoom) which prints its own GPS position, altitude, true azimuth,
    roll and pitch onto the picture.  Everything needed is therefore in the
    frame; the metadata is transcribed into a `Frame` by hand because it is
    burned-in text, not EXIF.

    WHAT COMES OUT

      * `extract_skyline` -- the terrain/sky boundary, per column, sub-pixel.
      * `fit_elevation`   -- the elevation the frame centre was ACTUALLY pointing
        at, from the DEM horizon; minus the app's own pitch readout, that is the
        pitch bias.

    WHAT DOES NOT.  Position.  See `RESULTS.md`: an 8x zoom sees ~9.5 degrees of
    azimuth, and over that window this horizon has 0.5 degrees of relief against
    6 degrees available in the full sector, so a free azimuth search lands 20
    degrees from the truth with a BETTER residual than the truth gets.  A narrow
    field is the wrong instrument for terrain resection, however good the DEM is.

    (c) 2026.  MIT License (see LICENSE file).
'''

import numpy as np

from ..terrain_resection import DemTiles, render_skyline


class Frame:
    ''' One AR-app frame: the picture plus the numbers burned into it. '''

    def __init__(self, path, lat, lon, alt_ft, az_mils, roll_deg, pitch_deg,
                 width=2622, height=1206, px_per_deg=276.0, time=None):
        self.path, self.lat, self.lon = path, lat, lon
        self.alt_m = alt_ft * 0.3048
        # Theodolite prints mils alongside a WHOLE-degree azimuth, and the mils
        # are just that whole degree converted -- 179 deg -> 3182 mils -- so the
        # bearing is known to +/-0.5 deg, not to the 0.056 deg a mil suggests.
        self.az_deg = az_mils * 0.05625
        self.roll_deg, self.pitch_deg = roll_deg, pitch_deg
        self.w, self.h = width, height
        self.xc, self.yc = (width - 1) / 2.0, (height - 1) / 2.0
        self.f = px_per_deg
        self.time = time


def extract_skyline(path, y_top=330, y_bot=660, drop=13.0, smooth=9,
                    sky_min=105.0, jump=6.0, win=45):
    ''' Sub-pixel terrain/sky boundary row for every column, NaN where unfound.

        Uses the RED channel.  Through summer haze a distant ridge is only ~30 DN
        darker than the sky in red but barely 8 DN in blue, because the haze that
        washes out the ridge is itself blue -- so red carries nearly four times
        the contrast of a luminance image.

        The sky reference is taken PER COLUMN from the top of the band.  That is
        what lets the extractor work straight through the app's translucent
        overlay: a uniform grey wash lifts the reference and the ridge together
        and cancels.  Masking every column carrying overlay furniture instead
        threw away two thirds of the frame and, with it, the azimuth baseline.

        Only genuinely opaque obstructions are rejected, and they are found in
        the data rather than hand-drawn: a dark pole drags the sky reference
        below `sky_min`, and anything else makes the detected row jump more than
        `jump` pixels away from its neighbours.
    '''
    from PIL import Image
    R = np.asarray(Image.open(path).convert("RGB"), float)[:, :, 0]
    n = R.shape[1]
    half = smooth // 2
    band = R[y_top:y_bot, :]
    # EDGE-REPLICATE.  np.convolve(mode='same') zero-pads, so the first samples
    # of every smoothed column dive toward zero and mimic a skyline edge at
    # index 0 -- which any "not at the very top" guard then rejects, silently
    # discarding the whole column.  It looks exactly like undetectable terrain.
    pad = np.vstack([np.repeat(band[:1], half, 0), band,
                     np.repeat(band[-1:], half, 0)])
    ker = np.ones(smooth) / smooth
    y = np.full(n, np.nan)
    for x in range(n):
        col = np.convolve(pad[:, x], ker, mode="valid")
        sky = np.median(col[:22])
        if sky < sky_min:
            continue
        below = np.where(col < sky - drop)[0]
        if len(below) and below[0] > 3:
            i = below[0]
            a, b = col[i - 1], col[i]
            frac = 0.0 if b == a else (sky - drop - a) / (b - a)
            y[x] = y_top + (i - 1) + frac
    for _ in range(2):
        med = np.full(n, np.nan)
        for x in range(n):
            w = y[max(0, x - win):min(n, x + win + 1)]
            w = w[np.isfinite(w)]
            if len(w) >= 8:
                med[x] = np.median(w)
        y[np.isfinite(y) & np.isfinite(med) & (np.abs(y - med) > jump)] = np.nan
    return _drop_straight_lines(y)


def _drop_straight_lines(y, tol=0.6, min_run=45):
    ''' Discard detections lying on a perfectly horizontal line.

        The app's translucent panels have sharp horizontal top edges, and to a
        detector looking for "the row where it first gets darker" a panel edge is
        indistinguishable from a ridge -- same sign, comparable size.  The give-
        away is not contrast but SHAPE: a panel edge returns the identical
        sub-pixel row for every one of the two or three hundred columns it spans,
        and real terrain 50 km away never does that, because at 276 px/deg a
        third of a pixel is a tenth of an arcminute of relief.

        Left in, these segments are worse than missing data: they are long,
        confident, and perfectly smooth, so a least-squares fit weights them
        heavily and quietly drags the horizon toward the app's layout.
    '''
    out = y.copy()
    good = np.isfinite(out)
    if good.sum() < min_run:
        return out
    vals = out[good]
    order = np.argsort(vals)
    sv = vals[order]
    idx = np.where(good)[0][order]
    i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] - sv[i] <= tol:
            j += 1
        if j - i + 1 >= min_run:
            out[idx[i:j + 1]] = np.nan
        i = j + 1
    return out


def dem_profile(frame, dem, span=9.0, step=0.01, d_max_km=190.0):
    return render_skyline(dem, frame.lat, frame.lon, frame.alt_m,
                          az_start=frame.az_deg - span,
                          az_end=frame.az_deg + span, az_step=step,
                          d_max_km=d_max_km, d_step_km=0.25)


def fit_elevation(frame, y, azs, els, az_slack=0.5, az_step=0.01):
    ''' Elevation the frame centre actually pointed at, from the DEM horizon.

        Roll is held at the app's reading (a levelled-device roll is its most
        trustworthy output) and azimuth is allowed to move only within the
        +/-0.5 deg the whole-degree readout leaves open.  The elevation offset is
        then the only real freedom, so what comes back is a measurement of where
        the camera was looking -- not a curve fit with enough slack to say
        anything.

        Returns dict(el_centre, bias, az, rms_px, rms_deg, n, corr).  `corr` is
        the correlation between the observed and predicted PROFILE SHAPES after
        removing the offset: it says whether the terrain pattern was actually
        recognised, or only its average height matched.
    '''
    x = np.arange(frame.w, dtype=float)
    ok = np.isfinite(y)
    x, yy = x[ok], y[ok]
    tilt = np.tan(np.radians(frame.roll_deg)) * (x - frame.xc)
    obs_el = -((yy - tilt) - frame.yc) / frame.f
    best = None
    for az0 in np.arange(frame.az_deg - az_slack, frame.az_deg + az_slack, az_step):
        e = np.interp(az0 + (x - frame.xc) / frame.f, azs, els)
        r = obs_el - e
        s = float(r.std())
        if best is None or s < best[0]:
            best = (s, az0, float(r.mean()), e)
    s, az0, off, e = best
    o = obs_el - obs_el.mean()
    p = e - e.mean()
    corr = float(np.corrcoef(o, p)[0, 1]) if p.std() > 1e-9 else float("nan")
    # SIGN.  `obs_el` is measured DOWNWARD from the frame centre outward, so a
    # feature at absolute elevation e sits at obs_el = e - el_centre.  The mean
    # residual is therefore MINUS the centre's elevation, not plus it.  Flipping
    # this turns a +1.5 deg bias into +0.15 deg and makes the four frames
    # disagree, which reads like noise rather than like an error.
    el_centre = -off
    return dict(el_centre=el_centre, bias=el_centre - frame.pitch_deg, az=az0,
                rms_deg=s, rms_px=s * frame.f, n=int(ok.sum()), corr=corr,
                obs_sd=float(obs_el.std()))


def pitch_bias(frames, dem_dir="imu_fusion/dem"):
    ''' Pitch bias per frame and pooled.  Independent frames, so the scatter is
        an honest error bar rather than a formal one. '''
    dem = DemTiles(dem_dir)
    out = []
    for fr in frames:
        azs, els = dem_profile(fr, dem)
        r = fit_elevation(fr, extract_skyline(fr.path), azs, els)
        r["frame"] = fr
        out.append(r)
    b = np.array([r["bias"] for r in out])
    return out, dict(mean=float(b.mean()),
                     sd=float(b.std(ddof=1)) if len(b) > 1 else float("nan"))
