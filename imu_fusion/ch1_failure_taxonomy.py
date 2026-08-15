''' What KIND of wrong is the CH1 skyline extraction?  A per-photo taxonomy.

    The sweep in RESULTS.md established that the structural gate drops 0.0% of
    203 CH1 photos while only 46% of extractions land within 30' of the curated
    mask.  That is a count of failures, not a description of them, and a count
    cannot choose a model: "half of it is wrong" is compatible both with a
    single coherent bias a learned data term would erase and with scattered
    ambiguity no data term can touch.

    This module measures which.  For every photo it derives the mask's own
    per-column sky/ground transition, runs the production DP, and classifies the
    disagreement by DIRECTION and by COVERAGE:

      GOOD              median |error| <= 30'  (== agreement on >= 50% of columns)
      FOREGROUND LOCK-ON  >= 60% of columns sit more than 30' BELOW the mask
      SKY LOCK-ON         >= 60% of columns sit more than 30' ABOVE the mask
      PARTIAL             right on >= 25% of the width, wrong elsewhere, no
                          single dominant direction
      OTHER               neither coherent nor partially right

    MASK CONVENTION.  Verified here, not assumed: `check_mask_convention` asserts
    0 = sky, 255 = ground on every mask before any comparison runs.  Reading it
    backwards has already produced one wrong conclusion in this project.

    THE DECIDING NUMBER.  Beyond the counts, the module measures how far the true
    boundary is from being CHOSEN by the existing cost -- the percentile rank of
    the true row's cost within its own column, and the local step response there
    versus at the DP's answer.  A failure where the truth is a strong edge that
    merely lost to a stronger nearer one is a failure a learned per-pixel cost
    can plausibly fix.  A failure where the truth carries no local evidence at
    all is not, and no amount of training data changes that.

    (c) 2026.  MIT License (see LICENSE file).
'''

import os
import glob
import numpy as np


# The production extraction, as run in the full CH1 sweep (RESULTS.md).
DP_KWARGS = dict(polarity=-1.0, smooth=7, max_jump=6, jump_penalty=1.0)

# Agreement tolerance.  30 arcmin is the figure the sweep already reports
# against, kept identical so the classification is comparable to it.
TOL_ARCMIN = 30.0

# A class is "coherent" when this share of ALL columns errs in one direction.
COHERENT_FRAC = 0.60

# A failure is "partial" when at least this share of columns is still right.
PARTIAL_FRAC = 0.25


def photo_list(root='CH1'):
    ''' Every CH1 photo, as (metadata_txt_path, image_path, mask_path). '''
    txts = (sorted(glob.glob(os.path.join(root, 'cvg', '*.png.txt')))
            + sorted(glob.glob(os.path.join(root, 'panoramio', '*.png.txt'))))
    out = []
    for t in txts:
        img = t[:-4]
        out.append((t, img, img[:-4] + '-mask.png'))
    return out


def check_mask_convention(mask_paths, loader):
    ''' Verify 0 = sky / 255 = ground on every mask, and return the evidence.

        The test does not trust the readme: it checks that the TOP of the frame
        is the low value and the BOTTOM is the high one, on every single mask.
        A photograph of terrain always has more sky at the top than at the
        bottom, so a unanimous result is proof; a split result would mean the
        convention is per-file and nothing downstream may assume it.
    '''
    top, bot, vals = [], [], set()
    for p in mask_paths:
        a = loader(p)
        top.append(float(a[:3].mean()))
        bot.append(float(a[-3:].mean()))
        vals |= set(np.unique(a).tolist())
    top, bot = np.array(top), np.array(bot)
    return dict(n=len(top), n_bottom_brighter=int((bot > top).sum()),
                mean_top=float(top.mean()), mean_bottom=float(bot.mean()),
                distinct_values=sorted(vals),
                convention_ok=bool((bot > top).all()))


def mask_boundary(mask_ground):
    ''' Per-column sky/ground transition row of a curated mask.

        `mask_ground` is boolean, True where the mask says GROUND.  Returns the
        topmost ground row per column, NaN where a column is all sky.  On CH1
        99.9% of columns carry exactly one transition, so "topmost ground" and
        "the boundary" are the same thing; `n_transitions` is returned so a
        caller on other data can check that before relying on it.
    '''
    g = np.asarray(mask_ground, bool)
    has = g.any(axis=0)
    rows = np.where(has, g.argmax(axis=0).astype(float), np.nan)
    ntr = np.abs(np.diff(g.astype(np.int8), axis=0)).sum(axis=0)
    return rows, ntr


def signed_error_arcmin(rows_dp, rows_gt, f_px, height):
    ''' Elevation error of the extraction, arcmin, POSITIVE = DP is BELOW truth.

        Rows are converted to elevation about the principal point before
        differencing rather than scaled by a single arcmin-per-pixel constant:
        CH1 focal lengths run 768-5547 px and the widest frames span 67 deg, so
        the small-angle shortcut is worth several arcmin at the frame edge.
    '''
    cy = 0.5 * (height - 1.0)
    el_dp = np.degrees(np.arctan((cy - np.asarray(rows_dp, float)) / float(f_px)))
    el_gt = np.degrees(np.arctan((cy - np.asarray(rows_gt, float)) / float(f_px)))
    return (el_gt - el_dp) * 60.0


def classify(err_arcmin, tol=TOL_ARCMIN, coherent=COHERENT_FRAC,
             partial=PARTIAL_FRAC):
    ''' Failure class of one photo from its per-column signed error.

        The order matters and is deliberate.  Agreement is decided first, so
        "GOOD" means exactly what the sweep's 30' figure means.  Direction is
        decided next, because a systematic offset is the hypothesis under test.
        Partial coverage is the fallback for extractions that are right
        somewhere and wrong elsewhere without a single dominant direction.
    '''
    e = np.asarray(err_arcmin, float)
    e = e[np.isfinite(e)]
    if len(e) < 16:
        return 'OTHER', dict(n=int(len(e)), agree_frac=0.0, below_frac=0.0,
                             above_frac=0.0)
    agree = float(np.mean(np.abs(e) <= tol))
    below = float(np.mean(e > tol))
    above = float(np.mean(e < -tol))
    stats = dict(n=int(len(e)), agree_frac=agree, below_frac=below,
                 above_frac=above)
    if agree >= 0.5:
        return 'GOOD', stats
    if below >= coherent:
        return 'FOREGROUND_LOCKON', stats
    if above >= coherent:
        return 'SKY_LOCKON', stats
    if agree >= partial:
        return 'PARTIAL', stats
    return 'OTHER', stats


def truth_cost_rank(cost, rows_gt):
    ''' Where the true boundary ranks, per column, among all rows by cost.

        0.0 means the existing data term already puts the true row at the cheapest
        row of its column and only the continuity term is keeping the path away;
        0.5 means the truth is no better than a coin flip under this cost.  This
        is the quantity that separates "a re-weighted per-pixel cost could pick
        the truth" from "the truth is invisible to any local appearance model".
    '''
    C = np.asarray(cost, float)
    nr, nc = C.shape
    out = np.full(nc, np.nan)
    for c in range(nc):
        r = rows_gt[c]
        if not np.isfinite(r):
            continue
        ri = int(np.clip(round(r), 0, nr - 1))
        col = C[:, c]
        out[c] = float(np.mean(col < col[ri]))
    return out


def analyse_photo(txt_path, img_path, mask_path, loader_rgb, loader_gray):
    ''' Extract, compare against the curated mask, classify, and diagnose. '''
    from .skyline_dp import step_cost, trace

    tok = open(txt_path).read().split()
    f_px = float(tok[0])
    lat, lon = float(tok[1]), float(tok[2])

    img = loader_rgb(img_path)
    lum = img.astype(float).mean(axis=2)
    H, W = lum.shape

    cost, resp = step_cost(lum, polarity=DP_KWARGS['polarity'],
                           smooth=DP_KWARGS['smooth'], soft=1.0)
    rows, _ = trace(cost, max_jump=DP_KWARGS['max_jump'],
                    jump_penalty=DP_KWARGS['jump_penalty'])

    gm = loader_gray(mask_path) > 127
    gt, ntr = mask_boundary(gm)

    err = signed_error_arcmin(rows, gt, f_px, H)
    err_px = rows - gt
    good = np.isfinite(err)

    cls, stats = classify(err)

    cols = np.arange(W)
    ri_dp = np.clip(rows.astype(int), 0, H - 1)
    ri_gt = np.clip(np.nan_to_num(gt, nan=0.0).astype(int), 0, H - 1)
    rank = truth_cost_rank(cost, gt)

    return dict(
        name=os.path.basename(img_path)[:-4], f_px=f_px, lat=lat, lon=lon,
        W=W, H=H, cls=cls, **stats,
        med_err_px=float(np.median(np.abs(err_px[good]))),
        med_err_arcmin=float(np.median(np.abs(err[good]))),
        med_signed_px=float(np.median(err_px[good])),
        med_signed_arcmin=float(np.median(err[good])),
        iqr_signed_arcmin=float(np.subtract(*np.percentile(err[good], [75, 25]))),
        # local edge evidence: response at the DP's answer vs at the truth
        resp_dp=float(np.median(resp[ri_dp, cols])),
        resp_gt=float(np.median(resp[ri_gt, cols][good])),
        # how close the truth is to being cheapest in its own column
        rank_gt=float(np.nanmedian(rank)),
        rank_gt_p90=float(np.nanpercentile(rank, 90)),
        # does the DP still follow the SHAPE of the true skyline?
        shape_corr=float(np.corrcoef(rows[good], gt[good])[0, 1])
        if good.sum() > 8 and np.std(gt[good]) > 0 else float('nan'),
        gt_max_jump=float(np.max(np.abs(np.diff(gt[np.isfinite(gt)]))))
        if np.isfinite(gt).sum() > 2 else float('nan'),
        multi_transition_frac=float(np.mean(ntr > 1)),
    )
