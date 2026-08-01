''' Does this extracted skyline deserve to reach the solver?

    THE LESSON THAT FORCED THIS MODULE.  On CH1, `skyline_region`'s guarded mode
    collapsed to a CONSTANT row on 7 of 10 scenes, and the collapse was invisible
    to median pixel error -- a flat line near the truth's central tendency scores
    better than a structured-but-biased track.  Worse, fed to the solver the flat
    line produced the BEST-looking residual (6.58 arcmin) at the WORST position
    (11.3 km), because a flat elevation curve is trivially consistent with any
    wide DEM plateau.

    So the gate has to run BEFORE the solver, and it must not use ground truth --
    in production there is no mask, only the photograph.  Every check here is
    computable from the image and the extracted row vector alone.

    WHAT IT CANNOT DO.  It cannot tell you the skyline is in the RIGHT PLACE; it
    can only tell you the extraction is structurally trustworthy (has shape, sits
    on real image evidence, is not pinned or degenerate).  A confident wrong
    boundary on a real edge -- our DP locking onto a nearer bright ridge -- passes
    this gate and is caught later, by the solver's separation metric.  The two
    gates are complementary and neither replaces the other.

    EMBEDDED NOTE.  All checks are O(W) or O(W*k) over the extracted row and a
    few image columns; no DEM, no search, no allocation larger than the image.
    They are meant to run on-device before spending the expensive render budget.

    (c) 2026.  MIT License (see LICENSE file).
'''

import numpy as np


# Thresholds are deliberately loose: the gate exists to reject STRUCTURAL
# failure, not to grade quality.  Anything it passes still faces the solver.
DEFAULT_LIMITS = dict(
    min_row_std_px=2.0,        # flat-line collapse (CH1: 7/10 scenes, std<1)
    max_edge_frac=0.5,         # pinning to the search band / frame edge
    min_unique_frac=0.02,      # degenerate quantisation of the path
    min_contrast=1.5,          # boundary must sit on real image evidence
    max_jump_frac=0.04,        # runaway discontinuities (occlusion, HUD, trees)
    min_valid_frac=0.6,        # enough columns actually carry a boundary
)


def _step_response(gray, rows, half=6):
    ''' Mean brightness ABOVE minus BELOW the boundary, per column.

        A real sky/terrain boundary has a consistent signed step; a path drawn
        across flat sky or flat rock does not.  Uses the same matched-filter
        idea as `skyline_dp.step_cost`, evaluated only on the chosen path.
    '''
    g = np.asarray(gray, float)
    H, W = g.shape
    r = np.clip(np.asarray(rows, float), 0, H - 1).astype(int)
    cols = np.arange(W)
    out = np.zeros(W, float)
    for d in range(1, half + 1):
        up = np.clip(r - d, 0, H - 1)
        dn = np.clip(r + d, 0, H - 1)
        out += g[up, cols] - g[dn, cols]
    return out / half


def assess(image, rows, y_band=None, limits=None, half=6):
    ''' Structural verdict on one extracted skyline.

        `image`  : HxWx3 or HxW array (only used for contrast evidence).
        `rows`   : per-column boundary row; NaN where no boundary was found.
        `y_band` : optional (lo, hi) arrays/scalars the extractor searched in;
                   edge-pinning is measured against these when given, else
                   against the frame.

        Returns a dict of metrics plus `ok` (bool) and `reasons` (list of str).
        Metrics are reported even when they pass, so a sweep can histogram them
        and the thresholds can be re-tuned against evidence rather than taste.
    '''
    lim = dict(DEFAULT_LIMITS)
    if limits:
        lim.update(limits)
    img = np.asarray(image)
    gray = img.astype(float).mean(axis=2) if img.ndim == 3 else img.astype(float)
    H, W = gray.shape
    r = np.asarray(rows, float)
    valid = np.isfinite(r)
    m = dict(width=int(W), valid_frac=float(valid.mean()))
    reasons = []

    if m["valid_frac"] < lim["min_valid_frac"]:
        reasons.append(f"only {m['valid_frac']:.0%} of columns have a boundary")
        m.update(row_std_px=0.0, unique_frac=0.0, edge_frac=1.0,
                 contrast=0.0, jump_frac=1.0, ok=False, reasons=reasons)
        return m

    rv = r[valid]
    m["row_std_px"] = float(rv.std())
    m["unique_frac"] = float(len(np.unique(np.round(rv))) / max(1, len(rv)))

    if y_band is None:
        lo = np.zeros(W); hi = np.full(W, H - 1.0)
    else:
        lo = np.broadcast_to(np.asarray(y_band[0], float), (W,))
        hi = np.broadcast_to(np.asarray(y_band[1], float), (W,))
    near = (np.abs(r - lo) <= 4.0) | (np.abs(r - hi) <= 4.0)
    m["edge_frac"] = float(np.mean(near[valid]))

    resp = _step_response(gray, np.where(valid, r, 0.0), half=half)
    m["contrast"] = float(np.median(np.abs(resp[valid])))

    d = np.abs(np.diff(rv))
    m["jump_frac"] = float(np.mean(d > max(6.0, 0.02 * H))) if len(d) else 1.0

    if m["row_std_px"] < lim["min_row_std_px"]:
        reasons.append(f"flat/degenerate path (row std {m['row_std_px']:.2f} px)")
    if m["unique_frac"] < lim["min_unique_frac"]:
        reasons.append(f"quantised path ({m['unique_frac']:.1%} unique rows)")
    if m["edge_frac"] > lim["max_edge_frac"]:
        reasons.append(f"pinned to search-band edge ({m['edge_frac']:.0%} of columns)")
    if m["contrast"] < lim["min_contrast"]:
        reasons.append(f"no image evidence at the boundary (contrast {m['contrast']:.2f})")
    if m["jump_frac"] > lim["max_jump_frac"]:
        reasons.append(f"discontinuous ({m['jump_frac']:.0%} of steps are jumps)")

    m["ok"] = not reasons
    m["reasons"] = reasons
    return m


def usable_span_deg(rows, f_px, valid=None):
    ''' Angular span of vertical structure in the extraction, degrees.

        A skyline with 0.1 deg of relief carries almost no position information
        however clean it is -- this is the extraction-side twin of
        `resection_geometry.sensitivity`, and it is what separates "extracted
        fine, worth solving" from "extracted fine, pointless".
    '''
    r = np.asarray(rows, float)
    v = np.isfinite(r) if valid is None else np.asarray(valid, bool)
    if v.sum() < 8:
        return 0.0
    rv = r[v]
    lo, hi = np.percentile(rv, [2.0, 98.0])
    return float(np.degrees(np.arctan(abs(hi - lo) / float(f_px))))
