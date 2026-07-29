''' Continuity-constrained boundary extraction: one optimal path, not N guesses.

    `skyline_extract` decides each column independently and then repairs the
    damage afterwards (`drop_straight_runs`, `drop_outliers`).  That works when
    the boundary is high-contrast and fails when it is not: on the Bodrum dusk
    frames the land/water step is only ~15 units of blueness, per-column
    thresholding returned 55% coverage that scattered into open water, and in one
    frame it jumped to the ridge top entirely.

    The information being discarded is that a coastline or a skyline is
    CONTINUOUS.  The classical horizon-detection literature encodes that as a
    minimum-cost path across the columns rather than a per-column decision, and
    that is what this module does: cost = how unlike a boundary this pixel is,
    plus a penalty for moving vertically between neighbouring columns.  A weak
    but coherent edge then beats a strong but isolated one, which is exactly the
    trade the Bodrum frames need.

    The path is global, so it always returns a full-width answer -- including
    across occlusions, where it interpolates rather than admits ignorance.  Use
    the returned per-column cost to mask those; a boundary hidden behind a tree
    is a guess, and the caller should be able to tell.

    See `related_work.md`: shortest-path / dynamic-programming boundary tracing
    is the standard baseline this project should have started from.
'''

import numpy as np


def step_cost(band, polarity=1.0, smooth=9, soft=1.0):
    ''' Per-pixel cost of "the boundary is here", from a signed vertical step.

        `band` is any 2-D field in which the boundary is a step DOWN the frame:
        luminance for a bright sky over dark land, or blue-minus-red for hazy
        land over water.  `polarity` +1 means the value INCREASES below the
        boundary, -1 that it decreases.

        The response is the difference of the means over a short window below and
        above each row -- a matched filter for a step, which is far less noisy
        than a bare gradient on a hazy edge.  Cost is the negated, softened
        response, so a strong step is cheap.

        `soft` sets how sharply cost falls with response; it is a scale, not a
        threshold, which is the point -- nothing is rejected here, the path
        decides.
    '''
    b = np.asarray(band, float)
    if b.ndim != 2:
        raise ValueError('band must be 2-D (rows, cols)')
    n = int(smooth)
    if n < 1:
        raise ValueError('smooth must be >= 1')
    k = np.ones(n) / n
    pad = np.vstack([np.repeat(b[:1], n, axis=0), b, np.repeat(b[-1:], n, axis=0)])
    m = np.apply_along_axis(lambda v: np.convolve(v, k, 'same'), 0, pad)[n:-n]
    below = np.vstack([m[n:], np.repeat(m[-1:], n, axis=0)])
    above = np.vstack([np.repeat(m[:1], n, axis=0), m[:-n]])
    resp = polarity * (below - above)
    return -resp / float(soft), resp


def trace(cost, max_jump=6, jump_penalty=1.0, band=None):
    ''' Minimum-cost continuous path across columns, one row per column.

        `cost[r, c]`   cost of putting the boundary at row r of column c
        `max_jump`     largest vertical move between neighbouring columns, px
        `jump_penalty` cost per pixel of vertical movement
        `band`         optional (lo, hi) row limits, scalars or per-column arrays

        Returns (rows, per_column_cost).  `rows` is float and always full width.

        A large `jump_penalty` buys smoothness at the price of following real
        terrain; the default is deliberately mild, because the failure this
        module exists to fix is incoherent jumping, not gentle wandering.
    '''
    C = np.array(cost, dtype=float)
    if C.ndim != 2:
        raise ValueError('cost must be 2-D (rows, cols)')
    nr, nc = C.shape
    J = int(max_jump)
    if J < 1:
        raise ValueError('max_jump must be >= 1')
    if band is not None:
        lo, hi = band
        lo = np.broadcast_to(np.asarray(lo, float), (nc,))
        hi = np.broadcast_to(np.asarray(hi, float), (nc,))
        rows = np.arange(nr)[:, None]
        C = np.where((rows >= lo[None, :]) & (rows <= hi[None, :]), C, np.inf)
    if not np.isfinite(C).any():
        raise ValueError('every pixel is masked out; check the band limits')

    shifts = np.arange(-J, J + 1)
    pen = jump_penalty * np.abs(shifts).astype(float)
    D = C[:, 0].copy()
    back = np.zeros((nr, nc), dtype=np.int32)
    big = np.inf
    for c in range(1, nc):
        stack = np.full((len(shifts), nr), big)
        for i, s in enumerate(shifts):
            # candidate: previous column's row is (r - s)
            if s > 0:
                stack[i, s:] = D[:nr - s] + pen[i]
            elif s < 0:
                stack[i, :nr + s] = D[-s:] + pen[i]
            else:
                stack[i] = D + pen[i]
        best = np.argmin(stack, axis=0)
        D = stack[best, np.arange(nr)] + C[:, c]
        back[:, c] = np.arange(nr) - shifts[best]

    out = np.empty(nc)
    r = int(np.argmin(D))
    for c in range(nc - 1, -1, -1):
        out[c] = r
        r = int(back[r, c])
    per_col = np.array([cost[int(out[c]), c] for c in range(nc)], dtype=float)
    return out, per_col


def extract(band, polarity=1.0, smooth=9, soft=1.0, max_jump=6,
            jump_penalty=1.0, y_band=None, reject=None):
    ''' Convenience: step cost, then trace, then optionally mask weak columns.

        `reject` masks (returns NaN for) columns whose step response is below
        this many units.  Leave it None to keep the full path -- an unmasked path
        is still the best available guess, but it should not be mistaken for a
        measurement where the boundary was invisible.
    '''
    cost, resp = step_cost(band, polarity=polarity, smooth=smooth, soft=soft)
    rows, _ = trace(cost, max_jump=max_jump, jump_penalty=jump_penalty,
                    band=y_band)
    if reject is None:
        return rows
    idx = np.clip(rows.astype(int), 0, resp.shape[0] - 1)
    strength = resp[idx, np.arange(resp.shape[1])]
    out = rows.astype(float).copy()
    out[strength < float(reject)] = np.nan
    return out
