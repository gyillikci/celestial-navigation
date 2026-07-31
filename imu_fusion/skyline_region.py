''' Region-integrating skyline extraction: the Baatz/Saurer data term.

    WHY THIS EXISTS.  Benchmarking on CH1 (the ETH dataset behind "Large Scale
    Visual Geo-Localization of Images in Mountainous Terrain") showed our
    `skyline_dp` extractor disagreeing with the curated masks by 67-120 px on
    layered alpine frames, while the matcher fed with those masks solved 6/10
    scenes to a median 383 m.  Extraction, not matching, was the weak stage.

    Reading their paper (bundled with the dataset) explains it in one sentence:

        "To obtain the data term for a candidate height in a column we sum all
         foreground costs below the candidate contour and sky costs above the
         contour, where we have trained a pixel's color and gradient likelihoods
         for ground and sky."

    That is a REGION-INTEGRATING cost: every pixel of the column votes on where
    the boundary is.  `skyline_dp` uses a LOCAL STEP cost -- a matched filter a
    few pixels tall -- so only pixels adjacent to the candidate boundary vote.
    The difference is exactly the CH1 failure mode: a strong internal edge (a
    snow line, a shaded ridge below the true crest, lens flare) beats the true
    boundary locally, and a local detector has no way to know better.  A region
    cost cannot be fooled the same way, because placing the boundary at an
    internal edge means declaring a large block of mountain to be sky and paying
    for that over hundreds of pixels.

    THE ALGEBRA THAT MAKES IT CHEAP.  Writing c_sky and c_gnd for the negative
    log-likelihoods of a pixel under each class, the cost of a boundary at row r
    in column c is

        sum_{y<r} c_sky(y,c) + sum_{y>=r} c_gnd(y,c)
      = sum_{y<r} [c_sky - c_gnd](y,c)  +  sum_{all y} c_gnd(y,c)

    and the second term is constant per column, so it cannot change which path
    wins.  The whole data term is therefore ONE cumulative sum of the
    log-likelihood ratio, evaluated in O(1) per candidate -- no slower than the
    local cost it replaces.

    WHAT IS NOT COPIED.  Their smoothness term is steered by a dehazing-derived
    depth gradient, and they allow corrective user strokes.  Neither is
    implemented here; this module reproduces the data term only, which is the
    part the CH1 evidence says matters.  Their code is not used -- the method is
    reimplemented from the paper's description.

    (c) 2026.  MIT License (see LICENSE file).
'''

import numpy as np

from .skyline_dp import trace


class ColorGradientModel:
    ''' Pixel likelihoods for sky and ground, from colour and gradient.

        Deliberately simple and non-parametric: a quantised RGB histogram plus a
        gradient-magnitude histogram, treated as independent and Laplace
        smoothed.  The paper says only "trained a pixel's color and gradient
        likelihoods"; this is the least-assuming way to honour that.
    '''

    def __init__(self, color_bins=10, grad_bins=12, grad_cap=60.0):
        self.cb, self.gb, self.gcap = int(color_bins), int(grad_bins), float(grad_cap)
        self.log_c = None      # (2, cb, cb, cb)
        self.log_g = None      # (2, gb)

    # -- features ---------------------------------------------------------- #
    def _index(self, rgb):
        q = np.clip((np.asarray(rgb, float) / 256.0 * self.cb).astype(int),
                    0, self.cb - 1)
        return q[..., 0], q[..., 1], q[..., 2]

    def _grad_index(self, rgb):
        g = np.asarray(rgb, float).mean(axis=2)
        gy, gx = np.gradient(g)
        mag = np.hypot(gx, gy)
        return np.clip((mag / self.gcap * self.gb).astype(int), 0, self.gb - 1)

    # -- training ---------------------------------------------------------- #
    def fit(self, images, sky_masks, alpha=1.0):
        ''' `sky_masks[i]` is True where the pixel is SKY. '''
        hc = np.full((2, self.cb, self.cb, self.cb), float(alpha))
        hg = np.full((2, self.gb), float(alpha))
        for img, sky in zip(images, sky_masks):
            img = np.asarray(img, float)
            r, g, b = self._index(img)
            gi = self._grad_index(img)
            sky = np.asarray(sky, bool)
            for cls, sel in ((0, sky), (1, ~sky)):
                np.add.at(hc[cls], (r[sel], g[sel], b[sel]), 1.0)
                np.add.at(hg[cls], (gi[sel],), 1.0)
        self.log_c = np.log(hc / hc.sum(axis=(1, 2, 3), keepdims=True))
        self.log_g = np.log(hg / hg.sum(axis=1, keepdims=True))
        return self

    # -- inference --------------------------------------------------------- #
    def log_likelihood_ratio(self, image):
        ''' log P(pixel | ground) - log P(pixel | sky), per pixel.

            Positive where the pixel looks like ground.  This single array IS
            the data term (see the module docstring).
        '''
        if self.log_c is None:
            raise RuntimeError("model not fitted")
        img = np.asarray(image, float)
        r, g, b = self._index(img)
        gi = self._grad_index(img)
        return ((self.log_c[1][r, g, b] - self.log_c[0][r, g, b])
                + (self.log_g[1][gi] - self.log_g[0][gi]))


def region_cost(llr, jump_scale=1.0):
    ''' Cumulative region cost: cost[r, c] for a boundary at row r of column c.

        `llr` is `ColorGradientModel.log_likelihood_ratio`.  Rows above the
        boundary are claimed to be sky, so each ground-looking pixel above it
        costs; the constant per-column term is dropped because it cannot change
        which path wins.
    '''
    llr = np.asarray(llr, float)
    cost = np.cumsum(llr, axis=0)
    # normalise per column so the jump penalty means the same thing everywhere
    span = np.ptp(cost, axis=0)
    span = np.where(span > 1e-9, span, 1.0)
    return (cost - cost.min(axis=0)) / span * float(jump_scale)


def extract(image, model, max_jump=6, jump_penalty=1.0, y_band=None):
    ''' Skyline row per column, by minimum-cost path over the region cost.

        Same path machinery as `skyline_dp.extract`; only the data term differs,
        which is the whole point of the comparison.
    '''
    llr = model.log_likelihood_ratio(image)
    cost = region_cost(llr)
    if y_band is not None:
        lo, hi = y_band
        rows = np.arange(cost.shape[0])[:, None]
        blocked = (rows < np.asarray(lo)[None, :]) | (rows > np.asarray(hi)[None, :])
        cost = np.where(blocked, cost.max() + 10.0, cost)
    return trace(cost, max_jump=max_jump, jump_penalty=jump_penalty)[0].astype(float)


def extract_hybrid(image, model, edge_weight=0.5, max_jump=6, jump_penalty=1.0,
                   y_band=None, smooth=9):
    ''' Region cost PLUS a local edge term -- the fix for the region cost's
        own failure mode.

        A crude colour/gradient model (this is a research-grade reimplementation
        of a description the source paper only "sketches roughly", not their
        trained system) can be so noisy on hazy or overcast sky that the
        cumulative sum has no clean minimum and the path collapses to a frame
        edge: measured on one CH1 scene, total failure (every column at row 0)
        against curated-mask truth.  A pure local edge cost has the opposite
        weakness -- it is confidently wrong at a strong internal edge, which is
        the CH1 failure mode in the first place.

        Blending both means the path must satisfy two independent signals at
        once: a real boundary is where GLOBAL region statistics change AND a
        LOCAL edge exists.  Neither alone is enough on hard alpine frames.
    '''
    from .skyline_dp import step_cost
    llr = model.log_likelihood_ratio(image)
    rc = region_cost(llr)
    gray = np.asarray(image, float).mean(axis=2)
    ec, _ = step_cost(gray, polarity=-1.0, smooth=smooth)
    ec = (ec - ec.min()) / max(ec.ptp(), 1e-9)
    cost = (1.0 - edge_weight) * rc + edge_weight * ec
    if y_band is not None:
        lo, hi = y_band
        rows = np.arange(cost.shape[0])[:, None]
        blocked = (rows < np.asarray(lo)[None, :]) | (rows > np.asarray(hi)[None, :])
        cost = np.where(blocked, cost.max() + 10.0, cost)
    return trace(cost, max_jump=max_jump, jump_penalty=jump_penalty)[0].astype(float)


def extract_guarded(image, model, max_jump=6, jump_penalty=1.0, y_band=None,
                    collapse_frac=0.6, edge_margin_px=4, smooth=9):
    ''' extract_hybrid, with a fallback for the region cost's total-failure mode.

        Measured on CH1: a hazy/overcast sky can make the colour+gradient model
        so noisy that the cumulative sum has NO interior minimum and the path
        pins to a frame edge for nearly every column -- one scene, complete
        failure (max error 296 px, the worst of any extractor tried), and
        `edge_weight` alone could not rescue it because a pinned region term
        still outvotes the edge term almost everywhere.

        The fix is a cheap, checkable SYMPTOM rather than a smarter model:
        count how many columns of the hybrid path sit within `edge_margin_px`
        of the search band's top or bottom.  A real skyline is essentially
        never systematically hugging a frame edge across the width; if more
        than `collapse_frac` of columns do, the region term is untrustworthy
        for this image and the pure local-edge extractor is used instead.
    '''
    band = y_band
    hyb = extract_hybrid(image, model, edge_weight=0.5, max_jump=max_jump,
                         jump_penalty=jump_penalty, y_band=band, smooth=smooth)
    nr = np.asarray(image).shape[0]
    lo = 0.0 if band is None else np.broadcast_to(np.asarray(band[0], float),
                                                   hyb.shape)
    hi = float(nr - 1) if band is None else np.broadcast_to(
        np.asarray(band[1], float), hyb.shape)
    pinned = (np.abs(hyb - lo) <= edge_margin_px) | (np.abs(hyb - hi) <= edge_margin_px)
    if float(np.mean(pinned)) > collapse_frac:
        from .skyline_dp import extract as local_extract
        gray = np.asarray(image, float).mean(axis=2)
        return local_extract(gray, polarity=-1.0, smooth=smooth,
                             max_jump=max_jump, jump_penalty=jump_penalty,
                             y_band=band)
    return hyb
