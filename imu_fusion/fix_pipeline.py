''' Layer 5 of ARCHITECTURE.md: the coarse-to-fine position search.

    Everything here is arranged around one measured asymmetry: rendering a
    candidate costs ~70x more than scoring one (273 ms against 4 ms on the
    Istanbul ultrawide search), because the DEM ray-march runs at 4-5 M
    lookups/s while the scoring interpolation runs at 50 M/s.  Three design
    rules follow, and they ARE the module:

      * NUISANCES LIVE INSIDE THE RENDER.  Heading, focal length and a free
        pitch are evaluated against an already-rendered profile, so their grids
        can be dense (101 headings x 4 focals costs 4 ms).  Position is the only
        parameter that costs a render; it is the only one searched coarsely.

      * COARSE APPROXIMATES THE PRUNING ORDER, NEVER THE ANSWER.  A pass at
        az 0.2 deg / range 0.4 km / every 40th column is 35x cheaper and still
        ranked the true winner #3 of 1492 real candidates; refining the top 50
        returned the identical rank-1 in 26 s against 414 s flat.  Because that
        claim is empirical, `solve_fix` reports `coarse_rank_of_winner` so every
        run audits it -- a winner that coarse ranks near `top_k` is a warning
        that the margin is thin, caught in the result rather than lost.

      * PLAUSIBILITY IS AN INPUT, NOT A HOPE.  The heading slack is set by the
        caller from measured magnetometer behaviour.  On the real frame the
        unconstrained search preferred rms 7.62 arcmin at 650 m by claiming a
        +7.4 deg compass error; slack held to the +0.0..+2.5 deg the
        magnetometer actually exhibited 30 s later gave 359 m at separation
        1.68x.  A better residual bought with an implausible nuisance is a
        worse answer, so the gate is a constructor argument, not a footnote.

    Pitch enters in one of three MODES, decided by extraction (layer 1):
      * measured  -- `horizon_row` given: the sea horizon was verified COLLINEAR,
        so pitch is computed from the horizon row per focal length (it depends
        on f) and never fitted;
      * fixed     -- `pitch_deg` given: an external attitude source;
      * free      -- neither given: pitch enters linearly, so its optimum is the
        residual mean, eliminated in closed form per hypothesis -- never gridded.

    (c) 2026.  MIT License (see LICENSE file).
'''

import time
from dataclasses import dataclass, field, replace
from math import atan, cos, degrees, radians, sqrt

import numpy as np

from .landfall import K_REFRACTION, effective_radius_km
from .resection_geometry import distance_km, image_ray_angles
from .terrain_resection import render_skyline


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #

@dataclass
class SkylineObservation:
    ''' One frame's extracted skyline, in pixels, plus what is known about it.

        `columns`/`rows` are the boundary samples (already cleaned by layer 1).
        Exactly one pitch mode applies -- see the module docstring; supplying
        both `horizon_row` and `pitch_deg` is refused rather than resolved.
    '''
    columns: np.ndarray
    rows: np.ndarray
    f_px: float
    cx: float
    cy: float
    roll_deg: float = 0.0
    pitch_deg: float = None            # fixed pitch, degrees
    horizon_row: float = None          # sea-horizon row at cx -> pitch measured
    eye_above_water_m: float = 2.0
    k: float = K_REFRACTION

    def __post_init__(self):
        self.columns = np.asarray(self.columns, dtype=float)
        self.rows = np.asarray(self.rows, dtype=float)
        if self.horizon_row is not None and self.pitch_deg is not None:
            raise ValueError("give horizon_row (measured) OR pitch_deg (fixed),"
                             " not both")

    @property
    def pitch_mode(self):
        if self.horizon_row is not None:
            return "measured"
        return "free" if self.pitch_deg is None else "fixed"

    def dip_deg(self):
        h = max(float(self.eye_above_water_m), 0.0)
        return degrees(sqrt(2.0 * h / 1000.0 / effective_radius_km(self.k)))

    def pitch_for(self, f_px):
        ''' Pitch at a given focal length.

            Measured mode recomputes per f because the horizon-row-to-angle
            conversion depends on it -- freezing pitch at the nominal f and then
            gridding f would quietly decouple the two.
        '''
        if self.horizon_row is not None:
            return degrees(atan((float(self.horizon_row) - self.cy) / f_px)) \
                - self.dip_deg()
        return 0.0 if self.pitch_deg is None else float(self.pitch_deg)

    def rays(self, f_px, column_step=1):
        ''' (azimuth offset, elevation) for every column_step-th sample. '''
        sel = slice(None, None, int(column_step))
        return image_ray_angles(
            self.columns[sel], self.rows[sel], f_px,
            pitch_deg=self.pitch_for(f_px), roll_deg=self.roll_deg,
            cx=self.cx, cy=self.cy)


@dataclass
class SearchPrior:
    ''' What is believed before the terrain says anything.

        The GPS coordinate CENTRES the box and never scores a candidate.  The
        heading slack is a plausibility gate: set it from measured magnetometer
        behaviour (drift across nearby frames), not from optimism -- an
        unconstrained heading lets the search buy residual with impossible
        compass errors.
    '''
    lat: float
    lon: float
    heading_deg: float
    heading_slack_deg: tuple = (-1.0, 4.0)
    heading_step_deg: float = 0.1
    focal_factors: tuple = (0.98, 1.00, 1.02, 1.04)
    box_lat_km: float = 8.0            # half-extent of the search box
    box_lon_km: float = 8.0
    cell_km: float = 0.25
    ground_range_m: tuple = (0.0, 25.0)

    def headings(self):
        lo, hi = self.heading_slack_deg
        return np.arange(self.heading_deg + lo,
                         self.heading_deg + hi + 1e-9, self.heading_step_deg)


@dataclass
class RenderGrid:
    ''' Resolution of the forward model, plus the observation subsampling. '''
    az_step_deg: float = 0.05
    d_min_km: float = 1.0
    d_max_km: float = 25.0
    d_step_km: float = 0.05
    column_step: int = 8
    water_level_m: float = 0.0
    visibility_km: float = None

    def coarsened(self, az=4.0, d=8.0, col=5.0):
        ''' The benchmarked coarse pass: 35x cheaper, winner still ranked #3. '''
        return replace(self, az_step_deg=self.az_step_deg * az,
                       d_step_km=self.d_step_km * d,
                       column_step=int(self.column_step * col))


# --------------------------------------------------------------------------- #
# Candidates
# --------------------------------------------------------------------------- #

def coastal_candidates(dem, prior: SearchPrior):
    ''' Grid cells inside the prior box whose ground height fits the observer.

        The height filter is the eye-height constraint doing real work: fixing
        the eye ABSOLUTELY (observer says "5 m above the water") and admitting
        only shore-adjacent ground improved the telephoto solve from 250 m to
        188 m -- it removes the freedom to stand a candidate on 40 m of ground
        and rescale every predicted elevation to suit.
    '''
    dlat = prior.cell_km / 111.19
    dlon = prior.cell_km / (111.19 * cos(radians(prior.lat)))
    lats = np.arange(prior.lat - prior.box_lat_km / 111.19,
                     prior.lat + prior.box_lat_km / 111.19 + 1e-12, dlat)
    lons = np.arange(
        prior.lon - prior.box_lon_km / (111.19 * cos(radians(prior.lat))),
        prior.lon + prior.box_lon_km / (111.19 * cos(radians(prior.lat)))
        + 1e-12, dlon)
    lo_g, hi_g = prior.ground_range_m
    cells = []
    for la in lats:
        hs = dem.elevation(np.full(len(lons), la), lons)
        for lo, hh in zip(lons, hs):
            if lo_g <= hh <= hi_g:
                cells.append((float(la), float(lo), float(hh)))
    return cells


# --------------------------------------------------------------------------- #
# Scoring one candidate (cheap part: runs against an existing render)
# --------------------------------------------------------------------------- #

def _render_window(obs: SkylineObservation, prior: SearchPrior, margin_deg=2.0):
    ''' Azimuth span the render must cover: every heading x every ray angle. '''
    f_min = obs.f_px * min(prior.focal_factors)
    az_off, _ = obs.rays(f_min)
    hd = prior.headings()
    return (float(hd.min() + az_off.min() - margin_deg),
            float(hd.max() + az_off.max() + margin_deg))


def score_candidate(dem, lat, lon, obs: SkylineObservation, prior: SearchPrior,
                    grid: RenderGrid, az_window):
    ''' Render once, then sweep every (heading, focal) hypothesis against it.

        Returns (rms_deg, heading_deg, focal_factor).  A free pitch is
        eliminated as the per-hypothesis residual mean -- closed form, so it
        never appears as a grid axis.
    '''
    azs, prof = render_skyline(
        dem, lat, lon, obs.eye_above_water_m + 0.0,
        az_start=az_window[0], az_end=az_window[1],
        az_step=grid.az_step_deg, d_min_km=grid.d_min_km,
        d_max_km=grid.d_max_km, d_step_km=grid.d_step_km,
        k=obs.k, water_level_m=grid.water_level_m,
        visibility_km=grid.visibility_km)
    hd = prior.headings()
    free_pitch = obs.pitch_mode == "free"
    best = None
    for factor in prior.focal_factors:
        az_off, el = obs.rays(obs.f_px * factor, grid.column_step)
        pred = np.interp((hd[:, None] + az_off[None, :]) % 360.0, azs, prof)
        resid = el[None, :] - pred
        if free_pitch:
            resid = resid - resid.mean(axis=1, keepdims=True)
        rms = np.sqrt((resid * resid).mean(axis=1))
        j = int(np.argmin(rms))
        if best is None or rms[j] < best[0]:
            best = (float(rms[j]), float(hd[j]), float(factor))
    return best


# --------------------------------------------------------------------------- #
# The scheduler
# --------------------------------------------------------------------------- #

def solve_fix(dem, obs: SkylineObservation, prior: SearchPrior,
              fine: RenderGrid = None, coarse: RenderGrid = None,
              top_k: int = 50, separation_km: float = 0.5, cells=None,
              progress=None):
    ''' Coarse pass over every candidate, fine pass over the top_k survivors.

        Returns a dict:
          fix                 (lat, lon) of the fine rank-1
          rms_arcmin          its residual
          heading_deg, focal_factor
          separation          rms of the best candidate > separation_km away,
                              divided by the winner's.  THE quality number: a
                              genuine fix reads well above 1 (0.5x ultrawide:
                              1.68x); every failed search this project ran sat
                              at 1.00-1.02x while its rms looked plausible.
          coarse_rank_of_winner   audit of the pruning -- near top_k means the
                              coarse margin was thin and top_k should grow
          results             fine list, sorted, as (rms_arcmin, lat, lon,
                              heading, focal_factor)
          n_cells, seconds_coarse, seconds_fine

        `cells` overrides candidate generation (pass the previous run's grid to
        re-solve after changing only the observation).  `top_k >= len(cells)`
        degenerates to a flat fine search, which is the correctness reference
        the tests compare against.
    '''
    fine = fine or RenderGrid()
    coarse = coarse or fine.coarsened()
    if cells is None:
        cells = coastal_candidates(dem, prior)
    if not cells:
        raise ValueError("no candidate cells: box, cell size and ground range "
                         "admit nothing -- the prior and the terrain disagree")
    window = _render_window(obs, prior)

    t0 = time.perf_counter()
    ranked = []
    for i, (la, lo, _h) in enumerate(cells):
        s, _hd, _f = score_candidate(dem, la, lo, obs, prior, coarse, window)
        ranked.append((s, i))
        if progress and i % 200 == 0:
            progress(f"coarse {i}/{len(cells)} best "
                     f"{min(v[0] for v in ranked) * 60:.1f}'")
    ranked.sort()
    t_coarse = time.perf_counter() - t0

    survivors = [i for _s, i in ranked[:max(1, int(top_k))]]
    t0 = time.perf_counter()
    out = []
    for n, i in enumerate(survivors):
        la, lo, _h = cells[i]
        s, hd, f = score_candidate(dem, la, lo, obs, prior, fine, window)
        out.append((s * 60.0, la, lo, hd, f, i))
        if progress and n % 20 == 0:
            progress(f"fine {n}/{len(survivors)} best "
                     f"{min(v[0] for v in out):.2f}'")
    out.sort()
    t_fine = time.perf_counter() - t0

    s0, la0, lo0, hd0, f0, i0 = out[0]
    far = [v for v in out
           if distance_km(la0, lo0, v[1], v[2]) > separation_km]
    coarse_pos = {i: r for r, (_s, i) in enumerate(ranked)}
    return dict(
        fix=(la0, lo0), rms_arcmin=s0, heading_deg=hd0, focal_factor=f0,
        separation=(far[0][0] / s0) if far and s0 > 0 else float("inf"),
        coarse_rank_of_winner=coarse_pos[i0] + 1,
        results=[v[:5] for v in out], n_cells=len(cells),
        seconds_coarse=t_coarse, seconds_fine=t_fine)
