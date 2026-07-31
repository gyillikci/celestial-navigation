''' Synthetic scenes over real terrain: ground truth the field can never give.

    The real-photo results (182-359 m coastal, no-fix through a windscreen)
    validate the pipeline but cannot ISOLATE anything: every frame carries
    extraction noise, compass error, canopy bias and haze at once, in unknown
    amounts.  This module generates the observation a camera WOULD record from
    an exactly known position over a real DEM, with each error injected alone
    and in known quantity -- so the error budget becomes measurable instead of
    inferred.

    The generator must be exactly consistent with the solver's forward model,
    or the experiment measures the inconsistency instead of the error.  Two
    pieces guarantee that:

      * `pixels_from_angles` is the algebraic INVERSE of
        `resection_geometry.image_ray_angles` (round-trip tested to 1e-9 deg);
      * `synthesize_rows` projects the same `render_skyline` profile the solver
        will compare against, so a zero-noise scene must solve back to its own
        cell exactly -- the pipeline's end-to-end identity, now a unit test.

    CANOPY is the error the real data said matters most (the joint all-lens
    solve pinned a ~200 m floor on correlated systematics, and SRTM's
    half-seen treetops are the leading suspect).  `CanopyDem` grows trees on
    the TERRAIN THE SCENE IS GENERATED FROM while the solver keeps the bare
    DEM -- exactly the mismatch a real photograph has.

    (c) 2026.  MIT License (see LICENSE file).
'''

import math
from dataclasses import dataclass, field

import numpy as np

from .fix_pipeline import SkylineObservation
from .landfall import K_REFRACTION
from .resection_geometry import image_ray_angles
from .terrain_resection import render_skyline


def focal_px(f35_mm, width_px, height_px):
    ''' Pixel focal length for a 35 mm-equivalent focal and a sensor raster. '''
    half_diag_deg = math.degrees(math.atan(math.hypot(36.0, 24.0) / 2.0 / f35_mm))
    return (math.hypot(width_px, height_px) / 2.0) / math.tan(
        math.radians(half_diag_deg))


def pixels_from_angles(az_off_deg, el_deg, f_px, pitch_deg=0.0, roll_deg=0.0,
                       cx=0.0, cy=0.0):
    ''' (azimuth offset, elevation) -> pixel coordinates.

        Exact inverse of `image_ray_angles`: build the world ray, rotate INTO
        the pitched camera, project, then apply roll.  Returns (x_px, y_px)
        and a validity mask (rays behind the camera project nowhere).
    '''
    az = np.radians(np.asarray(az_off_deg, dtype=float))
    el = np.radians(np.asarray(el_deg, dtype=float))
    wx = np.cos(el) * np.sin(az)
    wy = np.sin(el)
    wz = np.cos(el) * np.cos(az)
    cp, sp = math.cos(math.radians(pitch_deg)), math.sin(math.radians(pitch_deg))
    # inverse of  wz = vz*cp - vy*sp ; wy = vz*sp + vy*cp
    vz = wz * cp + wy * sp
    vy = -wz * sp + wy * cp
    vx = wx
    ok = vz > 1e-9
    vz_safe = np.where(ok, vz, 1.0)
    xp = float(f_px) * vx / vz_safe
    yp = -float(f_px) * vy / vz_safe
    cr, sr = math.cos(math.radians(roll_deg)), math.sin(math.radians(roll_deg))
    # inverse of  xp = x*cr + y*sr ; yp = -x*sr + y*cr
    x = xp * cr - yp * sr
    y = xp * sr + yp * cr
    return x + float(cx), y + float(cy), ok


class CanopyDem:
    ''' A DEM wearing trees: heights raised by `canopy_m` wherever the bare
        ground exceeds `above_m`.  Generate the scene from this, solve against
        the bare DEM, and the mismatch is exactly a forest SRTM half-sees.
    '''

    def __init__(self, dem, canopy_m=10.0, above_m=2.0):
        self.dem, self.canopy_m, self.above_m = dem, float(canopy_m), float(above_m)

    def elevation(self, lat, lon):
        h = self.dem.elevation(lat, lon)
        return np.where(h > self.above_m, h + self.canopy_m, h)


@dataclass
class SceneSpec:
    ''' Everything the synthetic camera is and does, plus the exact truth. '''
    lat: float
    lon: float
    heading_deg: float
    f35_mm: float
    width_px: int = 4032
    height_px: int = 3024
    pitch_deg: float = 2.0
    roll_deg: float = 0.5
    eye_above_water_m: float = None    # coastal: absolute eye
    eye_agl_m: float = None            # inland: eye above local ground
    water_level_m: float = 0.0
    d_min_km: float = 0.5
    d_max_km: float = 30.0
    visibility_km: float = None
    k: float = K_REFRACTION


def synthesize_rows(dem, spec: SceneSpec, az_step=0.02, d_step_km=0.05):
    ''' The skyline rows this camera records from the truth position.

        Renders the profile (from `dem`, which may be a CanopyDem), projects it
        through the exact camera model, and interpolates the boundary row at
        every integer column.  NaN where the skyline leaves the frame.
    '''
    if spec.eye_agl_m is not None:
        ground = float(dem.elevation(np.array([spec.lat]), np.array([spec.lon]))[0])
        cam = ground + spec.eye_agl_m
        water = None
    else:
        cam = spec.eye_above_water_m
        water = spec.water_level_m
    half = math.degrees(math.atan(
        spec.width_px / 2.0 / focal_px(spec.f35_mm, spec.width_px, spec.height_px)))
    azs, prof = render_skyline(
        dem, spec.lat, spec.lon, cam,
        az_start=spec.heading_deg - half - 3.0,
        az_end=spec.heading_deg + half + 3.0,
        az_step=az_step, d_min_km=spec.d_min_km, d_max_km=spec.d_max_km,
        d_step_km=d_step_km, k=spec.k, water_level_m=water,
        visibility_km=spec.visibility_km)
    f = focal_px(spec.f35_mm, spec.width_px, spec.height_px)
    cx, cy = (spec.width_px - 1) / 2.0, (spec.height_px - 1) / 2.0
    px, py, ok = pixels_from_angles(azs - spec.heading_deg, prof, f,
                                    pitch_deg=spec.pitch_deg,
                                    roll_deg=spec.roll_deg, cx=cx, cy=cy)
    cols = np.arange(spec.width_px, dtype=float)
    order = np.argsort(px[ok])
    rows = np.interp(cols, px[ok][order], py[ok][order],
                     left=np.nan, right=np.nan)
    return rows


def corrupt_rows(rows, rng, sigma_px=0.0, corr_len_px=150.0):
    ''' Correlated extraction noise: white noise smoothed to `corr_len_px`,
        rescaled to `sigma_px`.  Correlated, because that is what real
        extraction error is -- white noise flatters the solver.
    '''
    if sigma_px <= 0.0:
        return rows.copy()
    n = len(rows)
    w = rng.normal(size=n + int(6 * corr_len_px))
    k = np.exp(-0.5 * (np.arange(-3 * corr_len_px, 3 * corr_len_px + 1)
                       / corr_len_px) ** 2)
    z = np.convolve(w, k / np.sqrt((k * k).sum()), mode='same')[:n]
    z = z / (z.std() or 1.0)
    return rows + sigma_px * z


def observation_from_scene(dem, spec: SceneSpec, rng=None, sigma_px=0.0,
                           corr_len_px=150.0, focal_error=0.0,
                           pitch_known=False):
    ''' Build the SkylineObservation the solver will see.

        `focal_error` corrupts the focal length the SOLVER is told (the truth
        stays in the generator).  Compass error is injected at the prior, not
        here, because that is where a real compass error lives.
    '''
    rows = synthesize_rows(dem, spec)
    good = np.isfinite(rows)
    if rng is not None:
        rows = corrupt_rows(rows, rng, sigma_px, corr_len_px)
    kw = dict(
        columns=np.arange(spec.width_px, dtype=float)[good],
        rows=rows[good],
        f_px=focal_px(spec.f35_mm, spec.width_px, spec.height_px)
             * (1.0 + focal_error),
        cx=(spec.width_px - 1) / 2.0, cy=(spec.height_px - 1) / 2.0,
        roll_deg=spec.roll_deg, k=spec.k)
    if spec.eye_agl_m is not None:
        kw['eye_agl_m'] = spec.eye_agl_m
        kw['eye_above_water_m'] = 999.0        # unused in AGL mode
    else:
        kw['eye_above_water_m'] = spec.eye_above_water_m
    if pitch_known:
        kw['pitch_deg'] = spec.pitch_deg
    return SkylineObservation(**kw)
