''' Terrain resection: fixing position by matching a photographed skyline to a DEM.

    This is the "landfall" fix in its practical form.  Where `landfall.py` models
    the classical observables analytically (dipping range, vertical-angle range,
    horizontal-angle circle fix), this module implements the modern one: render
    the horizon that a digital elevation model predicts from a candidate
    position, compare it with the skyline extracted from a photograph or a
    panned video, and rank the candidates.

    THE LESSON THIS MODULE ENCODES.  The scoring function MUST use both the
    AZIMUTHS of the skyline features and their ELEVATION angles.  Validated on a
    real case (a viewpoint above Bitez, Bodrum, looking across the strait to Kos,
    with the true position known independently):

        scorer                                   blind rank-1 error
        dense skyline, single still frame              1985 m
        summit azimuths only, free scale               2059 m
        summit azimuths only, constrained scale        2059 m
        summit azimuths + ELEVATION                     297 m     <-- 7x better

    With an azimuth-only score the solution drifts to higher, inland terrain --
    nothing in the cost function knows how high the horizon should stand, so a
    hilltop several km away fits as well as the true coastal viewpoint.  Adding
    the elevation residual breaks that degeneracy.  `score_match` therefore has
    no option to disable it, and `test_imu_fusion.TestTerrainResection` asserts
    that an azimuth-only variant is measurably worse.

    HONEST LIMITS.  297 m is one sample at one site with strong landmarks, and
    the winning margin was thin (231.7 vs 230.7).  Published skyline-matching
    systems reach ~50 m using cylindrical registration and a rendered skyline
    database; this implementation ray-marches per candidate and leaves the
    angular scale (focal length) free, which remains the weakest link.  Treat the
    output as a few-hundred-metre fix, and prefer near-field landmarks
    (`landfall.two_landmark_circle_fix`) when they are available -- those give
    tens of metres.

    DEM DATA.  Tiles are 1-arc-second SRTM `.hgt` (3601x3601, big-endian int16),
    not committed to the repository.  Fetch the ones you need, e.g.

        curl -o N37E027.hgt.gz \\
          https://s3.amazonaws.com/elevation-tiles-prod/skadi/N37/N37E027.hgt.gz
        gunzip N37E027.hgt.gz

    and point `DemTiles` at the directory (or set `IMU_FUSION_DEM_DIR`).

    (c) 2026.  MIT License (see LICENSE file).
'''

import os
from math import radians, degrees, cos, sqrt

import numpy as np

from .landfall import K_REFRACTION, effective_radius_km

_DEG_PER_KM_LAT = 1.0 / 111.32


# --------------------------------------------------------------------------- #
# DEM access
# --------------------------------------------------------------------------- #

class DemTiles:
    ''' Lazy reader for a directory of 1-arc-second SRTM `.hgt` tiles.

        `elevation(lat, lon)` takes arrays and returns metres, bilinearly
        interpolated.  Missing tiles and voids read as 0 m (sea level), so a view
        that runs off the available data degrades gracefully instead of raising.
    '''

    def __init__(self, directory: str = None):
        self.directory = directory or os.environ.get(
            "IMU_FUSION_DEM_DIR",
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "dem"))
        self._cache: dict = {}

    def available(self) -> bool:
        ''' True if the directory holds at least one .hgt tile. '''
        return (os.path.isdir(self.directory)
                and any(f.endswith(".hgt") for f in os.listdir(self.directory)))

    def _tile(self, lat_i: int, lon_i: int):
        key = (lat_i, lon_i)
        if key in self._cache:
            return self._cache[key]
        name = (f"{'N' if lat_i >= 0 else 'S'}{abs(lat_i):02d}"
                f"{'E' if lon_i >= 0 else 'W'}{abs(lon_i):03d}.hgt")
        path = os.path.join(self.directory, name)
        tile = None
        if os.path.exists(path):
            side = int(sqrt(os.path.getsize(path) / 2))
            arr = np.fromfile(path, ">i2").reshape(side, side).astype(np.float32)
            arr[arr < -1000] = 0.0                       # voids -> sea level
            tile = (arr, side)
        self._cache[key] = tile
        return tile

    def elevation(self, lat, lon):
        ''' Bilinear terrain elevation (m) for array-like lat/lon. '''
        lat = np.asarray(lat, float)
        lon = np.asarray(lon, float)
        out = np.zeros(lat.shape, np.float32)
        li = np.floor(lat).astype(int)
        oi = np.floor(lon).astype(int)
        for a_i, o_i in set(zip(li.ravel().tolist(), oi.ravel().tolist())):
            tile = self._tile(a_i, o_i)
            if tile is None:
                continue
            arr, side = tile
            m = (li == a_i) & (oi == o_i)
            fy = np.clip((a_i + 1 - lat[m]) * (side - 1), 0, side - 1.001)
            fx = np.clip((lon[m] - o_i) * (side - 1), 0, side - 1.001)
            y0 = fy.astype(int); x0 = fx.astype(int)
            dy = fy - y0; dx = fx - x0
            out[m] = (arr[y0, x0] * (1 - dx) * (1 - dy)
                      + arr[y0, x0 + 1] * dx * (1 - dy)
                      + arr[y0 + 1, x0] * (1 - dx) * dy
                      + arr[y0 + 1, x0 + 1] * dx * dy)
        return out


class SyntheticDem:
    ''' Analytic stand-in for `DemTiles`: a sum of Gaussian hills.  Used by the
        tests so the resection logic is exercised without downloading tiles. '''

    def __init__(self, hills):
        # hills: iterable of (lat, lon, height_m, sigma_deg)
        self.hills = list(hills)

    def available(self) -> bool:
        return True

    def elevation(self, lat, lon):
        lat = np.asarray(lat, float); lon = np.asarray(lon, float)
        out = np.zeros(np.broadcast(lat, lon).shape, float)
        for hlat, hlon, h, s in self.hills:
            out = out + h * np.exp(-(((lat - hlat) ** 2 + (lon - hlon) ** 2)
                                     / (2.0 * s * s)))
        return out


# --------------------------------------------------------------------------- #
# Skyline rendering
# --------------------------------------------------------------------------- #

def render_skyline(dem, lat: float, lon: float, cam_height_m: float,
                   az_start: float = 0.0, az_end: float = 360.0,
                   az_step: float = 0.05, d_min_km: float = 0.10,
                   d_max_km: float = 45.0, d_step_km: float = 0.05,
                   k: float = K_REFRACTION):
    ''' Horizon elevation angle (degrees) versus azimuth, as seen from
        (lat, lon) at `cam_height_m` above sea level.

        Rays are marched outward and the maximum elevation angle along each ray
        is the horizon.  Earth curvature and refraction enter through the
        effective radius, exactly as in `landfall.vertical_angle_deg`:

            alpha = (H - h)/d - d / (2 R_eff)
    '''
    azs = np.arange(az_start, az_end, az_step)
    ds = np.arange(d_min_km, d_max_km, d_step_km)
    dlon_per_km = 1.0 / (111.32 * cos(radians(lat)))
    A, D = np.meshgrid(np.radians(azs), ds, indexing="ij")
    H = dem.elevation(lat + D * np.cos(A) * _DEG_PER_KM_LAT,
                      lon + D * np.sin(A) * dlon_per_km)
    alpha = np.degrees((H - cam_height_m) / 1000.0 / D
                       - D / (2.0 * effective_radius_km(k)))
    return azs, alpha.max(axis=1)


def skyline_peaks(azs, profile, window: int = 20, min_prominence: float = 0.10):
    ''' Prominent summits of a rendered skyline: array of (azimuth, elevation). '''
    n = len(profile)
    picked = []
    for i in range(n):
        lo = max(0, i - window); hi = min(n, i + window + 1)
        w = profile[lo:hi]
        if profile[i] >= w.max() and (profile[i] - w.min()) > min_prominence:
            picked.append([azs[i], profile[i]])
    merged = []
    for a, p in picked:
        if merged and abs(a - merged[-1][0]) < 0.8:
            if p > merged[-1][1]:
                merged[-1] = [a, p]
        else:
            merged.append([a, p])
    return np.array(merged) if merged else np.zeros((0, 2))


# --------------------------------------------------------------------------- #
# Observation + matching
# --------------------------------------------------------------------------- #

class SkylineObservation:
    ''' Summits measured in an image or a stitched panorama.

        x        : horizontal pixel position of each summit
        row      : vertical pixel position of the same summit (elevation!)
        weight   : per-summit confidence (e.g. topographic prominence in pixels)

        `row` is REQUIRED.  The elevation information it carries is what breaks
        the inland/high-ground degeneracy -- see the module docstring.
    '''

    def __init__(self, x, row, weight=None):
        self.x = np.asarray(x, float)
        self.row = np.asarray(row, float)
        if self.x.shape != self.row.shape:
            raise ValueError("x and row must have the same length")
        self.weight = (np.ones_like(self.x) if weight is None
                       else np.asarray(weight, float))
        self.xc = self.x - self.x.mean()

    def __len__(self):
        return len(self.x)


def score_match(obs: SkylineObservation, azs, profile, peaks,
                f_px_per_deg: float, az0: float, sign: int,
                tol_deg: float = 0.7, elev_scale_px: float = 20.0,
                min_inliers: int = 4):
    ''' Score one (scale, heading, direction) hypothesis.

        Returns None if too few summits align, else a dict with the azimuth
        inlier count/residual, the ELEVATION residual, and the combined score.

        The score deliberately combines both terms:

            score = sum(weight of azimuth inliers) / (1 + elev_resid / scale)

        so a hypothesis that lines the summits up in bearing but puts the horizon
        at the wrong height is penalised.  There is no azimuth-only mode.
    '''
    if len(peaks) == 0 or len(obs) == 0:
        return None
    model_az = peaks[:, 0]
    a = (az0 + sign * obs.xc / f_px_per_deg) % 360.0
    diff = np.abs((a[:, None] - model_az[None, :] + 180.0) % 360.0 - 180.0)
    nearest = diff.min(axis=1)
    inl = nearest < tol_deg
    if inl.sum() < min_inliers:
        return None
    el = np.interp(a[inl], azs, profile, period=360)
    predicted = -f_px_per_deg * el
    cy = np.mean(obs.row[inl] - predicted)
    elev_resid = float(np.sqrt(np.mean((obs.row[inl] - (cy + predicted)) ** 2)))
    az_resid = float(np.sqrt(np.mean(nearest[inl] ** 2)))
    score = float(obs.weight[inl].sum() / (1.0 + elev_resid / elev_scale_px))
    return dict(score=score, n_inliers=int(inl.sum()), az_resid_deg=az_resid,
                elev_resid_px=elev_resid, az0_deg=az0,
                f_px_per_deg=f_px_per_deg, sign=int(sign), cy=float(cy))


def best_match(obs: SkylineObservation, azs, profile, peaks,
               f_list, az_step: float = 0.5, tol_deg: float = 0.7,
               max_sweep_deg: float = 340.0, min_inliers: int = 4):
    ''' Best hypothesis over scale, heading and pan direction. '''
    best = None
    for f in f_list:
        span = (obs.xc.max() - obs.xc.min()) / f
        if span > max_sweep_deg:
            continue
        for sign in (1, -1):
            for az0 in np.arange(0.0, 360.0, az_step):
                m = score_match(obs, azs, profile, peaks, f, az0, sign,
                                tol_deg=tol_deg, min_inliers=min_inliers)
                if m is None:
                    continue
                if best is None or (m["score"], -m["elev_resid_px"]) > \
                        (best["score"], -best["elev_resid_px"]):
                    best = m
    return best


def resect(dem, obs: SkylineObservation, candidates, cam_above_ground_m: float = 2.0,
           f_list=None, az_step: float = 0.5, render_kw=None, peak_kw=None,
           min_inliers: int = 4):
    ''' Rank candidate positions by how well the DEM horizon explains `obs`.

        candidates : iterable of (lat, lon)
        Returns a list of dicts sorted best-first, each carrying the match
        diagnostics plus the candidate position and its terrain elevation.
    '''
    f_list = np.arange(30.0, 140.0, 2.0) if f_list is None else f_list
    render_kw = render_kw or {}
    peak_kw = peak_kw or {}
    out = []
    for lat, lon in candidates:
        ground = float(dem.elevation(np.array([lat]), np.array([lon]))[0])
        azs, prof = render_skyline(dem, lat, lon, ground + cam_above_ground_m,
                                   **render_kw)
        peaks = skyline_peaks(azs, prof, **peak_kw)
        m = best_match(obs, azs, prof, peaks, f_list, az_step=az_step,
                       min_inliers=min_inliers)
        if m is None:
            continue
        m.update(lat=lat, lon=lon, ground_elev_m=ground)
        out.append(m)
    out.sort(key=lambda r: -r["score"])
    return out
