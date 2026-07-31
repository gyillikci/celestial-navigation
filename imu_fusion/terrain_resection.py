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


def tile_name(lat: float, lon: float) -> str:
    ''' SRTM tile filename covering a coordinate, e.g. 37.03,27.36 -> N37E027. '''
    la, lo = int(np.floor(lat)), int(np.floor(lon))
    return (f"{'N' if la >= 0 else 'S'}{abs(la):02d}"
            f"{'E' if lo >= 0 else 'W'}{abs(lo):03d}")


def fetch_tiles(lat_min: float, lat_max: float, lon_min: float, lon_max: float,
                directory: str = None, base_url: str = None, quiet: bool = False):
    ''' Download the 1-arc-second SRTM tiles covering a bounding box.

        Tiles come from the public AWS "elevation-tiles-prod" skadi mirror and
        are ~25 MB each (gzipped ~6 MB), so they are NOT committed -- fetch the
        ones your area needs:

            from imu_fusion.terrain_resection import fetch_tiles
            fetch_tiles(36.9, 37.2, 27.1, 27.6)      # the Bodrum peninsula

        Returns the list of local .hgt paths that now exist.  Already-present
        tiles are skipped, and a tile that cannot be fetched is reported and
        skipped rather than raising (the reader treats absent data as sea level).
    '''
    import gzip
    import urllib.request

    directory = directory or DemTiles().directory
    base_url = base_url or "https://s3.amazonaws.com/elevation-tiles-prod/skadi"
    os.makedirs(directory, exist_ok=True)
    got = []
    for la in range(int(np.floor(lat_min)), int(np.floor(lat_max)) + 1):
        for lo in range(int(np.floor(lon_min)), int(np.floor(lon_max)) + 1):
            name = tile_name(la + 0.5, lo + 0.5)
            path = os.path.join(directory, name + ".hgt")
            if os.path.exists(path):
                got.append(path)
                continue
            url = f"{base_url}/{name[:3]}/{name}.hgt.gz"
            try:
                with urllib.request.urlopen(url, timeout=180) as r:
                    data = gzip.decompress(r.read())
                with open(path, "wb") as f:
                    f.write(data)
                got.append(path)
                if not quiet:
                    print(f"fetched {name} ({len(data) // (1 << 20)} MB)")
            except Exception as exc:                      # pragma: no cover
                if not quiet:
                    print(f"could not fetch {name}: {exc}")
    return got


class SyntheticDem:
    ''' Analytic stand-in for `DemTiles`: a sum of Gaussian hills.  Used by the
        tests so the resection logic is exercised without downloading tiles.

        NOTE: Gaussian hills are much smoother than real rocky terrain and give
        an OPTIMISTIC picture — on real SRTM data over the Bodrum peninsula the
        same experiment gives a ~2.5x larger median error (362 m vs 141 m) and a
        ~20% outright failure rate.  Use `DemTiles` + `fetch_tiles` for any
        accuracy claim. '''

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


def _dip_deg(eye_above_water_m, k=K_REFRACTION):
    ''' Depression of the sea horizon for an eye this far above the surface. '''
    h = max(float(eye_above_water_m), 0.0)
    return degrees(sqrt(2.0 * h / 1000.0 / effective_radius_km(k)))


def render_skyline(dem, lat: float, lon: float, cam_height_m: float,
                   az_start: float = 0.0, az_end: float = 360.0,
                   az_step: float = 0.05, d_min_km: float = 0.10,
                   d_max_km: float = 45.0, d_step_km: float = 0.05,
                   k: float = K_REFRACTION, water_level_m: float = None,
                   clamp_water_horizon: bool = True,
                   visibility_km: float = None,
                   scale_height_m: float = None,
                   contrast_threshold: float = None):
    ''' Horizon elevation angle (degrees) versus azimuth, as seen from
        (lat, lon) at `cam_height_m` above sea level.

        Rays are marched outward and the maximum elevation angle along each ray
        is the horizon.  Earth curvature and refraction enter through the
        effective radius, exactly as in `landfall.vertical_angle_deg`:

            alpha = (H - h)/d - d / (2 R_eff)

        `water_level_m` clamps terrain heights UP to that level, and on any
        coastal or lake scene you almost certainly want it.  The tiles this
        project fetches (AWS `elevation-tiles-prod`, the Mapzen terrain product)
        are SRTM on land MERGED WITH BATHYMETRY: measured on N40E028, 4.1 million
        cells lie below -100 m and the minimum is -1308 m.  What a camera sees
        over water is the water SURFACE, not the sea floor, so pass 0.0 at the
        coast and the lake level inland (1897.0 for Lake Tahoe).

        Leaving it None keeps raw heights, which is right only where genuine
        below-datum LAND is in view (Dead Sea, Death Valley) and wrong wherever
        there is water.  It stays the default so existing callers are unchanged.

        A second, subtler trap has no parameter because it needs judgement:
        SRTM's ~30 m posting SMEARS COASTLINES, so a camera standing at the
        water's edge sees a few metres of spurious "land" a few hundred metres
        out along seaward azimuths, which this function faithfully reports as a
        horizon floor of half a degree or more at EVERY bearing.  Measured at
        Istanbul: a phantom +0.65 deg from 6 m of DEM shoreline at 300 m.  Raise
        `d_min_km` past the smear (1.2 km sufficed there) whenever the camera is
        on a shore.

        `clamp_water_horizon` guards a THIRD trap, and the interesting thing about
        it is how narrow it turns out to be.  Along an azimuth with no land in
        range the profile is the maximum of -h/d - d/(2 R_eff), and that maximum
        is attained at d = sqrt(2 h R_eff) -- the blind range -- where its value
        is exactly -sqrt(2h/R_eff), the true dip.  So the march finds the sea
        horizon by itself PROVIDED `d_max_km` reaches the blind range: 8.6 km for
        a 5 m eye, 24 km for 40 m.  Short-range renders truncate before that and
        report the water too low (5 arcmin too low at d_max = 2 km, 0.6 at 5 km).
        Nothing can appear below the water horizon, so the profile is clamped.

        This was first written to fix a supposed 10.6 arcmin bias on an ultrawide
        frame whose render used d_max = 40 km.  That diagnosis was WRONG -- the
        march was already returning -4.016 arcmin, the correct dip -- and a unit
        test caught it.  The clamp stays because short-range renders are real, but
        it is a guard, not a repair.

        The clamp only acts when `water_level_m` is set and the camera is above
        it, so land-only callers are unaffected; pass False to disable.

        `visibility_km` is the FOURTH trap and the one that cost the most.  This
        function answers "what is GEOMETRICALLY visible" -- it marches until the
        earth curves away.  A photograph answers "what is ATMOSPHERICALLY
        visible".  Measured on the Istanbul ultrawide: the render placed the
        Samanli mountains on the far shore of the Marmara, 680 m at 39.5 km, in
        the skyline at +50 arcmin, correctly, because such a peak clears a 5 m
        eye's horizon out to 100 km.  The photograph shows only sea there, so the
        extractor traced the water horizon and every one of those samples entered
        the fit as a 50 arcmin error.  The residual was 32.8 arcmin and NO
        position could lower it, because the fault was in the forward model.

        Passing a visual range applies Koschmieder extinction along the actual
        SLANT path (see `visibility.py`), so terrain whose apparent contrast falls
        below `contrast_threshold` is excluded from the horizon.  A flat `d_max`
        cap cannot do this: aerosol sits near the ground, so a summit and a
        sea-level islet at the same range are not equally visible.  Leaving it
        None keeps the pure-geometry behaviour, which is right for a synthetic
        study and wrong for a photograph.
    '''
    azs = np.arange(az_start, az_end, az_step)
    ds = np.arange(d_min_km, d_max_km, d_step_km)
    dlon_per_km = 1.0 / (111.32 * cos(radians(lat)))
    A, D = np.meshgrid(np.radians(azs), ds, indexing="ij")
    H = dem.elevation(lat + D * np.cos(A) * _DEG_PER_KM_LAT,
                      lon + D * np.sin(A) * dlon_per_km)
    if water_level_m is not None:
        H = np.maximum(H, float(water_level_m))
    alpha = np.degrees((H - cam_height_m) / 1000.0 / D
                       - D / (2.0 * effective_radius_km(k)))
    if visibility_km is not None:
        from .visibility import (is_detectable, AEROSOL_SCALE_HEIGHT_M,
                                 CONTRAST_THRESHOLD)
        seen = is_detectable(
            D, cam_height_m, H, visibility_km,
            scale_height_m=(AEROSOL_SCALE_HEIGHT_M if scale_height_m is None
                            else scale_height_m),
            threshold=(CONTRAST_THRESHOLD if contrast_threshold is None
                       else contrast_threshold))
        # -inf, not 0: an invisible sample must never win the max, and terrain
        # below the observer legitimately has negative elevation angles.
        alpha = np.where(seen, alpha, -np.inf)
    profile = alpha.max(axis=1)
    if visibility_km is not None:
        # an azimuth with NOTHING detectable still has a sky: the water horizon
        # if we know the water level, otherwise the geometric horizon.
        profile = np.where(np.isfinite(profile), profile, -_dip_deg(
            cam_height_m - (0.0 if water_level_m is None else float(water_level_m)), k))
    if water_level_m is not None and clamp_water_horizon:
        eye = cam_height_m - float(water_level_m)
        if eye > 0.0:
            dip = degrees(sqrt(2.0 * eye / 1000.0 / effective_radius_km(k)))
            profile = np.maximum(profile, -dip)
    return azs, profile


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


# --------------------------------------------------------------------------- #
# Priors: the realistic on-device entry point
#
# Blind resection is not what a phone does.  It has a dead-reckoned position
# (metres to kilometres), a magnetometer heading (degrees), and it knows its own
# lens.  Each prunes a dimension of the search, and the pruning is what makes the
# problem real time.  Measured on one core:
#
#     blind             55 scales x 720 headings x 2 = 79200 hyp    1294 ms
#     + magnetometer +/-20 deg                       =  8800 hyp     150 ms
#     + known focal length                           =   160 hyp     3.1 ms
#     + magnetometer +/-2 deg                        =    32 hyp     0.5 ms
#
# The magnetometer is a poor MEASUREMENT (see `terrain_factors`: ~360 m from
# bearings at 1.5 deg, versus ~18 m from compass-free angles) but an excellent
# SEARCH PRIOR -- worth about 2600x here, and it is what makes landmark
# IDENTIFICATION unambiguous in the first place.
# --------------------------------------------------------------------------- #

def candidate_grid(lat: float, lon: float, radius_km: float, step_m: float):
    ''' Square grid of candidate positions centred on a rough position. '''
    n = max(1, int(radius_km * 1000.0 / step_m))
    dlat = step_m * _DEG_PER_KM_LAT / 1000.0
    dlon = step_m / (111320.0 * cos(radians(lat)))
    return [(lat + i * dlat, lon + j * dlon)
            for i in range(-n, n + 1) for j in range(-n, n + 1)]


def resect_with_priors(dem, obs: SkylineObservation, rough_lat: float,
                       rough_lon: float, prior_radius_km: float = 1.0,
                       grid_step_m: float = 200.0,
                       mag_heading_deg: float = None,
                       mag_sigma_deg: float = 2.0, mag_n_sigma: float = 3.0,
                       f_px_per_deg: float = None, f_list=None,
                       az_step: float = 0.25, cam_above_ground_m: float = 2.0,
                       min_ground_elev_m: float = None,
                       render_pad_deg: float = 15.0, render_kw=None,
                       peak_kw=None, min_inliers: int = 4):
    ''' Resection from a rough position plus, optionally, a magnetometer heading
        and a known angular scale.

        mag_heading_deg is the bearing of the observation's MEAN x column (frame
        or panorama centre), matching score_match's az0.  Only headings within
        mag_n_sigma * mag_sigma_deg are tried, and only the arc the camera could
        have seen is rendered.

        Returns the same ranked list of dicts as `resect`.

        ACCURACY IS GRID-LIMITED, NOT NOISE-LIMITED.  On synthetic data the
        rank-1 error tracks `grid_step_m` and is almost flat in measurement
        noise (mean error 141 m at 0 px and 141 m at 6 px of pointing noise;
        254 / 149 / 101 m for 200 / 100 / 50 m grids).  This search is therefore
        an IDENTIFIER, not the final estimator: use it to decide which DEM
        summits were photographed, then hand those to `terrain_factors` for a
        continuous least-squares fix with a covariance.

        ON REAL TERRAIN (SRTM, Bodrum peninsula, 24 viewpoints, exact synthetic
        observations, 200 m grid) it is markedly harder than smooth synthetic
        hills: 5 of 24 produced no match at all, and of the 19 that matched the
        median error was 362 m.

        THE ELEVATION RESIDUAL IS A SELF-ASSESSMENT.  `elev_resid_px` correlates
        with the true error (r = +0.58) and makes a usable accept/reject gate:

            elev_resid <  3 px   n= 6   median  80 m   max 279 m
            elev_resid >= 3 px   n=13   median 519 m

        So gate on it -- report a fix when the residual is small and decline
        when it is not, rather than always returning rank 1.
    '''
    render_kw = dict(render_kw or {})
    peak_kw = peak_kw or {}
    if f_px_per_deg is not None:
        f_list = np.array([float(f_px_per_deg)])
    elif f_list is None:
        f_list = np.arange(30.0, 140.0, 2.0)

    if mag_heading_deg is None:
        az_candidates = np.arange(0.0, 360.0, az_step)
    else:
        half = mag_n_sigma * mag_sigma_deg
        az_candidates = np.arange(mag_heading_deg - half,
                                  mag_heading_deg + half + 1e-9, az_step)
        sweep = (obs.xc.max() - obs.xc.min()) / float(np.min(f_list))
        az_lo = mag_heading_deg - half - sweep / 2.0 - render_pad_deg
        az_hi = mag_heading_deg + half + sweep / 2.0 + render_pad_deg
        if az_hi - az_lo < 360.0:
            render_kw.setdefault("az_start", az_lo)
            render_kw.setdefault("az_end", az_hi)

    out = []
    for lat, lon in candidate_grid(rough_lat, rough_lon, prior_radius_km,
                                   grid_step_m):
        ground = float(dem.elevation(np.array([lat]), np.array([lon]))[0])
        if min_ground_elev_m is not None and ground < min_ground_elev_m:
            continue
        azs, prof = render_skyline(dem, lat, lon, ground + cam_above_ground_m,
                                   **render_kw)
        peaks = skyline_peaks(azs, prof, **peak_kw)
        best = None
        for f in f_list:
            for sign in (1, -1):
                for az0 in az_candidates:
                    m = score_match(obs, azs, prof, peaks, f, az0 % 360.0, sign,
                                    min_inliers=min_inliers)
                    if m is None:
                        continue
                    if best is None or (m["score"], -m["elev_resid_px"]) > \
                            (best["score"], -best["elev_resid_px"]):
                        best = m
        if best is None:
            continue
        best.update(lat=lat, lon=lon, ground_elev_m=ground)
        out.append(best)
    out.sort(key=lambda r: -r["score"])
    return out


# --------------------------------------------------------------------------- #
# Synthetic observations (for tests and for error-budget experiments)
# --------------------------------------------------------------------------- #

def synth_skyline_observation(dem, lat: float, lon: float,
                              az_center_deg: float, fov_deg: float,
                              width_px: int = 2000, cam_above_ground_m: float = 2.0,
                              noise_px: float = 0.0, rng=None,
                              cy_px: float = 500.0, render_kw=None,
                              peak_kw=None):
    ''' Photograph the DEM: render the horizon from a known position and project
        its summits into image coordinates, exactly as `score_match` inverts.

            x   = (azimuth - az_center) * f
            row = cy - f * elevation_deg,     f = width_px / fov_deg

        The mapping is linear in angle (panorama / cylindrical convention), which
        is what the matcher assumes.  Returns (observation, truth) where truth
        carries the generating parameters and the summits used.
    '''
    import random as _random
    rng = rng or _random.Random(0)
    f = width_px / float(fov_deg)
    render_kw = dict(render_kw or {})
    render_kw.setdefault("az_start", az_center_deg - fov_deg)
    render_kw.setdefault("az_end", az_center_deg + fov_deg)
    ground = float(dem.elevation(np.array([lat]), np.array([lon]))[0])
    azs, prof = render_skyline(dem, lat, lon, ground + cam_above_ground_m,
                               **render_kw)
    peaks = skyline_peaks(azs, prof, **(peak_kw or {}))
    keep = [p for p in peaks
            if abs(((p[0] - az_center_deg + 180.0) % 360.0) - 180.0) <= fov_deg / 2.0]
    if not keep:
        return None, dict(reason="no summits in field of view")
    keep = np.array(keep)
    x = ((keep[:, 0] - az_center_deg + 180.0) % 360.0 - 180.0) * f
    row = cy_px - f * keep[:, 1]
    if noise_px:
        x = x + np.array([rng.gauss(0.0, noise_px) for _ in x])
        row = row + np.array([rng.gauss(0.0, noise_px) for _ in row])
    obs = SkylineObservation(x, row, weight=np.ones(len(x)) * 10.0)
    truth = dict(lat=lat, lon=lon, ground_elev_m=ground, az_center_deg=az_center_deg,
                 fov_deg=fov_deg, f_px_per_deg=f, cy_px=cy_px,
                 summits=keep, n_summits=len(keep))
    return obs, truth
