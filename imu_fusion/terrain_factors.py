''' Terrain landmarks as GTSAM factors — landfall observables inside the graph.

    `terrain_resection.py` answers "where am I?" by scoring a grid of candidate
    positions against a DEM-rendered horizon.  That is a good bootstrap, but it
    is a search, not an estimator: it has no covariance, it cannot be fused with
    the celestial sights or the IMU, and it cannot be run incrementally.

    This module closes that gap.  Once skyline summits have been IDENTIFIED --
    matched to DEM peaks with known (lat, lon, height) -- each one becomes an
    ordinary measurement of a surveyed point, and the whole apparatus of the
    celestial factor graph applies.  Three factors, in increasing order of how
    much they trust the compass:

      * `landmark_bearing_factor`      -- absolute bearing to one landmark.
        Needs a heading reference; the magnetometer's error enters the sigma
        directly, which is why this is the weakest of the three.
      * `landmark_elevation_factor`    -- vertical angle to a summit of known
        height: a RANGE, hence a distance circle.  Compass-free, but limited by
        refraction (see `landfall.range_bias_from_refraction_km`).
      * `landmark_horizontal_angle_factor` -- the angle SUBTENDED by two
        landmarks.  Compass-free by construction: a heading bias is common to
        both bearings and cancels in the difference.  This is the accurate one.

    MAGNETOMETER.  A phone magnetometer is good to roughly 1-2 deg, which at 10
    km is a 175-350 m cross-range error -- worse than the fix we are trying to
    make.  `bearing_sigma_deg` composes it with the pixel-level pointing error so
    the bearing factor is weighted honestly.  The recommended use of a
    magnetometer heading is therefore NOT as a precision observable but as the
    disambiguator that makes landmark IDENTIFICATION possible; the geometry is
    then carried by the compass-free factors.  `test_imu_fusion` asserts this:
    a heading bias wrecks the bearing-only fix and leaves the horizontal-angle
    fix untouched.

    All factors are unary on a Pose3 whose translation is the local ENU offset
    from (lat0, lon0), matching `celestial_factor_graph`, and they reuse its
    reduced (east, north) Jacobian.

    (c) 2026.  MIT License (see LICENSE file).
'''

from math import sin, cos, atan2, radians, degrees, sqrt

import numpy as np
import gtsam

from .astro import enu_to_latlon, great_circle_km
from .celestial_factor_graph import _reduced_en_jacobian
from .landfall import vertical_angle_deg, K_REFRACTION

_DEG = 1.0


# --------------------------------------------------------------------------- #
# Geometry helpers
# --------------------------------------------------------------------------- #

def bearing_deg(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    ''' Initial great-circle bearing from (lat1,lon1) to (lat2,lon2), degrees
        true, 0 = north. '''
    p1, p2 = radians(lat1), radians(lat2)
    dl = radians(lon2 - lon1)
    y = sin(dl) * cos(p2)
    x = cos(p1) * sin(p2) - sin(p1) * cos(p2) * cos(dl)
    return degrees(atan2(y, x)) % 360.0


def wrap180(x: float) -> float:
    return (x + 180.0) % 360.0 - 180.0


def bearing_sigma_deg(pixel_sigma_deg: float, mag_sigma_deg: float) -> float:
    ''' 1-sigma of an ABSOLUTE bearing: the pointing error within the image
        combined with the heading reference (magnetometer) error.  The
        magnetometer term normally dominates. '''
    return sqrt(pixel_sigma_deg ** 2 + mag_sigma_deg ** 2)


class Landmark:
    ''' A surveyed terrain point identified in the image. '''

    def __init__(self, name: str, lat: float, lon: float, height_m: float = 0.0):
        self.name = name
        self.lat = lat
        self.lon = lon
        self.height_m = height_m

    def __repr__(self):
        return f"Landmark({self.name!r}, {self.lat:.5f}, {self.lon:.5f}, {self.height_m:.0f} m)"


# --------------------------------------------------------------------------- #
# Factors
# --------------------------------------------------------------------------- #

def landmark_bearing_factor(key, landmark: Landmark, meas_bearing_deg: float,
                            sigma_deg: float, lat0: float, lon0: float,
                            bias_key=None):
    ''' Absolute bearing to a surveyed landmark -> a line of position through it.

        With `bias_key`, the factor becomes BINARY on (pose, compass bias) and a
        single bias variable is shared by every bearing in the frame.  That is
        the honest model and it changes the answer, not just the bookkeeping:

          * a compass error is COMMON to all bearings taken at one moment, not
            independent per landmark, so inflating each sigma by the
            magnetometer error both over-counts the noise and lets the fit
            pretend the errors average down;
          * moving across the line of sight swings every bearing by nearly the
            same amount, which is exactly what a compass bias looks like -- so
            with the bias estimated, the graph's own covariance reports the
            resulting line-of-position elongation instead of hiding it.

        Without `bias_key` the factor stays unary and `sigma_deg` must already
        carry the heading error (see `bearing_sigma_deg`).
    '''
    noise = gtsam.noiseModel.Isotropic.Sigma(1, sigma_deg * _DEG)
    keys = [key] if bias_key is None else [key, bias_key]

    def predict(pose):
        t = pose.translation()
        lat, lon = enu_to_latlon(t[0], t[1], lat0, lon0)
        return bearing_deg(lat, lon, landmark.lat, landmark.lon)

    def error(this, values, H):
        pose = values.atPose3(this.keys()[0])
        bias = (float(values.atVector(this.keys()[1])[0])
                if bias_key is not None else 0.0)
        base = predict(pose)
        resid = wrap180(base + bias - meas_bearing_deg)
        if H is not None:
            # wrap-safe reduced difference on (east, north)
            jac = np.zeros((1, 6), order="F")
            eps = 1e-6
            for i in (3, 4):
                d = np.zeros(6); d[i] = eps
                jac[0, i] = wrap180(predict(pose.retract(d)) - base) / eps
            H[0] = jac
            if bias_key is not None:
                H[1] = np.ones((1, 1), order="F")
        return np.array([resid])

    return gtsam.CustomFactor(noise, keys, error)


def landmark_elevation_factor(key, landmark: Landmark, meas_elev_deg: float,
                              sigma_deg: float, lat0: float, lon0: float,
                              cam_height_m: float = 2.0,
                              k: float = K_REFRACTION, bias_key=None):
    ''' Vertical angle to a summit of known height -> a RANGE (distance circle).

        Compass-free.  Uses the same curvature+refraction model as
        `landfall.vertical_angle_deg`.

        With `bias_key`, a single PITCH bias is shared by every elevation in the
        frame.  Measured on real Theodolite frames this bias was +1.50 deg from a
        moving car and +0.87 deg standing still -- large enough that treating
        elevations as unbiased is not an option, and common-mode enough that
        inflating each sigma separately would be the wrong correction. '''
    noise = gtsam.noiseModel.Isotropic.Sigma(1, sigma_deg * _DEG)
    keys = [key] if bias_key is None else [key, bias_key]

    def predict(pose):
        t = pose.translation()
        lat, lon = enu_to_latlon(t[0], t[1], lat0, lon0)
        d = great_circle_km(lat, lon, landmark.lat, landmark.lon)
        d = max(d, 1e-3)
        return vertical_angle_deg(d, landmark.height_m, cam_height_m, k)

    def error(this, values, H):
        pose = values.atPose3(this.keys()[0])
        bias = (float(values.atVector(this.keys()[1])[0])
                if bias_key is not None else 0.0)
        base = predict(pose)
        if H is not None:
            H[0] = _reduced_en_jacobian(pose, predict, base)
            if bias_key is not None:
                H[1] = np.ones((1, 1), order="F")
        return np.array([base + bias - meas_elev_deg])

    return gtsam.CustomFactor(noise, keys, error)


def landmark_horizontal_angle_factor(key, lm_a: Landmark, lm_b: Landmark,
                                     meas_angle_deg: float, sigma_deg: float,
                                     lat0: float, lon0: float):
    ''' The horizontal angle SUBTENDED by two surveyed landmarks -> a circle of
        position (the classical horizontal-sextant-angle fix).

        COMPASS-FREE: any heading bias is common to both bearings and cancels,
        so `sigma_deg` carries only the in-image pointing error -- typically a
        few arc-minutes, versus degrees for the magnetometer. '''
    noise = gtsam.noiseModel.Isotropic.Sigma(1, sigma_deg * _DEG)

    def predict(pose):
        t = pose.translation()
        lat, lon = enu_to_latlon(t[0], t[1], lat0, lon0)
        ba = bearing_deg(lat, lon, lm_a.lat, lm_a.lon)
        bb = bearing_deg(lat, lon, lm_b.lat, lm_b.lon)
        return wrap180(bb - ba)

    def error(this, values, H):
        pose = values.atPose3(this.keys()[0])
        base = predict(pose)
        resid = wrap180(base - meas_angle_deg)
        if H is not None:
            jac = np.zeros((1, 6), order="F")
            eps = 1e-6
            for i in (3, 4):
                d = np.zeros(6); d[i] = eps
                jac[0, i] = wrap180(predict(pose.retract(d)) - base) / eps
            H[0] = jac
        return np.array([resid])

    return gtsam.CustomFactor(noise, [key], error)


# --------------------------------------------------------------------------- #
# A small standalone landmark fix (the terrain analogue of `solve`)
# --------------------------------------------------------------------------- #

def calibrated_sigma_deg(azimuths, residuals_deg, floor_deg=0.005):
    ''' Per-observation sigma that accounts for CORRELATED residual.

        A dense skyline match returns hundreds of samples, and treating them as
        independent is how an uncertainty comes to look an order of magnitude
        better than it is.  On the Istanbul frames, 839 samples at 0.14 deg would
        have given ~15 m of lateral position if the noise were white; they gave
        nothing, because the residual is smooth in azimuth -- extraction bias,
        haze-dependent edge placement, DEM error, a pitch bias drifting between
        frames -- and correlated error is precisely what a position shift looks
        like.

        So the sigma handed to a factor is the raw scatter INFLATED by
        sqrt(n / n_eff), which is equivalent to admitting you have n_eff
        independent looks rather than n.  Returns dict(sigma_deg, raw_deg,
        n, n_eff, inflation).
    '''
    from .resection_geometry import effective_samples
    r = np.asarray(residuals_deg, float)
    ok = np.isfinite(r)
    raw = float(np.std(r[ok])) if ok.sum() > 1 else float("nan")
    e = effective_samples(np.asarray(azimuths, float)[ok], r[ok])
    return dict(sigma_deg=max(raw * e["sigma_inflation"], floor_deg),
                raw_deg=raw, n=e["n"], n_eff=e["n_eff"],
                inflation=e["sigma_inflation"])


def solve_landmark_fix(bearings=(), elevations=(), horizontal_angles=(),
                       lat0: float = 0.0, lon0: float = 0.0,
                       prior_en_m=(0.0, 0.0), prior_sigma_km: float = 30.0,
                       cam_height_m: float = 2.0, k: float = K_REFRACTION,
                       compass_bias_sigma_deg: float = 0.0,
                       pitch_bias_sigma_deg: float = 0.0):
    ''' Fuse terrain-landmark observations into one position fix.

        bearings          : [(Landmark, measured_bearing_deg, sigma_deg), ...]
        elevations        : [(Landmark, measured_elev_deg,   sigma_deg), ...]
        horizontal_angles : [(LandmarkA, LandmarkB, measured_angle_deg, sigma_deg), ...]

        NUISANCE BIASES.  Give `compass_bias_sigma_deg` and/or
        `pitch_bias_sigma_deg` and the corresponding bias becomes an ESTIMATED
        VARIABLE shared by every observation of that type, with a prior of that
        width, rather than being smeared into each measurement's sigma.  This is
        the difference between a covariance you can act on and one that flatters
        itself: the two biases absorb exactly the two first-order signals a
        position shift produces, so with them in the graph the reported ellipse
        shows the real line-of-position elongation.  `resection_geometry` predicts
        the same ellipse from geometry alone, and the two agree.

        Set them to 0 (the default) to keep the old unary behaviour, in which
        case `sigma_deg` must already carry the heading error.

        Returns dict(lat, lon, east_m, north_m, cov_en, sigma_km, n_factors,
        compass_bias_deg, pitch_bias_deg, semi_major_km, semi_minor_km,
        major_bearing_deg).
    '''
    X = gtsam.symbol_shorthand.X
    B = gtsam.symbol_shorthand.B
    graph = gtsam.NonlinearFactorGraph()

    prior_pose = gtsam.Pose3(gtsam.Rot3(),
                             gtsam.Point3(prior_en_m[0], prior_en_m[1], 0.0))
    sig = np.array([1e-4, 1e-4, 1e-4,
                    prior_sigma_km * 1000.0, prior_sigma_km * 1000.0, 1e-3])
    graph.add(gtsam.PriorFactorPose3(
        X(0), prior_pose, gtsam.noiseModel.Diagonal.Sigmas(sig)))

    ck = B(0) if compass_bias_sigma_deg > 0 else None
    pk = B(1) if pitch_bias_sigma_deg > 0 else None
    if ck is not None:
        graph.add(gtsam.PriorFactorVector(
            ck, np.array([0.0]),
            gtsam.noiseModel.Isotropic.Sigma(1, compass_bias_sigma_deg)))
    if pk is not None:
        graph.add(gtsam.PriorFactorVector(
            pk, np.array([0.0]),
            gtsam.noiseModel.Isotropic.Sigma(1, pitch_bias_sigma_deg)))

    n = 0
    for lm, meas, s_ in bearings:
        graph.add(landmark_bearing_factor(X(0), lm, meas, s_, lat0, lon0,
                                          bias_key=ck)); n += 1
    for lm, meas, s_ in elevations:
        graph.add(landmark_elevation_factor(X(0), lm, meas, s_, lat0, lon0,
                                            cam_height_m=cam_height_m, k=k,
                                            bias_key=pk)); n += 1
    for lm_a, lm_b, meas, s_ in horizontal_angles:
        graph.add(landmark_horizontal_angle_factor(X(0), lm_a, lm_b, meas, s_,
                                                   lat0, lon0)); n += 1

    initial = gtsam.Values()
    initial.insert(X(0), prior_pose)
    if ck is not None:
        initial.insert(ck, np.array([0.0]))
    if pk is not None:
        initial.insert(pk, np.array([0.0]))
    result = gtsam.LevenbergMarquardtOptimizer(
        graph, initial, gtsam.LevenbergMarquardtParams()).optimize()

    pose = result.atPose3(X(0))
    t = pose.translation()
    lat, lon = enu_to_latlon(t[0], t[1], lat0, lon0)
    out = dict(lat=lat, lon=lon, east_m=float(t[0]), north_m=float(t[1]),
               cov_en=None, sigma_km=float("nan"), n_factors=n,
               compass_bias_deg=(float(result.atVector(ck)[0]) if ck else 0.0),
               pitch_bias_deg=(float(result.atVector(pk)[0]) if pk else 0.0),
               semi_major_km=float("nan"), semi_minor_km=float("nan"),
               major_bearing_deg=float("nan"))
    try:
        marg = gtsam.Marginals(graph, result)
        cov = marg.marginalCovariance(X(0))[3:5, 3:5]
        out["cov_en"] = cov
        out["sigma_km"] = float(np.sqrt(np.trace(cov)) / 1000.0)
        # the ellipse itself, because "one sigma" hides a line of position
        w, v = np.linalg.eigh(cov)
        out["semi_major_km"] = float(np.sqrt(max(w[1], 0.0)) / 1000.0)
        out["semi_minor_km"] = float(np.sqrt(max(w[0], 0.0)) / 1000.0)
        out["major_bearing_deg"] = float(
            degrees(atan2(v[0, 1], v[1, 1])) % 180.0)
    except Exception:                                    # pragma: no cover
        pass
    return out


def synthesize_measurements(lat: float, lon: float, landmarks,
                            cam_height_m: float = 2.0,
                            k: float = K_REFRACTION,
                            heading_bias_deg: float = 0.0):
    ''' Exact (noise-free) observations of `landmarks` from (lat, lon).

        `heading_bias_deg` is added to every BEARING (and therefore cancels in
        every horizontal angle) — the lever the tests use to show which factors
        depend on the compass.
    '''
    bearings, elevations, angles = [], [], []
    for lm in landmarks:
        b = bearing_deg(lat, lon, lm.lat, lm.lon) + heading_bias_deg
        bearings.append((lm, b % 360.0))
        d = max(great_circle_km(lat, lon, lm.lat, lm.lon), 1e-3)
        elevations.append((lm, vertical_angle_deg(d, lm.height_m,
                                                  cam_height_m, k)))
    for i in range(len(landmarks) - 1):
        a, b = landmarks[i], landmarks[i + 1]
        angles.append((a, b, wrap180(bearing_deg(lat, lon, b.lat, b.lon)
                                     - bearing_deg(lat, lon, a.lat, a.lon))))
    return bearings, elevations, angles
