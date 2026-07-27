''' GTSAM factor-graph fusion of Sun+Moon altitude sights with phone IMU.

    Estimates the observer trajectory (horizontal position, and for moving
    regimes velocity) from a sequence of daytime Sun/Moon altitude measurements,
    optionally linked by IMU preintegration between shots.

    Graph structure
    ---------------
      Variables
        X(i) : Pose3   -- observer pose at keyframe i.  Translation is a local
                          ENU offset (metres) from the anchor (lat0, lon0);
                          rotation is the (level) body attitude.
        V(i) : Vector3 -- ENU velocity (only when use_imu=True).
        B    : imuBias.ConstantBias -- shared IMU bias (use_imu=True).

      Factors
        * Celestial altitude  : a gtsam.CustomFactor per (keyframe, body).
          residual = predicted_altitude(pos(X_i), body_GP) - measured_altitude.
          This is the classic celestial line of position; two bodies with
          separated azimuths fix the horizontal position at each epoch.
        * Attitude+height anchor : a PriorFactorPose3 per keyframe with TIGHT
          rotation & Up sigmas (the phone always knows attitude from gravity and
          height from baro/known elevation) but LOOSE East/North (the celestial
          factors, not the prior, determine horizontal position).
        * IMU (optional) : gtsam.ImuFactor between consecutive keyframes from
          preintegrated phone IMU, tying the epochs through the velocity state
          so many noisy fixes are smoothed into one trajectory.
        * Velocity & bias priors to fix the remaining gauge.

    Solve with Levenberg-Marquardt; report per-keyframe position and its
    marginal covariance (the error ellipse).

    (c) 2026.  MIT License (see LICENSE file).
'''

import numpy as np
import gtsam
from gtsam import Pose3, Rot3, Point3
from gtsam.symbol_shorthand import X, V, B

from .astro import predicted_altitude, predicted_azimuth, enu_to_latlon, \
    great_circle_km
from .iphone_model import G0

_ARCMIN = 1.0 / 60.0


# --------------------------------------------------------------------------- #
# Measurement factors (Python CustomFactor closures)
# --------------------------------------------------------------------------- #

def _numeric_pose_jacobian(pose, scalar_fn, base_val, eps=1e-6):
    ''' Finite-difference Jacobian (1x6, Fortran order) of a scalar function of
        a Pose3, taken on the Pose3 retract tangent [rx,ry,rz,tx,ty,tz]. '''
    jac = np.zeros((1, 6), order="F")
    for i in range(6):
        d = np.zeros(6)
        d[i] = eps
        jac[0, i] = (scalar_fn(pose.retract(d)) - base_val) / eps
    return jac


def celestial_altitude_factor(key, gp, meas_alt_deg, sigma_arcmin,
                              lat0, lon0):
    ''' A unary altitude line-of-position factor on a Pose3. '''
    noise = gtsam.noiseModel.Isotropic.Sigma(1, sigma_arcmin * _ARCMIN)

    def predict(pose):
        t = pose.translation()
        lat, lon = enu_to_latlon(t[0], t[1], lat0, lon0)
        return predicted_altitude(lat, lon, gp)

    def error(this, values, H):
        pose = values.atPose3(this.keys()[0])
        resid = predict(pose) - meas_alt_deg
        if H is not None:
            H[0] = _numeric_pose_jacobian(pose, predict, predict(pose))
        return np.array([resid])

    return gtsam.CustomFactor(noise, [key], error)


def celestial_azimuth_factor(key, gp, meas_az_deg, sigma_arcmin, lat0, lon0):
    ''' A unary azimuth factor on a Pose3 (heading-dependent, weaker). '''
    noise = gtsam.noiseModel.Isotropic.Sigma(1, sigma_arcmin * _ARCMIN)

    def predict(pose):
        t = pose.translation()
        lat, lon = enu_to_latlon(t[0], t[1], lat0, lon0)
        return predicted_azimuth(lat, lon, gp)

    def error(this, values, H):
        pose = values.atPose3(this.keys()[0])
        # Wrap the azimuth residual into [-180, 180].
        resid = (predict(pose) - meas_az_deg + 180.0) % 360.0 - 180.0
        if H is not None:
            H[0] = _numeric_pose_jacobian(pose, predict, predict(pose))
        return np.array([resid])

    return gtsam.CustomFactor(noise, [key], error)


# --------------------------------------------------------------------------- #
# IMU preintegration
# --------------------------------------------------------------------------- #

def _imu_params(imu):
    ''' ENU ("U") preintegration params from the phone IMU spec. '''
    params = gtsam.PreintegrationParams.MakeSharedU(G0)
    a2 = imu.accel_sigma() ** 2
    g2 = imu.gyro_sigma() ** 2
    params.setAccelerometerCovariance(np.eye(3) * a2)
    params.setGyroscopeCovariance(np.eye(3) * g2)
    params.setIntegrationCovariance(np.eye(3) * 1e-8)
    return params


def _preintegrate(samples, params, bias):
    ''' Integrate a list of (accel, gyro, dt) into a PIM. '''
    pim = gtsam.PreintegratedImuMeasurements(params, bias)
    for acc, gyr, dt in samples:
        pim.integrateMeasurement(np.array(acc), np.array(gyr), dt)
    return pim


# --------------------------------------------------------------------------- #
# Graph construction & solve
# --------------------------------------------------------------------------- #

def build_graph(scenario, use_imu=True, use_azimuth=False,
                pos_prior_km=1000.0):
    ''' Build the factor graph and initial values for a scenario. '''
    graph = gtsam.NonlinearFactorGraph()
    initial = gtsam.Values()
    kfs = scenario.keyframes
    lat0, lon0 = scenario.lat0, scenario.lon0
    imu = scenario.imu

    # Attitude(tight) + Up(tight) + horizontal(loose) anchor per keyframe.
    pos_prior_m = pos_prior_km * 1000.0
    pose_prior_sigmas = np.array([1e-3, 1e-3, 1e-3,          # rotation (rad)
                                  pos_prior_m, pos_prior_m,  # E, N (loose)
                                  2.0])                      # U (tight)
    pose_prior_noise = gtsam.noiseModel.Diagonal.Sigmas(pose_prior_sigmas)
    anchor_pose = Pose3(Rot3(), Point3(0.0, 0.0, 0.0))

    for kf in kfs:
        graph.add(gtsam.PriorFactorPose3(X(kf.index), anchor_pose,
                                         pose_prior_noise))
        # Dead-reckoned initial guess: anchor origin advanced by the (roughly
        # known) velocity.  Horizontal position is ~dr_error_km from truth until
        # the celestial factors correct it.
        dr = Pose3(Rot3(), Point3(scenario.vel_east * kf.time_s,
                                  scenario.vel_north * kf.time_s, 0.0))
        initial.insert(X(kf.index), dr)
        for o in kf.observations:
            graph.add(celestial_altitude_factor(
                X(kf.index), o.gp, o.meas_alt, o.alt_sigma_arcmin, lat0, lon0))
            if use_azimuth:
                graph.add(celestial_azimuth_factor(
                    X(kf.index), o.gp, o.meas_az, o.az_sigma_arcmin,
                    lat0, lon0))

    if use_imu:
        params = _imu_params(imu)
        bias0 = gtsam.imuBias.ConstantBias()
        v_init = np.array([scenario.vel_east, scenario.vel_north, 0.0])
        for kf in kfs:
            initial.insert(V(kf.index), v_init)
        initial.insert(B(0), bias0)
        # Priors to fix velocity/bias gauge.
        graph.add(gtsam.PriorFactorVector(
            V(0), v_init, gtsam.noiseModel.Isotropic.Sigma(3, 5.0)))
        graph.add(gtsam.PriorFactorConstantBias(
            B(0), bias0, gtsam.noiseModel.Isotropic.Sigma(6, 0.1)))
        for k in range(len(kfs) - 1):
            pim = _preintegrate(kfs[k].imu_to_next, params, bias0)
            graph.add(gtsam.ImuFactor(X(k), V(k), X(k + 1), V(k + 1),
                                      B(0), pim))
    return graph, initial


def solve(scenario, use_imu=True, use_azimuth=False, pos_prior_km=1000.0):
    ''' Optimise and return a result dict with per-keyframe estimates,
        covariances and errors vs. ground truth. '''
    graph, initial = build_graph(scenario, use_imu, use_azimuth, pos_prior_km)
    params = gtsam.LevenbergMarquardtParams()
    opt = gtsam.LevenbergMarquardtOptimizer(graph, initial, params)
    result = opt.optimize()

    try:
        marg = gtsam.Marginals(graph, result)
    except Exception:
        marg = None

    lat0, lon0 = scenario.lat0, scenario.lon0
    per_kf = []
    errs = []
    for kf in scenario.keyframes:
        pose = result.atPose3(X(kf.index))
        t = pose.translation()
        lat, lon = enu_to_latlon(t[0], t[1], lat0, lon0)
        err_km = great_circle_km(lat, lon, kf.true_lat, kf.true_lon)
        errs.append(err_km)
        cov_en = None
        if marg is not None:
            try:
                c = marg.marginalCovariance(X(kf.index))
                cov_en = c[3:5, 3:5].copy()      # E,N block (m^2)
            except Exception:
                cov_en = None
        per_kf.append(dict(index=kf.index, est_lat=lat, est_lon=lon,
                           true_lat=kf.true_lat, true_lon=kf.true_lon,
                           east=t[0], north=t[1], err_km=err_km,
                           cov_en=cov_en))

    errs = np.array(errs)
    # 1-sigma horizontal uncertainty (mean over keyframes), km.
    sig_km = None
    covs = [d["cov_en"] for d in per_kf if d["cov_en"] is not None]
    if covs:
        drms = [np.sqrt(np.trace(c)) / 1000.0 for c in covs]
        sig_km = float(np.mean(drms))

    return dict(regime=scenario.regime, use_imu=use_imu,
                use_azimuth=use_azimuth, per_kf=per_kf,
                rms_err_km=float(np.sqrt(np.mean(errs ** 2))),
                mean_err_km=float(np.mean(errs)),
                max_err_km=float(np.max(errs)),
                final_err_km=float(errs[-1]),
                sigma_km=sig_km,
                error_ellipse_cov=covs[-1] if covs else None)


if __name__ == "__main__":
    import random
    from .scenario import build_scenario
    for r in ("land", "sea", "air"):
        sc = build_scenario(r, random.Random(7), n_shots=12)
        res = solve(sc, use_imu=True)
        print(f"{r:5} FG(IMU)  RMS={res['rms_err_km']:.2f} km  "
              f"final={res['final_err_km']:.2f} km  "
              f"sigma={res['sigma_km']}")
