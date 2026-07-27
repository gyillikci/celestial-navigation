''' Streaming (near-real-time) estimator for on-device use.

    The batch `celestial_factor_graph.solve` re-optimises the whole trajectory on
    every call, so its cost grows with trip length — fine for an offline study,
    wrong for a phone. This module runs the same fusion INCREMENTALLY with a
    fixed-lag smoother: each gated Sun+Moon shot is one `update()` over only the
    new factors, and keyframes older than a time window are marginalised out, so
    the per-update cost and memory are BOUNDED regardless of voyage length.

    Celestial fixes are inherently low-rate (a clean shot every few seconds); the
    continuous position between shots comes from IMU dead-reckoning (the
    preintegration accumulated here). So "real-time" means each per-shot update
    finishes well within the inter-shot interval — tens of ms.

    Uses `gtsam_unstable.IncrementalFixedLagSmoother` (the fixed-lag smoothers
    live in gtsam_unstable in the 4.2.1 wheel) with the analytic-Jacobian
    factors from `celestial_factor_graph`.

    (c) 2026.  MIT License (see LICENSE file).
'''

import numpy as np
import gtsam
import gtsam_unstable
from gtsam import Pose3, Rot3, Point3
from gtsam.symbol_shorthand import X, V, B

from .astro import enu_to_latlon, great_circle_km
from .celestial_factor_graph import (celestial_altitude_factor,
                                     celestial_azimuth_factor,
                                     parallactic_angle_factor,
                                     _imu_params)


class StreamingEstimator:
    ''' Incremental fixed-lag smoother fed one keyframe (shot) at a time. '''

    def __init__(self, scenario, lag_s=90.0, use_imu=True, use_azimuth=True,
                 use_parallactic=True, pos_prior_km=1000.0,
                 relinearize_threshold=0.01):
        self.sc = scenario
        self.lat0, self.lon0 = scenario.lat0, scenario.lon0
        self.use_imu = use_imu
        self.use_azimuth = use_azimuth
        self.use_parallactic = use_parallactic

        params = gtsam.ISAM2Params()
        params.setRelinearizeThreshold(relinearize_threshold)
        params.relinearizeSkip = 1
        self.smoother = gtsam_unstable.IncrementalFixedLagSmoother(lag_s, params)

        pm = pos_prior_km * 1000.0
        self.pose_prior_noise = gtsam.noiseModel.Diagonal.Sigmas(
            np.array([1e-3, 1e-3, 1e-3, pm, pm, 2.0]))
        self.anchor = Pose3(Rot3(), Point3(0.0, 0.0, 0.0))
        self.v_init = np.array([scenario.vel_east, scenario.vel_north, 0.0])

        self.imu_params = _imu_params(scenario.imu)
        self.bias0 = gtsam.imuBias.ConstantBias()
        self.pim = gtsam.PreintegratedImuMeasurements(self.imu_params,
                                                      self.bias0)
        self.results = []            # per-keyframe estimate dicts
        self._last_time = None

    # -- online IMU ------------------------------------------------------- #
    def integrate_imu(self, samples):
        ''' Accumulate a leg of IMU samples between keyframes as they arrive. '''
        for acc, gyr, dt in samples:
            self.pim.integrateMeasurement(np.array(acc), np.array(gyr), dt)

    # -- one shot --------------------------------------------------------- #
    def add_keyframe(self, kf):
        ''' Add one gated Sun+Moon shot and return its updated estimate. '''
        graph = gtsam.NonlinearFactorGraph()
        values = gtsam.Values()
        ts = gtsam_unstable.FixedLagSmootherKeyTimestampMap()
        i = kf.index
        t = kf.time_s

        graph.add(gtsam.PriorFactorPose3(X(i), self.anchor,
                                         self.pose_prior_noise))
        dr = Pose3(Rot3(), Point3(self.sc.vel_east * t,
                                  self.sc.vel_north * t, 0.0))
        values.insert(X(i), dr)
        ts.insert((X(i), t))

        for o in kf.observations:
            graph.add(celestial_altitude_factor(
                X(i), o.gp, o.meas_alt, o.alt_sigma_arcmin, self.lat0, self.lon0))
            if self.use_azimuth:
                graph.add(celestial_azimuth_factor(
                    X(i), o.gp, o.meas_az, o.az_sigma_arcmin,
                    self.lat0, self.lon0))
            if self.use_parallactic and o.par_valid:
                graph.add(parallactic_angle_factor(
                    X(i), o.gp, o.par_meas, o.par_sigma_deg,
                    self.lat0, self.lon0))

        if self.use_imu:
            values.insert(V(i), self.v_init)
            ts.insert((V(i), t))
            if i == 0:
                values.insert(B(0), self.bias0)
                graph.add(gtsam.PriorFactorVector(
                    V(0), self.v_init, gtsam.noiseModel.Isotropic.Sigma(3, 5.0)))
                graph.add(gtsam.PriorFactorConstantBias(
                    B(0), self.bias0,
                    gtsam.noiseModel.Isotropic.Sigma(6, 0.1)))
            else:
                graph.add(gtsam.ImuFactor(X(i - 1), V(i - 1), X(i), V(i),
                                          B(0), self.pim))
                self.pim.resetIntegration()
            # Keep the shared bias inside the lag window (re-timestamp to now).
            ts.insert((B(0), t))

        self.smoother.update(graph, values, ts)
        est = self.smoother.calculateEstimate()
        pose = est.atPose3(X(i))
        tr = pose.translation()
        lat, lon = enu_to_latlon(tr[0], tr[1], self.lat0, self.lon0)
        err = great_circle_km(lat, lon, kf.true_lat, kf.true_lon)
        cov = None
        try:
            cov = self.smoother.getISAM2().marginalCovariance(X(i))[3:5, 3:5]
        except Exception:
            cov = None
        out = dict(index=i, est_lat=lat, est_lon=lon, err_km=err,
                   true_lat=kf.true_lat, true_lon=kf.true_lon, cov_en=cov)
        self.results.append(out)
        self._last_time = t
        return out


def solve_streaming(scenario, lag_s=90.0, use_imu=True, use_azimuth=True,
                    use_parallactic=True, pos_prior_km=1000.0):
    ''' Run a whole scenario through the streaming estimator (feeding IMU legs
        online) and return a result dict matching the batch `solve`. '''
    est = StreamingEstimator(scenario, lag_s, use_imu, use_azimuth,
                             use_parallactic, pos_prior_km)
    kfs = scenario.keyframes
    for k, kf in enumerate(kfs):
        est.add_keyframe(kf)
        if use_imu and k < len(kfs) - 1:
            est.integrate_imu(kf.imu_to_next)

    errs = np.array([r["err_km"] for r in est.results])
    covs = [r["cov_en"] for r in est.results if r["cov_en"] is not None]
    sig_km = (float(np.mean([np.sqrt(np.trace(c)) / 1000.0 for c in covs]))
              if covs else None)
    return dict(regime=scenario.regime, per_kf=est.results,
                rms_err_km=float(np.sqrt(np.mean(errs ** 2))),
                mean_err_km=float(np.mean(errs)),
                final_err_km=float(errs[-1]),
                sigma_km=sig_km)


if __name__ == "__main__":
    import random
    from .scenario import build_scenario
    for r in ("land", "sea", "air"):
        sc = build_scenario(r, random.Random(4), n_shots=14, horizon_mode="fused",
                            use_azimuth=True, heading_source="optical",
                            use_parallactic=True, sun_spots=True)
        res = solve_streaming(sc)
        print(f"{r:5} streaming RMS={res['rms_err_km']:.2f} km  "
              f"final={res['final_err_km']:.2f} km  sigma={res['sigma_km']}")
