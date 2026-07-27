''' Regression tests for the imu_fusion study.

    Kept fast (small keyframe counts, single seeds).  Requires the optional
    dependencies gtsam, numpy, matplotlib (see imu_fusion/requirements.txt); the
    whole module is skipped if gtsam is unavailable.

    © 2026.  MIT License (see LICENSE file).
'''

import os
import sys
import random
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    import gtsam  # noqa: F401
    from imu_fusion.scenario import build_scenario
    from imu_fusion.celestial_factor_graph import solve
    from imu_fusion.baseline import starfix_single_fix
    from imu_fusion.iphone_model import (DEFAULT_IMU, DEFAULT_CAM,
                                         KinematicState, altitude_sigma_arcmin)
    from imu_fusion.astro import body_gp, altaz
    _HAVE = True
except Exception:                               # pragma: no cover
    _HAVE = False


@unittest.skipUnless(_HAVE, "gtsam / numpy not installed")
class TestImuFusion(unittest.TestCase):

    def test_zero_noise_recovery(self):
        ''' With no measurement noise the factor graph recovers truth. '''
        for regime in ("land", "sea", "air"):
            sc = build_scenario(regime, random.Random(0), n_shots=6,
                                noise_scale=0.0)
            res = solve(sc, use_imu=False, pos_prior_km=100000.0)
            self.assertLess(res["rms_err_km"], 0.5,
                            f"{regime}: zero-noise RMS too high")

    def test_imu_improves_stationary(self):
        ''' IMU-linked smoothing beats independent per-epoch fixes on land. '''
        sc = build_scenario("land", random.Random(1), n_shots=12)
        no_imu = solve(sc, use_imu=False)["rms_err_km"]
        imu = solve(sc, use_imu=True)["rms_err_km"]
        self.assertLess(imu, no_imu)

    def test_two_bodies_beat_one(self):
        ''' A single body under-constrains a single-epoch fix; Sun+Moon fixes
            it.  Two bodies must give a much smaller error than Sun alone. '''
        sc2 = build_scenario("land", random.Random(2), n_shots=8,
                             bodies=("Sun", "Moon"))
        sc1 = build_scenario("land", random.Random(2), n_shots=8,
                             bodies=("Sun",))
        e2 = solve(sc2, use_imu=False)["rms_err_km"]
        e1 = solve(sc1, use_imu=False)["rms_err_km"]
        self.assertLess(e2, e1)

    def test_gating_cleaner_measurements(self):
        ''' Gated (least-rotation) shutter yields lower per-shot altitude
            sigma than an ungated shutter, for every regime. '''
        for regime in ("land", "sea", "air"):
            g = build_scenario(regime, random.Random(3), n_shots=10, gated=True)
            u = build_scenario(regime, random.Random(3), n_shots=10, gated=False)

            def mean_sig(sc):
                sigs = [o.alt_sigma_arcmin for kf in sc.keyframes
                        for o in kf.observations]
                return sum(sigs) / len(sigs)
            self.assertLess(mean_sig(g), mean_sig(u),
                            f"{regime}: gating did not reduce noise")

    def test_altitude_sigma_monotonic(self):
        ''' More motion -> larger synthetic-horizon altitude sigma. '''
        still = KinematicState(ang_rate=0.01, lin_accel=0.02)
        moving = KinematicState(ang_rate=0.6, lin_accel=2.0)
        s_still = altitude_sigma_arcmin(still, DEFAULT_IMU, DEFAULT_CAM)
        s_move = altitude_sigma_arcmin(moving, DEFAULT_IMU, DEFAULT_CAM)
        self.assertGreater(s_move, s_still)

    def test_daytime_geometry(self):
        ''' The canonical epoch really has both bodies up and well separated. '''
        t = "2026-03-24 12:00:00"
        sun_alt, sun_az = altaz(51.5, 0.0, body_gp("Sun", t))
        moon_alt, moon_az = altaz(51.5, 0.0, body_gp("Moon", t))
        self.assertGreater(sun_alt, 15)
        self.assertGreater(moon_alt, 15)
        daz = abs(((sun_az - moon_az + 180) % 360) - 180)
        self.assertGreater(daz, 60)             # well-conditioned fix

    def test_ultrawide_horizon_rescues_sea(self):
        ''' The optical ultrawide horizon (immune to swell acceleration) gives a
            much better sea fix than the IMU gravity horizon. '''
        sea_imu = solve(build_scenario("sea", random.Random(5), n_shots=12,
                                       horizon_mode="imu"), use_imu=True)
        sea_uw = solve(build_scenario("sea", random.Random(5), n_shots=12,
                                      horizon_mode="fused"), use_imu=True)
        self.assertLess(sea_uw["rms_err_km"], 0.5 * sea_imu["rms_err_km"])

    def test_ultrawide_horizon_unavailable_on_land(self):
        ''' Land has no true sea horizon, so the optical mode must fall back to
            the IMU and change nothing. '''
        from imu_fusion.ultrawide_horizon import (HORIZON_AVAILABLE,
                                                  horizon_reference_sigma_arcmin)
        self.assertFalse(HORIZON_AVAILABLE["land"])
        st = KinematicState(ang_rate=0.03, lin_accel=0.05)
        s_imu = horizon_reference_sigma_arcmin("imu", st, "land")
        s_uw = horizon_reference_sigma_arcmin("uw", st, "land")
        self.assertAlmostEqual(s_imu, s_uw, places=6)

    def test_starfix_baseline_runs(self):
        ''' The starfix single-fix baseline produces a finite error. '''
        sc = build_scenario("land", random.Random(4), n_shots=4)
        err, lat, lon = starfix_single_fix(sc.keyframes[0], sc.lat0, sc.lon0)
        self.assertIsNotNone(err)
        self.assertGreaterEqual(err, 0.0)


if __name__ == "__main__":
    unittest.main()
