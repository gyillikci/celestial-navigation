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

    def test_optical_disk_geometry(self):
        ''' Bright-limb PA, illuminated fraction and parallactic angle are
            sane at the canonical epoch, and q depends on latitude. '''
        from imu_fusion.optical_attitude import (moon_illuminated_fraction,
                                                 bright_limb_pa_deg,
                                                 parallactic_angle_deg,
                                                 moon_limb_available)
        from imu_fusion.astro import body_gp
        t = "2026-03-24 12:00:00"
        k = moon_illuminated_fraction(t)
        self.assertTrue(0.4 < k < 0.9)                 # waxing gibbous
        self.assertTrue(moon_limb_available(t))
        self.assertTrue(0.0 <= bright_limb_pa_deg(t) < 360.0)
        gm = body_gp("Moon", t)
        q1 = parallactic_angle_deg(50.0, 0.0, gm)
        q2 = parallactic_angle_deg(52.0, 0.0, gm)
        self.assertGreater(abs(q1 - q2), 1.0)          # q tracks latitude

    def test_optical_heading_beats_magnetometer(self):
        ''' The optical disk heading is tighter than the phone magnetometer. '''
        from imu_fusion.optical_attitude import optical_heading_sigma_deg
        from imu_fusion.iphone_model import heading_sigma_arcmin
        st = KinematicState(ang_rate=0.05, lin_accel=0.2)
        opt = optical_heading_sigma_deg("Moon", st)
        mag = heading_sigma_arcmin(st, DEFAULT_IMU) / 60.0
        self.assertLess(opt, mag)

    def test_optical_observables_help_sea(self):
        ''' Optical azimuth + parallactic line improve the sea fix on the (weak)
            IMU horizon, and zero-noise still recovers truth. '''
        base = build_scenario("sea", random.Random(6), n_shots=12,
                              horizon_mode="imu")
        opt = build_scenario("sea", random.Random(6), n_shots=12,
                             horizon_mode="imu", use_azimuth=True,
                             heading_source="optical", use_parallactic=True)
        e_base = solve(base, use_imu=True)["rms_err_km"]
        e_opt = solve(opt, use_imu=True, use_azimuth=True,
                      use_parallactic=True)["rms_err_km"]
        self.assertLess(e_opt, e_base)
        zn = build_scenario("sea", random.Random(6), n_shots=6, noise_scale=0.0,
                            use_azimuth=True, heading_source="optical",
                            use_parallactic=True)
        r = solve(zn, use_imu=False, use_azimuth=True, use_parallactic=True,
                  pos_prior_km=100000.0)
        self.assertLess(r["rms_err_km"], 0.5)

    def test_full_fusion_all_observables(self):
        ''' The unified graph with every observable on runs, recovers zero-noise
            truth, and beats the ultrawide-horizon-only baseline is not required
            but it must be well under 10 km at sea. '''
        full_sc = dict(horizon_mode="fused", use_azimuth=True,
                       heading_source="optical", use_parallactic=True,
                       sun_spots=True, gated=True)
        full_sv = dict(use_imu=True, use_azimuth=True, use_parallactic=True)
        sc = build_scenario("sea", random.Random(9), n_shots=12, **full_sc)
        self.assertLess(solve(sc, **full_sv)["rms_err_km"], 10.0)
        zn = build_scenario("sea", random.Random(9), n_shots=6, noise_scale=0.0,
                            **full_sc)
        self.assertLess(solve(zn, pos_prior_km=100000.0,
                              **full_sv)["rms_err_km"], 0.5)

    def test_ablation_ultrawide_matters_at_sea(self):
        ''' Removing the ultrawide horizon from the full fusion must hurt the
            sea fix (it is the sea's dominant observable). '''
        full_sc = dict(horizon_mode="fused", use_azimuth=True,
                       heading_source="optical", use_parallactic=True,
                       sun_spots=True)
        full_sv = dict(use_imu=True, use_azimuth=True, use_parallactic=True)
        full = solve(build_scenario("sea", random.Random(11), n_shots=12,
                                    **full_sc), **full_sv)["rms_err_km"]
        no_uw = solve(build_scenario("sea", random.Random(11), n_shots=12,
                                     **{**full_sc, "horizon_mode": "imu"}),
                      **full_sv)["rms_err_km"]
        self.assertLess(full, no_uw)

    def test_analytic_altitude_jacobian(self):
        ''' The closed-form altitude Jacobian matches central finite diff. '''
        import numpy as np
        import gtsam
        from gtsam import Pose3, Rot3, Point3
        from imu_fusion.astro import (body_gp, enu_to_latlon,
                                      predicted_altitude, predicted_azimuth)
        from imu_fusion.celestial_factor_graph import _altitude_analytic_jacobian
        from imu_fusion.numerical_derivative import numericalDerivative11
        gp = body_gp("Sun", "2026-03-24 12:00:00")
        pose = Pose3(Rot3(), Point3(2500.0, -1500.0, 30.0))

        def predict(p):
            t = p.translation()
            lat, lon = enu_to_latlon(t[0], t[1], 51.5, 0.0)
            return predicted_altitude(lat, lon, gp)
        num = numericalDerivative11(predict, pose, 1e-4).flatten()
        lat, lon = enu_to_latlon(pose.translation()[0], pose.translation()[1],
                                 51.5, 0.0)
        ana = np.asarray(_altitude_analytic_jacobian(
            pose, predicted_azimuth(lat, lon, gp))).flatten()
        np.testing.assert_allclose(ana, num, atol=1e-6)

    def test_streaming_matches_batch(self):
        ''' The streaming smoother's current-position error matches batch, and
            recovers zero-noise truth. '''
        from imu_fusion.realtime import solve_streaming
        full = dict(horizon_mode="fused", use_azimuth=True,
                    heading_source="optical", use_parallactic=True,
                    sun_spots=True)
        sc = build_scenario("sea", random.Random(21), n_shots=12, **full)
        batch = solve(sc, use_imu=True, use_azimuth=True,
                      use_parallactic=True)["final_err_km"]
        stream = solve_streaming(sc)["final_err_km"]
        self.assertLess(abs(batch - stream), 1.0)     # within noise
        zn = build_scenario("sea", random.Random(21), n_shots=8,
                            noise_scale=0.0, **full)
        self.assertLess(solve_streaming(zn, pos_prior_km=100000.0)["rms_err_km"],
                        0.6)

    def test_ephemeris_cache(self):
        ''' The GP cache returns an identical GP and is populated. '''
        from imu_fusion import astro
        astro.clear_gp_cache()
        g1 = astro.body_gp("Moon", "2026-03-24 12:00:00")
        g2 = astro.body_gp("Moon", "2026-03-24 12:00:00")
        self.assertEqual((g1.get_lat(), g1.get_lon()),
                         (g2.get_lat(), g2.get_lon()))
        self.assertGreaterEqual(len(astro._GP_CACHE), 1)

    def test_teleconverter_does_not_move_fix(self):
        ''' A 3x external optic sharpens pointing but the horizon-limited fix is
            unchanged; the camera contribution stays far below the horizon. '''
        from dataclasses import replace
        from imu_fusion.iphone_model import DEFAULT_CAM
        full = dict(horizon_mode="fused", use_azimuth=True,
                    heading_source="optical", use_parallactic=True,
                    sun_spots=True)
        sv = dict(use_imu=True, use_azimuth=True, use_parallactic=True)
        e1 = solve(build_scenario("sea", random.Random(31), n_shots=12,
                                  cam=DEFAULT_CAM, **full), **sv)["rms_err_km"]
        cam3 = replace(DEFAULT_CAM, teleconverter=3.0)
        e3 = solve(build_scenario("sea", random.Random(31), n_shots=12,
                                  cam=cam3, **full), **sv)["rms_err_km"]
        self.assertAlmostEqual(e1, e3, delta=0.2)          # fix unchanged
        # pointing improves 3x but is negligible vs the horizon reference
        self.assertLess(cam3.pointing_sigma_arcmin() * 3, 0.2)

    def test_horizon_lens_altitude_dependence(self):
        ''' The wide lens frames the horizon only for low bodies; at the
            canonical epoch (bodies ~32-40 deg) forcing it loses the horizon and
            wrecks the fix, while adaptive matches ultrawide. '''
        from imu_fusion.ultrawide_horizon import (lens_sees_horizon, WIDE_LENS,
                                                  ULTRAWIDE_LENS,
                                                  best_horizon_lens)
        self.assertTrue(lens_sees_horizon(WIDE_LENS, 20))
        self.assertFalse(lens_sees_horizon(WIDE_LENS, 40))
        self.assertTrue(lens_sees_horizon(ULTRAWIDE_LENS, 40))
        self.assertFalse(lens_sees_horizon(ULTRAWIDE_LENS, 60))
        self.assertEqual(best_horizon_lens(20, "sea").name, "wide")
        self.assertEqual(best_horizon_lens(40, "sea").name, "ultrawide")
        self.assertIsNone(best_horizon_lens(60, "sea"))       # IMU only
        full = dict(horizon_mode="fused", use_azimuth=True,
                    heading_source="optical", use_parallactic=True,
                    sun_spots=True)
        sv = dict(use_imu=True, use_azimuth=True, use_parallactic=True)

        def mean_rms(lens):
            return sum(solve(build_scenario("sea", random.Random(40 + s),
                                            n_shots=12, horizon_lens=lens,
                                            **full), **sv)["rms_err_km"]
                       for s in range(4)) / 4.0
        wide, adap, uw = mean_rms("wide"), mean_rms("adaptive"), mean_rms("ultrawide")
        self.assertGreater(wide, 2 * adap)      # wide-only loses the horizon
        self.assertAlmostEqual(adap, uw, delta=0.3)

    def test_starfix_baseline_runs(self):
        ''' The starfix single-fix baseline produces a finite error. '''
        sc = build_scenario("land", random.Random(4), n_shots=4)
        err, lat, lon = starfix_single_fix(sc.keyframes[0], sc.lat0, sc.lon0)
        self.assertIsNotNone(err)
        self.assertGreaterEqual(err, 0.0)


if __name__ == "__main__":
    unittest.main()
