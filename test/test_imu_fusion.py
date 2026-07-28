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
        ''' Without a horizon-free surface-feature anchor, the ultrawide optical
            horizon is the sea's dominant vertical reference, so removing it must
            hurt.  (WITH resolved sunspots/craters the horizon-free parallactic
            substitutes and the ultrawide becomes marginal -- exercised in
            test_sunspot_anchor_makes_horizon_redundant.) '''
        base = dict(use_azimuth=True, heading_source="magnetometer",
                    use_parallactic=False, sun_spots=False)
        sv = dict(use_imu=True, use_azimuth=True, use_parallactic=False)
        full = sum(solve(build_scenario("sea", random.Random(11 + s), n_shots=12,
                                        horizon_mode="fused", **base),
                         **sv)["rms_err_km"] for s in range(4)) / 4
        no_uw = sum(solve(build_scenario("sea", random.Random(11 + s), n_shots=12,
                                         horizon_mode="imu", **base),
                          **sv)["rms_err_km"] for s in range(4)) / 4
        self.assertLess(full, no_uw)

    def test_sunspot_anchor_makes_horizon_redundant(self):
        ''' With both disks resolved, the horizon-free differential Sun-Moon
            parallactic (plus the anchored IMU vertical) carries the fix, so
            losing the ultrawide optical horizon at sea barely hurts. '''
        sc = dict(use_azimuth=True, heading_source="optical",
                  use_parallactic=True, sun_spots=True, imu_anchor=True)
        sv = dict(use_imu=True, use_azimuth=True, use_parallactic=True)
        full = sum(solve(build_scenario("sea", random.Random(100 + s), n_shots=12,
                                        horizon_mode="fused", **sc),
                         **sv)["rms_err_km"] for s in range(4)) / 4
        no_uw = sum(solve(build_scenario("sea", random.Random(100 + s), n_shots=12,
                                         horizon_mode="imu", **sc),
                          **sv)["rms_err_km"] for s in range(4)) / 4
        self.assertLess(abs(no_uw - full), 0.6)        # ultrawide now marginal

    def test_differential_parallactic_is_horizon_free(self):
        ''' The differential Sun-Moon parallactic observable's sigma does not
            depend on the horizon: it is identical with a good optical horizon
            and with none, because the platform roll cancels between the two
            disks (no vertical reference enters). '''
        import random as _r
        fused = build_scenario("sea", _r.Random(7), n_shots=8,
                               horizon_mode="fused", use_parallactic=True,
                               sun_spots=True)
        imu = build_scenario("sea", _r.Random(7), n_shots=8,
                             horizon_mode="imu", use_parallactic=True,
                             sun_spots=True)
        sig_f = [kf.diff_q_sigma for kf in fused.keyframes if kf.diff_valid]
        sig_i = [kf.diff_q_sigma for kf in imu.keyframes if kf.diff_valid]
        self.assertTrue(sig_f and len(sig_f) == len(sig_i))
        for a, b in zip(sig_f, sig_i):
            self.assertAlmostEqual(a, b, places=9)     # horizon plays no part
        self.assertLess(max(sig_f), 0.3)               # ~0.13 deg, sub-degree

    def test_differential_rescues_horizon_denied(self):
        ''' With NO optical horizon and NO IMU anchor, the horizon-free
            differential Sun-Moon parallactic line still sharply improves the sea
            fix; removing just the differential factor (use_differential=False)
            loses that gain. '''
        def mean(par, diff):
            tot = 0.0
            for s in range(4):
                sc = build_scenario("sea", random.Random(300 + s), n_shots=12,
                                    horizon_mode="imu", use_azimuth=True,
                                    heading_source="optical", use_parallactic=par,
                                    sun_spots=True, imu_anchor=False)
                tot += solve(sc, use_imu=True, use_azimuth=True,
                             use_parallactic=par,
                             use_differential=diff)["rms_err_km"]
            return tot / 4
        no_par = mean(False, False)
        with_diff = mean(True, True)
        par_no_diff = mean(True, False)
        self.assertLess(with_diff, 0.6 * no_par)       # big horizon-free rescue
        self.assertLess(with_diff, par_no_diff)        # the differential does it

    def test_sequential_gap_grows_sigma_but_dr_removes_bias(self):
        ''' One phone = sequential shots.  The differential sigma grows with the
            inter-shot slew gap (gyro roll carried across the slew), while the
            dead-reckoned Moon offset removes the v*gap translation, so even a big
            gap in the AIR recovers truth under zero noise (no bias). '''
        sigmas = []
        for gap in (0.0, 5.0, 30.0):
            sc = build_scenario("air", random.Random(9), n_shots=8,
                                horizon_mode="imu", use_parallactic=True,
                                sun_spots=True, intershot_gap_s=gap)
            s = [kf.diff_q_sigma for kf in sc.keyframes if kf.diff_valid]
            off = [kf.diff_moon_off_en for kf in sc.keyframes if kf.diff_valid]
            sigmas.append(s[0])
            if gap > 0:
                self.assertGreater(abs(off[0][0]) + abs(off[0][1]), 1.0)  # DR moved
        self.assertLess(sigmas[0], sigmas[1])          # sigma grows with gap
        self.assertLess(sigmas[1], sigmas[2])
        # zero-noise air with a large gap: DR offset removes the translation bias
        sc = build_scenario("air", random.Random(10), n_shots=6, noise_scale=0.0,
                            horizon_mode="imu", use_azimuth=True,
                            heading_source="optical", use_parallactic=True,
                            sun_spots=True, intershot_gap_s=20.0)
        rms = solve(sc, use_imu=True, use_azimuth=True,
                    use_parallactic=True)["rms_err_km"]
        self.assertLess(rms, 0.5)                       # no translation bias

    def test_graph_export_matches_real_graph(self):
        ''' The exported graph structure (for the interactive viewer) has exactly
            the same factor count and variable ids as the real GTSAM graph, so the
            visualisation cannot drift from what is actually solved. '''
        from imu_fusion.graph_export import graph_structure, write_graph_viewer
        sc = build_scenario("sea", random.Random(0), n_shots=5,
                            horizon_mode="fused", use_azimuth=True,
                            heading_source="optical", use_parallactic=True,
                            sun_spots=True)
        st = graph_structure(sc, lag_s=30.0, validate=True)   # asserts internally
        ids = {n["id"] for n in st["nodes"]}
        for f in st["factors"]:                                # edges are valid
            for v in f["vars"]:
                self.assertIn(v, ids)
        # the horizon-free differential appears once per resolved keyframe
        ndiff = sum(1 for f in st["factors"] if f["type"] == "diff")
        self.assertEqual(ndiff, sum(1 for kf in sc.keyframes if kf.diff_valid))
        # marginalisation is causal (a factor leaves the window after it enters)
        for f in st["factors"]:
            if f["marg"] is not None:
                self.assertGreater(f["marg"], f["step"])
        # the HTML renders self-contained (no external hosts, data injected)
        import tempfile, os
        p = os.path.join(tempfile.mkdtemp(), "g.html")
        write_graph_viewer(st, p)
        html = open(p).read()
        self.assertNotIn("__DATA__", html)
        self.assertNotIn("http://", html.replace("http://www.w3.org/2000/svg", ""))

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
        # Isolate the horizon-LENS dependence with no horizon-free surface-feature
        # rescue (no sunspots/parallactic), so the optical horizon is the vertical.
        full = dict(horizon_mode="fused", use_azimuth=True,
                    heading_source="magnetometer", use_parallactic=False,
                    sun_spots=False)
        sv = dict(use_imu=True, use_azimuth=True, use_parallactic=False)

        def mean_rms(lens):
            return sum(solve(build_scenario("sea", random.Random(40 + s),
                                            n_shots=12, horizon_lens=lens,
                                            **full), **sv)["rms_err_km"]
                       for s in range(4)) / 4.0
        wide, adap, uw = mean_rms("wide"), mean_rms("adaptive"), mean_rms("ultrawide")
        self.assertGreater(wide, 1.4 * adap)    # wide-only loses the horizon
        self.assertAlmostEqual(adap, uw, delta=0.3)

    def test_visual_anchor_bounds_drift_and_is_motion_immune(self):
        ''' Gyro-only attitude diverges with time; the anchor stays bounded and
            is nearly identical moving vs stationary (acceleration-immune). '''
        from imu_fusion.visual_anchor import (gyro_only_attitude_arcmin,
                                              anchor_attitude_arcmin)
        from imu_fusion.iphone_model import KinematicState
        self.assertGreater(gyro_only_attitude_arcmin(300),
                           5 * gyro_only_attitude_arcmin(10))     # diverges
        still = KinematicState(0.02, 0.03)
        moving = KinematicState(0.30, 1.5)
        a_still = anchor_attitude_arcmin(30, still)
        a_move = anchor_attitude_arcmin(30, moving)
        self.assertLess(a_still, 6.0)                             # bounded
        self.assertAlmostEqual(a_still, a_move, delta=0.5)        # motion-immune

    def test_visual_anchor_rescues_lost_horizon(self):
        ''' With no optical horizon (IMU horizon) the anchor sharply improves the
            fix; zero-noise still recovers truth. '''
        full = dict(horizon_mode="imu", use_azimuth=True,
                    heading_source="optical", use_parallactic=True,
                    sun_spots=True)
        sv = dict(use_imu=True, use_azimuth=True, use_parallactic=True)

        def mean_rms(anchor):
            return sum(solve(build_scenario("sea", random.Random(50 + s),
                                            n_shots=12, imu_anchor=anchor,
                                            **full), **sv)["rms_err_km"]
                       for s in range(4)) / 4.0
        self.assertLess(mean_rms(True), 0.5 * mean_rms(False))
        zn = build_scenario("sea", random.Random(3), n_shots=8, noise_scale=0.0,
                            imu_anchor=True, **full)
        self.assertLess(solve(zn, pos_prior_km=100000.0, **sv)["rms_err_km"], 0.6)

    def test_cloud_degrades_gracefully(self):
        ''' Cloud drops obscured sights; the fix degrades smoothly and never
            crashes even with sparse/empty keyframes; cloud=None is unchanged. '''
        from imu_fusion.cloud import CloudSpec
        full = dict(horizon_mode="fused", use_azimuth=True,
                    heading_source="optical", use_parallactic=True,
                    sun_spots=True, imu_anchor=True)
        sv = dict(use_imu=True, use_azimuth=True, use_parallactic=True)

        def mean_rms(cf):
            cloud = None if cf >= 1.0 else CloudSpec(clear_fraction=cf)
            return sum(solve(build_scenario("sea", random.Random(60 + s),
                                            n_shots=14, cloud=cloud, **full),
                             **sv)["rms_err_km"] for s in range(4)) / 4.0
        clear = mean_rms(1.0)
        heavy = mean_rms(0.3)
        self.assertGreater(heavy, clear)              # worse under cloud
        self.assertLess(heavy, 15.0)                  # but still bounded

    def test_coast_budget(self):
        ''' Position dead-reckoning stays sub-km for a minute or two and grows;
            the calibrated gyro coasts no worse than the uncalibrated. '''
        from imu_fusion.visual_anchor import (deadreckon_position_km,
                                              coast_attitude_arcmin)
        a60 = coast_attitude_arcmin(60)
        self.assertLess(deadreckon_position_km(60, 0.5, a60), 0.5)     # <500 m/min
        self.assertGreater(deadreckon_position_km(300, 0.5,
                           coast_attitude_arcmin(300)), 2.0)           # blows up
        self.assertLessEqual(coast_attitude_arcmin(120, calibrated=True),
                             coast_attitude_arcmin(120, calibrated=False))

    def test_starfix_baseline_runs(self):
        ''' The starfix single-fix baseline produces a finite error. '''
        sc = build_scenario("land", random.Random(4), n_shots=4)
        err, lat, lon = starfix_single_fix(sc.keyframes[0], sc.lat0, sc.lon0)
        self.assertIsNotNone(err)
        self.assertGreaterEqual(err, 0.0)

    def test_parse_angle_negative_zero_degree(self):
        ''' A leading "-" on a zero-degree field must make the whole angle
            negative (regression for the Moon-Dec sign bug found by the
            independent ground-truth check: "-00:47.6" is -0.793°, not +0.793°). '''
        from starfix import parse_angle_string
        self.assertAlmostEqual(parse_angle_string("-00:47.6"), -47.6 / 60.0, places=4)
        self.assertAlmostEqual(parse_angle_string("-00:30:00"), -0.5, places=6)
        # Sign of a nonzero degree field is unchanged.
        self.assertAlmostEqual(parse_angle_string("-1:30"), -1.5, places=6)
        self.assertAlmostEqual(parse_angle_string("12:30"), 12.5, places=6)

    def test_subpixel_limb_recovers_synthetic_disk(self):
        ''' The NCC sub-pixel limb fit recovers a known synthetic disk centre
            and radius to well under 0.1 px. '''
        import numpy as np
        from math import erf
        from imu_fusion.disk_metrology import subpixel_limb, plate_scale_arcsec_px
        H, W = 700, 900
        cx0, cy0, R0 = 451.3, 352.8, 180.6
        yy, xx = np.mgrid[0:H, 0:W]
        r = np.hypot(xx - cx0, yy - cy0)
        # bright disk on dark sky with an erf-profile limb (~2 px), mild noise
        verf = np.vectorize(lambda t: 0.5 * (1 - erf(t / 2.0)))
        g = 12 + 200 * verf(r - R0)
        rng = np.random.RandomState(0)
        g = np.clip(g + rng.normal(0, 1.5, g.shape), 0, 255)
        fit = subpixel_limb(g)
        self.assertLess(abs(fit["cx"] - cx0), 0.1)
        self.assertLess(abs(fit["cy"] - cy0), 0.1)
        self.assertLess(abs(fit["R"] - R0), 0.15)     # unbiased to <0.15 px
        self.assertLess(fit["rmse"], 0.5)              # sub-pixel scatter
        self.assertGreater(plate_scale_arcsec_px(fit["R"]), 0.0)

    def test_bright_limb_sigma_is_phase_limited_and_degenerate_at_full(self):
        ''' Bright-limb heading sigma is ~2 deg near half phase, minimal there,
            and diverges toward full Moon; and it is always looser than the Sun's
            sharp-disk heading (so daytime leans on the Sun). '''
        from imu_fusion.optical_attitude import (
            bright_limb_sigma_deg, moon_bright_limb_heading_sigma_deg,
            optical_heading_sigma_deg)
        from imu_fusion.iphone_model import KinematicState
        half = bright_limb_sigma_deg(0.5)
        gibbous = bright_limb_sigma_deg(0.85)
        near_full = bright_limb_sigma_deg(0.99)
        self.assertLess(half, 2.5)                     # ~1.8 deg at half phase
        self.assertGreater(gibbous, half)              # worse toward full
        self.assertGreater(near_full, 8.0)             # diverging near full
        st = KinematicState(ang_rate=0.01, lin_accel=0.05)
        moon_head = moon_bright_limb_heading_sigma_deg(0.5, st)
        sun_head = optical_heading_sigma_deg("Sun", st)
        self.assertGreater(moon_head, sun_head)        # Sun disk is the better heading

    def test_subpixel_limb_seed_robust_on_partial_phase(self):
        ''' The limb fit must survive a DIM, half-lit (quarter-phase) disk, where
            a brightness-threshold seed latches onto a sliver -- the gradient
            sky-limb seed recovers the full disk. '''
        import numpy as np
        from math import erf
        from imu_fusion.disk_metrology import subpixel_limb
        H, W = 700, 900
        cx0, cy0, R0 = 455.0, 350.0, 175.0
        yy, xx = np.mgrid[0:H, 0:W]
        r = np.hypot(xx - cx0, yy - cy0)
        verf = np.vectorize(lambda t: 0.5 * (1 - erf(t / 2.0)))
        disk = 60 * verf(r - R0)                        # dim disk (max L~60)
        disk[xx < cx0] *= 0.03                          # left half in shadow
        rng = np.random.RandomState(1)
        g = np.clip(disk + rng.normal(0, 1.2, disk.shape), 0, 255)
        fit = subpixel_limb(g)
        self.assertLess(abs(fit["R"] - R0), 1.0)       # full disk, not a sliver
        self.assertLess(abs(fit["cx"] - cx0), 1.0)     # TRUE centre from the arc
        self.assertLess(abs(fit["cy"] - cy0), 1.0)
        # The full-circle centre must be the true centre, NOT the lit-blob
        # centroid -- the unlit half is still there even though it is not seen.
        lit = g > 0.5 * g.max()
        ly, lx = np.nonzero(lit)
        lit_centroid_x = lx.mean()
        self.assertGreater(abs(lit_centroid_x - cx0), 30)   # centroid is biased
        self.assertLess(abs(fit["cx"] - cx0),
                        abs(lit_centroid_x - cx0))          # fit is far better

    def test_crater_ncc_recovers_roll_sub_degree(self):
        ''' Rotational NCC of the resolved Moon disk against a reference recovers
            the in-plane rotation (camera roll) to well under a degree -- the
            crater/mare pattern is a far stronger orientation reference than the
            ~2 deg bright limb, and it works at full Moon. '''
        from imu_fusion.lunar_orientation import render_moon, recover_roll, _rotate
        ref, (cx, cy, r) = render_moon(size=221, libration=(3.0, -2.0))
        for true in (6.4, -15.2, 28.0):
            target = _rotate(ref, true, cx, cy)        # rotate real-ish texture
            out = recover_roll(target, ref, (cx, cy), r, coarse=1.0)
            self.assertLess(abs(out["roll"] - true), 0.5)
            self.assertGreater(out["ncc"], 0.95)
        # Works at QUARTER phase too (left half shadowed) -- terminator crater
        # relief is a strong fiducial, unlike the bright limb which fails at full.
        import numpy as np
        yy, xx = np.mgrid[0:221, 0:221]
        halflit = ref.copy(); halflit[xx < cx] *= 0.05
        for true in (7.0, -20.0):
            out = recover_roll(_rotate(halflit, true, cx, cy), halflit,
                               (cx, cy), r, coarse=1.0)
            self.assertLess(abs(out["roll"] - true), 0.6)

    def test_crater_roll_beats_bright_limb_for_parallactic(self):
        ''' The crater-roll parallactic reference (horizon-free) is tighter than
            the bright-limb heading and does not need a horizon. '''
        from imu_fusion.optical_attitude import (parallactic_sigma_deg,
                                                 CRATER_ROLL_SIGMA_DEG,
                                                 SUN_SPOT_ROLL_SIGMA_DEG)
        from imu_fusion.iphone_model import KinematicState
        st = KinematicState(ang_rate=0.01, lin_accel=0.05)
        # no horizon available (horizon sigma huge) but pattern roll rescues it
        s_pat = parallactic_sigma_deg("Moon", st, 600.0, pattern_roll=True)
        s_nohorizon = parallactic_sigma_deg("Moon", st, 600.0, pattern_roll=False)
        self.assertLess(s_pat, s_nohorizon)
        self.assertLess(s_pat, 0.5)                    # sub-degree, horizon-free
        self.assertLess(CRATER_ROLL_SIGMA_DEG, 0.15)   # measured ~0.06 deg
        # Sun with a resolved sunspot pattern also gives a horizon-free roll
        s_sun = parallactic_sigma_deg("Sun", st, 600.0, pattern_roll=True)
        self.assertLess(s_sun, 0.5)
        self.assertLess(SUN_SPOT_ROLL_SIGMA_DEG, 0.15)  # measured ~0.09 deg

    def test_elongation_budget(self):
        ''' Elongation is a strong TIME/longitude observable (dE/dt ~0.5 deg/hr)
            but a negligible DIRECT position line (parallax-only); the two-body
            altitude fix is a few km with good crossing geometry. '''
        from imu_fusion.elongation import position_budget
        b = position_budget("2026-03-24 12:00:00", 51.5, 0.0,
                            sigma_alt_arcmin=2.0, sigma_sep_arcmin=2.0)
        self.assertAlmostEqual(b["alt_sun_km"], 2.0 * 1.852, places=2)
        self.assertTrue(0.3 < abs(b["dE_dt_deg_per_hr"]) < 0.7)   # Moon motion
        self.assertGreater(b["delta_az_deg"], 60)                 # good geometry
        self.assertLess(b["two_lop_fix_km"], 8.0)                 # few-km fix
        self.assertLess(b["elong_longitude_km"], 120)            # time->longitude
        self.assertTrue(b["elong_direct_negligible"])            # weak direct LOP

    def test_groundtruth_matches_independent_ephemeris(self):
        ''' If an independent engine is installed, the almanac's Sun/Moon GHA and
            Dec must agree with it to well under an arc-minute over a time grid.
            Skipped when no engine (astropy/skyfield) is available. '''
        from imu_fusion import validate_ephemeris as V
        if V.ENGINE is None:
            self.skipTest("no independent ephemeris engine installed")
        res = V.compare(V.default_grid(n_days=20, step_hours=12),
                        locations=((51.5, 0.0),))
        for body in ("Sun", "Moon"):
            gha_max = max(abs(x) for x in res[body]["gha_as"])
            dec_max = max(abs(x) for x in res[body]["dec_as"])
            self.assertLess(gha_max, 60.0, f"{body} GHA residual {gha_max:.1f}\"")
            self.assertLess(dec_max, 60.0, f"{body} Dec residual {dec_max:.1f}\"")


try:
    from imu_fusion import corrections as _C
    from imu_fusion.astro import body_gp as _body_gp, body_distance_km as _dist, \
        predicted_altitude as _palt, gp_dec_gha as _gpdg
    _HAVE_CORR = True
except Exception:                               # pragma: no cover
    _HAVE_CORR = False


@unittest.skipUnless(_HAVE_CORR, "starfix / corrections not importable")
class TestCorrections(unittest.TestCase):
    ''' The apparent-altitude corrections (refraction + topocentric parallax)
        that bridge the study's geometric geocentric altitudes to what a phone
        actually measures.  These mirror the Swift port and are validated against
        IAU ERFA / astropy when those are installed. '''

    MEAN_MOON_KM = 384400.0

    def test_parallax_zero_at_zenith_max_at_horizon(self):
        p_zen = _C.parallax_in_altitude_deg(90.0, self.MEAN_MOON_KM)
        p_hor = _C.parallax_in_altitude_deg(0.0, self.MEAN_MOON_KM)
        hp = _C.horizontal_parallax_deg(self.MEAN_MOON_KM)
        self.assertAlmostEqual(p_zen, 0.0, places=6)
        self.assertAlmostEqual(p_hor, hp, places=6)          # horizon parallax == HP
        self.assertTrue(0.9 < hp < 1.0)                      # Moon HP ~57'

    def test_refraction_positive_decreasing_zeroing(self):
        r_hor = _C.refraction_deg(0.5)
        r_mid = _C.refraction_deg(30.0)
        r_zen = _C.refraction_deg(90.0)
        self.assertGreater(r_hor, r_mid)                     # bigger low down
        self.assertGreater(r_mid, r_zen)
        self.assertAlmostEqual(r_zen, 0.0, places=3)         # ~0 at the zenith
        self.assertTrue(0.4 < r_hor < 0.7)                   # ~30-40' near horizon

    def test_roundtrip_apparent_geometric(self):
        for dist in (356500.0, 384400.0, 406700.0):
            for geo in (8.0, 20.0, 45.0, 70.0, 88.0):
                app = _C.apparent_from_geometric(geo, dist)
                back = _C.geometric_from_apparent(app, dist)
                self.assertAlmostEqual(back, geo, places=6,
                                       msg=f"roundtrip d={dist} geo={geo}")

    def test_moon_net_negative_sun_net_small(self):
        # Moon at low altitude: parallax dominates -> apparent well BELOW geometric.
        net_moon = _C.total_correction_deg(15.0, self.MEAN_MOON_KM)
        self.assertLess(net_moon, -0.7)
        # Sun: negligible parallax, small positive refraction lift.
        net_sun = _C.total_correction_deg(15.0, 1.495978707e8)
        self.assertTrue(0.0 < net_sun < 0.1)

    def test_refraction_matches_erfa(self):
        from imu_fusion import validate_ephemeris as V
        rows = V.compare_refraction()
        if not rows:
            self.skipTest("pyerfa not available")
        worst = max(abs(r["resid_am"]) for r in rows)
        self.assertLess(worst, 0.15, f"Bennett vs ERFA max {worst:.3f}'")

    def test_apparent_matches_astropy(self):
        from imu_fusion import validate_ephemeris as V
        if V.ENGINE is None:
            self.skipTest("no independent ephemeris engine")
        app = V.compare_apparent(V.default_grid(n_days=20, step_hours=12))
        if not app or not any(app[b]["app_am"] for b in app):
            self.skipTest("astropy AltAz witness unavailable")
        for b, d in app.items():
            if not d["app_am"]:
                continue
            mx = max(abs(x) for x in d["app_am"])
            self.assertLess(mx, 0.5, f"{b} apparent-alt residual {mx:.3f}'")

    def test_golden_vectors_self_consistent(self):
        ''' The exported goldens must be reproducible by the same formulas the
            Swift port implements (guards the generator and the Swift target). '''
        from imu_fusion.export_golden import build_records
        recs = build_records()
        self.assertGreaterEqual(len(recs), 4)
        for r in recs:
            app = _C.apparent_from_geometric(r["geometricAltDeg"], r["distanceKm"])
            self.assertAlmostEqual(app, r["apparentAltDeg"], places=5)
            back = _C.geometric_from_apparent(r["apparentAltDeg"], r["distanceKm"])
            self.assertAlmostEqual(back, r["geometricAltDeg"], places=5)


try:
    from imu_fusion import stellarium_source as _SS
    _HAVE_SS = True
except Exception:                               # pragma: no cover
    _HAVE_SS = False


@unittest.skipUnless(_HAVE_SS, "stellarium_source not importable")
class TestStellariumSource(unittest.TestCase):
    ''' Ingestion of a Stellarium export as the authoritative astronomical
        source: the tolerant CSV loader, linear interpolation (with RA wrap), the
        engine-free sidereal time, and the guarantee that with NO export present
        the study transparently falls back to the starfix almanac. '''

    _CSV = (
        "#STELLARIUM_EXPORT v1\n"
        "#SCHEMA (Moon @ ...): ra = ...\n"
        "utc,body,ra_deg,dec_deg,dist_au,alt_deg,az_deg,elong_deg,phase,size_arcsec\n"
        "2026-03-24T00:00:00,Moon,359.0,10.0,0.00250,,,,,\n"
        "2026-03-24T01:00:00,Moon,1.0,11.0,0.00260,,,,,\n"   # RA wraps 359->1
        "2026-03-24T00:00:00,Sun,100.0,-5.0,1.0,,,,,\n"
        "2026-03-24T01:00:00,Sun,100.5,-4.8,1.0,,,,,\n"
        "\n#END\n"
    )

    def _write(self, text):
        import tempfile
        fd, path = tempfile.mkstemp(suffix=".csv"); os.close(fd)
        with open(path, "w") as f:
            f.write(text)
        self.addCleanup(os.remove, path)
        return path

    def test_loader_skips_comments_and_optional_blanks(self):
        rows = _SS.load_csv(self._write(self._CSV))
        self.assertEqual(len(rows), 4)                      # 2 Moon + 2 Sun
        moon = [r for r in rows if r["body"] == "moon"]
        self.assertEqual(moon[0]["dist_au"], 0.0025)
        self.assertIsNone(moon[0]["alt_deg"])               # blank optional -> None

    def test_interpolation_and_ra_wrap(self):
        table = _SS.Table(_SS.load_csv(self._write(self._CSV)))
        from datetime import datetime, timezone
        mid = datetime(2026, 3, 24, 0, 30, tzinfo=timezone.utc)
        dec, gha = table.gp_dec_gha("Moon", mid)
        self.assertAlmostEqual(dec, 10.5, places=6)         # 10 -> 11 midpoint
        # RA 359 -> 1 must interpolate through 0 (=> 360), not backwards to 180.
        gast = _SS.gast_deg(mid)
        self.assertAlmostEqual(gha, (gast - 0.0) % 360.0, places=4)
        dist = table.distance_km("Moon", mid)
        self.assertAlmostEqual(dist, 0.00255 * 149_597_870.7, delta=1.0)

    def test_interpolation_clamps_outside_range(self):
        table = _SS.Table(_SS.load_csv(self._write(self._CSV)))
        from datetime import datetime, timezone
        before = datetime(2026, 3, 23, tzinfo=timezone.utc)
        dec, _ = table.gp_dec_gha("Moon", before)
        self.assertAlmostEqual(dec, 10.0, places=6)         # clamped to first row

    def test_gast_matches_astropy_where_available(self):
        from datetime import datetime, timezone
        try:
            from astropy.time import Time
        except Exception:
            self.skipTest("astropy not available")
        import warnings
        worst = 0.0
        for iso in ("2026-03-24 12:00:00", "2026-06-15 03:00:00"):
            dt = datetime.fromisoformat(iso).replace(tzinfo=timezone.utc)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ap = Time(dt.replace(tzinfo=None), scale="utc").sidereal_time(
                    "apparent", "greenwich").deg % 360.0
            d = abs(((_SS.gast_deg(dt) - ap + 180) % 360) - 180) * 3600.0
            worst = max(worst, d)
        self.assertLess(worst, 2.0, f"engine-free GAST vs astropy {worst:.2f}\"")

    def test_no_export_means_none_and_starfix_fallback(self):
        # No stellarium_ephemeris.csv is shipped, so the authoritative table is
        # absent and the study keeps using the starfix almanac unchanged.
        _SS.clear_cache()
        self.assertIsNone(_SS.get_table())
        if _HAVE:
            gp = body_gp("Moon", "2026-03-24 12:00:00")
            self.assertTrue(-29 < gp.get_lat() < 29)        # sane Dec from starfix


@unittest.skipUnless(_HAVE, "starfix / numpy not installed")
class TestMissionPlan(unittest.TestCase):
    ''' Route sky-forecast used for mission planning (Istanbul -> Ankara). '''

    def test_leg_geometry(self):
        from imu_fusion.mission_plan import ISTANBUL_ANKARA as leg
        self.assertTrue(300 < leg.distance_km < 400)          # ~350 km
        a = leg.position_at(0.0)
        b = leg.position_at(1.0)
        self.assertAlmostEqual(a[0], leg.lat_a, places=6)
        self.assertAlmostEqual(b[1], leg.lon_b, places=6)
        mid = leg.position_at(0.5)                            # between the ends
        self.assertTrue(min(leg.lat_a, leg.lat_b) <= mid[0] <= max(leg.lat_a, leg.lat_b))
        self.assertTrue(min(leg.lon_a, leg.lon_b) <= mid[1] <= max(leg.lon_a, leg.lon_b))

    def test_two_lop_sigma_degenerates_at_small_crossing_angle(self):
        from imu_fusion.mission_plan import two_lop_sigma_km
        good = two_lop_sigma_km(90.0, 2.0)
        poor = two_lop_sigma_km(10.0, 2.0)
        self.assertLess(good, poor)                           # 90 deg is best
        self.assertAlmostEqual(good, 2.0 * 1.852 * (2 ** 0.5), places=6)
        self.assertEqual(two_lop_sigma_km(0.0, 2.0), float("inf"))   # parallel LOPs

    def test_forecast_reports_both_bodies_and_corrections(self):
        from datetime import datetime, timezone
        from imu_fusion.mission_plan import ISTANBUL_ANKARA as leg, forecast_leg
        rows = forecast_leg(leg, datetime(2026, 8, 18, 12, tzinfo=timezone.utc),
                            n_points=3)
        self.assertEqual(len(rows), 3)
        for r in rows:
            self.assertTrue(-90 <= r["sun_alt"] <= 90)
            self.assertTrue(0 <= r["sun_az"] < 360)
            self.assertTrue(-180 <= r["daz"] <= 180)
            # the Moon's parallax is the dominant correction when it is up
            if r["moon_alt"] > 5:
                self.assertGreater(r["moon_parallax_deg"], 0.0)

    # --- the planner product: availability gate, windows, brief -------------

    def test_leg_from_names_and_unknown_place(self):
        from imu_fusion.mission_plan import leg_from_names
        leg = leg_from_names("Istanbul", "Ankara")
        self.assertTrue(300 < leg.distance_km < 400)
        with self.assertRaises(KeyError):
            leg_from_names("Istanbul", "Atlantis")

    def test_plan_flags_moon_unavailable_near_full(self):
        ''' 2026-07-28: the Moon is at ~167 deg elongation and never shares the
            daytime sky with the Sun -> the brief must say SUN-ONLY and warn. '''
        from datetime import datetime, timezone
        from imu_fusion.mission_plan import ISTANBUL_ANKARA as leg, plan_leg, brief_text
        b = plan_leg(leg, datetime(2026, 7, 28, tzinfo=timezone.utc), step_min=20)
        self.assertFalse(b.moon_available)
        self.assertEqual(b.fix_mode, "sun-only")
        self.assertEqual(b.windows, [])
        self.assertIsNone(b.recommended)
        self.assertTrue(any("SUN-ONLY" in w for w in b.warnings))
        self.assertIn("MOON AVAILABLE (daytime, with the Sun): NO", brief_text(b))

    def test_plan_finds_window_and_recommends_the_longest(self):
        ''' 2026-08-18: a wide daytime Sun+Moon window over the leg. '''
        from datetime import datetime, timezone
        from imu_fusion.mission_plan import ISTANBUL_ANKARA as leg, plan_leg
        b = plan_leg(leg, datetime(2026, 8, 18, tzinfo=timezone.utc), step_min=20)
        self.assertTrue(b.moon_available)
        self.assertEqual(b.fix_mode, "sun+moon")
        self.assertTrue(b.windows)
        longest = max(b.windows,
                      key=lambda w: (w["end"] - w["start"]).total_seconds())
        self.assertIs(b.recommended, longest)     # longest, NOT best crossing

    def test_find_next_opportunity_after_a_blank_day(self):
        from datetime import datetime, timezone
        from imu_fusion.mission_plan import ISTANBUL_ANKARA as leg, find_next_opportunity
        found = find_next_opportunity(leg, datetime(2026, 7, 28, tzinfo=timezone.utc),
                                      max_days=14, min_window_min=60)
        self.assertTrue(found, "expected usable days within two weeks")
        for day, brief in found:
            self.assertTrue(brief.moon_available)


if __name__ == "__main__":
    unittest.main()
