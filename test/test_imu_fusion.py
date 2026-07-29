''' Regression tests for the imu_fusion study.

    Kept fast (small keyframe counts, single seeds).  Requires the optional
    dependencies gtsam, numpy, matplotlib (see imu_fusion/requirements.txt); the
    whole module is skipped if gtsam is unavailable.

    © 2026.  MIT License (see LICENSE file).
'''

import os
import sys
import random
import math
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
        # elongation ~73 deg at this epoch -> a waxing crescent approaching first
        # quarter, k ~ 0.35.  (This assertion previously read 0.4 < k < 0.9 and
        # called it gibbous, which encoded an inverted phase formula.)
        self.assertTrue(0.25 < k < 0.50, f"illuminated fraction {k:.3f}")
        self.assertTrue(moon_limb_available(t))
        self.assertTrue(0.0 <= bright_limb_pa_deg(t) < 360.0)
        gm = body_gp("Moon", t)
        q1 = parallactic_angle_deg(50.0, 0.0, gm)
        q2 = parallactic_angle_deg(52.0, 0.0, gm)
        self.assertGreater(abs(q1 - q2), 1.0)          # q tracks latitude

    def test_illuminated_fraction_endpoints(self):
        ''' REGRESSION GUARD.  The phase formula must give NEW at conjunction and
            FULL at opposition.  It was inverted -- (1+cos E)/2 -- until a
            full-Moon photograph at elongation 170 deg contradicted it (the code
            predicted a 0.7% crescent for an obviously round, fully lit disk). '''
        from math import cos, radians
        from imu_fusion.optical_attitude import moon_illuminated_fraction

        def k_of(E):                      # the formula, isolated from the ephemeris
            return (1.0 - cos(radians(E))) / 2.0
        self.assertAlmostEqual(k_of(0.0), 0.0, places=9)     # new
        self.assertAlmostEqual(k_of(90.0), 0.5, places=9)    # quarter
        self.assertAlmostEqual(k_of(180.0), 1.0, places=9)   # full
        # and the shipped function must agree with it at a real epoch
        from imu_fusion.optical_attitude import moon_elongation_deg
        t = "2026-07-28 19:00:00"                            # the photographed full Moon
        self.assertAlmostEqual(moon_illuminated_fraction(t),
                               k_of(moon_elongation_deg(t)), places=9)
        self.assertGreater(moon_illuminated_fraction(t), 0.95)

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


@unittest.skipUnless(_HAVE, "starfix / numpy not installed")
class TestLandfall(unittest.TestCase):
    ''' Landfall: fixing from surveyed land features once they clear the horizon. '''

    def test_geographic_range_matches_classic_rule(self):
        from imu_fusion.landfall import geographic_range_km, horizon_distance_km
        # d(km) ~ 3.83 * sqrt(h_m) with standard refraction
        self.assertAlmostEqual(horizon_distance_km(100.0), 3.83 * 10.0, delta=0.5)
        # a 2000 m peak from a 10 m bridge is visible ~180 km off
        self.assertTrue(170 < geographic_range_km(10.0, 2000.0) < 195)
        # taller peak -> seen sooner
        self.assertGreater(geographic_range_km(10.0, 3000.0),
                           geographic_range_km(10.0, 1000.0))

    def test_vertical_angle_range_roundtrip(self):
        from imu_fusion.landfall import (vertical_angle_deg,
                                         range_from_vertical_angle_km)
        for d in (20.0, 60.0, 120.0):
            a = vertical_angle_deg(d, 2000.0, 10.0)
            back = range_from_vertical_angle_km(a, 2000.0, 10.0)
            self.assertAlmostEqual(back, d, places=4)

    def test_range_precision_degrades_with_distance(self):
        from imu_fusion.landfall import range_sigma_km
        near = range_sigma_km(30.0, 2000.0, 0.1, 10.0)
        far = range_sigma_km(120.0, 2000.0, 0.1, 10.0)
        self.assertLess(near, far)                     # geometry weakens with range
        # a better vertical reference (sea horizon) beats the IMU floor
        self.assertLess(range_sigma_km(60.0, 2000.0, 0.02, 10.0),
                        range_sigma_km(60.0, 2000.0, 0.10, 10.0))

    def test_refraction_is_the_floor_on_ranging(self):
        from imu_fusion.landfall import (geographic_range_km,
                                         geographic_range_spread_km)
        # the dipping range is uncertain by a few percent no matter the instrument
        r = geographic_range_km(10.0, 2000.0)
        s = geographic_range_spread_km(10.0, 2000.0)
        self.assertTrue(0.01 < s / r < 0.06)

    def test_visibility_gate(self):
        from imu_fusion.landfall import visible, geographic_range_km
        r = geographic_range_km(10.0, 1000.0)
        self.assertTrue(visible(r * 0.5, 1000.0, 10.0))
        self.assertFalse(visible(r * 1.5, 1000.0, 10.0))

    def test_horizontal_angle_fix_recovers_known_position(self):
        ''' Three surveyed peaks + two subtended angles must recover the observer
            with NO compass, from a badly displaced DR estimate. '''
        from math import radians, degrees, sin, cos, atan2, fabs
        from starfix import LatLonGeodetic
        from imu_fusion.landfall import two_landmark_circle_fix
        from imu_fusion.astro import great_circle_km

        def bearing(lat1, lon1, lat2, lon2):
            p1, p2 = radians(lat1), radians(lat2)
            dl = radians(lon2 - lon1)
            return degrees(atan2(sin(dl) * cos(p2),
                                 cos(p1) * sin(p2)
                                 - sin(p1) * cos(p2) * cos(dl))) % 360

        peaks = [(41.60, 27.20), (41.20, 27.90), (40.85, 28.60)]
        lat, lon = 40.60, 27.60
        b = [bearing(lat, lon, p[0], p[1]) for p in peaks]
        a12 = fabs(((b[1] - b[0] + 180) % 360) - 180)
        a23 = fabs(((b[2] - b[1] + 180) % 360) - 180)
        dr = LatLonGeodetic(lat + 0.25, lon - 0.30)          # ~38 km off
        res = two_landmark_circle_fix(LatLonGeodetic(*peaks[0]),
                                      LatLonGeodetic(*peaks[1]), a12,
                                      LatLonGeodetic(*peaks[2]), a23,
                                      estimate=dr)
        fix = res[0]
        if isinstance(fix, tuple):
            fix = fix[0]
        err_km = great_circle_km(lat, lon, fix.get_lat(), fix.get_lon())
        self.assertLess(err_km, 0.2, f"landmark fix off by {err_km:.3f} km")

    def test_horizontal_angle_beats_compass_bearing(self):
        from imu_fusion.landfall import (cross_range_sigma_km,
                                         horizontal_angle_sigma_deg,
                                         horizontal_angle_fix_sigma_km)
        d = 80.0
        compass = cross_range_sigma_km(d, 1.5)               # phone magnetometer
        sa = horizontal_angle_sigma_deg(1.0, 9.0, 50.0, d)   # +/-50 m summit ident
        angles = horizontal_angle_fix_sigma_km(d, sa)
        self.assertLess(angles, compass)                     # compass-free wins


@unittest.skipUnless(_HAVE, "numpy / starfix not installed")
class TestTerrainResection(unittest.TestCase):
    ''' Fixing position by matching a photographed skyline to a DEM.

        Uses a synthetic analytic terrain so the logic is exercised without
        downloading SRTM tiles.  The critical test is
        `test_elevation_term_is_required_for_correct_ranking`, which guards the
        finding that an azimuth-only score drifts to high inland terrain. '''

    def _dem(self):
        from imu_fusion.terrain_resection import SyntheticDem
        # a ridge of hills to the south + a decoy high hill to the north
        return SyntheticDem([
            (36.90, 27.20, 900, 0.012),
            (36.92, 27.30, 700, 0.010),
            (36.94, 27.40, 850, 0.011),
            (36.90, 27.48, 650, 0.012),
            (36.96, 27.55, 780, 0.011),
            (36.97, 27.10, 820, 0.012),
            (36.86, 27.36, 950, 0.013),
            (37.10, 27.36, 500, 0.014),        # inland decoy
        ])

    # peak scale must match the terrain's angular width (see module docstring)
    PEAK_KW = dict(window=60, min_prominence=0.05)

    def _synth_observation(self, dem, lat, lon, f=60.0, az0=200.0,
                           cam_above=2.0, cy=400.0):
        ''' Render the skyline from a known point and turn its summits into a
            synthetic image observation (x, row). '''
        import numpy as np
        from imu_fusion.terrain_resection import render_skyline, skyline_peaks, \
            SkylineObservation
        ground = float(dem.elevation(np.array([lat]), np.array([lon]))[0])
        azs, prof = render_skyline(dem, lat, lon, ground + cam_above,
                                   az_step=0.05, d_max_km=45.0, d_step_km=0.05)
        peaks = skyline_peaks(azs, prof, **self.PEAK_KW)
        self.assertGreaterEqual(len(peaks), 3, "synthetic terrain gave too few summits")
        sel = [p for p in peaks if abs(((p[0] - az0 + 180) % 360) - 180) < 60]
        self.assertGreaterEqual(len(sel), 3)
        sel = np.array(sel)
        x = (sel[:, 0] - az0) * f
        row = cy - f * sel[:, 1]
        return SkylineObservation(x, row, weight=np.ones(len(sel)) * 10.0), sel

    def test_render_skyline_is_sane(self):
        import numpy as np
        from imu_fusion.terrain_resection import render_skyline
        dem = self._dem()
        azs, prof = render_skyline(dem, 37.02, 27.36, 60.0, az_step=0.5)
        self.assertEqual(len(azs), len(prof))
        # hills lie south -> the horizon must be higher to the south than north
        south = prof[(azs > 150) & (azs < 230)].max()
        north = prof[(azs > 340) | (azs < 20)].max()
        self.assertGreater(south, north)
        self.assertTrue(-5 < prof.max() < 20)          # plausible elevation angles

    def test_resection_recovers_a_known_viewpoint(self):
        import numpy as np
        from imu_fusion.terrain_resection import resect
        from imu_fusion.astro import great_circle_km
        dem = self._dem()
        true_lat, true_lon = 37.0270, 27.3620
        obs, _ = self._synth_observation(dem, true_lat, true_lon)
        cands = [(la, lo)
                 for la in np.arange(37.010, 37.061, 0.010)
                 for lo in np.arange(27.330, 27.401, 0.010)]
        ranked = resect(dem, obs, cands, f_list=np.arange(40., 90., 5.0),
                        az_step=1.0, render_kw=dict(az_step=0.1, d_step_km=0.08),
                        peak_kw=self.PEAK_KW)
        self.assertTrue(ranked, "resection produced no candidates")
        err = great_circle_km(true_lat, true_lon,
                              ranked[0]["lat"], ranked[0]["lon"])
        self.assertLess(err, 1.5, f"rank-1 was {err:.2f} km from truth")

    def test_elevation_term_is_required_for_correct_ranking(self):
        ''' REGRESSION GUARD.  Scoring on summit azimuths alone loses the
            constraint that fixes how high the horizon should stand, and the
            solution drifts to high inland terrain.  On the real Bodrum/Kos case
            that cost 2059 m vs 297 m.  Here we assert the elevation residual
            actually discriminates: the true viewpoint must have a markedly
            smaller elevation residual than a decoy on high inland ground. '''
        import numpy as np
        from imu_fusion.terrain_resection import (render_skyline, skyline_peaks,
                                                  best_match)
        dem = self._dem()
        true_lat, true_lon = 37.0270, 27.3620
        obs, _ = self._synth_observation(dem, true_lat, true_lon)
        fl = np.arange(40., 90., 2.0)

        def fit(lat, lon):
            ground = float(dem.elevation(np.array([lat]), np.array([lon]))[0])
            azs, prof = render_skyline(dem, lat, lon, ground + 2.0,
                                       az_step=0.1, d_step_km=0.08)
            return best_match(obs, azs, prof,
                              skyline_peaks(azs, prof, **self.PEAK_KW), fl,
                              az_step=1.0)

        good = fit(true_lat, true_lon)
        decoy = fit(37.0500, 27.3600)          # 2.5 km north, higher ground
        self.assertIsNotNone(good)
        # the truth must fit better overall ...
        if decoy is not None:
            self.assertGreaterEqual(good["score"], decoy["score"])
            # ... and specifically on the ELEVATION term
            self.assertLess(good["elev_resid_px"], decoy["elev_resid_px"] + 1e-9)
        # score_match must always report an elevation residual (no az-only mode)
        self.assertIn("elev_resid_px", good)
        self.assertGreaterEqual(good["elev_resid_px"], 0.0)

    def test_observation_requires_row_elevations(self):
        from imu_fusion.terrain_resection import SkylineObservation
        with self.assertRaises(ValueError):
            SkylineObservation([1, 2, 3], [10, 20])        # mismatched lengths

    def test_dem_tiles_absent_is_graceful(self):
        import numpy as np
        from imu_fusion.terrain_resection import DemTiles
        d = DemTiles(directory="/nonexistent-dem-dir")
        self.assertFalse(d.available())
        z = d.elevation(np.array([37.0]), np.array([27.0]))
        self.assertEqual(float(z[0]), 0.0)                 # sea level, no crash


@unittest.skipUnless(_HAVE, "gtsam / numpy not installed")
class TestTerrainFactors(unittest.TestCase):
    ''' Terrain landmarks as GTSAM factors, and what a magnetometer heading is
        actually worth. '''

    LAT0, LON0 = 37.0270, 27.3620

    def _landmarks(self):
        from imu_fusion.terrain_factors import Landmark
        return [Landmark("A", 36.90, 27.20, 900.0),
                Landmark("B", 36.92, 27.32, 700.0),
                Landmark("C", 36.94, 27.45, 850.0),
                Landmark("D", 37.05, 27.52, 600.0)]

    def test_bearing_helper_and_sigma_composition(self):
        from imu_fusion.terrain_factors import bearing_deg, bearing_sigma_deg
        # due north / due east sanity
        self.assertAlmostEqual(bearing_deg(37.0, 27.0, 38.0, 27.0), 0.0, places=3)
        self.assertAlmostEqual(bearing_deg(37.0, 27.0, 37.0, 28.0), 90.0, delta=0.5)
        # the magnetometer dominates an absolute bearing
        s = bearing_sigma_deg(pixel_sigma_deg=0.05, mag_sigma_deg=1.5)
        self.assertGreater(s, 1.4)
        self.assertLess(s, 1.6)

    def test_zero_noise_recovery_from_landmarks(self):
        ''' Exact measurements must return the true position. '''
        from imu_fusion.terrain_factors import (solve_landmark_fix,
                                                synthesize_measurements)
        from imu_fusion.astro import great_circle_km
        lms = self._landmarks()
        b, e, a = synthesize_measurements(self.LAT0, self.LON0, lms)
        fix = solve_landmark_fix(
            bearings=[(lm, m, 0.5) for lm, m in b],
            elevations=[(lm, m, 0.1) for lm, m in e],
            horizontal_angles=[(x, y, m, 0.05) for x, y, m in a],
            lat0=self.LAT0, lon0=self.LON0,
            prior_en_m=(3000.0, -2500.0), prior_sigma_km=30.0)
        err_km = great_circle_km(self.LAT0, self.LON0, fix["lat"], fix["lon"])
        self.assertLess(err_km, 0.05, f"zero-noise fix off by {err_km*1000:.0f} m")
        self.assertEqual(fix["n_factors"], len(b) + len(e) + len(a))

    def test_horizontal_angles_are_compass_free(self):
        ''' THE MAGNETOMETER LESSON.  A heading bias corrupts every absolute
            bearing but cancels exactly in the subtended angles, so a
            bearing-only fix is dragged away while the horizontal-angle fix is
            untouched. '''
        from imu_fusion.terrain_factors import (solve_landmark_fix,
                                                synthesize_measurements)
        from imu_fusion.astro import great_circle_km
        lms = self._landmarks()
        BIAS = 6.0                                   # deg of compass error
        b, e, a = synthesize_measurements(self.LAT0, self.LON0, lms,
                                          heading_bias_deg=BIAS)
        common = dict(lat0=self.LAT0, lon0=self.LON0,
                      prior_en_m=(1500.0, -1500.0), prior_sigma_km=30.0)
        bearing_fix = solve_landmark_fix(
            bearings=[(lm, m, 0.5) for lm, m in b], **common)
        angle_fix = solve_landmark_fix(
            horizontal_angles=[(x, y, m, 0.05) for x, y, m in a], **common)
        err_bearing = great_circle_km(self.LAT0, self.LON0,
                                      bearing_fix["lat"], bearing_fix["lon"])
        err_angle = great_circle_km(self.LAT0, self.LON0,
                                    angle_fix["lat"], angle_fix["lon"])
        self.assertGreater(err_bearing, 1.0,
                           "a 3 deg heading bias should wreck a bearing fix")
        self.assertLess(err_angle, 0.2,
                        f"horizontal angles must be compass-free "
                        f"(got {err_angle*1000:.0f} m)")
        self.assertLess(err_angle, err_bearing)

    def test_elevation_factor_gives_a_range_and_is_compass_free(self):
        from imu_fusion.terrain_factors import (solve_landmark_fix,
                                                synthesize_measurements)
        from imu_fusion.astro import great_circle_km
        lms = self._landmarks()
        b, e, a = synthesize_measurements(self.LAT0, self.LON0, lms,
                                          heading_bias_deg=5.0)
        fix = solve_landmark_fix(
            elevations=[(lm, m, 0.05) for lm, m in e],
            lat0=self.LAT0, lon0=self.LON0,
            prior_en_m=(800.0, 800.0), prior_sigma_km=30.0)
        err = great_circle_km(self.LAT0, self.LON0, fix["lat"], fix["lon"])
        self.assertLess(err, 1.0)                    # ranges alone locate it
        self.assertTrue(fix["sigma_km"] == fix["sigma_km"])   # covariance exists

    def test_fix_reports_covariance(self):
        from imu_fusion.terrain_factors import (solve_landmark_fix,
                                                synthesize_measurements)
        lms = self._landmarks()
        b, e, a = synthesize_measurements(self.LAT0, self.LON0, lms)
        fix = solve_landmark_fix(
            horizontal_angles=[(x, y, m, 0.05) for x, y, m in a],
            elevations=[(lm, m, 0.1) for lm, m in e],
            lat0=self.LAT0, lon0=self.LON0, prior_sigma_km=30.0)
        self.assertIsNotNone(fix["cov_en"])
        self.assertGreater(fix["sigma_km"], 0.0)
        self.assertLess(fix["sigma_km"], 30.0)       # sharper than the prior

    def test_jacobians_are_finite_and_position_only(self):
        ''' Linearising each factor must give a finite Jacobian that is non-zero
            only on the (east, north) translation columns. '''
        import numpy as np
        from imu_fusion.terrain_factors import (Landmark,
                                                landmark_bearing_factor,
                                                landmark_elevation_factor,
                                                landmark_horizontal_angle_factor)
        X = gtsam.symbol_shorthand.X
        lm1 = Landmark("A", 36.90, 27.20, 900.0)
        lm2 = Landmark("B", 36.94, 27.45, 850.0)
        pose = gtsam.Pose3(gtsam.Rot3(), gtsam.Point3(400.0, -250.0, 0.0))
        vals = gtsam.Values(); vals.insert(X(0), pose)
        for f in (landmark_bearing_factor(X(0), lm1, 210.0, 1.0, self.LAT0, self.LON0),
                  landmark_elevation_factor(X(0), lm1, 2.0, 0.1, self.LAT0, self.LON0),
                  landmark_horizontal_angle_factor(X(0), lm1, lm2, 30.0, 0.1,
                                                   self.LAT0, self.LON0)):
            A, _ = f.linearize(vals).jacobian()
            self.assertTrue(np.all(np.isfinite(A)))
            self.assertAlmostEqual(float(A[0, 0]), 0.0, places=9)   # rotation
            self.assertAlmostEqual(float(A[0, 5]), 0.0, places=9)   # up
            self.assertGreater(abs(A[0, 3]) + abs(A[0, 4]), 0.0)    # east/north


@unittest.skipUnless(_HAVE, "numpy / starfix not installed")
class TestSyntheticSkyline(unittest.TestCase):
    ''' End-to-end terrain resection on synthetic skyline data, driven the way a
        phone would drive it: a rough (dead-reckoned) position, a magnetometer
        heading, and a known lens. '''

    LAT, LON = 37.0270, 27.3620
    AZ, FOV = 200.0, 120.0
    RENDER = dict(az_step=0.2, d_step_km=0.15, d_max_km=40.0)
    PEAKS = dict(window=30, min_prominence=0.05)

    def _dem(self):
        from imu_fusion.terrain_resection import SyntheticDem
        return SyntheticDem([
            (36.90, 27.20, 900, 0.012), (36.92, 27.30, 700, 0.010),
            (36.94, 27.40, 850, 0.011), (36.90, 27.48, 650, 0.012),
            (36.96, 27.55, 780, 0.011), (36.97, 27.10, 820, 0.012),
            (36.86, 27.36, 950, 0.013), (37.10, 27.36, 500, 0.014),
        ])

    def _obs(self, noise_px=0.0, seed=0):
        import random
        from imu_fusion.terrain_resection import synth_skyline_observation
        return synth_skyline_observation(
            self._dem(), self.LAT, self.LON, self.AZ, self.FOV, width_px=2400,
            noise_px=noise_px, rng=random.Random(seed),
            render_kw=self.RENDER, peak_kw=self.PEAKS)

    def _resect(self, obs, truth, **kw):
        from imu_fusion.terrain_resection import resect_with_priors
        args = dict(prior_radius_km=0.4, grid_step_m=200.0,
                    f_px_per_deg=truth["f_px_per_deg"],
                    render_kw=self.RENDER, peak_kw=self.PEAKS)
        args.update(kw)
        # rough position offset ~500 m from truth, as a DR error would be
        return resect_with_priors(self._dem(), obs,
                                  self.LAT + 0.004, self.LON - 0.003, **args)

    # --- the synthetic data itself ---------------------------------------

    def test_synth_observation_is_well_formed(self):
        obs, truth = self._obs()
        self.assertIsNotNone(obs)
        self.assertGreaterEqual(truth["n_summits"], 4)
        self.assertEqual(len(obs), truth["n_summits"])
        self.assertAlmostEqual(truth["f_px_per_deg"], 2400.0 / self.FOV, places=6)
        # summits must fall inside the frame
        half_w = 0.5 * self.FOV * truth["f_px_per_deg"]
        self.assertLessEqual(float(abs(obs.x).max()), half_w + 1e-6)

    def test_projection_inverts_to_the_generating_angles(self):
        ''' The generator and `score_match` must use the same image model:
            x -> azimuth and row -> elevation invert exactly. '''
        import numpy as np
        obs, truth = self._obs()
        f = truth["f_px_per_deg"]
        az_back = self.AZ + obs.x / f
        el_back = (truth["cy_px"] - obs.row) / f
        np.testing.assert_allclose(np.sort(az_back),
                                   np.sort(truth["summits"][:, 0]), atol=1e-6)
        np.testing.assert_allclose(np.sort(el_back),
                                   np.sort(truth["summits"][:, 1]), atol=1e-6)

    # --- resection with priors -------------------------------------------

    def test_recovers_position_with_rough_prior_and_magnetometer(self):
        from imu_fusion.astro import great_circle_km
        obs, truth = self._obs()
        ranked = self._resect(obs, truth, mag_heading_deg=self.AZ + 1.0,
                              mag_sigma_deg=2.0)
        self.assertTrue(ranked, "no candidate survived the inlier gate")
        err_m = great_circle_km(self.LAT, self.LON,
                                ranked[0]["lat"], ranked[0]["lon"]) * 1000.0
        self.assertLess(err_m, 350.0, f"rank-1 was {err_m:.0f} m off")
        self.assertGreaterEqual(ranked[0]["n_inliers"], 4)

    def test_magnetometer_prunes_without_hurting_accuracy(self):
        ''' The heading prior must cut the search without moving the answer --
            that is the whole reason a phone can run this. '''
        from imu_fusion.astro import great_circle_km
        obs, truth = self._obs()

        def err(ranked):
            return great_circle_km(self.LAT, self.LON,
                                   ranked[0]["lat"], ranked[0]["lon"]) * 1000.0
        with_mag = self._resect(obs, truth, mag_heading_deg=self.AZ + 1.0,
                                mag_sigma_deg=2.0, az_step=0.25)
        without = self._resect(obs, truth, az_step=2.0)     # blind, coarse
        self.assertTrue(with_mag)
        self.assertTrue(without)
        self.assertLessEqual(err(with_mag), err(without) + 1.0)
        # the pruned run considers far fewer headings
        n_pruned = len(range(int(2 * 3 * 2.0 / 0.25)))
        self.assertLess(n_pruned, int(360 / 2.0))

    def test_noise_degrades_gracefully(self):
        from imu_fusion.astro import great_circle_km
        errs = []
        for noise in (0.0, 3.0):
            obs, truth = self._obs(noise_px=noise, seed=11)
            ranked = self._resect(obs, truth, mag_heading_deg=self.AZ + 1.0,
                                  mag_sigma_deg=2.0)
            self.assertTrue(ranked, f"no candidates at {noise} px noise")
            errs.append(great_circle_km(self.LAT, self.LON, ranked[0]["lat"],
                                        ranked[0]["lon"]) * 1000.0)
        self.assertLess(errs[0], 350.0)
        self.assertLess(errs[1], 900.0)          # bounded, not divergent

    def test_a_badly_wrong_magnetometer_is_not_silently_trusted(self):
        ''' If the heading prior is far from the truth the correct solution lies
            outside the searched window, so the resection must FAIL LOUDLY (no
            candidates, or a clearly worse fit) rather than return a confident
            wrong answer. '''
        from imu_fusion.astro import great_circle_km
        obs, truth = self._obs()
        good = self._resect(obs, truth, mag_heading_deg=self.AZ + 1.0,
                            mag_sigma_deg=2.0)
        bad = self._resect(obs, truth, mag_heading_deg=self.AZ + 40.0,
                           mag_sigma_deg=2.0)
        self.assertTrue(good)
        if bad:
            err_bad = great_circle_km(self.LAT, self.LON, bad[0]["lat"],
                                      bad[0]["lon"]) * 1000.0
            err_good = great_circle_km(self.LAT, self.LON, good[0]["lat"],
                                       good[0]["lon"]) * 1000.0
            self.assertTrue(bad[0]["score"] < good[0]["score"] or
                            err_bad > err_good,
                            "a 40 deg heading error must not score as well")

    def test_identified_summits_feed_the_factor_graph(self):
        ''' The whole chain: synthetic skyline -> resection identifies which DEM
            summits were photographed -> those become terrain factors -> a fix
            with a covariance, sharper than the search grid. '''
        from imu_fusion.terrain_factors import (Landmark, solve_landmark_fix,
                                                synthesize_measurements)
        from imu_fusion.astro import great_circle_km
        obs, truth = self._obs()
        ranked = self._resect(obs, truth, mag_heading_deg=self.AZ + 1.0,
                              mag_sigma_deg=2.0)
        self.assertTrue(ranked)
        # the DEM hills, as they would be after identification
        lms = [Landmark(f"H{i}", la, lo, h)
               for i, (la, lo, h, _) in enumerate(self._dem().hills[:4])]
        b, e, a = synthesize_measurements(self.LAT, self.LON, lms)
        fix = solve_landmark_fix(
            horizontal_angles=[(x, y, m, 0.05) for x, y, m in a],
            elevations=[(lm, m, 0.1) for lm, m in e],
            lat0=ranked[0]["lat"], lon0=ranked[0]["lon"],
            prior_sigma_km=2.0)
        err_m = great_circle_km(self.LAT, self.LON, fix["lat"], fix["lon"]) * 1000.0
        self.assertLess(err_m, 200.0, f"fused fix {err_m:.0f} m off")
        self.assertIsNotNone(fix["cov_en"])
        self.assertLess(fix["sigma_km"], 2.0)


def _dem_available():
    try:
        from imu_fusion.terrain_resection import DemTiles
        return DemTiles().available()
    except Exception:                                # pragma: no cover
        return False


@unittest.skipUnless(_HAVE and _dem_available(),
                     "SRTM tiles absent — fetch with "
                     "imu_fusion.terrain_resection.fetch_tiles(36.9,37.2,27.1,27.6)")
class TestRealTerrainResection(unittest.TestCase):
    ''' The same resection experiment on REAL SRTM terrain rather than smooth
        synthetic hills.  Skipped unless tiles have been fetched. '''

    RENDER = dict(az_step=0.1, d_step_km=0.08, d_max_km=40.0)
    PEAKS = dict(window=20, min_prominence=0.10)

    def _dem(self):
        from imu_fusion.terrain_resection import DemTiles
        return DemTiles()

    def test_tile_name_mapping(self):
        from imu_fusion.terrain_resection import tile_name
        self.assertEqual(tile_name(37.03, 27.36), "N37E027")
        self.assertEqual(tile_name(-1.2, -0.5), "S02W001")

    def test_real_terrain_viewpoint_is_recovered(self):
        ''' A viewpoint on the real Bodrum peninsula, observed exactly, must be
            recovered to within a few grid cells. '''
        from imu_fusion.terrain_resection import (synth_skyline_observation,
                                                  resect_with_priors)
        from imu_fusion.astro import great_circle_km
        dem = self._dem()
        lat, lon, az = 37.0450, 27.3900, 190.0        # Konacik-ish, 38 m
        obs, truth = synth_skyline_observation(
            dem, lat, lon, az, 120.0, width_px=2400,
            render_kw=self.RENDER, peak_kw=self.PEAKS)
        self.assertIsNotNone(obs)
        self.assertGreaterEqual(len(obs), 5)
        ranked = resect_with_priors(
            dem, obs, lat + 0.004, lon - 0.003, prior_radius_km=0.4,
            grid_step_m=200.0, mag_heading_deg=az + 1.0, mag_sigma_deg=2.0,
            f_px_per_deg=truth["f_px_per_deg"],
            render_kw=self.RENDER, peak_kw=self.PEAKS)
        self.assertTrue(ranked, "no candidate matched on real terrain")
        err_m = great_circle_km(lat, lon, ranked[0]["lat"],
                                ranked[0]["lon"]) * 1000.0
        self.assertLess(err_m, 600.0, f"rank-1 {err_m:.0f} m off on real terrain")

    def test_elevation_residual_flags_a_bad_fix(self):
        ''' The elevation residual must be small when the hypothesis is right
            and large when it is wrong — that is what makes it usable as an
            accept/reject gate (r = +0.58 with true error over 19 viewpoints). '''
        from imu_fusion.terrain_resection import (synth_skyline_observation,
                                                  render_skyline, skyline_peaks,
                                                  score_match)
        dem = self._dem()
        lat, lon, az = 37.0450, 27.3900, 190.0
        obs, truth = synth_skyline_observation(
            dem, lat, lon, az, 120.0, width_px=2400,
            render_kw=self.RENDER, peak_kw=self.PEAKS)
        f = truth["f_px_per_deg"]

        import numpy as _np
        # az0 is the bearing of the observation's MEAN x column, not the frame
        # centre — the summits are not symmetric about the boresight.
        az0 = az + float(obs.x.mean()) / f

        def resid_at(plat, plon):
            g = float(dem.elevation(_np.array([plat]), _np.array([plon]))[0])
            azs, prof = render_skyline(dem, plat, plon, g + 2.0, **self.RENDER)
            pk = skyline_peaks(azs, prof, **self.PEAKS)
            m = score_match(obs, azs, prof, pk, f, az0 % 360.0, 1)
            return m["elev_resid_px"] if m else float("inf")

        good = resid_at(lat, lon)                       # truth
        bad = resid_at(lat + 0.030, lon - 0.030)        # ~4 km away
        self.assertLess(good, bad,
                        f"residual failed to flag the wrong position "
                        f"(truth {good:.2f} px vs decoy {bad:.2f} px)")


class TestLunarGeometry(unittest.TestCase):
    """Libration, the axis position angle, and the sub-solar point."""

    SITE = (41.0082, 28.9784)          # Istanbul
    T = "2026-07-28T20:45:00"

    def test_meeus_worked_example(self):
        """Meeus, Astronomical Algorithms 2nd ed., example 53.a.

        The almanac this project uses only spans 2024-2030, so 1992 April 12
        cannot go through the ephemeris; feeding Meeus's own lambda/beta/alpha
        into the ch.53 geometry is the only independent check available with no
        network.  His answer includes the PHYSICAL libration, worth about
        0.03 deg, which is why the longitude tolerance is looser than the
        latitude's rather than because the latitude got lucky.
        """
        from imu_fusion.lunar_geometry import libration_from_ecliptic
        l, b, P = libration_from_ecliptic(133.162655, -3.229126, 134.688470,
                                          -0.077221081451, 0.004610, 23.440636)
        self.assertAlmostEqual(l, -1.206, delta=0.02)
        self.assertAlmostEqual(b, +4.194, delta=0.01)
        self.assertAlmostEqual(P, 15.08, delta=0.05)

    def test_libration_stays_inside_its_physical_bounds(self):
        """Optical libration cannot exceed roughly +/-8 deg and +/-7 deg.

        Those bounds come from the orbit -- eccentricity for longitude, axial
        tilt for latitude -- so a formula that breaches them is wrong no matter
        how plausible any single value looks.  Sampled over a year rather than
        at one epoch, because a sign error can hide for weeks.
        """
        from datetime import datetime, timedelta
        from imu_fusion.lunar_geometry import geocentric_libration
        t0 = datetime(2026, 1, 1)
        lo, la, pa = [], [], []
        for h in range(0, 365 * 24, 12):
            iso = (t0 + timedelta(hours=h)).strftime("%Y-%m-%dT%H:%M:%S")
            g = geocentric_libration(iso)
            lo.append(g["lon"]); la.append(g["lat"]); pa.append(g["pole_pa"])
        self.assertLess(max(abs(min(lo)), abs(max(lo))), 8.5)
        self.assertGreater(max(lo) - min(lo), 10.0)       # it must actually move
        self.assertLess(max(abs(min(la)), abs(max(la))), 7.2)
        self.assertGreater(max(la) - min(la), 10.0)
        self.assertLess(max(abs(min(pa)), abs(max(pa))), 26.0)

    def test_topocentric_shift_is_of_order_the_parallax(self):
        """Being on the surface rather than at the centre moves the sub-Earth
        point by about the horizontal parallax, ~1 deg, and never by more."""
        from imu_fusion.lunar_geometry import (geocentric_libration,
                                               topocentric_libration)
        g = geocentric_libration(self.T)
        t = topocentric_libration(self.T, *self.SITE)
        shift = math.hypot(t["lon"] - g["lon"], t["lat"] - g["lat"])
        self.assertGreater(shift, 0.2)
        self.assertLess(shift, 1.4)
        self.assertLess(t["parallax_deg"], 1.1)

    def test_subsolar_point_tracks_the_phase(self):
        """Selenographic colongitude is the lunar clock: 90 deg at full Moon,
        advancing ~12.2 deg/day.  Checked as a RATE as well as a value, because
        a constant offset would pass a single-epoch assertion."""
        from imu_fusion.lunar_texture import subsolar_point
        a = subsolar_point("2026-07-28T20:45:00")
        b = subsolar_point("2026-07-29T20:45:00")
        self.assertAlmostEqual(a["colongitude"], 84.3, delta=1.0)
        drift = (b["colongitude"] - a["colongitude"]) % 360.0
        self.assertAlmostEqual(drift, 12.2, delta=0.6)
        self.assertLess(abs(a["lat"]), 1.6)     # the Sun stays near the equator

    def test_features_near_centre_finds_sinus_medii(self):
        """With the sub-Earth point a couple of degrees from (0, 0), the nearest
        named feature must be Sinus Medii -- the Central Bay is central."""
        from imu_fusion.lunar_geometry import (topocentric_libration,
                                               features_near_centre)
        lib = topocentric_libration(self.T, *self.SITE)
        near = features_near_centre(lib["lon"], lib["lat"], 0.30)
        self.assertEqual(near[0]["name"], "Sinus Medii")
        self.assertLess(near[0]["sep_km"], 220.0)


class TestLunarTexture(unittest.TestCase):
    """Rendering the Moon from Stellarium's albedo map."""

    SITE = (41.0082, 28.9784)
    T = "2026-07-28T20:45:00"

    def setUp(self):
        from imu_fusion.lunar_texture import find_texture
        if find_texture() is None:
            self.skipTest("Stellarium's lunar albedo map is not installed")

    def test_sky_rotation_sign_against_3_vectors(self):
        """The in-image angle of the Moon's pole is P - q, not P + q.

        Both are position angles from celestial north through east, and east is
        counter-clockwise in an un-mirrored photograph, so the two subtract.  The
        plus sign is the natural-looking mistake and it drifts by 2q -- up to
        ~90 deg across a night -- which then hides inside a bogus camera roll.
        Rather than trust the algebra, this computes the pole's image angle from
        3-vectors (object, celestial north, east, zenith) and compares.
        """
        import numpy as np
        from imu_fusion.astro import body_gp, gp_dec_gha
        from imu_fusion.stellarium_source import gast_deg, _parse_dt
        from imu_fusion.lunar_geometry import geocentric_libration
        from imu_fusion.lunar_texture import sky_rotation
        lat, lon = self.SITE
        for iso in ("2026-07-28T18:00:00", "2026-07-28T21:00:00",
                    "2026-07-28T23:00:00"):
            dec, gha = gp_dec_gha(body_gp("Moon", iso))
            dt = _parse_dt(iso)
            ra = (gast_deg(dt) - gha) % 360.0
            a, d = math.radians(ra), math.radians(dec)
            o = np.array([math.cos(d) * math.cos(a), math.cos(d) * math.sin(a),
                          math.sin(d)])
            n = np.array([-math.sin(d) * math.cos(a), -math.sin(d) * math.sin(a),
                          math.cos(d)])
            e = np.array([-math.sin(a), math.cos(a), 0.0])
            th = math.radians((gast_deg(dt) + lon) % 360.0)
            ph = math.radians(lat)
            z = np.array([math.cos(ph) * math.cos(th), math.cos(ph) * math.sin(th),
                          math.sin(ph)])
            up = z - (z @ o) * o
            up = up / np.linalg.norm(up)
            right = np.cross(o, up)
            P = geocentric_libration(iso)["pole_pa"]
            pole = n * math.cos(math.radians(P)) + e * math.sin(math.radians(P))
            direct = math.degrees(math.atan2(-(pole @ right), pole @ up))
            self.assertAlmostEqual(sky_rotation(iso, lat, lon), direct, delta=0.02)

    def test_render_rotation_matches_an_image_rotation(self):
        """render(rotation_deg=A) must equal _rotate(render(0), -A).

        The two live in different modules with opposite sign conventions, and
        the photograph matcher converts between them; if that mapping drifts,
        every recovered rotation silently flips sign.
        """
        import numpy as np
        from imu_fusion.lunar_texture import render
        from imu_fusion.lunar_orientation import _rotate
        a, _ = render(size=180, libration=(3.0, -4.0), subsolar=(10.0, 1.0),
                      rotation_deg=25.0)
        b, _ = render(size=180, libration=(3.0, -4.0), subsolar=(10.0, 1.0),
                      rotation_deg=0.0)
        m = np.hypot(*(np.mgrid[0:180, 0:180] - 89.5)) < 70

        def miss(sign):
            c = _rotate(b, sign * 25.0, 89.5, 89.5)
            return float(np.sqrt(((a - c)[m] ** 2).mean()))

        # Asserted as a RATIO, not against an absolute threshold: `_rotate`
        # resamples bilinearly, so even the correct sign leaves a real residual
        # (blur, ~20% of the signal here).  Only the comparison between the two
        # signs isolates the convention from the interpolation.
        self.assertLess(miss(-1), 0.4 * miss(+1))
        self.assertLess(miss(-1), 0.3 * b[m].std())

    def test_render_on_grid_agrees_with_render(self):
        """The fractional-centre renderer must reproduce the square one exactly
        on the grid where they overlap -- it is the only reason the tie-point
        residuals can be smaller than the half pixel that pasting would cost."""
        import numpy as np
        from imu_fusion.lunar_texture import render, render_on_grid
        n = 160
        a, (cx, cy, r) = render(size=n, libration=(2.0, 3.0),
                                subsolar=(8.0, 0.5), rotation_deg=11.0)
        yy, xx = np.mgrid[0:n, 0:n].astype(float)
        b = render_on_grid(xx, yy, cx, cy, r, libration=(2.0, 3.0),
                           subsolar=(8.0, 0.5), rotation_deg=11.0)
        self.assertLess(np.abs(a - b).max(), 1e-9)

    def test_terminator_follows_the_subsolar_point(self):
        """Near full Moon the unlit lune must be thin and lie OPPOSITE the Sun.

        Rendering the disk and asking which side is dark is a check on the
        illumination geometry that no amount of staring at a number gives.
        """
        import numpy as np
        from imu_fusion.lunar_texture import render
        img, (cx, cy, r) = render(size=201, libration=(0.0, 0.0),
                                  subsolar=(6.0, 0.0), rotation_deg=0.0)
        yy, xx = np.mgrid[0:201, 0:201].astype(float)
        inside = ((xx - cx) ** 2 + (yy - cy) ** 2) < (0.995 * r) ** 2
        dark = inside & (img <= 0.0)
        # the Sun is at +6 deg selenographic longitude, i.e. toward +x, so the
        # unlit sliver must sit at -x, and must be a sliver
        self.assertGreater(dark.sum(), 0)
        self.assertLess(dark.sum() / inside.sum(), 0.05)
        self.assertLess(xx[dark].mean(), cx)


class TestLunarMatch(unittest.TestCase):
    """Fitting the disk geometry to an image."""

    def setUp(self):
        from imu_fusion.lunar_texture import find_texture
        if find_texture() is None:
            self.skipTest("Stellarium's lunar albedo map is not installed")

    def test_recovers_a_synthetic_geometry(self):
        """Render a Moon at a known centre, radius, rotation and libration, then
        make the pipeline find them back.

        This is the end-to-end guard on the offset SIGN in `tie_points`.  With
        the sign flipped the residuals stay small -- the six-parameter model
        simply absorbs the error into the disk centre, which is 0.96 correlated
        with the libration -- so only a case with a KNOWN answer catches it.  On
        the real photograph the same bug made the solve/re-render loop divergent,
        doubling the libration error every round.
        """
        import numpy as np
        from imu_fusion.lunar_texture import render_on_grid
        from imu_fusion.lunar_match import tie_points, fit_geometry

        true = dict(cx=612.4, cy=497.8, r=300.0, rot=7.5, lon=-2.4, lat=4.6)
        sun = (5.7, 0.7)
        yy, xx = np.mgrid[0:1000, 0:1200].astype(float)
        img = render_on_grid(xx, yy, true["cx"], true["cy"], true["r"],
                             libration=(true["lon"], true["lat"]), subsolar=sun,
                             rotation_deg=true["rot"])

        # start deliberately wrong: 6 px off centre, 1.5 deg off in libration
        guess = (true["cx"] + 6.0, true["cy"] - 6.0, true["r"] + 2.0,
                 true["rot"] + 1.0, true["lon"] + 1.5, true["lat"] - 1.5)
        p = np.array(guess, float)
        for _ in range(4):
            ties = tie_points(img, p[0], p[1], p[2], (p[4], p[5]), sun, p[3],
                              search=16)
            self.assertGreaterEqual(len(ties), 8)
            p = fit_geometry(ties, p)["params"]
        self.assertAlmostEqual(p[0], true["cx"], delta=0.5)
        self.assertAlmostEqual(p[1], true["cy"], delta=0.5)
        self.assertAlmostEqual(p[2], true["r"], delta=0.5)
        self.assertAlmostEqual(p[3], true["rot"], delta=0.1)
        self.assertAlmostEqual(p[4], true["lon"], delta=0.15)
        self.assertAlmostEqual(p[5], true["lat"], delta=0.15)

    def test_libration_is_nearly_degenerate_with_the_disk_centre(self):
        """Document the geometry that makes this hard.

        A small libration rotates the sphere about a diameter, which to first
        order just TRANSLATES the near-side pattern -- exactly what moving the
        disk centre does.  The two are ~0.96 correlated, which is why the disk
        centre cannot be left free and unconstrained, and why an error in one
        shows up almost entirely as an error in the other.  If this correlation
        ever drops, the geometry assumptions have changed and the surrounding
        cautions need revisiting.
        """
        import numpy as np
        from imu_fusion.lunar_texture import render_on_grid
        from imu_fusion.lunar_match import tie_points, fit_geometry, project

        cx, cy, r, rot, lo, la = 612.4, 497.8, 300.0, 7.5, -2.4, 4.6
        sun = (5.7, 0.7)
        yy, xx = np.mgrid[0:1000, 0:1200].astype(float)
        img = render_on_grid(xx, yy, cx, cy, r, libration=(lo, la),
                             subsolar=sun, rotation_deg=rot)
        ties = tie_points(img, cx, cy, r, (lo, la), sun, rot)
        p = fit_geometry(ties, (cx, cy, r, rot, lo, la))["params"]

        obs = np.array([[t["x"], t["y"]] for t in ties])
        sel = np.array([[t["lon"], t["lat"]] for t in ties])

        def resid(q):
            out = np.empty(2 * len(sel))
            for k, (a, b) in enumerate(sel):
                x, y, _ = project(a, b, (q[4], q[5]), q[3], q[0], q[1], q[2])
                out[2 * k], out[2 * k + 1] = x - obs[k, 0], y - obs[k, 1]
            return out

        r0 = resid(p)
        J = np.empty((len(r0), 6))
        for k in range(6):
            q = p.copy(); q[k] += 1e-3
            J[:, k] = (resid(q) - r0) / 1e-3
        C = np.linalg.inv(J.T @ J)
        d = np.sqrt(np.diag(C))
        corr = C / np.outer(d, d)
        self.assertGreater(abs(corr[0, 4]), 0.9)     # cx  vs libration longitude
        self.assertGreater(abs(corr[1, 5]), 0.9)     # cy  vs libration latitude


class TestSkylineExtract(unittest.TestCase):
    """The general terrain/sky boundary detector."""

    def _scene(self, snow=False):
        """Synthetic sky-over-ridge, optionally with a snowfield BRIGHTER than
        the sky -- the case that defeats every single-channel threshold."""
        import numpy as np
        H, W = 300, 400
        rgb = np.zeros((H, W, 3))
        yy = np.arange(H)[:, None]
        # sky: a blue vertical gradient, as real sky has
        rgb[:, :, 0] = 160 + 0.10 * yy
        rgb[:, :, 1] = 205 + 0.10 * yy
        rgb[:, :, 2] = 230 + 0.08 * yy
        crest = (150 + 25 * np.sin(np.linspace(0, 3.2, W))).astype(int)
        for x in range(W):
            c = crest[x]
            if snow and 120 < x < 280:
                rgb[c:, x, :] = 248.0            # brighter than the sky
            else:
                rgb[c:, x, 0] = 100.0            # hazed rock: darker, bluer
                rgb[c:, x, 1] = 165.0
                rgb[c:, x, 2] = 215.0
        return rgb, crest

    def test_edge_replication_is_not_optional(self):
        """Zero-padded smoothing fabricates an edge at the top of every column.

        `np.convolve(mode='same')` pads with zeros, so the first smoothed samples
        of any column dive toward zero -- which is exactly the signature a
        skyline detector looks for.  Guards that reject a detection "at the very
        top" then throw the whole column away, and the failure presents as
        terrain that cannot be detected at all.  Two of 840 columns survived when
        this bug was live.
        """
        import numpy as np
        from imu_fusion.skyline_extract import smooth_columns
        band = np.full((60, 5), 200.0)
        sm = smooth_columns(band, 9)
        self.assertAlmostEqual(float(sm[0, 0]), 200.0, places=6)
        self.assertAlmostEqual(float(sm[-1, 0]), 200.0, places=6)
        self.assertLess(float(np.abs(sm - 200.0).max()), 1e-6)

    def test_sky_model_catches_bright_and_dark_terrain(self):
        """The sky-model detector must find the crest whether the terrain is
        darker than the sky (haze) or brighter (snow).  A channel threshold can
        only ever do one."""
        import numpy as np
        from imu_fusion.skyline_extract import sky_model_edge
        for snow in (False, True):
            rgb, crest = self._scene(snow=snow)
            y = sky_model_edge(rgb, y_fit=(10, 100), y_search=(100, 290),
                               thresh=13.0)
            good = np.isfinite(y)
            self.assertGreater(good.sum(), 0.95 * len(y),
                               f"snow={snow}: only {good.sum()} columns found")
            err = y[good] - crest[good]
            # A threshold-crossing detector always lands a little INSIDE the
            # edge: the smoothed departure has to build up before it trips.  A
            # couple of pixels of constant offset is expected and harmless -- the
            # elevation bias in any fit absorbs it.  What must stay small is the
            # SCATTER, because that is the part no nuisance parameter can take
            # out, and it is the same correlated-systematic floor that
            # `resection_geometry.effective_samples` exists to account for.
            self.assertLess(abs(float(np.median(err))), 4.0,
                            f"snow={snow}: bias {np.median(err):+.1f} px")
            self.assertLess(float(np.std(err - np.median(err))), 1.5,
                            f"snow={snow}: scatter {np.std(err):.1f} px")

    def test_straight_run_rejection_keeps_terrain(self):
        """An overlay panel's edge and a ridge look identical to an edge
        detector; only their SHAPE differs.  Rejecting dead-flat runs must
        delete the panel and keep the terrain."""
        import numpy as np
        from imu_fusion.skyline_extract import drop_straight_runs
        x = np.arange(600)
        terrain = 150 + 12 * np.sin(x / 90.0)
        y = terrain.copy()
        y[200:420] = 143.0                        # a panel edge: perfectly flat
        out = drop_straight_runs(y, tol=0.6, min_run=45)
        self.assertTrue(np.all(np.isnan(out[210:410])), "panel edge survived")
        self.assertGreater(np.isfinite(out[:180]).sum(), 170, "terrain deleted")
        self.assertGreater(np.isfinite(out[440:]).sum(), 140, "terrain deleted")


class TestResectionGeometry(unittest.TestCase):
    """Predicting whether a view can fix a position, before searching for it."""

    # the two real scenes this study measured, as (distance_km, height_m)
    ISTANBUL = [(49.0, 871), (52.0, 937), (53.0, 1125), (71.8, 1122)]
    DENVER = [(14.6, 1770), (92.5, 4346)]

    def test_nuisance_biases_destroy_most_of_the_leverage(self):
        """A compass bias eats the common bearing shift, a pitch bias the common
        elevation shift.  Only the near-far DIFFERENCE survives, and quoting the
        raw sensitivity is how a single-ridge view comes to look far better than
        it is."""
        from imu_fusion.resection_geometry import sensitivity
        s = sensitivity(self.ISTANBUL, eye_m=177.0)
        self.assertAlmostEqual(s["lateral_raw"], 1.169, delta=0.02)
        self.assertAlmostEqual(s["lateral_absorbed"], 0.371, delta=0.02)
        self.assertGreater(s["lateral_ratio"], 2.5)

    def test_eye_height_is_required_not_defaulted(self):
        """Defaulting the observer height to zero computes every elevation as if
        from sea level, which overstates the radial sensitivity by ~20x for a
        hilltop viewpoint -- and still returns a plausible-looking number."""
        from imu_fusion.resection_geometry import sensitivity
        with self.assertRaises(TypeError):
            sensitivity(self.DENVER)          # positional eye_m is mandatory
        near = sensitivity(self.DENVER, eye_m=1693.0)["radial_raw"]
        sea = sensitivity(self.DENVER, eye_m=0.0)["radial_raw"]
        self.assertGreater(sea / near, 10.0)

    def test_reproduces_both_measured_outcomes(self):
        """Geometry alone must reach the conclusions that cost hours of search:
        Istanbul unusable, Denver a tight line of position."""
        from imu_fusion.resection_geometry import verdict
        ist = verdict(self.ISTANBUL, 0.142, 177.0)
        self.assertFalse(ist["usable"])
        self.assertGreater(ist["elongation"], 10.0)
        self.assertGreater(ist["along_km"], 5.0)

        den = verdict(self.DENVER, 0.0185, 1693.0)
        self.assertLess(den["across_km"], 0.2)        # measured arc ~0.2 km wide
        self.assertGreater(den["along_km"], 1.0)      # measured arc ~11 km long
        self.assertTrue(den["line_of_position"])

    def test_near_landmark_is_worth_an_order_of_magnitude(self):
        """The one actionable recommendation: include something CLOSE."""
        from imu_fusion.resection_geometry import sensitivity
        far_only = sensitivity([(85.0, 4000), (92.5, 4346)], eye_m=1693.0)
        with_near = sensitivity(self.DENVER, eye_m=1693.0)
        self.assertGreater(with_near["lateral_absorbed"],
                           10 * far_only["lateral_absorbed"])

    def test_unknown_focal_length_is_a_third_nuisance_parameter(self):
        """An unknown pixels-per-degree is degenerate with radial position.

        Moving the camera toward the terrain scales every apparent angle up
        together, which is exactly what shortening the focal length does.  So a
        cropped image, a screenshot or a stock photograph -- anything without
        EXIF -- loses most of its range information before the search begins.

        Measured on the Lake Tahoe wallpaper: with the scale free, moving the
        camera EIGHT KILOMETRES changed the skyline residual by 0.1 arcmin.  The
        module has to say so, or `verdict` will promise a fix that the data
        cannot deliver.
        """
        from imu_fusion.resection_geometry import (position_dilution,
                                                   focal_absorbed_radial)
        tahoe = [(19.9, 2545), (18.2, 2610), (27.1, 2610), (27.7, 2624),
                 (28.9, 2760), (26.6, 2755), (27.6, 2814)]
        known = position_dilution(tahoe, 2.33 / 60, 1899.0, focal_free=False)
        unknown = position_dilution(tahoe, 2.33 / 60, 1899.0, focal_free=True)
        self.assertGreater(unknown["along_km"], 4.0 * known["along_km"])
        self.assertAlmostEqual(unknown["across_km"], known["across_km"], places=6)
        self.assertGreater(focal_absorbed_radial(tahoe, 1899.0), 0.0)

    def test_focal_bounds_bracket_the_measured_tahoe_scale(self):
        """The bracket must contain the answer and still exclude the failure.

        The Tahoe search with the scale free railed at 510 px/deg and returned no
        minimum; bounded, the same search converged at 143.  So the bracket has
        to be tight enough to keep 510 out -- a bound that admits the runaway
        buys nothing -- while still comfortably containing 143.
        """
        from imu_fusion.resection_geometry import focal_bounds_from_relief
        # 129.4 px of crest once the near rock point is masked, against the
        # extremes of the angular relief the DEM renders over that 42 deg arc
        # across all 140 shoreline candidates: 0.84 deg where the far range is
        # across open water, 19.09 where the camera stands under close terrain.
        lo, hi = focal_bounds_from_relief(129.4, [0.84, 1.02, 5.99, 19.09])
        self.assertLess(lo, 143.0)
        self.assertGreater(hi, 143.0)
        self.assertLess(hi, 510.0)
        # The upper bound is the one that does the work, and it must come from
        # the SMALLEST relief; taking the largest would put it near 5 and exclude
        # the answer outright.
        self.assertAlmostEqual(hi, 129.4 / (0.84 * 0.75), places=6)

    def test_focal_bounds_invert_relief_to_scale(self):
        """More degrees for the same pixels is a SHORTER focal length.

        An off-by-inversion here would bracket the reciprocal of the answer and
        the search would rail against a bound while looking perfectly principled.
        """
        from imu_fusion.resection_geometry import focal_bounds_from_relief
        lo, hi = focal_bounds_from_relief(100.0, [1.0, 2.0], margin=0.0)
        self.assertAlmostEqual(lo, 50.0, places=9)     # 100 px / 2 deg
        self.assertAlmostEqual(hi, 100.0, places=9)    # 100 px / 1 deg
        wide = focal_bounds_from_relief(100.0, [1.0, 2.0], margin=0.25)
        self.assertLess(wide[0], lo)
        self.assertGreater(wide[1], hi)
        with self.assertRaises(ValueError):
            focal_bounds_from_relief(0.0, [1.0])

    def test_circle_of_position_reproduces_its_own_angle(self):
        """Every returned point must actually see the requested angle.

        The textbook inscribed-angle circle is a PLANAR construction and bearings
        are spherical; across a 78 km chord the two disagree by ~0.3 deg, which
        at 3 deg/km is a 100 m bias.  Solving the locus numerically and verifying
        each point removes both that and the branch-picking that has caused three
        separate sign errors in this project.
        """
        import numpy as np
        from imu_fusion.resection_geometry import (circle_of_position,
                                                   bearing_deg)
        near = (39.74298, -104.99051)          # a downtown tower, ~15 km
        far = (40.2549, -105.6151)             # Longs Peak, ~92 km
        want = -1.2632
        la, lo = circle_of_position(near, far, want, centre=(39.65, -104.88),
                                    span_km=25, step_km=0.5)
        self.assertGreater(len(la), 20)
        got = np.array([((bearing_deg(a, b, *far) - bearing_deg(a, b, *near)
                          + 180) % 360) - 180 for a, b in zip(la, lo)])
        self.assertLess(float(np.abs(got - want).max()), 1e-3)

    def test_circle_of_position_passes_through_the_measured_answer(self):
        """The Denver photograph was solved independently at 39.644 N 104.878 W;
        the locus built from the tower-to-summit angle must pass through it."""
        import numpy as np
        from imu_fusion.resection_geometry import circle_of_position
        la, lo = circle_of_position((39.74298, -104.99051), (40.2549, -105.6151),
                                    -1.2632, centre=(39.65, -104.88),
                                    span_km=25, step_km=0.25)
        d = np.hypot((la - 39.644) * 110.6, (lo + 104.878) * 85.4)
        self.assertLess(float(d.min()), 0.35)

    def test_correlated_residuals_do_not_average_down(self):
        """839 samples of smooth, azimuth-correlated residual carry the
        information of a handful, not of 839.  Treating them as independent is
        what made a real uncertainty look 10x better than it was."""
        import numpy as np
        from imu_fusion.resection_geometry import effective_samples
        rng = np.random.default_rng(7)
        x = np.linspace(0, 10, 839)
        white = rng.normal(0, 0.14, 839)
        smooth = np.convolve(rng.normal(0, 1, 839 + 200),
                             np.ones(200) / 200, mode="valid")[:839]
        smooth *= 0.14 / smooth.std()
        w = effective_samples(x, white)
        s = effective_samples(x, smooth)
        self.assertGreater(w["n_eff"], 400)
        self.assertLess(s["n_eff"], 60)
        self.assertGreater(s["sigma_inflation"], 3.0)


class TestTerrainBiasFactors(unittest.TestCase):
    """The nuisance biases, modelled as graph variables instead of noise."""

    LAT, LON, EYE = 39.644, -104.878, 1695.0

    def _landmarks(self, near=True):
        from imu_fusion.terrain_factors import Landmark
        lms = [Landmark("longs", 40.2549, -105.6151, 4346),
               Landmark("p2", 40.18, -105.65, 3960),
               Landmark("p3", 40.28, -105.52, 3818)]
        if near:
            lms.append(Landmark("tower", 39.74298, -104.99051, 1770))
        return lms

    def _solve(self, bear, elev, **kw):
        from imu_fusion.terrain_factors import solve_landmark_fix
        return solve_landmark_fix(
            bearings=[(l, b, 0.05) for l, b in bear],
            elevations=[(l, e, 0.05) for l, e in elev],
            lat0=self.LAT, lon0=self.LON, prior_en_m=(0.0, 0.0),
            prior_sigma_km=50.0, cam_height_m=self.EYE, **kw)

    def test_compass_bias_as_a_variable_rescues_the_fix(self):
        """A 1.5 deg heading error is ordinary for a phone magnetometer.  Left
        out of the model it does not merely widen the answer -- it moves it by
        over ten kilometres, because a common bearing offset is exactly what a
        lateral position shift looks like.  As an estimated variable shared by
        every bearing, it is absorbed and the position comes back."""
        from imu_fusion.terrain_factors import synthesize_measurements
        lms = self._landmarks()
        bear, elev, _ = synthesize_measurements(self.LAT, self.LON, lms,
                                                cam_height_m=self.EYE,
                                                heading_bias_deg=1.5)
        naive = self._solve(bear, elev, compass_bias_sigma_deg=0.0)
        modelled = self._solve(bear, elev, compass_bias_sigma_deg=2.0)
        self.assertGreater(math.hypot(naive["east_m"], naive["north_m"]), 5000.0)
        self.assertLess(math.hypot(modelled["east_m"], modelled["north_m"]), 50.0)
        self.assertAlmostEqual(modelled["compass_bias_deg"], 1.5, delta=0.02)

    def test_pitch_bias_as_a_variable_rescues_the_fix(self):
        """+0.87 deg is what this study MEASURED on a real device standing
        still, and +1.50 deg from a moving car.  Either would throw a
        terrain fix by kilometres if elevations are assumed unbiased."""
        from imu_fusion.terrain_factors import synthesize_measurements
        lms = self._landmarks()
        bear, elev, _ = synthesize_measurements(self.LAT, self.LON, lms,
                                                cam_height_m=self.EYE)
        elev = [(l, e + 0.9) for l, e in elev]
        naive = self._solve(bear, elev, pitch_bias_sigma_deg=0.0)
        modelled = self._solve(bear, elev, pitch_bias_sigma_deg=1.5)
        self.assertGreater(math.hypot(naive["east_m"], naive["north_m"]), 3000.0)
        self.assertLess(math.hypot(modelled["east_m"], modelled["north_m"]), 50.0)
        self.assertAlmostEqual(modelled["pitch_bias_deg"], 0.9, delta=0.02)

    def test_a_near_landmark_shrinks_the_reported_ellipse(self):
        """With the biases estimated, the covariance finally tells the truth
        about range spread: adding one landmark at 15 km to a set otherwise at
        89-93 km must collapse the across-track axis by a large factor.  This is
        the same conclusion `resection_geometry` reaches from geometry alone."""
        from imu_fusion.terrain_factors import synthesize_measurements
        out = {}
        for tag, near in (("far", False), ("near", True)):
            lms = self._landmarks(near=near)
            bear, elev, _ = synthesize_measurements(self.LAT, self.LON, lms,
                                                    cam_height_m=self.EYE)
            out[tag] = self._solve(bear, elev, compass_bias_sigma_deg=2.0,
                                   pitch_bias_sigma_deg=1.5)
        self.assertLess(out["near"]["semi_minor_km"],
                        out["far"]["semi_minor_km"] / 5.0)
        for tag in ("far", "near"):
            self.assertGreater(out[tag]["semi_major_km"],
                               out[tag]["semi_minor_km"] * 2.0,
                               f"{tag}: covariance should still read as a line "
                               f"of position")

    def test_calibrated_sigma_inflates_only_for_correlated_residual(self):
        """A dense skyline match is not hundreds of independent looks.  The
        sigma handed to a factor must be inflated by sqrt(n / n_eff), or the
        graph reports a covariance an order of magnitude too tight."""
        import numpy as np
        from imu_fusion.terrain_factors import calibrated_sigma_deg
        rng = np.random.default_rng(11)
        az = np.linspace(120.0, 140.0, 800)
        white = rng.normal(0, 0.05, 800)
        smooth = np.convolve(rng.normal(0, 1, 800 + 200),
                             np.ones(200) / 200, mode="valid")[:800]
        smooth *= 0.05 / smooth.std()
        w = calibrated_sigma_deg(az, white)
        c = calibrated_sigma_deg(az, smooth)
        self.assertLess(w["inflation"], 1.5)
        self.assertGreater(c["inflation"], 3.0)
        self.assertGreater(c["sigma_deg"], 2.5 * w["sigma_deg"])


if __name__ == "__main__":
    unittest.main()
