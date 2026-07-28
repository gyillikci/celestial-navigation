''' Run the full Sun+Moon / IMU factor-graph study and emit results.

    Produces (in imu_fusion/results/):
      * results.json        -- all aggregated numbers
      * fig_main.png        -- RMS error by regime: single-fix vs FG(no IMU) vs FG(IMU)
      * fig_gating.png      -- gated (least-rotation) vs ungated shutter
      * fig_convergence.png -- error vs number of shots
      * fig_trigger.png     -- stillness trace with gated shutter instants
      * fig_ellipse.png     -- FG estimate scatter + 1-sigma covariance ellipse
    and RESULTS.md + dashboard.html summarising them.

    Run:  python -m imu_fusion.run_study
    All randomness is seeded for reproducibility.

    (c) 2026.  MIT License (see LICENSE file).
'''

import os
import json
import statistics as st
import random

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .scenario import build_scenario, EPOCH
from .celestial_factor_graph import solve
from .baseline import starfix_single_fix
from .capture_trigger import PROFILES, simulate_trace, find_shutter_instants, \
    stillness
from .iphone_model import DEFAULT_IMU, DEFAULT_CAM, KinematicState, \
    altitude_sigma_arcmin

REGIMES = ("land", "sea", "air")
REGIME_LABEL = {"land": "Land (stationary)",
                "sea": "Sea (vessel + swell)",
                "air": "Air (aircraft)"}
N_SEEDS = 8
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
COLORS = {"single": "#d1495b", "noimu": "#edae49", "imu": "#00798c",
          "gated": "#00798c", "ungated": "#d1495b"}


def _mean_std(xs):
    xs = [x for x in xs if x is not None]
    if not xs:
        return None, None
    return float(st.mean(xs)), float(st.pstdev(xs)) if len(xs) > 1 else 0.0


# --------------------------------------------------------------------------- #
# Experiments
# --------------------------------------------------------------------------- #

def exp_main(n_shots=12):
    ''' Single-fix (starfix) vs FG(no IMU) vs FG(IMU), gated, per regime. '''
    out = {}
    for r in REGIMES:
        single, noimu, imu = [], [], []
        for s in range(N_SEEDS):
            sc = build_scenario(r, random.Random(1000 + s), n_shots=n_shots)
            fixes = [starfix_single_fix(kf, sc.lat0, sc.lon0)[0]
                     for kf in sc.keyframes]
            fixes = [f for f in fixes if f is not None]
            single.append(float(np.sqrt(np.mean(np.square(fixes)))) if fixes
                          else None)
            noimu.append(solve(sc, use_imu=False)["rms_err_km"])
            imu.append(solve(sc, use_imu=True)["rms_err_km"])
        out[r] = dict(single=_mean_std(single), noimu=_mean_std(noimu),
                      imu=_mean_std(imu))
    return out


def exp_gating(n_shots=12):
    ''' Gated (least-rotation) vs ungated (periodic) shutter, FG(IMU). '''
    out = {}
    for r in REGIMES:
        g, u = [], []
        for s in range(N_SEEDS):
            scg = build_scenario(r, random.Random(2000 + s),
                                 n_shots=n_shots, gated=True)
            scu = build_scenario(r, random.Random(2000 + s),
                                 n_shots=n_shots, gated=False)
            g.append(solve(scg, use_imu=True)["rms_err_km"])
            u.append(solve(scu, use_imu=True)["rms_err_km"])
        out[r] = dict(gated=_mean_std(g), ungated=_mean_std(u))
    return out


def exp_horizon(n_shots=12):
    ''' Horizon reference: IMU gravity vs optical ultrawide vs fused, FG+IMU. '''
    out = {}
    for r in REGIMES:
        row = {}
        for mode in ("imu", "uw", "fused"):
            vals = [solve(build_scenario(r, random.Random(5000 + s),
                                         n_shots=n_shots, horizon_mode=mode),
                          use_imu=True)["rms_err_km"] for s in range(N_SEEDS)]
            row[mode] = _mean_std(vals)
        out[r] = row
    return out


def exp_horizon_budget():
    ''' Horizon-reference sigma (arcmin) per regime at gated motion. '''
    from .ultrawide_horizon import (horizon_reference_sigma_arcmin,
                                    HORIZON_AVAILABLE, dip_arcmin)
    rep = {"land": KinematicState(0.02, 0.03),
           "sea": KinematicState(0.10, 0.20),
           "air": KinematicState(0.05, 0.30)}
    out = {}
    for r in REGIMES:
        st_ = rep[r]
        out[r] = dict(
            imu=horizon_reference_sigma_arcmin("imu", st_, r),
            uw=horizon_reference_sigma_arcmin("uw", st_, r),
            fused=horizon_reference_sigma_arcmin("fused", st_, r),
            dip=dip_arcmin(r), available=HORIZON_AVAILABLE[r])
    return out


def exp_optical(n_shots=12, horizon_mode="imu"):
    ''' Tele-disk orientation observables, on top of the (weak) IMU horizon so
        their contribution is visible: altitude only, + azimuth via magnetometer
        heading, + azimuth via OPTICAL heading, + optical azimuth AND the
        parallactic-angle position line. '''
    methods = {
        "alt": (dict(), dict()),
        "az_mag": (dict(use_azimuth=True, heading_source="mag"),
                   dict(use_azimuth=True)),
        "az_opt": (dict(use_azimuth=True, heading_source="optical"),
                   dict(use_azimuth=True)),
        "optical": (dict(use_azimuth=True, heading_source="optical",
                         use_parallactic=True),
                    dict(use_azimuth=True, use_parallactic=True)),
    }
    out = {}
    for r in REGIMES:
        row = {}
        for name, (sc_kw, sv_kw) in methods.items():
            vals = [solve(build_scenario(r, random.Random(6000 + s),
                                         n_shots=n_shots,
                                         horizon_mode=horizon_mode, **sc_kw),
                          use_imu=True, **sv_kw)["rms_err_km"]
                    for s in range(N_SEEDS)]
            row[name] = _mean_std(vals)
        out[r] = row
    return out


def exp_fullstack(n_shots=12):
    ''' The deployed configuration: best horizon (ultrawide fused) with ALL
        sensors stacked -- vs the horizon alone -- to show the optical disk adds
        on top of a good horizon (and to explain that the larger numbers in the
        optical section are the weak-IMU-horizon isolation baseline, not a
        regression). '''
    out = {}
    for r in REGIMES:
        horizon_only = [solve(build_scenario(r, random.Random(7000 + s),
                                             n_shots=n_shots,
                                             horizon_mode="fused"),
                              use_imu=True)["rms_err_km"] for s in range(N_SEEDS)]
        full = [solve(build_scenario(r, random.Random(7000 + s), n_shots=n_shots,
                                     horizon_mode="fused", use_azimuth=True,
                                     heading_source="optical",
                                     use_parallactic=True, sun_spots=True),
                      use_imu=True, use_azimuth=True,
                      use_parallactic=True)["rms_err_km"]
                for s in range(N_SEEDS)]
        out[r] = dict(horizon_only=_mean_std(horizon_only),
                      full_stack=_mean_std(full))
    return out


# The unified full-fusion configuration: every observable turned on together.
FULL_SC = dict(horizon_mode="fused", use_azimuth=True, heading_source="optical",
               use_parallactic=True, sun_spots=True, gated=True)
FULL_SV = dict(use_imu=True, use_azimuth=True, use_parallactic=True)


def exp_groundtruth(n_days=90, step_hours=6):
    ''' Validate the Sun/Moon ground truth (GHA/Dec) against an independent
        ephemeris engine over a time grid.  Returns {} if no engine is available.
    '''
    from .validate_ephemeris import ENGINE, compare, default_grid
    if ENGINE is None:
        return {"engine": None}
    times = default_grid(n_days, step_hours)
    res = compare(times, locations=((51.5, 0.0),))
    t0 = times[0]
    out = {"engine": ENGINE,
           "days": [(t - t0).total_seconds() / 86400.0 for t in times]}
    for b in ("Sun", "Moon"):
        g, d = res[b]["gha_as"], res[b]["dec_as"]
        out[b] = dict(gha_as=g, dec_as=d,
                      gha_rms=float(np.sqrt(np.mean(np.square(g)))),
                      dec_rms=float(np.sqrt(np.mean(np.square(d)))),
                      gha_max=float(np.max(np.abs(g))),
                      dec_max=float(np.max(np.abs(d))))
    return out


def exp_cloud(clear_fractions=(1.0, 0.85, 0.7, 0.55, 0.4, 0.25), n_shots=14):
    ''' Graceful degradation: fix RMS and shot availability vs cloud cover. '''
    from .cloud import CloudSpec
    full = dict(horizon_mode="fused", use_azimuth=True, heading_source="optical",
                use_parallactic=True, sun_spots=True, imu_anchor=True)
    sv = dict(use_imu=True, use_azimuth=True, use_parallactic=True)
    out = {}
    for r in REGIMES:
        out[r] = {}
        for cf in clear_fractions:
            cloud = None if cf >= 1.0 else CloudSpec(clear_fraction=cf,
                                                     mean_passage_s=25.0)
            rms, avail = [], []
            for s in range(N_SEEDS):
                sc = build_scenario(r, random.Random(9800 + s), n_shots=n_shots,
                                    cloud=cloud, **full)
                nobs = sum(len(kf.observations) for kf in sc.keyframes)
                avail.append(nobs / (n_shots * len(sc.bodies)))
                rms.append(solve(sc, **sv)["rms_err_km"])
            out[r][cf] = dict(rms=_mean_std(rms),
                              avail=float(np.mean(avail)))
    return out


def exp_coast(times=(5, 15, 30, 60, 120, 180, 300)):
    ''' Coast budget while both bodies are clouded: attitude and DR position
        drift vs outage duration, with vs without prior anchor calibration. '''
    from .visual_anchor import coast_curve
    return coast_curve(times)


def exp_anchor_drift(times=(0.5, 1, 2, 5, 10, 30, 60, 120, 300)):
    ''' Attitude error vs tracking time: gyro-only vs accelerometer-aided vs
        visual-anchor-aided, stationary and moving. '''
    from .visual_anchor import attitude_error_curve
    return {"stationary": attitude_error_curve(times, moving=False),
            "moving": attitude_error_curve(times, moving=True)}


def exp_anchor_fix(n_shots=12):
    ''' Fix RMS with/without the visual anchor, when the optical horizon is
        available (fused) and when it is NOT (IMU horizon — a high sight, land,
        or horizon out of frame). '''
    out = {}
    for r in REGIMES:
        out[r] = {}
        for cond, hmode in (("optical horizon", "fused"),
                            ("no optical horizon", "imu")):
            row = {}
            for anchor in (False, True):
                vals = [solve(build_scenario(r, random.Random(9700 + s),
                                             n_shots=n_shots, horizon_mode=hmode,
                                             use_azimuth=True,
                                             heading_source="optical",
                                             use_parallactic=True,
                                             sun_spots=True, imu_anchor=anchor),
                              use_imu=True, use_azimuth=True,
                              use_parallactic=True)["rms_err_km"]
                        for s in range(N_SEEDS)]
                row["anchor" if anchor else "none"] = _mean_std(vals)
            out[r][cond] = row
    return out


def exp_lens(alts=tuple(range(10, 71, 5)), n_shots=12):
    ''' Wide vs ultrawide horizon lens: reference sigma vs body altitude (the
        field-of-view constraint) plus the fix under each lens policy. '''
    from .ultrawide_horizon import horizon_reference_sigma_lens
    k = KinematicState(0.10, 0.20)                   # sea, gated
    curve = {"alt": list(alts), "wide": [], "ultrawide": [], "adaptive": []}
    for a in alts:
        for pol in ("wide", "ultrawide", "adaptive"):
            curve[pol].append(
                horizon_reference_sigma_lens("uw", k, "sea", a, pol))
    fix = {}
    for r in ("sea", "air"):
        fix[r] = {}
        for pol in ("wide", "ultrawide", "adaptive"):
            vals = [solve(build_scenario(r, random.Random(9600 + s),
                                         n_shots=n_shots, horizon_mode="fused",
                                         use_azimuth=True,
                                         heading_source="optical",
                                         use_parallactic=True, sun_spots=True,
                                         horizon_lens=pol),
                          use_imu=True, use_azimuth=True,
                          use_parallactic=True)["rms_err_km"]
                    for s in range(N_SEEDS)]
            fix[r][pol] = _mean_std(vals)
    return {"curve": curve, "fix": fix}


def exp_zoom(zooms=(1, 2, 3), n_shots=12):
    ''' Does a clip-on teleconverter help?  Full-fusion RMS vs zoom, plus the
        error budget that explains the answer. '''
    from dataclasses import replace as _replace
    from .ultrawide_horizon import horizon_reference_sigma_arcmin
    out = {"rms": {}, "budget": {}}
    for r in REGIMES:
        out["rms"][r] = {}
        for z in zooms:
            cam = _replace(DEFAULT_CAM, teleconverter=z)
            vals = [solve(build_scenario(r, random.Random(9500 + s),
                                         n_shots=n_shots, horizon_mode="fused",
                                         use_azimuth=True,
                                         heading_source="optical",
                                         use_parallactic=True, sun_spots=True,
                                         cam=cam),
                          use_imu=True, use_azimuth=True,
                          use_parallactic=True)["rms_err_km"]
                    for s in range(N_SEEDS)]
            out["rms"][r][z] = _mean_std(vals)
    # Altitude error budget at a representative gated sea instant (arc minutes).
    sea_kin = KinematicState(0.10, 0.20)
    out["budget"] = {
        "camera pointing 1x": DEFAULT_CAM.pointing_sigma_arcmin(),
        "camera pointing 3x":
            _replace(DEFAULT_CAM, teleconverter=3).pointing_sigma_arcmin(),
        "optical horizon (sea)":
            horizon_reference_sigma_arcmin("uw", sea_kin, "sea"),
        "IMU horizon (sea)":
            horizon_reference_sigma_arcmin("imu", sea_kin, "sea"),
    }
    return out


def exp_realtime():
    ''' Batch-vs-streaming latency and accuracy parity (see bench.py). '''
    from . import bench
    return bench.run(shot_counts=(2, 6, 10, 16, 22, 30))


def exp_ablation(n_shots=12):
    ''' Leave-one-out from the unified full fusion: how much each observable is
        worth *inside the combined graph*. '''
    variants = {
        "full fusion": (FULL_SC, FULL_SV),
        "- ultrawide horizon": ({**FULL_SC, "horizon_mode": "imu"}, FULL_SV),
        "- IMU link": (FULL_SC, {**FULL_SV, "use_imu": False}),
        "- optical azimuth": ({**FULL_SC, "use_azimuth": False},
                              {**FULL_SV, "use_azimuth": False}),
        "- parallactic line": ({**FULL_SC, "use_parallactic": False},
                               {**FULL_SV, "use_parallactic": False}),
        "- Δq (Sun-Moon)": (FULL_SC, {**FULL_SV, "use_differential": False}),
        "- gating": ({**FULL_SC, "gated": False}, FULL_SV),
        "- Moon (Sun only)": ({**FULL_SC, "bodies": ("Sun",)}, FULL_SV),
    }
    out = {}
    for r in REGIMES:
        row = {}
        for name, (sc_kw, sv_kw) in variants.items():
            vals = [solve(build_scenario(r, random.Random(8000 + s),
                                         n_shots=n_shots, **sc_kw),
                          **sv_kw)["rms_err_km"] for s in range(N_SEEDS)]
            row[name] = _mean_std(vals)
        out[r] = row
    return out


def exp_heading_budget():
    ''' Heading sigma (deg): magnetometer vs optical disk.  The Moon column is
        its realistic phase-limited BRIGHT-LIMB compass (~2 deg), not the
        geometric feature axis -- so it reads far looser than the Sun's disk,
        which is why daytime heading should lean on the Sun. '''
    from .optical_attitude import (optical_heading_sigma_deg,
                                   moon_bright_limb_heading_sigma_deg,
                                   moon_illuminated_fraction)
    from .iphone_model import heading_sigma_arcmin
    k = moon_illuminated_fraction(EPOCH.strftime("%Y-%m-%d %H:%M:%S"))
    rep = {"land": KinematicState(0.02, 0.03),
           "sea": KinematicState(0.10, 0.20),
           "air": KinematicState(0.05, 0.30)}
    out = {}
    for r in REGIMES:
        st_ = rep[r]
        out[r] = dict(mag=heading_sigma_arcmin(st_, DEFAULT_IMU) / 60.0,
                      moon=moon_bright_limb_heading_sigma_deg(k, st_),
                      sun=optical_heading_sigma_deg("Sun", st_))
    return out


def exp_convergence(shot_counts=(4, 8, 12, 16, 20)):
    ''' RMS error vs number of shots, FG(IMU), gated. '''
    out = {r: {} for r in REGIMES}
    for r in REGIMES:
        for n in shot_counts:
            vals = [solve(build_scenario(r, random.Random(3000 + s), n_shots=n),
                          use_imu=True)["rms_err_km"] for s in range(N_SEEDS)]
            out[r][n] = _mean_std(vals)
    return out


def exp_sensor_budget():
    ''' Modelled altitude sigma at representative gated / ungated motion. '''
    imu, cam = DEFAULT_IMU, DEFAULT_CAM
    rows = {}
    for r in REGIMES:
        prof = PROFILES[r]
        t, w, a = simulate_trace(prof, 120.0, imu.sample_rate_hz,
                                 random.Random(7))
        g = find_shutter_instants(t, w, a, 12, 6.0, gated=True)
        u = find_shutter_instants(t, w, a, 12, 6.0, gated=False)
        gs = st.mean(altitude_sigma_arcmin(s, imu, cam) for _, s in g)
        us = st.mean(altitude_sigma_arcmin(s, imu, cam) for _, s in u)
        rows[r] = dict(gated_arcmin=gs, ungated_arcmin=us)
    return rows


def collect_ellipse(regime="land", n_shots=12, seeds=40):
    ''' FG(IMU) horizontal error samples (final keyframe) + mean E,N cov. '''
    de, dn, covs = [], [], []
    for s in range(seeds):
        sc = build_scenario(regime, random.Random(4000 + s), n_shots=n_shots)
        res = solve(sc, use_imu=True)
        last = res["per_kf"][-1]
        # true east/north of the final keyframe in the anchor frame:
        from .astro import latlon_to_enu
        te, tn = latlon_to_enu(last["true_lat"], last["true_lon"],
                               sc.lat0, sc.lon0)
        de.append((last["east"] - te) / 1000.0)
        dn.append((last["north"] - tn) / 1000.0)
        if last["cov_en"] is not None:
            covs.append(last["cov_en"])
    cov = np.mean(covs, axis=0) / 1e6 if covs else None   # km^2
    return np.array(de), np.array(dn), cov


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #

def _grouped_bar(ax, groups, series, title, ylabel):
    labels = list(groups)
    n = len(series)
    width = 0.8 / n
    x = np.arange(len(labels))
    for i, (name, vals, color) in enumerate(series):
        means = [vals[g][0] for g in labels]
        errs = [vals[g][1] for g in labels]
        ax.bar(x + i * width - 0.4 + width / 2, means, width, yerr=errs,
               capsize=3, label=name, color=color)
    ax.set_xticks(x)
    ax.set_xticklabels([REGIME_LABEL[g] for g in labels], fontsize=9)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)


def plot_main(main, path):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    _grouped_bar(ax, REGIMES,
                 [("starfix single-fix (RMS/epoch)",
                   {r: main[r]["single"] for r in REGIMES}, COLORS["single"]),
                  ("factor graph, no IMU",
                   {r: main[r]["noimu"] for r in REGIMES}, COLORS["noimu"]),
                  ("factor graph + IMU",
                   {r: main[r]["imu"] for r in REGIMES}, COLORS["imu"])],
                 "Position error: incumbent single-fix vs factor graph",
                 "RMS position error (km)")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_gating(gating, path):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    _grouped_bar(ax, REGIMES,
                 [("gated (least-rotation shutter)",
                   {r: gating[r]["gated"] for r in REGIMES}, COLORS["gated"]),
                  ("ungated (periodic shutter)",
                   {r: gating[r]["ungated"] for r in REGIMES}, COLORS["ungated"])],
                 "Least-rotation capture trigger (factor graph + IMU)",
                 "RMS position error (km)")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_horizon(hz, path):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    _grouped_bar(ax, REGIMES,
                 [("IMU gravity horizon",
                   {r: hz[r]["imu"] for r in REGIMES}, "#edae49"),
                  ("ultrawide optical horizon",
                   {r: hz[r]["uw"] for r in REGIMES}, "#00798c"),
                  ("fused (IMU + ultrawide)",
                   {r: hz[r]["fused"] for r in REGIMES}, "#3ddc97")],
                 "Optical horizon from the ultrawide camera rescues sea & air",
                 "RMS position error (km)")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_ablation(ab, path):
    ''' Leave-one-out degradation from full fusion, per regime. '''
    order = ["full fusion", "- ultrawide horizon", "- IMU link",
             "- optical azimuth", "- parallactic line", "- Δq (Sun-Moon)",
             "- gating", "- Moon (Sun only)"]
    fig, ax = plt.subplots(figsize=(9, 4.8))
    x = np.arange(len(order))
    w = 0.26
    cols = {"land": "#edae49", "sea": "#00798c", "air": "#3ddc97"}
    for i, r in enumerate(REGIMES):
        means = [ab[r][k][0] for k in order]
        errs = [ab[r][k][1] for k in order]
        ax.bar(x + (i - 1) * w, means, w, yerr=errs, capsize=2,
               label=REGIME_LABEL[r], color=cols[r])
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=25, ha="right", fontsize=8.5)
    ax.set_ylabel("RMS position error (km)")
    ax.set_title("Leave-one-out: each observable's worth inside the unified graph")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_groundtruth(gt, path):
    ''' Sun/Moon GHA & Dec residual vs an independent ephemeris, over time. '''
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.4), sharex=True)
    days = gt["days"]
    cols = {"Sun": "#c8992e", "Moon": "#5b6b86"}
    for b in ("Sun", "Moon"):
        a1.plot(days, gt[b]["gha_as"], color=cols[b], lw=0.9, label=b)
        a2.plot(days, gt[b]["dec_as"], color=cols[b], lw=0.9, label=b)
    for ax, ttl in ((a1, "GHA residual"), (a2, "Declination residual")):
        ax.axhspan(-6, 6, color="#3ddc97", alpha=0.12)
        ax.axhline(0, color="#888", lw=0.6)
        ax.set_xlabel("days from 2026-03-24")
        ax.set_ylabel("residual (arc-seconds)")
        ax.set_title(ttl)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    a1.text(days[0], 6.3, "±0.1′ almanac quantization", fontsize=7,
            color="#0a7d54")
    fig.suptitle(f"Ground truth vs independent ephemeris ({gt['engine']}) — "
                 "both bodies agree to a few arc-seconds", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def exp_elongation_budget():
    ''' Position error budget at the canonical epoch, including the Sun-Moon
        elongation as a TIME (longitude) observable vs a direct position line. '''
    from .elongation import position_budget
    iso = EPOCH.strftime("%Y-%m-%d %H:%M:%S")
    out = {"iso": iso}
    for s in (0.5, 1.0, 2.0, 4.0):
        out[f"sep_{s}"] = position_budget(iso, 51.5, 0.0, sigma_alt_arcmin=2.0,
                                          sigma_sep_arcmin=s)
    out["base"] = out["sep_2.0"]
    return out


def plot_elongation(eb, path):
    ''' Position budget: altitude LOPs and the two-body fix vs the elongation
        time->longitude channel across separation-measurement precision. '''
    base = eb["base"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4))
    # left: per-observable position sigma (km)
    labels = ["Sun alt\nLOP", "Moon alt\nLOP", "two-body\nfix", "elong→\nlongitude"]
    vals = [base["alt_sun_km"], base["alt_moon_km"], base["two_lop_fix_km"],
            base["elong_longitude_km"]]
    cols = ["#c8992e", "#5b6b86", "#3ddc97", "#00b4d8"]
    a1.bar(labels, vals, color=cols)
    for i, v in enumerate(vals):
        a1.text(i, v + 1, f"{v:.0f}", ha="center", fontsize=9)
    a1.set_ylabel("position 1σ (km)")
    a1.set_title(f"Single-epoch budget  (σ_alt=2′, σ_sep=2′, ΔAz={base['delta_az_deg']:.0f}°)")
    a1.text(0.5, 0.92, "elongation direct position line: negligible (parallax-only)",
            transform=a1.transAxes, ha="center", fontsize=7.5, color="#b06a12")
    # right: elongation->longitude vs separation-measurement precision
    seps = [0.5, 1.0, 2.0, 4.0]
    lon = [eb[f"sep_{s}"]["elong_longitude_km"] for s in seps]
    tim = [eb[f"sep_{s}"]["elong_time_s"] for s in seps]
    a2.plot(seps, lon, "o-", color="#00b4d8", label="longitude (km)")
    a2.set_xlabel("Sun–Moon separation measurement σ (arc-min)")
    a2.set_ylabel("longitude 1σ (km)", color="#00b4d8")
    a2b = a2.twinx()
    a2b.plot(seps, tim, "s--", color="#9db0c0", label="clock (s)")
    a2b.set_ylabel("recovered clock 1σ (s)", color="#7a8a99")
    a2.set_title(f"Lunar-distance clock  (dE/dt={base['dE_dt_deg_per_hr']:.2f}°/hr)")
    a2.grid(alpha=0.3)
    fig.suptitle("Position error budget — altitude sights vs Sun–Moon elongation",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_cloud(cd, path):
    ''' Fix RMS (bars) and shot availability (line) vs cloud cover. '''
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    cfs = sorted(next(iter(cd.values())).keys(),
                 key=lambda c: float(c), reverse=True)
    covers = [round((1 - float(c)) * 100) for c in cfs]
    x = np.arange(len(cfs))
    w = 0.26
    cols = {"land": "#edae49", "sea": "#00798c", "air": "#3ddc97"}
    for i, r in enumerate(REGIMES):
        means = [cd[r][c]["rms"][0] for c in cfs]
        errs = [cd[r][c]["rms"][1] for c in cfs]
        ax.bar(x + (i - 1) * w, means, w, yerr=errs, capsize=2,
               color=cols[r], label=REGIME_LABEL[r])
    ax.set_xticks(x)
    ax.set_xticklabels([f"{c}%" for c in covers])
    ax.set_xlabel("cloud cover")
    ax.set_ylabel("RMS position error (km)")
    ax.set_title("Graceful degradation as cloud obscures the bodies")
    ax.legend(loc="upper left")
    ax.grid(axis="y", alpha=0.3)
    ax2 = ax.twinx()
    avail = [np.mean([cd[r][c]["avail"] for r in REGIMES]) * 100 for c in cfs]
    ax2.plot(x, avail, color="#666", marker="o", ls="--", label="shots usable")
    ax2.set_ylabel("shots usable (%)", color="#666")
    ax2.set_ylim(0, 105)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_coast(cc, path):
    ''' Coast budget: attitude and DR position drift vs outage duration. '''
    t = cc["time"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.4))
    a1.plot(t, cc["att_calibrated"], marker="o", color="#00798c",
            label="anchor-calibrated gyro")
    a1.plot(t, cc["att_uncalibrated"], marker="s", color="#d1495b",
            label="uncalibrated gyro")
    a1.set_xlabel("cloud outage (s)")
    a1.set_ylabel("attitude / horizon error (arc min)")
    a1.set_title("Attitude coast")
    a1.grid(alpha=0.3)
    a1.legend(fontsize=8)
    a2.plot(t, cc["pos_calibrated_km"], marker="o", color="#00798c",
            label="anchor-calibrated")
    a2.plot(t, cc["pos_uncalibrated_km"], marker="s", color="#d1495b",
            label="uncalibrated")
    a2.axhline(1.0, color="#666", ls=":", lw=1)
    a2.text(t[0], 1.05, "1 km", fontsize=7, color="#666")
    a2.set_xlabel("cloud outage (s)")
    a2.set_ylabel("dead-reckoned position drift (km)")
    a2.set_title("Position coast (both bodies clouded)")
    a2.grid(alpha=0.3)
    a2.legend(fontsize=8)
    fig.suptitle("Coast budget when the sky is lost", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_anchor_drift(ad, path):
    ''' Attitude error vs tracking time — stationary and moving panels. '''
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5), sharey=True)
    for ax, key, title in ((axes[0], "stationary", "Stationary"),
                           (axes[1], "moving", "Moving (maneuvering)")):
        c = ad[key]
        t = c["time"]
        ax.plot(t, c["gyro_only"], marker="o", color="#d1495b",
                label="gyro only (drifts)")
        ax.plot(t, c["accel_aided"], marker="s", color="#edae49",
                label="accelerometer-aided")
        ax.plot(t, c["anchor_aided"], marker="^", color="#00798c",
                label="visual-anchor-aided")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("tracking time (s, log)")
        ax.set_title(title)
        ax.grid(alpha=0.3, which="both")
    axes[0].set_ylabel("attitude / horizon error (arc min, log)")
    axes[1].legend(fontsize=8, loc="upper left")
    fig.suptitle("A tracked sunspot bounds gyro drift and is immune to motion",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_anchor_fix(af, path):
    ''' Fix RMS with/without the anchor, optical-horizon vs no-optical-horizon. '''
    fig, ax = plt.subplots(figsize=(8.8, 4.6))
    conds = ["optical horizon", "no optical horizon"]
    x = np.arange(len(REGIMES) * len(conds))
    labels, none_v, anch_v = [], [], []
    for r in REGIMES:
        for cond in conds:
            labels.append(f"{REGIME_LABEL[r].split(' ')[0]}\n{cond}")
            none_v.append(af[r][cond]["none"][0])
            anch_v.append(af[r][cond]["anchor"][0])
    w = 0.38
    ax.bar(x - w / 2, none_v, w, color="#d1495b", label="no anchor")
    ax.bar(x + w / 2, anch_v, w, color="#00798c", label="with sunspot anchor")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel("RMS position error (km)")
    ax.set_title("The anchor rescues the fix when the optical horizon is lost")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_lens(ld, path):
    ''' Horizon reference sigma vs body altitude for the wide vs ultrawide lens,
        with the field-of-view zones marked. '''
    c = ld["curve"]
    alts = c["alt"]
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.plot(alts, c["wide"], marker="o", color="#00798c",
            label="wide (main) lens — sharper, narrow field")
    ax.plot(alts, c["ultrawide"], marker="s", color="#edae49",
            label="ultrawide lens — wider field, softer")
    ax.plot(alts, c["adaptive"], color="#3ddc97", lw=3, alpha=0.6,
            label="best available (adaptive)")
    ax.set_yscale("log")
    ax.axvspan(0, 30, color="#00798c", alpha=0.07)
    ax.axvspan(30, 52, color="#edae49", alpha=0.08)
    ax.axvspan(52, 72, color="#d1495b", alpha=0.08)
    ax.text(15, ax.get_ylim()[1] * 0.5, "wide\nzone", ha="center", fontsize=8,
            color="#0a6")
    ax.text(41, ax.get_ylim()[1] * 0.5, "ultrawide\nzone", ha="center",
            fontsize=8, color="#a8791f")
    ax.text(62, ax.get_ylim()[1] * 0.5, "IMU-only\n(no horizon)", ha="center",
            fontsize=8, color="#a03040")
    for alt, name in ((32, "Moon"), (40, "Sun")):
        ax.axvline(alt, ls=":", color="#666", lw=1)
        ax.text(alt, ax.get_ylim()[0] * 1.3, name, rotation=90, fontsize=7,
                va="bottom", color="#666")
    ax.set_xlabel("body altitude (degrees)")
    ax.set_ylabel("horizon reference σ (arc minutes, log)")
    ax.set_title("Which horizon lens? It depends on the body's altitude")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_zoom(zd, path):
    ''' Left: full-fusion RMS vs teleconverter (flat). Right: altitude error
        budget (log) — the camera is already ~100x below the horizon, so more
        zoom is optically wasted. '''
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.4))
    zooms = sorted(next(iter(zd["rms"].values())).keys(),
                   key=lambda z: int(z) if isinstance(z, str) else z)
    zi = [int(z) if isinstance(z, str) else z for z in zooms]
    cols = {"land": "#edae49", "sea": "#00798c", "air": "#3ddc97"}
    for r in REGIMES:
        means = [zd["rms"][r][z][0] for z in zooms]
        errs = [zd["rms"][r][z][1] for z in zooms]
        a1.errorbar(zi, means, yerr=errs, marker="o", capsize=3,
                    color=cols[r], label=REGIME_LABEL[r])
    a1.set_xticks(zi)
    a1.set_xticklabels([f"{z}×" for z in zi])
    a1.set_xlabel("external teleconverter")
    a1.set_ylabel("RMS position error (km)")
    a1.set_ylim(bottom=0)
    a1.set_title("More zoom does not move the fix")
    a1.grid(alpha=0.3)
    a1.legend(fontsize=8)

    b = zd["budget"]
    labels = ["camera\npointing 3×", "camera\npointing 1×",
              "optical horizon\n(sea)", "IMU horizon\n(sea)"]
    keys = ["camera pointing 3x", "camera pointing 1x",
            "optical horizon (sea)", "IMU horizon (sea)"]
    vals = [b[k] for k in keys]
    barcols = ["#9bbfc7", "#00798c", "#3ddc97", "#d1495b"]
    a2.barh(labels, vals, color=barcols)
    a2.set_xscale("log")
    a2.set_xlabel("altitude error contribution (arc minutes, log)")
    a2.set_title("The camera is ~100× below the horizon")
    for i, v in enumerate(vals):
        a2.text(v * 1.1, i, f"{v:.2f}′", va="center", fontsize=8)
    a2.grid(axis="x", alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_realtime(rt, path):
    ''' Per-fix latency vs trip length: batch re-solve (grows) vs streaming
        fixed-lag update (flat), averaged over regimes. '''
    ns = rt["shot_counts"]

    def _get(d, n):
        return d[n] if n in d else d[str(n)]     # JSON stringifies int keys

    def avg(method):
        means, stds = [], []
        for n in ns:
            vals = [_get(rt["latency"][r][method], n)[0] for r in REGIMES]
            means.append(np.mean(vals))
            stds.append(np.std(vals))
        return np.array(means), np.array(stds)

    bm, bs = avg("batch")
    sm, ss = avg("stream")
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.errorbar(ns, bm, yerr=bs, marker="o", capsize=3, color="#d1495b",
                label="batch re-solve (whole trajectory)")
    ax.errorbar(ns, sm, yerr=ss, marker="s", capsize=3, color="#00798c",
                label="streaming fixed-lag update")
    ax.set_yscale("log")
    ax.set_xlabel("number of Sun+Moon shots so far")
    ax.set_ylabel("per-fix latency (ms, log scale)")
    ax.axhspan(0, 50, color="#3ddc97", alpha=0.12)
    ax.text(ns[0], 52, "≤ 50 ms real-time budget", fontsize=8, color="#0a7d54")
    ax.set_title("Streaming keeps per-fix latency flat as the voyage grows")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_factorgraph(path):
    ''' Schematic of the unified factor graph (3 keyframes). '''
    from matplotlib.patches import Circle, FancyBboxPatch
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(-1.5, 13.5)
    ax.set_ylim(-4.6, 4.2)
    ax.axis("off")
    acc = "#00798c"
    ink = "#23303a"

    def var(xy, label, r=0.36, color="#e8f1f3"):
        ax.add_patch(Circle(xy, r, fc=color, ec=ink, lw=1.6, zorder=3))
        ax.text(xy[0], xy[1], label, ha="center", va="center",
                fontsize=10, zorder=4)

    def fac(xy, s=0.16, color=acc):
        ax.add_patch(FancyBboxPatch((xy[0] - s, xy[1] - s), 2 * s, 2 * s,
                     boxstyle="square,pad=0", fc=color, ec=ink, lw=1.1, zorder=3))

    def link(a, b, color=ink, lw=1.0, ls="-"):
        ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=lw, ls=ls, zorder=1)

    xs = [1.0, 6.0, 11.0]
    X = [(x, 0.0) for x in xs]
    V = [(x, 2.4) for x in xs]
    Bpos = (6.0, -3.9)
    var(Bpos, "B", r=0.34, color="#f2e2c2")
    ax.text(Bpos[0], Bpos[1] - 0.62, "IMU bias", ha="center", fontsize=8,
            color=ink)

    cel = ["alt☉", "alt☾", "az☉", "az☾", "q☉", "q☾"]
    for i, (xp, vp) in enumerate(zip(X, V)):
        var(vp, f"V{i}", r=0.32, color="#e7efe8")
        var(xp, f"X{i}", r=0.4)
        # velocity <-> pose
        link(xp, vp, color="#9bb", lw=0.8)
        # prior / coarse GPS + attitude reference feeding this keyframe
        fp = (xp[0] - 1.15, xp[1] + 0.9)
        fac(fp, color="#d1495b")
        link(fp, xp, color="#d1495b")
        ax.text(fp[0] - 0.05, fp[1] + 0.32, "prior", ha="center", fontsize=7.5,
                color="#a03040")
        # celestial measurement factors fanned below the pose
        for j, name in enumerate(cel):
            fx = xp[0] - 1.5 + j * 0.6
            fpos = (fx, -1.9)
            fac(fpos, s=0.14)
            link(fpos, xp, color=acc, lw=0.8)
            ax.text(fx, -2.28, name, ha="center", fontsize=7.0, color=acc)
        # horizon-free DIFFERENTIAL Sun-Moon parallactic factor (roll cancels)
        dp = (xp[0] + 1.12, xp[1] + 0.9)
        fac(dp, color="#3ddc97")
        link(dp, xp, color="#2fae79")
        ax.text(dp[0] + 0.02, dp[1] + 0.32, "Δq ☉−☾", ha="center",
                fontsize=6.8, color="#0a7d54")
    # IMU factors between consecutive keyframes
    for i in range(len(X) - 1):
        mid = ((xs[i] + xs[i + 1]) / 2, 1.2)
        fac(mid, s=0.2, color="#edae49")
        for pt in (X[i], V[i], X[i + 1], V[i + 1], Bpos):
            link(mid, pt, color="#c9962f", lw=0.9)
        ax.text(mid[0], mid[1] + 0.34, "IMU", ha="center", fontsize=8,
                color="#a8791f")

    ax.text(6.0, 3.7, "Unified factor graph — one keyframe per Sun+Moon shot",
            ha="center", fontsize=12, weight="bold", color=ink)
    ax.text(6.0, -2.95,
            "alt = altitude line   az = azimuth line   q = per-body parallactic "
            "(needs the horizon)   (☉ Sun, ☾ Moon)\n"
            "Δq ☉−☾ = differential Sun−Moon parallactic — HORIZON-FREE (the shared "
            "platform roll cancels between the two resolved disks).\n"
            "Each celestial factor's covariance is built from the fused attitude "
            "references: IMU gravity + ultrawide horizon + tele-disk orientation.",
            ha="center", fontsize=8.5, color=ink)
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_optical(opt, path):
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    _grouped_bar(ax, REGIMES,
                 [("altitude only",
                   {r: opt[r]["alt"] for r in REGIMES}, "#d1495b"),
                  ("+ azimuth (magnetometer)",
                   {r: opt[r]["az_mag"] for r in REGIMES}, "#edae49"),
                  ("+ azimuth (optical heading)",
                   {r: opt[r]["az_opt"] for r in REGIMES}, "#00798c"),
                  ("+ optical azimuth & parallactic line",
                   {r: opt[r]["optical"] for r in REGIMES}, "#3ddc97")],
                 "Tele-disk orientation: magnetometer-free heading + "
                 "parallactic line (IMU horizon)",
                 "RMS position error (km)")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_convergence(conv, path):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for r in REGIMES:
        ns = sorted(conv[r])
        means = [conv[r][n][0] for n in ns]
        errs = [conv[r][n][1] for n in ns]
        ax.errorbar(ns, means, yerr=errs, marker="o", capsize=3,
                    label=REGIME_LABEL[r])
    ax.set_xlabel("number of Sun+Moon shots fused")
    ax.set_ylabel("RMS position error (km)")
    ax.set_title("Error decreases as more gated shots are fused (FG + IMU)")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_trigger(path, regime="sea"):
    prof = PROFILES[regime]
    t, w, a = simulate_trace(prof, 60.0, DEFAULT_IMU.sample_rate_hz,
                             random.Random(11))
    cost = [stillness(wi, ai) for wi, ai in zip(w, a)]
    g = find_shutter_instants(t, w, a, 6, 6.0, gated=True)
    u = find_shutter_instants(t, w, a, 6, 6.0, gated=False)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(t, cost, color="#333", lw=0.8, label="stillness cost |ω|+|a|/g")
    for tt, _ in g:
        ax.axvline(tt, color=COLORS["gated"], lw=1.4, alpha=0.9)
    ax.scatter([tt for tt, _ in g], [stillness(s.ang_rate, s.lin_accel)
               for _, s in g], color=COLORS["gated"], zorder=5,
               label="gated shutter (calm minima)")
    ax.scatter([tt for tt, _ in u], [stillness(s.ang_rate, s.lin_accel)
               for _, s in u], color=COLORS["ungated"], marker="x", zorder=5,
               label="ungated shutter (periodic)")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("stillness cost (lower = calmer)")
    ax.set_title(f"Least-rotation shutter picks calm instants ({regime})")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_ellipse(path, regime="land"):
    de, dn, cov = collect_ellipse(regime)
    fig, ax = plt.subplots(figsize=(5.6, 5.4))
    ax.scatter(de, dn, s=18, color=COLORS["imu"], alpha=0.7,
               label="FG+IMU estimates")
    ax.scatter([0], [0], color="k", marker="*", s=160, label="truth")
    if cov is not None:
        vals, vecs = np.linalg.eigh(cov)
        ang = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
        from matplotlib.patches import Ellipse
        for k, alpha in ((1, 0.35), (2, 0.18)):
            w, h = 2 * k * np.sqrt(np.maximum(vals, 0))
            ax.add_patch(Ellipse((0, 0), w, h, angle=ang, fill=True,
                                  color=COLORS["imu"], alpha=alpha,
                                  label=f"{k}σ covariance" if k == 1 else None))
    ax.set_aspect("equal")
    ax.set_xlabel("East error (km)")
    ax.set_ylabel("North error (km)")
    ax.set_title(f"Factor-graph fix + covariance ({regime})")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Report + dashboard
# --------------------------------------------------------------------------- #

def write_results_md(data, path):
    def cell(ms):
        return "n/a" if ms is None or ms[0] is None else f"{ms[0]:.1f} ± {ms[1]:.1f}"
    m, g, c, b = (data["main"], data["gating"], data["convergence"],
                  data["sensor_budget"])
    L = []
    L.append("# Results — iPhone 17 Pro Sun+Moon daytime sighting with "
             "IMU + factor graph\n")
    L.append(f"*Epoch {EPOCH:%Y-%m-%d %H:%M} UTC, observer near Greenwich "
             f"(51.5°N, 0°). Sun+Moon ~94° apart in azimuth. "
             f"{N_SEEDS} seeds per cell; values are mean ± std of RMS position "
             f"error in km.*\n")
    # Improvement factors for the narrative.
    fac = {r: (m[r]["single"][0] / m[r]["imu"][0]) for r in REGIMES}
    L.append("## Key findings\n")
    L.append(f"1. **Fusing many Sun+Moon shots with IMU cuts error "
             f"{min(fac.values()):.0f}–{max(fac.values()):.0f}×** versus a "
             f"single-epoch two-body fix — from tens/hundreds of km down to "
             f"{m['land']['imu'][0]:.0f} km (land), {m['air']['imu'][0]:.0f} km "
             f"(air) and {m['sea']['imu'][0]:.0f} km (sea).\n")
    L.append("2. **A single hand-held phone sight is a weak instrument.** The "
             "synthetic horizon's tilt error is ~6′ braced on land but grows to "
             "tens of arc-minutes at sea/air (vs ~1–2′ for a marine sextant), so "
             "the whole value is in fusing many shots, not in any one shot.\n")
    L.append("3. **The least-rotation shutter clearly helps on land and in the "
             "air** (≈2× lower error) by cutting the per-shot horizon noise "
             "3–6×. **At sea it is a wash** — even the calmest swell instant is "
             "too tilted, so the IMU gravity horizon alone is not enough there.\n")
    hz = data["horizon"]
    L.append(f"4. **The ultrawide camera fixes the sea (and air) problem.** "
             f"Shooting the ultrawide horizon at the same instant as the tele "
             f"body gives an *optical* horizon that is immune to acceleration. "
             f"It drops the sea fix from {hz['sea']['imu'][0]:.0f} km to "
             f"**{hz['sea']['fused'][0]:.1f} km** and the air fix from "
             f"{hz['air']['imu'][0]:.0f} km to **{hz['air']['fused'][0]:.1f} "
             f"km** — bringing the moving platforms to land-class accuracy. "
             f"On land there is no true sea horizon, so it falls back to the "
             f"IMU (no change).\n")
    hdb = data["heading_budget"]
    op = data["optical"]
    L.append(f"5. **The tele lens is more than a pointer.** Resolving the disk "
             f"gives a magnetometer-free heading, but the two bodies are not "
             f"equal: the **Sun's** sharp disk yields ~{hdb['sea']['sun']:.1f}° "
             f"(vs ~{hdb['sea']['mag']:.0f}° for the phone compass), while the "
             f"**Moon's bright limb** is only ~{hdb['sea']['moon']:.1f}° "
             f"(phase-limited, and degenerate near full) — *looser than the "
             f"magnetometer*. So in daytime the heading should come from the "
             f"Sun; the Moon's limb is a night / Sun-occluded backup. And the "
             f"**difference** of the two disks' orientations gives a genuinely "
             f"horizon-free position line (Δq, roll cancels); on a weak (IMU) "
             f"horizon at sea the optical stack cuts the fix from "
             f"{op['sea']['alt'][0]:.0f} km to {op['sea']['optical'][0]:.0f} km.\n")
    L.append("6. **Geometry matters:** the fix is well-conditioned only when "
             "the Sun and Moon are well separated in azimuth (~90° here, a "
             "first-quarter Moon); near-parallel lines of position degrade it.\n")
    # ------- ground-truth validation -------
    gt = data.get("groundtruth", {})
    if gt.get("engine"):
        L.append("## Ground-truth validation — are the Sun/Moon positions right?\n")
        L.append(f"The whole study rests on the Sun/Moon positions from "
                 f"`starfix`'s almanac. Those are cross-checked here against an "
                 f"**independent** ephemeris (`{gt['engine']}` — a separate "
                 f"implementation from the almanac's Skyfield/JPL source), over a "
                 f"90-day grid. Residuals (arc-seconds):\n")
        L.append("| Body | GHA rms / max | Dec rms / max | in km |")
        L.append("|---|---|---|---|")
        for _bd in ("Sun", "Moon"):
            gg = gt[_bd]
            km = max(gg["gha_max"], gg["dec_max"]) / 3600 * 111.2 * 0.6
            L.append(f"| {_bd} | {gg['gha_rms']:.1f}″ / {gg['gha_max']:.1f}″ | "
                     f"{gg['dec_rms']:.1f}″ / {gg['dec_max']:.1f}″ | ~{km:.2f} km |")
        L.append("\n![groundtruth](results/fig_groundtruth.png)\n")
        L.append("Both bodies agree to a **few arc-seconds** — far below the "
                 "study's arc-minute measurement noise — confirming the ground "
                 "truth. The residual floor is the almanac's 0.1′ table "
                 "quantization plus hourly linear interpolation.\n")
        L.append("> **This validation earned its keep.** The independent check "
                 "first flagged a Moon-declination error of up to ~1.6° near the "
                 "Dec≈0 crossings — a sign-handling bug in `starfix."
                 "parse_angle_string` for the almanac's `-00:MM.M` "
                 "negative-zero-degree format (`float(\"-00\")` is `-0.0`, and "
                 "`-0.0 < 0` is `False`). It is now fixed; the Moon residual "
                 "dropped from ~5700″ to a few arc-seconds. The study's canonical "
                 "epoch (Moon at Dec +28°) was never affected. A **Stellarium** "
                 "export can be added as a third witness — see "
                 "`stellarium_reference.md`.\n")

    # ------- elongation / position error budget -------
    eb = data.get("elongation")
    if eb:
        base = eb["base"]
        L.append("## Position error budget — and the Sun–Moon elongation\n")
        L.append("Where does the fix error come from, and what does measuring "
                 "the **Sun–Moon angular separation** (elongation) add? Single "
                 "epoch, σ = 2′ per sight:\n")
        L.append("| Observable | position 1σ | role |")
        L.append("|---|---|---|")
        L.append(f"| Sun altitude LOP | {base['alt_sun_km']:.1f} km | position "
                 f"line (1′ = 1 nmi) |")
        L.append(f"| Moon altitude LOP | {base['alt_moon_km']:.1f} km | position "
                 f"line |")
        L.append(f"| **Two-body fix** | **{base['two_lop_fix_km']:.1f} km** | LOPs "
                 f"cross at ΔAz={base['delta_az_deg']:.0f}° (good geometry) |")
        L.append(f"| Elongation → **time** | {base['elong_longitude_km']:.0f} km | "
                 f"dE/dt={base['dE_dt_deg_per_hr']:.2f}°/hr → clock to "
                 f"{base['elong_time_s']:.0f} s → longitude |")
        L.append(f"| Elongation → direct pos | negligible | parallax-only, "
                 f"observer-independent |")
        L.append("\n![elongation budget](results/fig_elongation.png)\n")
        L.append("**Reading it.** The two altitude sights are the workhorses — a "
                 f"~{base['two_lop_fix_km']:.0f} km single-epoch fix, driven down "
                 "to the headline few-km numbers by fusing many gated shots. The "
                 "**elongation is a poor _direct_ position line** (it barely "
                 "changes with where you stand — only lunar parallax moves it), "
                 "**but a strong _time_ observable**: the Moon slides ~0.5°/hr "
                 "against the Sun, so measuring the separation to a few arc-"
                 "minutes fixes UTC to minutes and hence longitude — the classic "
                 "*lunar-distance* method. That is exactly the lever when a "
                 "photo's timestamp is missing: the sky itself carries the "
                 "clock.\n")

    # ------- unified full-fusion section -------
    ab = data["ablation"]
    L.append("## The unified factor graph — everything fused\n")
    L.append("One graph fuses every observable at once. Per Sun+Moon shot there "
             "is a pose `X(i)` (position + attitude), a velocity `V(i)` and a "
             "shared IMU bias `B`; consecutive keyframes are tied by IMU "
             "preintegration, and each shot contributes six celestial factors "
             "(altitude, azimuth and parallactic line, for the Sun and the "
             "Moon).\n")
    L.append("![factor graph](results/fig_factorgraph.png)\n")
    L.append("**Observables fused** (and where each comes from):\n")
    L.append("- **Altitude** of Sun and Moon — tele pointing measured against "
             "the horizon reference. *Position lines.*")
    L.append("- **Horizon reference** = **IMU gravity** ⊕ **ultrawide optical "
             "horizon** (acceleration-immune, dip-corrected). *Sets the altitude "
             "covariance.*")
    L.append("- **Azimuth** of Sun and Moon — usable because the **tele-disk "
             "orientation** (Moon bright limb, Sun sunspot P-angle) gives a "
             "**magnetometer-free heading**. *Position lines.*")
    L.append("- **Parallactic angle q** of each body — the disk orientation vs. "
             "the vertical. *A per-body position line, but it needs the horizon* "
             "(a single disk gives θ = PA − q − roll; q and roll don't separate "
             "without a vertical).")
    L.append("- **Differential Δq (Sun−Moon)** — the difference of the two disks' "
             "orientations. The shared platform roll **cancels**, so this is the "
             "genuinely **horizon-free** position line. *Needs both disks resolved.*")
    L.append("- **IMU preintegration** between shots (+ bias state). *Links the "
             "trajectory, smooths many fixes.*")
    L.append("- **Least-rotation gating** — selects the calm shutter instants "
             "that feed the graph.")
    L.append("- **Coarse position prior** (a stale ~30 km dead-reckoning). "
             "*Disambiguates; does not drive accuracy.*\n")
    L.append("Deployed accuracy of the **full fusion** (RMS km, "
             f"{N_SEEDS} seeds):\n")
    L.append("| Regime | full fusion |")
    L.append("|---|---|")
    for r in REGIMES:
        L.append(f"| {REGIME_LABEL[r]} | **{cell(ab[r]['full fusion'])}** |")
    L.append("\n### What each observable is worth (leave-one-out)\n")
    L.append("Starting from the full fusion and removing one observable at a "
             "time:\n")
    keys = ["full fusion", "- ultrawide horizon", "- IMU link",
            "- optical azimuth", "- parallactic line", "- Δq (Sun-Moon)",
            "- gating", "- Moon (Sun only)"]
    L.append("| Regime | " + " | ".join(keys) + " |")
    L.append("|" + "---|" * (len(keys) + 1))
    for r in REGIMES:
        L.append(f"| {REGIME_LABEL[r]} | " +
                 " | ".join(cell(ab[r][k]) for k in keys) + " |")
    L.append("\n![ablation](results/fig_ablation.png)\n")
    L.append("### Deliberately left out (and why)\n")
    L.append("- **Lunar distance for time/longitude.** The Sun-Moon separation "
             "(or the Moon's phase) gives chronometer-free GMT — but a phone "
             "already has an accurate clock, so time is known and this is not "
             "needed. It would matter only for a long clock outage.")
    L.append("- **Atmospheric refraction and the Moon's ~1° horizontal "
             "parallax.** Treated as pre-corrected here (geometric altitudes); "
             "`starfix.Sight` already models both and would be layered in for a "
             "field build.")
    L.append("- **Per-shot device attitude** is marginalised into the celestial "
             "factor covariances (the fused gravity/horizon/disk references), "
             "rather than carried as a free state — valid because the tilt "
             "errors are independent shot to shot.\n")

    # ------- real-time section -------
    if "realtime" in data:
        rt = data["realtime"]
        nmax = rt["shot_counts"][-1]

        def _rt(d, n):
            return d[n] if n in d else d[str(n)]
        L.append("## Real-time on iPhone 17 Pro — streaming fixed-lag smoother\n")
        L.append("The batch solver re-optimises the whole trajectory each fix, "
                 "so latency grows with the voyage; the streaming estimator "
                 "(`realtime.py`, `gtsam_unstable.IncrementalFixedLagSmoother`) "
                 "marginalises keyframes older than a time window and "
                 "preintegrates the IMU once online, so each per-shot update is "
                 "**bounded and flat**. Fixes are low-rate (a gated shot every "
                 "few seconds); IMU dead-reckoning gives the continuous position "
                 "between them.\n")
        L.append(f"Per-fix latency at {nmax} shots (host CPU, single thread; "
                 "an A19 Pro is comparable):\n")
        L.append("| Regime | batch re-solve | **streaming update** | speedup | "
                 "final-error parity (batch / stream) |")
        L.append("|---|---|---|---|---|")
        for r in REGIMES:
            bat_ms = _rt(rt["latency"][r]["batch"], nmax)[0]
            str_ms = _rt(rt["latency"][r]["stream"], nmax)[0]
            bf, sf = rt["parity"][r]
            L.append(f"| {REGIME_LABEL[r]} | {bat_ms:.0f} ms | "
                     f"**{str_ms:.1f} ms** | {bat_ms / str_ms:.0f}× | "
                     f"{bf:.2f} / {sf:.2f} km |")
        L.append("\n![realtime](results/fig_realtime.png)\n")
        L.append("The streaming current-position estimate matches batch to "
                 "~0.02 km — same accuracy, bounded cost. (The dominant saving "
                 "is doing IMU preintegration once online instead of re-"
                 "preintegrating every leg on each batch solve; analytic "
                 "Jacobians and the reduced two-DOF finite difference — only "
                 "east/north affect a celestial factor — remove the rest.)\n")

    # ------- teleconverter section -------
    if "zoom" in data:
        zd = data["zoom"]
        zooms = sorted(next(iter(zd["rms"].values())).keys(),
                       key=lambda z: int(z))
        L.append("## Would an external 3× teleconverter help? No.\n")
        L.append("A clip-on afocal optic triples the tele focal length "
                 "(sharper pointing, bigger disk), but the fix is **unchanged** "
                 "— the system is limited by the horizon/attitude reference, not "
                 "the camera. Full-fusion RMS (km) vs teleconverter:\n")
        L.append("| Regime | " + " | ".join(f"{z}×" for z in zooms) + " |")
        L.append("|" + "---|" * (len(zooms) + 1))
        for r in REGIMES:
            L.append(f"| {REGIME_LABEL[r]} | " +
                     " | ".join(cell(zd["rms"][r][z]) for z in zooms) + " |")
        bud = zd["budget"]
        L.append(f"\nWhy: the altitude error budget is dominated by the horizon "
                 f"(optical ~{bud['optical horizon (sea)']:.1f}′, IMU "
                 f"~{bud['IMU horizon (sea)']:.0f}′ at sea), while camera pointing "
                 f"is ~{bud['camera pointing 1x']:.2f}′ — already ~100× smaller, and "
                 f"3× zoom only shrinks that already-negligible term. Heading is "
                 f"floored by astronomical/model residuals (libration, seeing, "
                 f"P-angle) that zoom cannot improve.\n")
        L.append("![zoom](results/fig_zoom.png)\n")
        L.append("The lever that *would* help is a better **horizon/attitude** "
                 "reference — a tripod/gimbal (kills the swell/tremor tilt), a "
                 "longer ultrawide baseline, or a better AHRS — not more zoom. "
                 "Downsides of a 3× optic: 3× narrower field (harder to acquire "
                 "the body), more motion/blur sensitivity, and added "
                 "weight/alignment/aberration.\n")

    # ------- horizon lens section -------
    if "lens" in data:
        ld = data["lens"]
        L.append("## Wide vs. ultrawide horizon lens — an altitude question\n")
        L.append("All the phone's lenses point the same way, so to sight a body "
                 "at altitude *h* the cluster aims at elevation *h* and the "
                 "horizon sits *h* below the boresight. A lens captures the "
                 "horizon only while *h* stays inside its field. The main "
                 "**wide** lens gives a sharper horizon (less distortion) but its "
                 "narrower field loses the horizon above ~30°; the **ultrawide** "
                 "holds it to ~52°; above that neither sees it and the fix falls "
                 "back to the IMU.\n")
        L.append("![lens](results/fig_lens.png)\n")
        L.append("So it is the reverse of *\"wide lens for high hours\"*: **the "
                 "wide lens is the LOW-altitude choice** (a sharper horizon for "
                 "bodies below ~30°), while **high sights need the ultrawide** — "
                 "and very high sights (>52°) lose the optical horizon entirely. "
                 "At the canonical epoch (Moon ~32°, Sun ~40°) both bodies are "
                 "already in the ultrawide zone, so forcing the wide lens loses "
                 "the horizon and wrecks the fix:\n")
        L.append("| Regime | wide-only | ultrawide | adaptive |")
        L.append("|---|---|---|---|")
        for r in ("sea", "air"):
            f = ld["fix"][r]
            L.append(f"| {REGIME_LABEL[r]} | {cell(f['wide'])} | "
                     f"{cell(f['ultrawide'])} | {cell(f['adaptive'])} |")
        L.append("\nPractical guidance: for this method prefer **moderate-to-low "
                 "body altitudes (~15–30°)** — the wide lens then delivers the "
                 "sharpest horizon and refraction is still manageable (below "
                 "~15° refraction/dip uncertainty grows; `starfix` models it). "
                 "High-noon sights are the worst case for the optical horizon.\n")

    # ------- visual-anchor section -------
    if "anchor_drift" in data:
        ad = data["anchor_drift"]
        af = data["anchor_fix"]

        def _at(curve, key, t):
            i = curve["time"].index(t)
            return curve[key][i]
        L.append("## Sunspots as a visual anchor — increasing IMU precision\n")
        L.append("A tracked disk feature (a sunspot through a solar filter, or a "
                 "Moon crater/limb) is a star-tracker landmark: its direction is "
                 "**translation-invariant** (so it behaves the same moving or "
                 "stationary) and **acceleration-immune**. Tracking it pins the "
                 "gyro bias and gives a drift-free attitude.\n")
        st_g = _at(ad["stationary"], "gyro_only", 120)
        st_a = _at(ad["stationary"], "anchor_aided", 120)
        mv_ac = _at(ad["moving"], "accel_aided", 30)
        mv_an = _at(ad["moving"], "anchor_aided", 30)
        L.append("Attitude/horizon error (arc-minutes):\n")
        L.append("| | gyro only | accelerometer-aided | **anchor-aided** |")
        L.append("|---|---|---|---|")
        L.append(f"| Stationary, after 120 s | {st_g:.0f}′ (drifting) | "
                 f"{_at(ad['stationary'],'accel_aided',120):.0f}′ | "
                 f"**{st_a:.1f}′** |")
        L.append(f"| Moving, after 30 s | {_at(ad['moving'],'gyro_only',30):.0f}′ | "
                 f"{mv_ac:.0f}′ (motion-corrupted) | **{mv_an:.1f}′** |")
        L.append("\n![anchor drift](results/fig_anchor_drift.png)\n")
        L.append("The gyro alone diverges (~0.17′/s); the accelerometer is "
                 "bounded but wrecked by motion (~500′ while maneuvering); the "
                 "anchor stays a few arc-minutes **whether moving or "
                 "stationary**. That acceleration-immune attitude, plus the "
                 "position estimate, gives a vertical/horizon that needs neither "
                 "the accelerometer nor the sea horizon — so it **rescues the "
                 "fix when the optical horizon is unavailable** (a high sight, a "
                 "land skyline, or the horizon out of frame):\n")
        L.append("| Regime | optical horizon: no anchor / anchor | "
                 "no optical horizon: no anchor / anchor |")
        L.append("|---|---|---|")
        for r in REGIMES:
            oh = af[r]["optical horizon"]
            nh = af[r]["no optical horizon"]
            L.append(f"| {REGIME_LABEL[r]} | {cell(oh['none'])} / "
                     f"**{cell(oh['anchor'])}** | {cell(nh['none'])} / "
                     f"**{cell(nh['anchor'])}** |")
        L.append("\n![anchor fix](results/fig_anchor_fix.png)\n")
        L.append("So the anchor is transformative exactly where the accelerometer "
                 "and the optical horizon fail, and a modest sharpener elsewhere. "
                 "Cost: it needs a body **continuously tracked** in the tele "
                 "field (the Sun through a filter, or the Moon's features), and "
                 "the anchored vertical is only as good as the position estimate "
                 "(~1–2′ at a ~2 km fix).\n")

    # ------- cloud section -------
    if "cloud" in data:
        cd = data["cloud"]
        cc = data["coast"]

        def _ct(key, t):
            return cc[key][cc["time"].index(t)]
        L.append("## What if cloud obscures the tracked body?\n")
        L.append("Cloud hits both the **sight** for that body (no line of "
                 "position) and the **visual anchor** (the gyro stops being "
                 "calibrated). The system degrades gracefully: obscured shots are "
                 "dropped, the fix continues on whatever is clear, and when both "
                 "bodies are lost it **coasts** on the (freshly calibrated) IMU "
                 "until the sky clears and it re-anchors.\n")
        L.append("**Graceful degradation** — fix RMS (km) vs cloud cover:\n")
        cfs = sorted(next(iter(cd.values())).keys(), key=lambda c: float(c),
                     reverse=True)
        L.append("| Regime | " +
                 " | ".join(f"{round((1-float(c))*100)}% cloud" for c in cfs) +
                 " |")
        L.append("|" + "---|" * (len(cfs) + 1))
        for r in REGIMES:
            L.append(f"| {REGIME_LABEL[r]} | " +
                     " | ".join(cell(cd[r][c]["rms"]) for c in cfs) + " |")
        L.append("\n![cloud](results/fig_cloud.png)\n")
        L.append("**Coast budget** — when both bodies are clouded there are no "
                 "new sights, so position dead-reckons on the IMU. The anchor's "
                 "parting gift is a calibrated gyro, but the coast is ultimately "
                 "limited by gyro random-walk:\n")
        L.append("| Outage | attitude error | DR position drift |")
        L.append("|---|---|---|")
        for t in (30, 60, 120, 300):
            L.append(f"| {t} s | {_ct('att_calibrated', t):.0f}′ | "
                     f"{_ct('pos_calibrated_km', t)*1000:.0f} m |")
        L.append("\n![coast](results/fig_coast.png)\n")
        L.append("So the practical coast budget is **~1–2 minutes** (sub-km) "
                 "before attitude drift leaks into the position quadratically. "
                 "One body clouded → carry on with the other; both clouded → "
                 "coast a minute or two, then it's dead-reckoning until a gap in "
                 "the cloud lets it re-anchor and snap back. Persistent overcast "
                 "→ celestial is unavailable, like any sextant.\n")

    L.append("## 1. Factor graph vs. the incumbent single-fix\n")
    L.append("| Regime | starfix single-fix (per-epoch RMS) | "
             "Factor graph, no IMU | **Factor graph + IMU** |")
    L.append("|---|---|---|---|")
    for r in REGIMES:
        L.append(f"| {REGIME_LABEL[r]} | {cell(m[r]['single'])} | "
                 f"{cell(m[r]['noimu'])} | **{cell(m[r]['imu'])}** |")
    L.append("\n![main](results/fig_main.png)\n")
    L.append("## 2. Least-rotation capture trigger\n")
    L.append("Modelled synthetic-horizon altitude noise at the chosen shutter "
             "instants (arc-minutes):\n")
    L.append("| Regime | gated (least-rotation) | ungated (periodic) |")
    L.append("|---|---|---|")
    for r in REGIMES:
        L.append(f"| {REGIME_LABEL[r]} | {b[r]['gated_arcmin']:.1f}′ | "
                 f"{b[r]['ungated_arcmin']:.1f}′ |")
    L.append("\nResulting position error (factor graph + IMU):\n")
    L.append("| Regime | gated | ungated |")
    L.append("|---|---|---|")
    for r in REGIMES:
        L.append(f"| {REGIME_LABEL[r]} | {cell(g[r]['gated'])} | "
                 f"{cell(g[r]['ungated'])} |")
    L.append("\n![gating](results/fig_gating.png)\n")
    L.append("## 3. Optical horizon from the ultrawide camera\n")
    L.append("The ultrawide and tele lenses fire together: the tele resolves the "
             "body, the ultrawide sees the visible horizon line. Fitting that "
             "line gives an *optical* local-vertical that — unlike the "
             "accelerometer — is not fooled by linear acceleration. Modelled "
             "horizon-reference noise at gated motion (arc-minutes):\n")
    hb = data["horizon_budget"]
    L.append("| Regime | IMU gravity | ultrawide optical | fused | dip corrected |")
    L.append("|---|---|---|---|---|")
    for r in REGIMES:
        avail = "" if hb[r]["available"] else " *(no true horizon → n/a)*"
        uw = f"{hb[r]['uw']:.1f}′{avail}"
        L.append(f"| {REGIME_LABEL[r]} | {hb[r]['imu']:.1f}′ | {uw} | "
                 f"{hb[r]['fused']:.1f}′ | {hb[r]['dip']:.0f}′ |")
    L.append("\nResulting position error (factor graph + IMU, gated):\n")
    L.append("| Regime | IMU horizon | ultrawide horizon | **fused** |")
    L.append("|---|---|---|---|")
    for r in REGIMES:
        L.append(f"| {REGIME_LABEL[r]} | {cell(hz[r]['imu'])} | "
                 f"{cell(hz[r]['uw'])} | **{cell(hz[r]['fused'])}** |")
    L.append("\n![horizon](results/fig_horizon.png)\n")
    L.append("## 4. Tele-resolved disk: magnetometer-free heading + "
             "parallactic line\n")
    L.append("The tele lens resolves the disk, not just a dot. The Moon's "
             "bright limb (and the Sun's sunspot P-angle) give an *absolute* "
             "celestial orientation in the image, which yields the parallactic "
             "angle *q(lat, lon)* — a magnetometer-free heading and a position "
             "line. A single disk's *q* still needs a vertical (θ = PA − q − "
             "roll), but the **difference** of the two disks' orientations, "
             "**Δq(Sun−Moon)**, cancels the shared platform roll and so is a "
             "genuinely **horizon-free** position line (it needs both disks "
             "resolved).\n")
    hd = data["heading_budget"]
    L.append("Heading sigma (degrees) — magnetometer vs. optical disk:\n")
    L.append("| Regime | magnetometer | optical (Moon limb) |")
    L.append("|---|---|---|")
    for r in REGIMES:
        L.append(f"| {REGIME_LABEL[r]} | {hd[r]['mag']:.1f}° | "
                 f"{hd[r]['moon']:.1f}° |")
    op = data["optical"]
    L.append("\nPosition error (factor graph + IMU, gated, **IMU horizon** so "
             "the optical gain is visible):\n")
    L.append("| Regime | altitude only | + az (magnetometer) | "
             "+ az (optical) | **+ optical az & parallactic** |")
    L.append("|---|---|---|---|---|")
    for r in REGIMES:
        L.append(f"| {REGIME_LABEL[r]} | {cell(op[r]['alt'])} | "
                 f"{cell(op[r]['az_mag'])} | {cell(op[r]['az_opt'])} | "
                 f"**{cell(op[r]['optical'])}** |")
    L.append("\nThe optical disk gives a heading several times better than the "
             "magnetometer, which is what makes the azimuth lines of position "
             "usable; combined with the parallactic line it materially improves "
             "the fix when the horizon is weak (e.g. sea on the IMU horizon). "
             "The Sun's P-angle needs a solar filter and visible spots (matched "
             "to an observatory reference taken just before the journey), so the "
             "Moon's bright limb is the workhorse.\n")
    L.append("\n![optical](results/fig_optical.png)\n")
    fs = data["fullstack"]
    L.append("> **Why the numbers above are larger than §3.** This section runs "
             "the optical disk on the *weak IMU horizon* on purpose, to isolate "
             "its contribution. It is not a regression versus the ~2 km ultrawide "
             "horizon — the two use different horizon baselines. Stacking "
             "*everything* (ultrawide fused horizon **and** the optical disk, "
             "Moon + Sun) gives the deployed accuracy below:\n")
    L.append("| Regime | ultrawide horizon only | **full stack (horizon + optical)** |")
    L.append("|---|---|---|")
    for r in REGIMES:
        L.append(f"| {REGIME_LABEL[r]} | {cell(fs[r]['horizon_only'])} | "
                 f"**{cell(fs[r]['full_stack'])}** |")
    L.append("")
    L.append("## 5. Error vs. number of fused shots\n")
    L.append("![convergence](results/fig_convergence.png)\n")
    L.append("## 6. The trigger in action\n")
    L.append("![trigger](results/fig_trigger.png)\n")
    L.append("## 7. Fix with covariance (error ellipse)\n")
    L.append("![ellipse](results/fig_ellipse.png)\n")
    L.append("## Model assumptions\n")
    L.append("- Representative phone-class MEMS IMU and periscope tele camera "
             "(see `iphone_model.py`); numbers are order-of-magnitude, not "
             "datasheet values.\n"
             "- Spherical Earth; observer positions geocentric; Sun/Moon "
             "ephemeris reused verbatim from `starfix` (real nautical almanac).\n"
             "- Altitudes are geometric (refraction/parallax/semidiameter "
             "treated as pre-corrected). A ~30 km dead-reckoning offset seeds "
             "the coarse prior, so reported accuracy comes from the sky, not "
             "the prior.\n"
             "- The synthetic horizon comes from the phone's gravity vector "
             "(no sea-horizon dip); its dominant error is IMU tilt, which the "
             "least-rotation shutter minimises.\n")
    with open(path, "w") as f:
        f.write("\n".join(L))


def write_dashboard(data, path):
    def cell(ms):
        return "n/a" if ms is None or ms[0] is None else f"{ms[0]:.1f} ± {ms[1]:.1f}"
    m, g, b = data["main"], data["gating"], data["sensor_budget"]

    def row(r):
        return (f"<tr><td>{REGIME_LABEL[r]}</td>"
                f"<td>{cell(m[r]['single'])}</td>"
                f"<td>{cell(m[r]['noimu'])}</td>"
                f"<td class='hi'>{cell(m[r]['imu'])}</td></tr>")

    def grow(r):
        return (f"<tr><td>{REGIME_LABEL[r]}</td>"
                f"<td>{b[r]['gated_arcmin']:.1f}′</td>"
                f"<td>{b[r]['ungated_arcmin']:.1f}′</td>"
                f"<td>{cell(g[r]['gated'])}</td>"
                f"<td>{cell(g[r]['ungated'])}</td></tr>")

    hz = data["horizon"]

    def hrow(r):
        return (f"<tr><td>{REGIME_LABEL[r]}</td>"
                f"<td>{cell(hz[r]['imu'])}</td>"
                f"<td>{cell(hz[r]['uw'])}</td>"
                f"<td class='hi'>{cell(hz[r]['fused'])}</td></tr>")

    op, hdb, fs = data["optical"], data["heading_budget"], data["fullstack"]
    ab = data["ablation"]
    rt = data.get("realtime")
    rt_nmax = rt["shot_counts"][-1] if rt else None

    def _rtget(d, n):
        return d[n] if n in d else d[str(n)]

    def rtrow(r):
        bat = _rtget(rt["latency"][r]["batch"], rt_nmax)[0]
        strm = _rtget(rt["latency"][r]["stream"], rt_nmax)[0]
        bf, sf = rt["parity"][r]
        return (f"<tr><td>{REGIME_LABEL[r]}</td><td class='mono dim'>{bat:.0f} ms</td>"
                f"<td class='mono hi'>{strm:.1f} ms</td>"
                f"<td class='mono'>{bat/strm:.0f}×</td>"
                f"<td class='mono'>{bf:.2f} / {sf:.2f} km</td></tr>")
    abkeys = ["- ultrawide horizon", "- IMU link", "- optical azimuth",
              "- parallactic line", "- Δq (Sun-Moon)", "- gating",
              "- Moon (Sun only)"]

    def ffrow(r):
        return (f"<tr><td>{REGIME_LABEL[r]}</td>"
                f"<td class='hi'>{cell(ab[r]['full fusion'])}</td>" +
                "".join(f"<td class='mono dim'>{cell(ab[r][k])}</td>"
                        for k in abkeys) + "</tr>")

    def fsrow(r):
        return (f"<tr><td>{REGIME_LABEL[r]}</td>"
                f"<td>{cell(fs[r]['horizon_only'])}</td>"
                f"<td class='hi'>{cell(fs[r]['full_stack'])}</td></tr>")

    def orow(r):
        return (f"<tr><td>{REGIME_LABEL[r]}</td>"
                f"<td>{hdb[r]['mag']:.1f}° / {hdb[r]['moon']:.1f}°</td>"
                f"<td>{cell(op[r]['alt'])}</td>"
                f"<td>{cell(op[r]['az_mag'])}</td>"
                f"<td>{cell(op[r]['az_opt'])}</td>"
                f"<td class='hi'>{cell(op[r]['optical'])}</td></tr>")

    gt = data.get("groundtruth", {})
    gt_card = ""
    if gt.get("engine"):
        def gtrow(bd):
            gg = gt[bd]
            km = max(gg["gha_max"], gg["dec_max"]) / 3600 * 111.2 * 0.6
            return (f"<tr><td>{bd}</td>"
                    f"<td class='mono'>{gg['gha_rms']:.1f}″ / {gg['gha_max']:.1f}″</td>"
                    f"<td class='mono'>{gg['dec_rms']:.1f}″ / {gg['dec_max']:.1f}″</td>"
                    f"<td class='mono dim'>~{km:.2f} km</td></tr>")
        gt_card = f"""
<div class='card'>
<h2>Are the Sun/Moon positions right? Independent ground-truth check</h2>
<p>Every fix rests on the Sun/Moon positions from <code>starfix</code>'s almanac.
Those are cross-checked here against an <b>independent</b> ephemeris
(<code>{gt['engine']}</code> — a separate implementation from the almanac's
Skyfield/JPL source) over a 90-day grid. Both bodies agree to a
<b>few arc-seconds</b>, far below the study's arc-minute measurement noise.</p>
<table><tr><th>Body</th><th>GHA rms / max</th><th>Dec rms / max</th><th>≈ position</th></tr>
{gtrow('Sun')}{gtrow('Moon')}</table>
<figure><img src='results/fig_groundtruth.png'><figcaption>GHA/Dec residual vs an
independent ephemeris across 90 days — a few arc-seconds, set by the almanac's
0.1′ quantization and hourly interpolation.</figcaption></figure>
<p class='note'><b>This check earned its keep:</b> it first flagged a Moon-declination
error up to ~1.6° near the Dec≈0 crossings — a sign bug in
<code>starfix.parse_angle_string</code> for the <code>-00:MM.M</code> format
(<code>float("-00")</code> is <code>-0.0</code>, and <code>-0.0&nbsp;&lt;&nbsp;0</code>
is false). Now fixed; the canonical epoch (Moon at Dec&nbsp;+28°) was never affected.
A <b>Stellarium</b> export can be dropped in as a third witness — see
<code>stellarium_reference.md</code>.</p>
</div>
"""
    eb = data.get("elongation")
    el_card = ""
    if eb:
        base = eb["base"]
        el_card = f"""
<div class='card'>
<h2>Position error budget — and what the Sun–Moon elongation buys</h2>
<p>Single epoch, σ = 2′ per sight. The two <b>altitude</b> sights are the
workhorses — a ~{base['two_lop_fix_km']:.0f} km fix at ΔAz={base['delta_az_deg']:.0f}°,
pushed to the headline few-km numbers by fusing many gated shots. The <b>Sun–Moon
elongation</b> is a <em>poor direct position line</em> (only lunar parallax moves
it), but a <em>strong time observable</em>: the Moon slides
{base['dE_dt_deg_per_hr']:.2f}°/hr against the Sun, so the separation fixes UTC to
~{base['elong_time_s']:.0f} s → longitude ~{base['elong_longitude_km']:.0f} km — the
classic <b>lunar-distance</b> method, and exactly the lever when a photo's
timestamp is stripped: the sky carries the clock.</p>
<table><tr><th>Observable</th><th>position 1σ</th><th>role</th></tr>
<tr><td>Sun altitude LOP</td><td class='mono'>{base['alt_sun_km']:.1f} km</td><td>position line (1′=1 nmi)</td></tr>
<tr><td>Moon altitude LOP</td><td class='mono'>{base['alt_moon_km']:.1f} km</td><td>position line</td></tr>
<tr><td>Two-body fix</td><td class='mono hi'>{base['two_lop_fix_km']:.1f} km</td><td>LOPs cross at ΔAz={base['delta_az_deg']:.0f}°</td></tr>
<tr><td>Elongation → time</td><td class='mono'>{base['elong_longitude_km']:.0f} km</td><td>clock → longitude</td></tr>
<tr><td>Elongation → direct</td><td class='mono dim'>negligible</td><td>parallax-only</td></tr></table>
<figure><img src='results/fig_elongation.png'><figcaption>Left: per-observable
position σ. Right: the elongation-as-clock precision vs how well the Sun–Moon
separation is measured.</figcaption></figure>
</div>
"""
    imgs = "".join(
        f"<figure><img src='results/{fn}'><figcaption>{cap}</figcaption></figure>"
        for fn, cap in [
            ("fig_main.png", "Factor graph vs incumbent single-fix"),
            ("fig_horizon.png", "Ultrawide optical horizon rescues sea & air"),
            ("fig_optical.png", "Tele-disk heading + parallactic line"),
            ("fig_gating.png", "Least-rotation shutter"),
            ("fig_convergence.png", "Error vs number of shots"),
            ("fig_trigger.png", "Trigger picking calm instants"),
            ("fig_ellipse.png", "Fix with covariance ellipse")])
    html = f"""<title>iPhone Sun+Moon Celestial Fix — IMU + Factor Graph</title>
<style>
:root{{--bg:#0f1419;--card:#1a222c;--fg:#e6edf3;--mut:#9db0c0;--acc:#00b4d8;--hi:#3ddc97}}
@media (prefers-color-scheme:light){{:root{{--bg:#f5f7fa;--card:#fff;--fg:#12212e;--mut:#4a5b6b;--acc:#00798c;--hi:#0a7d54}}}}
:root[data-theme=dark]{{--bg:#0f1419;--card:#1a222c;--fg:#e6edf3;--mut:#9db0c0}}
:root[data-theme=light]{{--bg:#f5f7fa;--card:#fff;--fg:#12212e;--mut:#4a5b6b}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--fg);
font:16px/1.55 -apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif}}
.wrap{{max-width:1000px;margin:0 auto;padding:2rem 1.2rem}}
h1{{font-size:1.7rem;margin:.2rem 0}}.sub{{color:var(--mut);margin-bottom:1.6rem}}
.card{{background:var(--card);border-radius:14px;padding:1.1rem 1.3rem;margin:1rem 0;
box-shadow:0 1px 3px rgba(0,0,0,.25)}}
h2{{font-size:1.15rem;border-left:4px solid var(--acc);padding-left:.6rem}}
table{{width:100%;border-collapse:collapse;margin:.4rem 0;font-size:.93rem}}
th,td{{padding:.5rem .6rem;text-align:left;border-bottom:1px solid rgba(128,128,128,.2)}}
th{{color:var(--mut);font-weight:600}}td.hi{{color:var(--hi);font-weight:700}}
.note{{color:var(--mut);font-size:.88rem;margin:.7rem 0 .3rem}}.note b{{color:var(--fg)}}
figure{{margin:1rem 0}}img{{width:100%;border-radius:10px;background:#fff}}
figcaption{{color:var(--mut);font-size:.85rem;margin-top:.3rem}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:1rem}}
@media(max-width:640px){{.grid{{grid-template-columns:1fr}}}}
.k{{display:inline-block;background:var(--acc);color:#fff;border-radius:6px;
padding:.05rem .5rem;font-size:.8rem;margin-right:.3rem}}
</style>
<div class='wrap'>
<h1>Daytime Sun + Moon Celestial Fix</h1>
<div class='sub'>iPhone 17 Pro (modelled) · IMU synthetic horizon · GTSAM factor graph ·
{EPOCH:%Y-%m-%d %H:%M} UTC near Greenwich · Sun+Moon ~94° apart</div>

<div class='card'>
<h2>Headline: fusing many gated shots with IMU cuts error 4–10×</h2>
<p><span class='k'>{N_SEEDS} seeds</span> RMS horizontal position error (km), mean ± std.</p>
<table><tr><th>Regime</th><th>starfix single-fix</th><th>Factor graph, no IMU</th>
<th>Factor graph + IMU</th></tr>
{row('land')}{row('sea')}{row('air')}</table>
</div>

<div class='card'>
<h2>The unified factor graph — everything fused</h2>
<p>One graph fuses every observable at once: per Sun+Moon shot a pose (position
+ attitude), a velocity and a shared IMU bias, tied across shots by IMU
preintegration, with six celestial factors per shot (altitude, azimuth and
parallactic line, for Sun and Moon). Each factor's covariance is built from the
fused attitude references — IMU gravity, ultrawide optical horizon and
tele-disk orientation.</p>
<figure><img src='results/fig_factorgraph.png'><figcaption>Variables
(circles) and factors (squares): priors, IMU links, and the celestial lines of
position.</figcaption></figure>
<table><tr><th>Regime</th><th>full fusion</th>
<th>− ultrawide</th><th>− IMU</th><th>− opt. az</th><th>− parallactic</th>
<th>− Δq(☉−☾)</th><th>− gating</th><th>− Moon</th></tr>
{ffrow('land')}{ffrow('sea')}{ffrow('air')}</table>
<p class='note'>Leave-one-out: the first column is the full fusion; each other
column removes one observable, so a bigger number is a more valuable observable.</p>
<figure><img src='results/fig_ablation.png'><figcaption>What each observable is
worth inside the unified graph.</figcaption></figure>
</div>

<div class='card'>
<h2>Real-time on iPhone 17 Pro</h2>
<p>Batch re-solves the whole trajectory each fix, so latency grows with the
voyage. The streaming fixed-lag smoother marginalises old keyframes and
preintegrates the IMU once online, so each per-shot update is <b>flat and
bounded</b> — with the same current-position accuracy. Fixes are low-rate; IMU
dead-reckoning covers the gaps.</p>
<table><tr><th>Regime</th><th>batch @ {rt_nmax} shots</th>
<th>streaming update</th><th>speedup</th><th>final-error (batch / stream)</th></tr>
{rtrow('land')}{rtrow('sea')}{rtrow('air')}</table>
<figure><img src='results/fig_realtime.png'><figcaption>Per-fix latency vs trip
length: batch grows, streaming stays flat under the real-time budget.</figcaption></figure>
</div>

<div class='card'>
<h2>Ultrawide optical horizon — the sea &amp; air fix</h2>
<p>Firing the ultrawide lens with the tele gives an <em>optical</em> horizon,
immune to the acceleration that corrupts the IMU gravity horizon. It brings the
moving platforms to land-class accuracy; on land (no true sea horizon) it falls
back to the IMU.</p>
<table><tr><th>Regime</th><th>IMU horizon</th><th>ultrawide horizon</th>
<th>fused</th></tr>
{hrow('land')}{hrow('sea')}{hrow('air')}</table>
</div>

<div class='card'>
<h2>The tele lens as an instrument, not a pointer</h2>
<p>Resolving the disk gives an absolute celestial orientation: a
<em>magnetometer-free heading</em> (which makes the azimuth lines of position
usable) and a <em>parallactic</em> position line. A single disk's parallactic
still needs a vertical (θ = PA − q − roll), so it is <em>not</em> horizon-free on
its own; the genuinely horizon-free line is the <b>difference</b> of the two
disks' orientations, <b>Δq(Sun−Moon)</b>, where the shared platform roll cancels.
As a compass the bodies are unequal — the <b>Sun's</b> sharp disk is precise,
the <b>Moon's bright limb</b> only ~2° (phase-limited, degenerate near full),
looser than a magnetometer — so in daytime take the heading from the <b>Sun</b>.
Shown on the weak IMU horizon so the gain is visible.</p>
<table><tr><th>Regime</th><th>heading σ (mag / optical)</th><th>alt only</th>
<th>+az (mag)</th><th>+az (optical)</th><th>+optical &amp; parallactic</th></tr>
{orow('land')}{orow('sea')}{orow('air')}</table>
<p class="note"><b>These numbers use the weak IMU horizon on purpose</b>, to
isolate the optical-disk contribution — not a regression versus the ~2&nbsp;km
ultrawide horizon. Stacking <em>everything</em> (ultrawide horizon + optical
disk, Moon&nbsp;+&nbsp;Sun) gives the deployed accuracy:</p>
<table><tr><th>Regime</th><th>ultrawide horizon only</th>
<th>full stack (horizon + optical)</th></tr>
{fsrow('land')}{fsrow('sea')}{fsrow('air')}</table>
</div>

<div class='card'>
<h2>Would a 3× external teleconverter help? No.</h2>
<p>A clip-on optic triples the tele focal length, but the fix is unchanged — the
system is limited by the horizon/attitude reference, not the camera (already
~100× sharper). More zoom only shrinks an already-negligible term; heading is
floored by astronomical residuals (libration, seeing, P-angle) that zoom can't
fix. The lever that helps is a better horizon — a gimbal, longer ultrawide
baseline, or better AHRS.</p>
<figure><img src='results/fig_zoom.png'><figcaption>Left: full-fusion error is
flat across 1×/2×/3×. Right: the camera sits ~100× below the horizon in the
altitude budget.</figcaption></figure>
</div>

<div class='card'>
<h2>Wide vs ultrawide horizon lens — an altitude question</h2>
<p>All lenses point the same way, so sighting a body at altitude <i>h</i> puts
the horizon <i>h</i> below the boresight — captured only while it stays in the
horizon lens's field. The sharper <b>wide</b> lens frames the horizon only below
~30°; the <b>ultrawide</b> holds it to ~52°; above that the fix falls back to the
IMU. So it is the reverse of "wide lens for high hours": the wide lens is the
<em>low</em>-altitude choice, high sights need the ultrawide, and very high
sights lose the optical horizon. Prefer moderate/low sights (~15–30°).</p>
<figure><img src='results/fig_lens.png'><figcaption>Horizon reference σ vs body
altitude: wide zone, ultrawide zone, then IMU-only.</figcaption></figure>
</div>

<div class='card'>
<h2>Sunspots as a visual anchor — sharper IMU, moving or still</h2>
<p>A tracked disk feature (a sunspot through a filter, or a Moon crater/limb) is
a star-tracker landmark: translation-invariant and acceleration-immune. Tracking
it pins the gyro bias, so attitude stays a few arc-minutes <b>whether moving or
stationary</b> — while the gyro alone drifts (~0.17′/s) and the accelerometer is
wrecked by motion (~500′ maneuvering). That acceleration-immune vertical rescues
the fix when the optical horizon is unavailable (high sights, land, out of
frame):</p>
<figure><img src='results/fig_anchor_drift.png'><figcaption>Attitude error vs
tracking time: the anchor curve is identical stationary and moving.</figcaption></figure>
<figure><img src='results/fig_anchor_fix.png'><figcaption>With no optical horizon,
the sunspot anchor cuts the fix error several-fold.</figcaption></figure>
</div>

<div class='card'>
<h2>What if cloud obscures the tracked body?</h2>
<p>Cloud hits both the sight (no line of position) and the anchor (gyro
calibration stops). It degrades gracefully: obscured shots are dropped, the fix
continues on whatever is clear, and when both bodies are lost it coasts on the
freshly-calibrated IMU until a gap lets it re-anchor. The practical coast budget
is ~1–2 minutes (sub-km) before attitude drift leaks into position; persistent
overcast makes celestial unavailable, like any sextant.</p>
<figure><img src='results/fig_cloud.png'><figcaption>Fix error rises smoothly
with cloud cover; shot availability tracks the clear sky.</figcaption></figure>
<figure><img src='results/fig_coast.png'><figcaption>Coast budget when both
bodies are lost: sub-km for ~2 minutes, then dead-reckoning runs away.</figcaption></figure>
</div>
{gt_card}
{el_card}

<div class='card'>
<h2>Least-rotation capture trigger</h2>
<p>Firing the shutter at the calmest instants lowers the synthetic-horizon
tilt error 3–6×. That converts to ≈2× lower fix error on land and in the air;
at sea even the calmest swell instant is too tilted, so the swell floor
dominates and a gimbal/mount is needed.</p>
<table><tr><th>Regime</th><th>gated σ</th><th>ungated σ</th><th>gated fix</th>
<th>ungated fix</th></tr>
{grow('land')}{grow('sea')}{grow('air')}</table>
</div>

<div class='card'><h2>Figures</h2>{imgs}</div>

<div class='card'><h2>How it works</h2>
<p>The phone is a digital theodolite: its IMU gives the local vertical (a
gravity-derived artificial horizon, no sea-horizon dip) and the camera gives the
direction to the Sun and Moon. Each shot yields two altitude lines of position;
GTSAM fuses them across shots, linked by IMU preintegration, into one trajectory
with a covariance. A single hand-held phone sight is weak (~6–70′ of horizon
tilt); the win comes from fusing many <em>calm</em> shots.</p></div>
</div>"""
    with open(path, "w") as f:
        f.write(html)


# --------------------------------------------------------------------------- #

def main():
    os.makedirs(OUT, exist_ok=True)
    print("Running experiments (this takes a few minutes)…")
    data = dict(main=exp_main(), gating=exp_gating(),
                horizon=exp_horizon(), horizon_budget=exp_horizon_budget(),
                optical=exp_optical(), heading_budget=exp_heading_budget(),
                fullstack=exp_fullstack(), ablation=exp_ablation(),
                realtime=exp_realtime(), zoom=exp_zoom(), lens=exp_lens(),
                anchor_drift=exp_anchor_drift(), anchor_fix=exp_anchor_fix(),
                cloud=exp_cloud(), coast=exp_coast(),
                groundtruth=exp_groundtruth(),
                elongation=exp_elongation_budget(),
                convergence=exp_convergence(), sensor_budget=exp_sensor_budget())
    with open(os.path.join(OUT, "results.json"), "w") as f:
        json.dump(data, f, indent=2, default=lambda o: list(o)
                  if isinstance(o, tuple) else o)
    print("Plotting…")
    plot_main(data["main"], os.path.join(OUT, "fig_main.png"))
    plot_gating(data["gating"], os.path.join(OUT, "fig_gating.png"))
    plot_horizon(data["horizon"], os.path.join(OUT, "fig_horizon.png"))
    plot_optical(data["optical"], os.path.join(OUT, "fig_optical.png"))
    plot_elongation(data["elongation"], os.path.join(OUT, "fig_elongation.png"))
    plot_factorgraph(os.path.join(OUT, "fig_factorgraph.png"))
    plot_ablation(data["ablation"], os.path.join(OUT, "fig_ablation.png"))
    plot_realtime(data["realtime"], os.path.join(OUT, "fig_realtime.png"))
    plot_zoom(data["zoom"], os.path.join(OUT, "fig_zoom.png"))
    plot_lens(data["lens"], os.path.join(OUT, "fig_lens.png"))
    plot_anchor_drift(data["anchor_drift"],
                      os.path.join(OUT, "fig_anchor_drift.png"))
    plot_anchor_fix(data["anchor_fix"], os.path.join(OUT, "fig_anchor_fix.png"))
    plot_cloud(data["cloud"], os.path.join(OUT, "fig_cloud.png"))
    plot_coast(data["coast"], os.path.join(OUT, "fig_coast.png"))
    if data.get("groundtruth", {}).get("engine"):
        plot_groundtruth(data["groundtruth"],
                         os.path.join(OUT, "fig_groundtruth.png"))
    plot_convergence(data["convergence"], os.path.join(OUT, "fig_convergence.png"))
    plot_trigger(os.path.join(OUT, "fig_trigger.png"))
    plot_ellipse(os.path.join(OUT, "fig_ellipse.png"))
    write_results_md(data, os.path.join(HERE, "RESULTS.md"))
    write_dashboard(data, os.path.join(HERE, "dashboard.html"))
    print("Done. See imu_fusion/RESULTS.md and imu_fusion/results/.")
    for r in REGIMES:
        mm = data["main"][r]
        print(f"  {r:5} single={mm['single'][0]:.1f}  noIMU={mm['noimu'][0]:.1f}"
              f"  IMU={mm['imu'][0]:.1f} km")


if __name__ == "__main__":
    main()
