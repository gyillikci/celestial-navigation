''' Optical horizon from the ultrawide camera, captured simultaneously with the
    tele shot of the celestial body.

    The iPhone fires the ultrawide and tele lenses together (fixed, factory
    calibrated relative orientation).  The tele resolves the Sun/Moon; the
    ultrawide sees the wide scene INCLUDING the visible horizon line.  Fitting
    that line gives the camera's tilt relative to true horizontal -- an
    *optical* artificial horizon that competes with (and can be fused with) the
    IMU gravity vector.

    Why it matters
    --------------
    The IMU horizon is corrupted by linear acceleration (a horizontal specific
    force tilts the accel-derived "down"); this was the dominant error at sea in
    the baseline study, and gating barely helped.  The OPTICAL horizon is
    geometric -- it does not care about acceleration -- so a visible sea horizon
    is a far better reference on a pitching boat or a vibrating aircraft.

    Error sources modelled (added in quadrature):
      * line-fit detection   -- tiny: the line is fit across thousands of pixels;
      * residual lens distortion of the ultrawide after calibration;
      * dip-of-horizon residual -- the dip is corrected from height of eye
        (reusing starfix.get_dip_of_horizon), leaving a small residual from
        height uncertainty and refraction variability;
      * surface / visibility  -- wave roughening of the sea horizon, or haze
        aloft;
      * rotation during the (short) exposure.

    Validity
    --------
    A true geometric horizon needs an unobstructed sea (or a known-height air
    horizon).  On LAND the "horizon" is a terrain/building skyline at unknown
    elevation -- not the astronomical horizon -- so the optical horizon is
    treated as UNAVAILABLE there and the estimator falls back to the IMU.

    Like gravity, the horizon constrains only tilt (roll+pitch), not heading.

    (c) 2026.  MIT License (see LICENSE file).
'''

from dataclasses import dataclass
from math import sqrt, radians

from starfix import get_dip_of_horizon

from .iphone_model import (IphoneImuSpec, KinematicState, DEFAULT_IMU,
                           gravity_tilt_sigma_arcmin, ARCMIN_PER_RAD)

# Where a true optical horizon is usable.
HORIZON_AVAILABLE = {"land": False, "sea": True, "air": True}

# Representative height of eye per regime [m] (drives the dip correction).
REGIME_HEIGHT_M = {"land": 2.0, "sea": 3.0, "air": 10000.0}


@dataclass(frozen=True)
class UltrawideHorizonSpec:
    ''' Representative ultrawide-camera horizon-sensing model. '''
    fov_deg: float = 120.0
    width_px: int = 4032
    detect_px_sigma: float = 1.0          # per-pixel horizon edge localisation
    distortion_arcmin: float = 2.0        # residual lens distortion after calib
    # Regime-dependent contributions [arc minutes].
    dip_residual_arcmin: dict = None      # height/refraction residual of the dip
    surface_arcmin: dict = None           # wave roughness / haze of the horizon

    def __post_init__(self):
        # Frozen dataclass: set mutable defaults via object.__setattr__.
        if self.dip_residual_arcmin is None:
            object.__setattr__(self, "dip_residual_arcmin",
                               {"sea": 0.5, "air": 2.0})
        if self.surface_arcmin is None:
            object.__setattr__(self, "surface_arcmin",
                               {"sea": 1.5, "air": 2.0})

    def arcmin_per_px(self) -> float:
        return self.fov_deg * 60.0 / self.width_px

    def detect_sigma_arcmin(self, n_fit_points: int = 1500) -> float:
        ''' Horizon-line orientation sigma from fitting many edge pixels. '''
        return self.detect_px_sigma * self.arcmin_per_px() / sqrt(n_fit_points)


DEFAULT_UW = UltrawideHorizonSpec()


# --------------------------------------------------------------------------- #
# Horizon-lens choice (wide main lens vs ultrawide).
#
# On a phone all lenses point the SAME way, so to sight a body at altitude h the
# cluster is aimed at elevation h and the horizon sits h below the boresight.
# A lens captures the horizon only if h (plus a margin to clear the distorted
# edge) stays inside its half field of view.  The main WIDE lens has less
# distortion and finer resolution than the ultrawide -> a SHARPER horizon -> a
# lower tilt sigma -> a better fix -- but its narrower field loses the horizon at
# a lower body altitude.  So the right lens depends on the body's altitude.
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class HorizonLensSpec:
    ''' A lens used to image the horizon. '''
    name: str
    half_fov_deg: float          # half of the full field of view
    distortion_arcmin: float     # residual distortion after calibration
    fov_margin_deg: float = 5.0  # keep the horizon off the distorted extreme edge


# Representative: the ultrawide sees the horizon up to a high body altitude but
# is distorted; the main wide lens is sharper but frames the horizon only for
# lower bodies.
ULTRAWIDE_LENS = HorizonLensSpec("ultrawide", 57.0, 2.0)
WIDE_LENS = HorizonLensSpec("wide", 35.0, 0.7)
HORIZON_LENSES = (WIDE_LENS, ULTRAWIDE_LENS)     # sharpest first


def lens_sees_horizon(lens: HorizonLensSpec, body_alt_deg: float) -> bool:
    ''' True if the horizon is within this lens's field when sighting a body at
        `body_alt_deg`. '''
    return body_alt_deg + lens.fov_margin_deg <= lens.half_fov_deg


def lens_horizon_sigma_arcmin(lens: HorizonLensSpec, state: KinematicState,
                              regime: str, uw: UltrawideHorizonSpec = DEFAULT_UW):
    ''' Optical horizon tilt sigma [arcmin] for a specific horizon lens. '''
    dip_res = uw.dip_residual_arcmin.get(regime, 1.0)
    surface = uw.surface_arcmin.get(regime, 1.5)
    detect = uw.detect_sigma_arcmin()
    rot = state.ang_rate * state.exposure_s * ARCMIN_PER_RAD
    return sqrt(detect ** 2 + lens.distortion_arcmin ** 2 + dip_res ** 2 +
                surface ** 2 + rot ** 2)


def best_horizon_lens(body_alt_deg: float, regime: str):
    ''' The sharpest horizon lens that still frames the horizon for this body
        altitude, or None if no optical horizon is available (too high, or no
        true horizon on land). '''
    if not HORIZON_AVAILABLE.get(regime, False):
        return None
    for lens in HORIZON_LENSES:            # sharpest first
        if lens_sees_horizon(lens, body_alt_deg):
            return lens
    return None


def horizon_reference_sigma_lens(mode: str, state: KinematicState, regime: str,
                                 body_alt_deg: float, lens_policy: str = "adaptive",
                                 imu: IphoneImuSpec = DEFAULT_IMU,
                                 uw: UltrawideHorizonSpec = DEFAULT_UW) -> float:
    ''' Horizon (local-vertical) reference sigma [arcmin] with an explicit lens
        policy and the body altitude that gates the field of view.

        lens_policy: "ultrawide" | "wide" | "adaptive" (wide when it frames the
        horizon, else ultrawide).  If the optical horizon is unavailable (body
        too high for the chosen lens, or land), the IMU gravity horizon is used.
    '''
    s_imu = gravity_tilt_sigma_arcmin(state, imu)
    if mode == "imu":
        return s_imu

    if lens_policy == "wide":
        lens = WIDE_LENS if lens_sees_horizon(WIDE_LENS, body_alt_deg) else None
    elif lens_policy == "ultrawide":
        lens = (ULTRAWIDE_LENS if lens_sees_horizon(ULTRAWIDE_LENS, body_alt_deg)
                else None)
    else:                                  # adaptive
        lens = best_horizon_lens(body_alt_deg, regime)
    if lens is None or not HORIZON_AVAILABLE.get(regime, False):
        return s_imu                       # fall back to gravity
    s_opt = lens_horizon_sigma_arcmin(lens, state, regime, uw)
    if mode == "uw":
        return s_opt
    # fused: inverse-variance combine gravity + optical
    return 1.0 / sqrt(1.0 / s_imu ** 2 + 1.0 / s_opt ** 2)


def dip_arcmin(regime: str, height_m: float = None,
               temperature: float = 10.0) -> float:
    ''' Nominal dip of the horizon for a regime (reuses starfix). '''
    h = REGIME_HEIGHT_M[regime] if height_m is None else height_m
    return get_dip_of_horizon(h, temperature)


def ultrawide_horizon_sigma_arcmin(state: KinematicState, regime: str,
                                   uw: UltrawideHorizonSpec = DEFAULT_UW):
    ''' 1-sigma of the OPTICAL horizon tilt reference [arc minutes], and whether
        it is available for this regime.

        Returns (sigma_arcmin, available).
    '''
    if not HORIZON_AVAILABLE.get(regime, False):
        return None, False
    detect = uw.detect_sigma_arcmin()
    distortion = uw.distortion_arcmin
    dip_res = uw.dip_residual_arcmin.get(regime, 1.0)
    surface = uw.surface_arcmin.get(regime, 1.5)
    # Rotation during the exposure smears the fitted line too.
    rot = state.ang_rate * state.exposure_s * ARCMIN_PER_RAD
    sigma = sqrt(detect ** 2 + distortion ** 2 + dip_res ** 2 +
                 surface ** 2 + rot ** 2)
    return sigma, True


def horizon_reference_sigma_arcmin(mode: str, state: KinematicState,
                                   regime: str,
                                   imu: IphoneImuSpec = DEFAULT_IMU,
                                   uw: UltrawideHorizonSpec = DEFAULT_UW) -> float:
    ''' 1-sigma of the fused horizon (local-vertical) reference [arc minutes].

        mode = "imu"   : IMU gravity vector only (baseline study).
        mode = "uw"    : optical ultrawide horizon only (falls back to IMU where
                         no true horizon exists, e.g. land).
        mode = "fused" : inverse-variance fusion of gravity + optical horizon.
    '''
    s_imu = gravity_tilt_sigma_arcmin(state, imu)
    s_uw, valid = ultrawide_horizon_sigma_arcmin(state, regime, uw)
    if mode == "imu" or not valid:
        return s_imu
    if mode == "uw":
        return s_uw
    if mode == "fused":
        return 1.0 / sqrt(1.0 / s_imu ** 2 + 1.0 / s_uw ** 2)
    raise ValueError(f"unknown horizon mode {mode!r}")


def summarise() -> str:
    ''' Compare horizon references at representative gated motion. '''
    still = {"land": KinematicState(0.02, 0.03),
             "sea": KinematicState(0.10, 0.20),
             "air": KinematicState(0.05, 0.30)}
    lines = ["Horizon reference sigma (arcmin) at gated motion:",
             f"  {'regime':5} {'IMU':>8} {'ultrawide':>10} {'fused':>8}  dip"]
    for r in ("land", "sea", "air"):
        st = still[r]
        si = horizon_reference_sigma_arcmin("imu", st, r)
        su = horizon_reference_sigma_arcmin("uw", st, r)
        sf = horizon_reference_sigma_arcmin("fused", st, r)
        dip = dip_arcmin(r)
        avail = "" if HORIZON_AVAILABLE[r] else " (no optical horizon)"
        lines.append(f"  {r:5} {si:8.1f} {su:10.1f} {sf:8.1f}  "
                     f"{dip:6.1f}'{avail}")
    return "\n".join(lines)


if __name__ == "__main__":
    print(summarise())
