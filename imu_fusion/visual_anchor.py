''' Sunspot / disk-feature VISUAL ANCHOR for aiding the IMU.

    A resolved celestial disk with identifiable features (sunspots through a
    solar filter, or the Moon's craters/limb) is a star-tracker landmark:

      * its DIRECTION is translation-invariant (a distant body shows no parallax
        over a short baseline), so tracking a feature is a pure ATTITUDE
        reference that behaves identically whether the phone moves or is
        stationary — only rotation moves it in the image;
      * it is ACCELERATION-IMMUNE, unlike the accelerometer's gravity vertical,
        which is corrupted by any linear acceleration.

    Continuously tracking the feature across frames therefore lets the absolute
    (ephemeris-known) celestial reference PIN the gyroscope bias — bounding the
    attitude drift that a free-running gyro accumulates (~0.17'/s here) — and
    yields a drift-free attitude.  Combined with the position estimate that
    attitude gives a local vertical (hence a horizon) that needs neither the
    accelerometer nor the sea horizon.

    This module models the effect (a complementary filter), not a full VIO: it
    returns attitude-error-vs-time curves for gyro-only / accelerometer-aided /
    anchor-aided operation, and an anchor-based horizon-reference sigma.

    (c) 2026.  MIT License (see LICENSE file).
'''

from dataclasses import dataclass
from math import sqrt, degrees

from .iphone_model import (IphoneImuSpec, DEFAULT_IMU, DEFAULT_CAM,
                           TeleCameraSpec, KinematicState, G0, ARCMIN_PER_RAD)
from .optical_attitude import orientation_sigma_deg, OpticalDiskSpec, DEFAULT_DISK

# Earth angle scale: 1 degree of local-vertical direction ~ 111 km on the ground.
KM_PER_DEG = 111.195


@dataclass(frozen=True)
class AnchorTrackSpec:
    ''' How the feature is tracked for attitude. '''
    rate_hz: float = 20.0          # visual tracking / attitude-update rate
    body: str = "Sun"              # which body's disk provides the anchor


DEFAULT_TRACK = AnchorTrackSpec()


def gyro_only_attitude_arcmin(t_s: float, imu: IphoneImuSpec = DEFAULT_IMU
                              ) -> float:
    ''' Attitude error after propagating the gyro alone for t_s seconds
        [arcmin].  Bias gives a linear ramp; angle-random-walk adds a sqrt(t)
        term.  This is what DIVERGES without an absolute reference. '''
    bias_ramp = imu.gyro_bias * t_s                       # rad
    arw = imu.gyro_noise_density * sqrt(max(t_s, 0.0))    # rad (angle RW)
    biasrw = 0.5 * imu.gyro_bias_rw * (max(t_s, 0.0) ** 1.5)  # rad, growing bias
    return sqrt(bias_ramp ** 2 + arw ** 2 + biasrw ** 2) * ARCMIN_PER_RAD


def _anchor_frame_sigma_arcmin(state: KinematicState,
                               imu: IphoneImuSpec,
                               cam: TeleCameraSpec,
                               disk: OpticalDiskSpec,
                               track: AnchorTrackSpec) -> float:
    ''' Single-frame attitude sigma from the tracked disk feature [arcmin]. '''
    return orientation_sigma_deg(track.body, state, cam, disk) * 60.0


def anchor_attitude_arcmin(track_time_s: float, state: KinematicState,
                           imu: IphoneImuSpec = DEFAULT_IMU,
                           cam: TeleCameraSpec = DEFAULT_CAM,
                           disk: OpticalDiskSpec = DEFAULT_DISK,
                           track: AnchorTrackSpec = DEFAULT_TRACK) -> float:
    ''' Attitude error with the visual anchor after tracking for track_time_s
        [arcmin].  A complementary filter: the gyro carries the short term while
        the absolute anchor bounds the long term.  Averaging N = rate * time
        frames drives the anchor's per-frame error down and estimates (removes)
        the gyro bias, so the result is BOUNDED and acceleration-immune — it does
        not depend on whether the platform is moving.
    '''
    frame = _anchor_frame_sigma_arcmin(state, imu, cam, disk, track)
    n = max(1.0, track.rate_hz * track_time_s)
    averaged = frame / sqrt(n)
    # A small floor: residual feature-model / tracking error that does not
    # average away (libration model, spot proper motion, centroid bias).
    floor = 1.5
    # Short-term gyro contribution between anchor frames (one frame interval).
    gyro_between = gyro_only_attitude_arcmin(1.0 / track.rate_hz, imu)
    return sqrt(averaged ** 2 + floor ** 2 + gyro_between ** 2)


def accel_aided_attitude_arcmin(state: KinematicState,
                                imu: IphoneImuSpec = DEFAULT_IMU) -> float:
    ''' Attitude from the accelerometer-aided AHRS (the current model): a static
        floor plus the motion-dependent tilt from linear acceleration and
        rotation — bounded but ACCELERATION-CORRUPTED. '''
    from .iphone_model import gravity_tilt_sigma_arcmin
    return gravity_tilt_sigma_arcmin(state, imu)


def anchor_horizon_sigma_arcmin(pos_km: float, track_time_s: float,
                                state: KinematicState,
                                imu: IphoneImuSpec = DEFAULT_IMU,
                                cam: TeleCameraSpec = DEFAULT_CAM,
                                disk: OpticalDiskSpec = DEFAULT_DISK,
                                track: AnchorTrackSpec = DEFAULT_TRACK) -> float:
    ''' Horizon (local-vertical) sigma from the anchor-calibrated attitude plus
        the position estimate [arcmin].  The vertical direction in inertial space
        depends on position (~1 deg per 111 km), so the anchored vertical is the
        anchor attitude error combined with the position-induced vertical error.
        Acceleration-immune and body-altitude independent — available even when
        the optical sea-horizon is not (high sights, land, out of frame).
    '''
    att = anchor_attitude_arcmin(track_time_s, state, imu, cam, disk, track)
    # Position-induced vertical error: pos_km/111 deg -> arcmin.
    pos_vertical = (pos_km / KM_PER_DEG) * 60.0
    return sqrt(att ** 2 + pos_vertical ** 2)


# --------------------------------------------------------------------------- #
# Coasting when the anchor is lost (cloud occlusion).
# --------------------------------------------------------------------------- #

# Fraction of the raw gyro bias that survives after the anchor has calibrated it.
CALIBRATED_BIAS_FRACTION = 0.2


def coast_attitude_arcmin(outage_s: float, imu: IphoneImuSpec = DEFAULT_IMU,
                          calibrated: bool = True,
                          start_arcmin: float = 2.0) -> float:
    ''' Attitude error after the anchor is lost for `outage_s` seconds.

        The gyro coasts from the last anchored attitude (`start_arcmin`).  If the
        anchor had just CALIBRATED the gyro, only a small residual bias remains,
        so the coast is slow; without calibration the full raw bias drifts fast.
    '''
    frac = CALIBRATED_BIAS_FRACTION if calibrated else 1.0
    bias_ramp = frac * imu.gyro_bias * outage_s
    arw = imu.gyro_noise_density * sqrt(max(outage_s, 0.0))
    biasrw = 0.5 * imu.gyro_bias_rw * (max(outage_s, 0.0) ** 1.5)
    drift = sqrt(bias_ramp ** 2 + arw ** 2 + biasrw ** 2) * ARCMIN_PER_RAD
    return sqrt(start_arcmin ** 2 + drift ** 2)


def deadreckon_position_km(outage_s: float, vel_sigma_ms: float,
                           attitude_arcmin: float,
                           imu: IphoneImuSpec = DEFAULT_IMU) -> float:
    ''' Inertial dead-reckoning position drift over a celestial outage [km].

        Two dominant terms: the current velocity error carried forward
        (vel_sigma * t), and a tilt-leak — an attitude error theta tips the
        gravity subtraction, injecting a horizontal acceleration ~ g*theta that
        integrates twice over time.
    '''
    theta = attitude_arcmin / ARCMIN_PER_RAD                # rad
    vel_term = vel_sigma_ms * outage_s
    tilt_term = 0.5 * G0 * theta * outage_s ** 2
    return sqrt(vel_term ** 2 + tilt_term ** 2) / 1000.0


def coast_curve(times, vel_sigma_ms: float = 0.5) -> dict:
    ''' Attitude and DR-position drift vs outage duration, with vs without the
        anchor having pre-calibrated the gyro. '''
    att_cal = [coast_attitude_arcmin(t, calibrated=True) for t in times]
    att_raw = [coast_attitude_arcmin(t, calibrated=False) for t in times]
    pos_cal = [deadreckon_position_km(t, vel_sigma_ms, a)
               for t, a in zip(times, att_cal)]
    pos_raw = [deadreckon_position_km(t, vel_sigma_ms, a)
               for t, a in zip(times, att_raw)]
    return {"time": list(times), "att_calibrated": att_cal,
            "att_uncalibrated": att_raw, "pos_calibrated_km": pos_cal,
            "pos_uncalibrated_km": pos_raw}


def attitude_error_curve(times, moving: bool,
                         imu: IphoneImuSpec = DEFAULT_IMU,
                         cam: TeleCameraSpec = DEFAULT_CAM,
                         track: AnchorTrackSpec = DEFAULT_TRACK) -> dict:
    ''' Attitude-error-vs-time (arcmin) for the three regimes.  `moving` selects
        a maneuvering kinematic state (large linear accel) vs a near-still one;
        the anchor curve is nearly identical for both — that is the point. '''
    if moving:
        state = KinematicState(ang_rate=0.30, lin_accel=1.5)
    else:
        state = KinematicState(ang_rate=0.02, lin_accel=0.03)
    gyro = [gyro_only_attitude_arcmin(t, imu) for t in times]
    accel = [accel_aided_attitude_arcmin(state, imu) for _ in times]
    anchor = [anchor_attitude_arcmin(t, state, imu, cam, DEFAULT_DISK, track)
              for t in times]
    return {"time": list(times), "gyro_only": gyro,
            "accel_aided": accel, "anchor_aided": anchor}


def summarise() -> str:
    still = KinematicState(0.02, 0.03)
    moving = KinematicState(0.30, 1.5)
    lines = ["Visual-anchor attitude aiding (arcmin):"]
    for t in (1, 5, 30, 120):
        lines.append(
            f"  t={t:4}s  gyro-only={gyro_only_attitude_arcmin(t):7.1f}  "
            f"anchor(still)={anchor_attitude_arcmin(t, still):5.1f}  "
            f"anchor(moving)={anchor_attitude_arcmin(t, moving):5.1f}")
    lines.append(f"  accel-aided still={accel_aided_attitude_arcmin(still):.1f}  "
                 f"moving={accel_aided_attitude_arcmin(moving):.1f}")
    lines.append(f"  anchor horizon (2 km fix, 10 s track, moving) = "
                 f"{anchor_horizon_sigma_arcmin(2.0, 10.0, moving):.1f} arcmin")
    return "\n".join(lines)


if __name__ == "__main__":
    print(summarise())
