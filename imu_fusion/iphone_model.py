''' iPhone 17 Pro sensor & camera model for daytime celestial sighting.

    Models the two things a phone contributes to a celestial sight:

      1. The IMU (accelerometer + gyroscope) defines the *local vertical* from
         gravity -- a synthetic "artificial horizon" that replaces the sea
         horizon.  No dip correction is needed, but the horizon is only as good
         as the accelerometer's estimate of "down", which is corrupted whenever
         the phone is accelerating or rotating.

      2. The camera measures the *direction to the body* in the phone body frame
         (the pixel where the Sun/Moon disk centre falls).  Combined with the
         IMU attitude this yields a measured altitude and azimuth of the body --
         the phone acts as a digital theodolite.

    All numeric specs below are REPRESENTATIVE values for a modern phone-class
    MEMS IMU and a periscope tele camera.  They are order-of-magnitude correct
    and clearly documented as modelling assumptions, not manufacturer data.

    (c) 2026.  MIT License (see LICENSE file).
'''

from dataclasses import dataclass
from math import sqrt, radians, degrees

G0 = 9.80665           # m/s^2, standard gravity
ARCMIN_PER_RAD = 60.0 * 180.0 / 3.141592653589793


@dataclass(frozen=True)
class IphoneImuSpec:
    ''' Representative phone-class MEMS IMU noise model (per axis).

        Noise *densities* are converted to per-sample sigmas given the sample
        rate.  Bias instabilities are modelled as slow random walks handled by
        the factor graph's bias state.
    '''
    sample_rate_hz: float = 100.0
    # Accelerometer white noise density  [ (m/s^2) / sqrt(Hz) ]
    accel_noise_density: float = 1.5e-3      # ~150 ug/rtHz
    # Gyroscope white noise density       [ (rad/s) / sqrt(Hz) ]
    gyro_noise_density: float = 1.4e-4       # ~0.008 deg/s/rtHz
    # Bias random-walk (continuous)       [ (m/s^2)/sqrt(s), (rad/s)/sqrt(s) ]
    accel_bias_rw: float = 3.0e-4
    gyro_bias_rw: float = 2.0e-5
    # Static bias magnitudes (1-sigma)
    accel_bias: float = 2.0e-2               # ~2 mg
    gyro_bias: float = 5.0e-5                # ~10 deg/hr
    # Residual tilt floor of the attitude (AHRS) filter after in-motion
    # gyro/accel-bias self-calibration [arc minutes].  This is the empirical
    # ~0.1-0.2 deg accuracy floor of a phone "artificial horizon" and is the
    # reason a single phone sight is ~10x worse than a marine sextant.
    static_tilt_arcmin: float = 6.0
    # Window over which the AHRS averages the accelerometer for gravity [s].
    attitude_filter_window_s: float = 0.5

    def accel_sigma(self) -> float:
        ''' Per-sample accelerometer sigma [m/s^2]. '''
        return self.accel_noise_density * sqrt(self.sample_rate_hz)

    def filtered_accel_sigma(self) -> float:
        ''' Accelerometer sigma after AHRS averaging over the filter window. '''
        n = max(1.0, self.attitude_filter_window_s * self.sample_rate_hz)
        return self.accel_sigma() / sqrt(n)

    def gyro_sigma(self) -> float:
        ''' Per-sample gyroscope sigma [rad/s]. '''
        return self.gyro_noise_density * sqrt(self.sample_rate_hz)


@dataclass(frozen=True)
class TeleCameraSpec:
    ''' Representative iPhone-class periscope tele camera used for the sights.

        The daytime Sun and Moon disks (~0.53 deg) are well resolved on the tele
        lens; centroiding the disk to a fraction of a pixel gives an extremely
        precise *pointing* measurement -- so precise that the IMU attitude, not
        the camera, dominates the altitude error.  That is the central finding
        the numbers below are chosen to expose.
    '''
    width_px: int = 8064
    height_px: int = 6048
    arcsec_per_px: float = 9.0               # ~ tele plate scale (bare lens)
    centroid_px_sigma: float = 0.2           # sub-pixel disk-centre precision
    teleconverter: float = 1.0               # external afocal optic (e.g. 3x)

    def eff_arcsec_per_px(self) -> float:
        ''' Effective plate scale with any clip-on teleconverter: a 3x optic
            triples the focal length, so each pixel spans 1/3 the sky. '''
        return self.arcsec_per_px / self.teleconverter

    def pointing_sigma_arcmin(self) -> float:
        ''' Camera-only pointing sigma from disk centroiding [arc minutes]. '''
        return self.centroid_px_sigma * self.eff_arcsec_per_px() / 60.0


@dataclass(frozen=True)
class KinematicState:
    ''' The instantaneous handheld disturbance at a candidate shutter instant. '''
    ang_rate: float          # |omega|, body angular rate magnitude [rad/s]
    lin_accel: float         # residual linear acceleration magnitude [m/s^2]
    exposure_s: float = 0.008  # rolling-shutter / exposure window [s]


def gravity_tilt_sigma_arcmin(state: KinematicState, imu: IphoneImuSpec) -> float:
    ''' 1-sigma error of the AHRS-derived local vertical (the synthetic
        horizon), in arc minutes.

        Four contributions, added in quadrature:

          * a static floor (residual accel bias / calibration) -- the phone
            inclinometer accuracy floor;
          * filtered accelerometer white noise -> tilt ~ filtered_sigma / g;
          * residual linear acceleration masquerading as gravity: a horizontal
            specific force a tilts "down" by ~ a / g  (THE dominant term while
            the hand/platform is moving -- what the least-rotation shutter
            avoids);
          * rotation during the exposure smears the attitude by ~ omega * t_exp.
    '''
    tilt_floor = radians(imu.static_tilt_arcmin / 60.0)
    tilt_noise = imu.filtered_accel_sigma() / G0
    tilt_accel = state.lin_accel / G0
    tilt_rot = state.ang_rate * state.exposure_s
    tilt_rad = sqrt(tilt_floor ** 2 + tilt_noise ** 2 +
                    tilt_accel ** 2 + tilt_rot ** 2)
    return tilt_rad * ARCMIN_PER_RAD


def altitude_sigma_arcmin(state: KinematicState,
                          imu: IphoneImuSpec,
                          cam: TeleCameraSpec) -> float:
    ''' Total 1-sigma of a measured body ALTITUDE from this phone at this
        instant [arc minutes] = camera pointing (small) + synthetic-horizon tilt
        (dominant while moving), in quadrature.
    '''
    return sqrt(cam.pointing_sigma_arcmin() ** 2 +
                gravity_tilt_sigma_arcmin(state, imu) ** 2)


def heading_sigma_arcmin(state: KinematicState, imu: IphoneImuSpec,
                         base_heading_deg: float = 1.0) -> float:
    ''' 1-sigma of a measured body AZIMUTH [arc minutes].

        Azimuth needs absolute heading (magnetometer + gyro), which on a phone
        is far less certain than the gravity-derived vertical.  Modelled as a
        base magnetometer heading error inflated by rotation during exposure.
    '''
    base = base_heading_deg * 60.0
    rot = state.ang_rate * state.exposure_s * ARCMIN_PER_RAD
    return sqrt(base ** 2 + rot ** 2)


# Default instances used across the study.
DEFAULT_IMU = IphoneImuSpec()
DEFAULT_CAM = TeleCameraSpec()


def summarise() -> str:
    ''' Human-readable summary of the modelled sensor budget. '''
    imu, cam = DEFAULT_IMU, DEFAULT_CAM
    still = KinematicState(ang_rate=0.01, lin_accel=0.05)
    moving = KinematicState(ang_rate=0.5, lin_accel=1.5)
    lines = [
        "iPhone 17 Pro (representative) sensor budget",
        f"  camera pointing sigma      : {cam.pointing_sigma_arcmin():.3f} arcmin",
        f"  accel per-sample sigma     : {imu.accel_sigma():.4f} m/s^2",
        f"  gyro  per-sample sigma     : {degrees(imu.gyro_sigma()):.4f} deg/s",
        f"  altitude sigma (near-still): "
        f"{altitude_sigma_arcmin(still, imu, cam):.2f} arcmin",
        f"  altitude sigma (moving)    : "
        f"{altitude_sigma_arcmin(moving, imu, cam):.2f} arcmin",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(summarise())
