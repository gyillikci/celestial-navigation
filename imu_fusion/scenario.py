''' Sea / land / air scenario generator for the Sun+Moon daytime sighting study.

    Builds ground truth (observer trajectory + true body geometry from the real
    ephemeris) and the noisy measurements a phone would actually record:

      * a sequence of KEYFRAMES (one per camera shot), each with the measured
        altitude (+ optional azimuth) of the Sun and Moon, corrupted by the
        synthetic-horizon error of `iphone_model` at that shutter's motion state;
      * a high-rate IMU stream (specific force + angular rate) between keyframes,
        for the factor graph's IMU preintegration.

    Three regimes:
      land : stationary braced observer (a coastal/handheld fix).
      sea  : vessel on a straight constant-speed leg, strong swell.
      air  : aircraft on a straight constant-speed leg, high altitude, vibration.

    All randomness flows through an injected random.Random for reproducibility
    (the codebase forbids wall-clock randomness in this context).

    (c) 2026.  MIT License (see LICENSE file).
'''

from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone

from starfix import LatLonGeocentric

from .astro import body_gp, altaz, enu_to_latlon, latlon_to_enu
from .iphone_model import (IphoneImuSpec, TeleCameraSpec, KinematicState,
                           DEFAULT_IMU, DEFAULT_CAM, G0,
                           heading_sigma_arcmin)
from .ultrawide_horizon import (UltrawideHorizonSpec, DEFAULT_UW,
                                horizon_reference_sigma_arcmin,
                                horizon_reference_sigma_lens)
from .optical_attitude import (parallactic_angle_deg, parallactic_sigma_deg,
                               optical_heading_sigma_deg, moon_limb_available)
from .visual_anchor import anchor_horizon_sigma_arcmin, coast_attitude_arcmin
from .cloud import CloudSpec, body_clear_flags
from .capture_trigger import PROFILES, simulate_trace, find_shutter_instants

# Canonical daytime epoch.  Greenwich (the home of longitude) near local noon,
# 2026-03-24 12:00 UTC, with a first-quarter Moon.  Both bodies are well up and,
# crucially, ~94 deg apart in azimuth (Sun alt 40/az 178, Moon alt 32/az 84), so
# the two altitude lines of position cross near-perpendicularly -- a
# well-conditioned Sun+Moon fix.  This matches the half/gibbous daytime Moon in
# a blue sky seen in the reference photos.
EPOCH = datetime(2026, 3, 24, 12, 0, 0, tzinfo=timezone.utc)


@dataclass(frozen=True)
class Observation:
    ''' One body measured at one keyframe. '''
    body: str
    time_iso: str
    gp: LatLonGeocentric
    true_alt: float
    true_az: float
    meas_alt: float
    meas_az: float
    alt_sigma_arcmin: float
    az_sigma_arcmin: float
    kin: KinematicState
    # Optical parallactic-angle observable from the tele-resolved disk.
    par_meas: float = None        # measured parallactic angle q (deg)
    par_sigma_deg: float = None
    par_valid: bool = False


@dataclass
class Keyframe:
    ''' One camera shot: an observer state plus the bodies seen. '''
    index: int
    time_s: float
    time_iso: str
    true_lat: float
    true_lon: float
    true_east: float
    true_north: float
    observations: list = field(default_factory=list)
    # IMU samples (accel_body[3], gyro_body[3], dt) integrated to reach the NEXT
    # keyframe.  Empty for the last keyframe.
    imu_to_next: list = field(default_factory=list)


@dataclass
class Scenario:
    ''' A complete synthetic data set for one regime/run. '''
    regime: str
    lat0: float               # ENU anchor (also the coarse-prior position)
    lon0: float
    base_lat: float           # true starting position
    base_lon: float
    vel_east: float           # true constant velocity, ENU [m/s]
    vel_north: float
    keyframes: list
    imu: IphoneImuSpec
    bodies: tuple


# Regime definitions: (start lat/lon, speed m/s, heading deg true, altitude note)
_REGIMES = {
    # stationary observer near Greenwich
    "land": dict(lat=51.5, lon=0.0, speed=0.0, heading=0.0),
    # ~19 kn vessel heading NE
    "sea":  dict(lat=51.5, lon=0.0, speed=10.0, heading=45.0),
    # ~390 kn aircraft heading E
    "air":  dict(lat=51.5, lon=0.0, speed=200.0, heading=90.0),
}


def _iso(t_s: float) -> str:
    ''' ISO-8601 UTC time string for `t_s` seconds after EPOCH. '''
    return (EPOCH + timedelta(seconds=t_s)).strftime("%Y-%m-%d %H:%M:%S+00:00")


def _velocity_enu(speed: float, heading_deg: float) -> tuple[float, float]:
    ''' Constant velocity in ENU [m/s] from speed and true heading. '''
    from math import sin, cos, radians
    hr = radians(heading_deg)
    return speed * sin(hr), speed * cos(hr)     # east, north


def build_scenario(regime: str,
                   rng,
                   n_shots: int = 12,
                   shot_interval_s: float = 10.0,
                   gated: bool = True,
                   bodies: tuple = ("Sun", "Moon"),
                   imu: IphoneImuSpec = DEFAULT_IMU,
                   cam: TeleCameraSpec = DEFAULT_CAM,
                   use_azimuth: bool = False,
                   noise_scale: float = 1.0,
                   dr_error_km: float = 30.0,
                   dr_bearing_deg: float = 60.0,
                   horizon_mode: str = "imu",
                   uw: UltrawideHorizonSpec = DEFAULT_UW,
                   heading_source: str = "mag",
                   use_parallactic: bool = False,
                   sun_spots: bool = False,
                   horizon_lens: str = None,
                   imu_anchor: bool = False,
                   anchor_track_s: float = 10.0,
                   anchor_pos_km: float = 3.0,
                   cloud: CloudSpec = None) -> Scenario:
    ''' Generate a full scenario: trajectory, true geometry, noisy measurements
        and IMU stream.

        gated=True selects calm shutter instants (least-rotation); gated=False
        uses a naive periodic shutter.
    '''
    cfg = _REGIMES[regime]
    base_lat, base_lon = cfg["lat"], cfg["lon"]
    ve, vn = _velocity_enu(cfg["speed"], cfg["heading"])

    # Coarse-prior / ENU anchor, offset from the TRUE start by a realistic
    # GPS-denied dead-reckoning error.  The estimator must use the celestial
    # measurements to correct this offset; simply trusting the prior yields a
    # large error, so the reported accuracy is honest.
    from math import sin, cos, radians
    br = radians(dr_bearing_deg)
    a_east = dr_error_km * 1000.0 * sin(br)
    a_north = dr_error_km * 1000.0 * cos(br)
    lat0, lon0 = enu_to_latlon(a_east, a_north, base_lat, base_lon)

    duration = n_shots * shot_interval_s + shot_interval_s
    prof = PROFILES[regime]
    t_trace, w_trace, a_trace = simulate_trace(prof, duration,
                                               imu.sample_rate_hz, rng)
    shutters = find_shutter_instants(t_trace, w_trace, a_trace,
                                     n_shots, shot_interval_s * 0.6, gated)

    # Cloud occlusion: per-body clear/obscured over the shot times, plus the
    # tracked-anchor body's clear runs (for coasting the attitude during an
    # outage).  The anchor tracks the first body that is up.
    shot_times = [t for t, _ in shutters]
    clear_flags = body_clear_flags(shot_times, bodies, cloud, rng)
    anchor_body = bodies[0]
    last_anchor_clear_t = shot_times[0] if shot_times else 0.0

    keyframes = []
    for k, (t_s, kin) in enumerate(shutters):
        # True observer position at this time.
        east = ve * t_s
        north = vn * t_s
        lat, lon = enu_to_latlon(east, north, base_lat, base_lon)
        iso = _iso(t_s)

        kf = Keyframe(index=k, time_s=t_s, time_iso=iso,
                      true_lat=lat, true_lon=lon,
                      true_east=east, true_north=north)

        # Anchor coast: if the tracked body is clear now, the anchor is live;
        # otherwise the attitude coasts on the gyro since it was last clear.
        if clear_flags[anchor_body][k]:
            last_anchor_clear_t = t_s
        anchor_coast_s = t_s - last_anchor_clear_t

        # Horizon reference sigma for THIS shot: IMU gravity, optical ultrawide
        # horizon, or their fusion.  The tele-camera pointing error adds in
        # quadrature to give the altitude measurement sigma.
        # Horizon reference sigma per shot; with a lens policy it depends on the
        # body altitude (the horizon must stay in the horizon lens's field).
        href_common = horizon_reference_sigma_arcmin(horizon_mode, kin, regime,
                                                     imu, uw)
        for body in bodies:
            # Cloud: a body obscured at this shot yields no sight (drop it).
            if not clear_flags[body][k]:
                continue
            gp = body_gp(body, iso)
            t_alt, t_az = altaz(lat, lon, gp)
            if horizon_lens is None:
                href = href_common
            else:
                href = horizon_reference_sigma_lens(horizon_mode, kin, regime,
                                                    t_alt, horizon_lens, imu, uw)
            # A tracked sunspot / disk-feature anchor is an extra, acceleration-
            # immune, altitude-independent attitude reference: fuse it in
            # (inverse-variance).  It rescues the cases where the optical horizon
            # is unavailable (high sights, land, out of frame -> href ~ IMU) and
            # modestly sharpens the rest.
            if imu_anchor:
                if anchor_coast_s <= 0.0:
                    s_anchor = anchor_horizon_sigma_arcmin(anchor_pos_km,
                                                           anchor_track_s, kin,
                                                           imu, cam)
                else:
                    # Tracked body clouded: the anchor coasts on the gyro.
                    s_anchor = coast_attitude_arcmin(anchor_coast_s, imu,
                                                     calibrated=True)
                href = 1.0 / ((1.0 / href ** 2 + 1.0 / s_anchor ** 2) ** 0.5)
            a_sig = (cam.pointing_sigma_arcmin() ** 2 + href ** 2) ** 0.5
            meas_alt = t_alt + noise_scale * rng.gauss(0.0, a_sig / 60.0)

            # Azimuth line of position: heading from magnetometer or from the
            # optical disk orientation (magnetometer-free).
            if heading_source == "optical":
                h_sig = optical_heading_sigma_deg(body, kin) * 60.0
            else:
                h_sig = heading_sigma_arcmin(kin, imu)
            meas_az = (t_az + noise_scale * rng.gauss(0.0, h_sig / 60.0)
                       if use_azimuth else t_az)

            # Optical parallactic-angle observable from the resolved disk.
            par_meas = par_sig = None
            par_valid = False
            if use_parallactic:
                avail = (moon_limb_available(iso) if body.lower() == "moon"
                         else sun_spots)
                if avail:
                    q_true = parallactic_angle_deg(lat, lon, gp)
                    par_sig = parallactic_sigma_deg(body, kin, href, cam)
                    par_meas = q_true + noise_scale * rng.gauss(0.0, par_sig)
                    par_valid = True

            kf.observations.append(Observation(
                body=body, time_iso=iso, gp=gp,
                true_alt=t_alt, true_az=t_az,
                meas_alt=meas_alt, meas_az=meas_az,
                alt_sigma_arcmin=a_sig, az_sigma_arcmin=h_sig, kin=kin,
                par_meas=par_meas, par_sigma_deg=par_sig, par_valid=par_valid))
        keyframes.append(kf)

    # IMU stream between consecutive keyframes.  With a constant, known level
    # attitude (body frame == nav ENU frame) and straight constant-speed motion,
    # the true specific force is just gravity reaction (0,0,+g) and the true
    # angular rate is zero; we add bias + white noise.  This ties consecutive
    # positions through the velocity state without over-modelling manoeuvres.
    ab = (rng.gauss(0, imu.accel_bias), rng.gauss(0, imu.accel_bias),
          rng.gauss(0, imu.accel_bias))
    gb = (rng.gauss(0, imu.gyro_bias), rng.gauss(0, imu.gyro_bias),
          rng.gauss(0, imu.gyro_bias))
    dt = 1.0 / imu.sample_rate_hz
    asig, gsig = imu.accel_sigma(), imu.gyro_sigma()
    for k in range(len(keyframes) - 1):
        t_a = keyframes[k].time_s
        t_b = keyframes[k + 1].time_s
        n = max(1, int(round((t_b - t_a) * imu.sample_rate_hz)))
        samples = []
        for _ in range(n):
            acc = (ab[0] + rng.gauss(0, asig),
                   ab[1] + rng.gauss(0, asig),
                   G0 + ab[2] + rng.gauss(0, asig))
            gyr = (gb[0] + rng.gauss(0, gsig),
                   gb[1] + rng.gauss(0, gsig),
                   gb[2] + rng.gauss(0, gsig))
            samples.append((acc, gyr, dt))
        keyframes[k].imu_to_next = samples

    return Scenario(regime=regime, lat0=lat0, lon0=lon0,
                    base_lat=base_lat, base_lon=base_lon,
                    vel_east=ve, vel_north=vn, keyframes=keyframes,
                    imu=imu, bodies=bodies)


def scenario_summary(sc: Scenario) -> str:
    ''' One-line-per-keyframe human summary (truth vs measured altitude). '''
    lines = [f"[{sc.regime}] {len(sc.keyframes)} keyframes, "
             f"bodies={sc.bodies}, v_enu=({sc.vel_east:.1f},{sc.vel_north:.1f}) m/s"]
    for kf in sc.keyframes:
        parts = [f"  kf{kf.index:2d} t={kf.time_s:6.1f}s "
                 f"pos=({kf.true_lat:.4f},{kf.true_lon:.4f})"]
        for o in kf.observations:
            parts.append(f"{o.body}: Ho={o.meas_alt:6.2f} "
                         f"(true {o.true_alt:6.2f}, sig {o.alt_sigma_arcmin:5.1f}')")
        lines.append("  ".join(parts))
    return "\n".join(lines)


if __name__ == "__main__":
    import random
    for r in ("land", "sea", "air"):
        sc = build_scenario(r, random.Random(1), n_shots=4)
        print(scenario_summary(sc))
        print()
