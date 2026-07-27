''' Optical attitude from the TELE-resolved disk of the Sun and Moon.

    The tele lens does not just centroid the body -- it resolves the disk.  The
    orientation of disk features carries an absolute celestial reference:

      * Moon  -- the bright limb points at the Sun; its position angle chi (from
        celestial north) is pure ephemeris.  The terminator / maria give the same
        axis even when the limb is faint.
      * Sun   -- a sunspot at known heliographic coordinates, or the disk axis via
        the tabulated P / B0 angles, fixes the solar rotation axis's position
        angle P (ephemeris).

    Sunspot reference (deployment).  In practice the Sun's disk orientation is
    read through a solar filter by matching the live sunspot pattern to a
    REFERENCE image taken from an observatory just before the journey.  Over a
    short trip with a fresh reference the pattern is a stable fiducial, so the
    Sun's differential rotation is ignored here (set `sun_spots=True` in the
    scenario to use it).  Its value is mainly redundancy: it carries the heading
    when the Moon is near full or not visible, and gives a stronger parallactic
    line when the Sun is well off the meridian.

    Measured against the gravity vertical in the image, the feature axis gives

        theta_image = PA_ephemeris - q                                   (1)

    where q is the PARALLACTIC ANGLE at the body,

        q = atan2( sin H , tan(phi) cos(delta) - sin(delta) cos H ),     H = GHA + lon.

    Because PA is ephemeris (independent of the observer), it cancels between the
    synthetic truth and the graph's prediction -- so the real observable is
    q(lat, lon).  That buys two things the altitude sights cannot:

      1. a magnetometer-free HEADING (the same reference that orients the disk
         orients the platform), robust on steel hulls and in avionics;
      2. an independent, horizon-free POSITION line through q.

    Sign / frame notes.  Right-ascension differences are hour-angle differences:
    (alpha_sun - alpha_moon) = (GHA_moon - GHA_sun).  Longitudes are east-positive
    (as in starfix GP: lon = -GHA), so H = GHA + lon, matching astro.py.

    (c) 2026.  MIT License (see LICENSE file).
'''

from dataclasses import dataclass
from math import sin, cos, tan, atan2, acos, radians, degrees, sqrt

from .astro import body_gp, gp_dec_gha
from .iphone_model import TeleCameraSpec, DEFAULT_CAM, KinematicState, \
    ARCMIN_PER_RAD

# Angular radii of the disks (degrees); good to ~2% for this study.
MOON_RADIUS_DEG = 0.259
SUN_RADIUS_DEG = 0.266


def _sun_moon_dec_gha(time_iso: str):
    ds, gs = gp_dec_gha(body_gp("Sun", time_iso))
    dm, gm = gp_dec_gha(body_gp("Moon", time_iso))
    return ds, gs, dm, gm


def moon_elongation_deg(time_iso: str) -> float:
    ''' Sun-Moon angular separation (elongation), degrees. '''
    ds, gs, dm, gm = _sun_moon_dec_gha(time_iso)
    dsr, dmr = radians(ds), radians(dm)
    dgha = radians(gm - gs)
    c = sin(dsr) * sin(dmr) + cos(dsr) * cos(dmr) * cos(dgha)
    return degrees(acos(max(-1.0, min(1.0, c))))


def moon_illuminated_fraction(time_iso: str) -> float:
    ''' Fraction of the Moon's disk that is lit (0..1). '''
    return (1.0 + cos(radians(moon_elongation_deg(time_iso)))) / 2.0


def bright_limb_pa_deg(time_iso: str) -> float:
    ''' Position angle of the Moon's bright limb, from celestial north (deg).
        Meeus, Astronomical Algorithms, ch. 48. '''
    ds, gs, dm, gm = _sun_moon_dec_gha(time_iso)
    dsr, dmr = radians(ds), radians(dm)
    dalpha = radians(gm - gs)                    # alpha_sun - alpha_moon
    y = cos(dsr) * sin(dalpha)
    x = sin(dsr) * cos(dmr) - cos(dsr) * sin(dmr) * cos(dalpha)
    return degrees(atan2(y, x)) % 360.0


def parallactic_angle_deg(lat: float, lon: float, gp) -> float:
    ''' Parallactic angle q at a body seen from (lat, lon), degrees.
        q is the angle at the body between the directions to the zenith and to
        the celestial north pole. '''
    dec, gha = gp_dec_gha(gp)
    h = radians(gha + lon)
    latr, decr = radians(lat), radians(dec)
    y = sin(h)
    x = tan(latr) * cos(decr) - sin(decr) * cos(h)
    return degrees(atan2(y, x))


def moon_limb_available(time_iso: str,
                        lo: float = 0.03, hi: float = 0.98) -> bool:
    ''' The bright limb / terminator gives a usable axis except very near new or
        full Moon. '''
    k = moon_illuminated_fraction(time_iso)
    return lo < k < hi


@dataclass(frozen=True)
class OpticalDiskSpec:
    ''' How well the tele image pins the disk-feature orientation. '''
    # Floor on feature-axis orientation from model/libration/limb-fit residual.
    moon_axis_floor_deg: float = 0.35
    sun_axis_floor_deg: float = 0.60        # sunspot centroiding, when present
    edge_px_sigma: float = 0.5              # per-edge-point localisation


DEFAULT_DISK = OpticalDiskSpec()


def orientation_sigma_deg(body: str, state: KinematicState,
                          cam: TeleCameraSpec = DEFAULT_CAM,
                          disk: OpticalDiskSpec = DEFAULT_DISK) -> float:
    ''' 1-sigma of the measured feature-axis orientation in the image [deg].

        Fit of an axis across the disk rim: the longer the lever arm (disk
        radius in pixels), the better the angle.  A per-body floor accounts for
        libration / P-angle model error and spot-position uncertainty; rotation
        during the exposure smears the axis.
    '''
    r_deg = MOON_RADIUS_DEG if body.lower() == "moon" else SUN_RADIUS_DEG
    r_px = r_deg * 3600.0 / cam.eff_arcsec_per_px()
    # Axis-fit sigma ~ edge sigma / lever arm, over ~ (rim points) samples.
    n = max(1.0, r_px)                        # ~1 sample per pixel of rim length
    fit_deg = degrees(disk.edge_px_sigma / r_px) / sqrt(n)
    floor = (disk.moon_axis_floor_deg if body.lower() == "moon"
             else disk.sun_axis_floor_deg)
    rot_deg = degrees(state.ang_rate * state.exposure_s)
    return sqrt(fit_deg ** 2 + floor ** 2 + rot_deg ** 2)


def parallactic_sigma_deg(body: str, state: KinematicState,
                          horizon_sigma_arcmin: float,
                          cam: TeleCameraSpec = DEFAULT_CAM,
                          disk: OpticalDiskSpec = DEFAULT_DISK) -> float:
    ''' 1-sigma of the parallactic angle recovered from one disk [deg].

        q = PA - theta_image - roll: the error is the orientation-fit error and
        the error of the vertical/roll reference (the horizon reference sigma),
        in quadrature.
    '''
    s_orient = orientation_sigma_deg(body, state, cam, disk)
    s_roll = horizon_sigma_arcmin / 60.0
    return sqrt(s_orient ** 2 + s_roll ** 2)


def optical_heading_sigma_deg(body: str, state: KinematicState,
                              cam: TeleCameraSpec = DEFAULT_CAM,
                              disk: OpticalDiskSpec = DEFAULT_DISK) -> float:
    ''' Heading (azimuth-reference) sigma from the disk orientation [deg].
        The absolute celestial axis in the image yields the platform heading
        without a magnetometer; precision ~ the orientation-fit sigma. '''
    return orientation_sigma_deg(body, state, cam, disk)


def summarise(time_iso: str = "2026-03-24 12:00:00") -> str:
    k = moon_illuminated_fraction(time_iso)
    chi = bright_limb_pa_deg(time_iso)
    st = KinematicState(0.05, 0.2)
    lines = [f"Optical attitude @ {time_iso}",
             f"  Moon illuminated fraction : {k:.2f}  "
             f"(limb usable: {moon_limb_available(time_iso)})",
             f"  Moon bright-limb PA       : {chi:.1f} deg from N",
             f"  Moon orientation sigma    : "
             f"{orientation_sigma_deg('Moon', st):.2f} deg",
             f"  Sun  orientation sigma    : "
             f"{orientation_sigma_deg('Sun', st):.2f} deg (needs spots)"]
    for lat, lon, name in [(51.5, 0.0, "Greenwich")]:
        qm = parallactic_angle_deg(lat, lon, body_gp("Moon", time_iso))
        qs = parallactic_angle_deg(lat, lon, body_gp("Sun", time_iso))
        lines.append(f"  Parallactic q @ {name}: Moon {qm:.1f} deg, "
                     f"Sun {qs:.1f} deg")
    return "\n".join(lines)


if __name__ == "__main__":
    print(summarise())
