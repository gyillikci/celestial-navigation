''' Landfall: fixing position from distant land features once they rise over the
    horizon.

    Context.  The celestial stack gives a few-km fix with no horizon and no GPS.
    When land appears, a MUCH stronger class of observable becomes available,
    because the landmarks are at KNOWN, surveyed positions -- a peak is a fixed
    point on the Earth, unlike a body whose sub-point moves at 0.25 km/s.  This
    module models the three classical landfall observables, their error budgets,
    and how they attach to the existing factor graph.

    THE THREE OBSERVABLES (in increasing order of accuracy):

      1. BEARING to one identified peak -> a line of position through it.
         Limited by the phone's heading (magnetometer ~1-2 deg), which is the
         weakest sensor in the whole system.  Cross-range error = d * sigma.

      2. VERTICAL ANGLE to a summit of known height -> a RANGE (distance circle).
         The angle is small (a 2000 m peak at 100 km subtends ~0.7 deg above the
         horizontal) and it is measured against either the IMU synthetic horizon
         or -- much better at sea -- the visible sea horizon.  Its accuracy is
         capped not by the sensor but by ATMOSPHERIC REFRACTION: the coefficient
         k varies ~0.10-0.20 with the air/sea temperature profile, and the range
         scales as 1/sqrt(1-k), so a "dipping" range carries a few-percent
         systematic error that no better instrument removes.

      3. HORIZONTAL ANGLE between TWO identified peaks -> a circle of position
         (the classical horizontal-sextant-angle fix; the observer lies on the
         circle through both landmarks subtending that angle).  Two such angles
         from three peaks give a fix.  This is the ACCURATE one, for two reasons:
         it needs no compass (the angle is differential, so heading error
         cancels), and lateral refraction is far smaller than vertical, so the
         measurement is essentially unbiased.  With the tele camera the angle is
         pixel-limited, i.e. arc-seconds.

    So the ranking is the opposite of the intuitive one: the horizontal angles
    between peaks -- not the range to them -- carry the position information.

    Geometry conventions: spherical Earth of radius `EARTH_RADIUS` (as elsewhere
    in the study), refraction folded in through an effective radius
    R_eff = R / (1 - k).

    (c) 2026.  MIT License (see LICENSE file).
'''

from math import sqrt, tan, sin, cos, asin, atan2, radians, degrees, fabs

from starfix import EARTH_RADIUS, LatLonGeodetic, get_terrestrial_position

# Standard terrestrial refraction coefficient and its realistic spread.
# k ~ 0.13 is the textbook mean; it ranges ~0.10-0.20 in ordinary conditions and
# can invert entirely over a cold sea (superior mirage / looming).
K_REFRACTION = 0.13
K_SPREAD = 0.05

_R_KM = EARTH_RADIUS                 # spherical Earth radius, km


def effective_radius_km(k: float = K_REFRACTION) -> float:
    ''' Refraction-inflated Earth radius: light bends toward the surface, which
        is equivalent to a larger, less-curved Earth. '''
    return _R_KM / (1.0 - k)


# --------------------------------------------------------------------------- #
# 1. Visibility -- when does a peak break the horizon?
# --------------------------------------------------------------------------- #

def horizon_distance_km(height_m: float, k: float = K_REFRACTION) -> float:
    ''' Distance from an object of height `height_m` to its own sea horizon. '''
    if height_m <= 0:
        return 0.0
    return sqrt(2.0 * effective_radius_km(k) * (height_m / 1000.0))


def geographic_range_km(eye_height_m: float, peak_height_m: float,
                        k: float = K_REFRACTION) -> float:
    ''' The "geographic range": the distance at which a summit of height
        `peak_height_m` is exactly ON the horizon for an observer at
        `eye_height_m`.  Beyond this it is hidden by the Earth's bulge.

        This is the classic d = 3.83*(sqrt(h) + sqrt(H)) km rule (h, H in metres).
    '''
    return (horizon_distance_km(eye_height_m, k)
            + horizon_distance_km(peak_height_m, k))


def geographic_range_spread_km(eye_height_m: float, peak_height_m: float,
                               k_spread: float = K_SPREAD) -> float:
    ''' Half-spread of the geographic range caused by refraction uncertainty.

        This is the FLOOR on a "dipping range" fix: the moment a peak rises over
        the horizon gives a range circle, but that circle's radius is only known
        to this accuracy, however precisely the event is timed. '''
    hi = geographic_range_km(eye_height_m, peak_height_m, K_REFRACTION + k_spread)
    lo = geographic_range_km(eye_height_m, peak_height_m, K_REFRACTION - k_spread)
    return (hi - lo) / 2.0


# --------------------------------------------------------------------------- #
# 2. Vertical angle -> range
# --------------------------------------------------------------------------- #

def vertical_angle_deg(distance_km: float, peak_height_m: float,
                       eye_height_m: float = 2.0,
                       k: float = K_REFRACTION) -> float:
    ''' Apparent elevation of a summit above the ASTRONOMICAL horizontal (the
        gravity level, i.e. what the IMU gives), in degrees:

            alpha ~ (H - h)/d  -  d / (2 R_eff)

        The first term is the height difference over the distance; the second is
        the curvature drop of the far point.  Negative means the summit is below
        the horizontal (but it may still be visible above the sea horizon, which
        is itself dipped). '''
    if distance_km <= 0:
        return 90.0
    dh_km = (peak_height_m - eye_height_m) / 1000.0
    return degrees(dh_km / distance_km - distance_km / (2.0 * effective_radius_km(k)))


def range_from_vertical_angle_km(alpha_deg: float, peak_height_m: float,
                                 eye_height_m: float = 2.0,
                                 k: float = K_REFRACTION) -> float:
    ''' Invert `vertical_angle_deg` for the distance (the positive root). '''
    a = radians(alpha_deg)
    R = effective_radius_km(k)
    dh_km = (peak_height_m - eye_height_m) / 1000.0
    disc = a * a + 2.0 * dh_km / R
    if disc < 0:
        return float("nan")
    return R * (-a + sqrt(disc))


def range_sigma_km(distance_km: float, peak_height_m: float,
                   sigma_alpha_deg: float, eye_height_m: float = 2.0,
                   k: float = K_REFRACTION) -> float:
    ''' Range uncertainty from a vertical-angle measurement of 1-sigma
        `sigma_alpha_deg`.  Propagates through d(alpha)/d(distance). '''
    R = effective_radius_km(k)
    dh_km = (peak_height_m - eye_height_m) / 1000.0
    dalpha_dd = -dh_km / (distance_km ** 2) - 1.0 / (2.0 * R)   # rad per km
    if dalpha_dd == 0:
        return float("inf")
    return fabs(radians(sigma_alpha_deg) / dalpha_dd)


def range_bias_from_refraction_km(distance_km: float, peak_height_m: float,
                                  eye_height_m: float = 2.0,
                                  k_spread: float = K_SPREAD) -> float:
    ''' Systematic range error from not knowing k, at fixed measured angle.
        This is a BIAS, not noise -- averaging more sights does not reduce it. '''
    alpha = vertical_angle_deg(distance_km, peak_height_m, eye_height_m)
    hi = range_from_vertical_angle_km(alpha, peak_height_m, eye_height_m,
                                      K_REFRACTION + k_spread)
    lo = range_from_vertical_angle_km(alpha, peak_height_m, eye_height_m,
                                      K_REFRACTION - k_spread)
    return fabs(hi - lo) / 2.0


# --------------------------------------------------------------------------- #
# 3. Bearing and horizontal angle
# --------------------------------------------------------------------------- #

def cross_range_sigma_km(distance_km: float, sigma_bearing_deg: float) -> float:
    ''' Position error across the line of sight from a bearing of 1-sigma
        `sigma_bearing_deg` to a landmark at `distance_km`. '''
    return distance_km * radians(sigma_bearing_deg)


def horizontal_angle_sigma_deg(pixel_sigma: float, arcsec_per_pixel: float,
                               identification_sigma_m: float = 0.0,
                               distance_km: float = 50.0) -> float:
    ''' 1-sigma of a horizontal angle between two peaks measured in one image.

        Two error sources: locating each summit in the frame (pixels, hence the
        sqrt(2) for two peaks) and knowing WHICH point of the massif the surveyed
        coordinate refers to (`identification_sigma_m` of lateral position at
        `distance_km`).  The latter usually dominates in practice. '''
    ang = sqrt(2.0) * pixel_sigma * arcsec_per_pixel / 3600.0
    ident = degrees(identification_sigma_m / 1000.0 / distance_km) if distance_km > 0 else 0.0
    return sqrt(ang ** 2 + ident ** 2)


def two_landmark_circle_fix(p1, p2, angle_12, p3, angle_23, estimate=None):
    ''' Classical horizontal-angle fix: two subtended angles between three
        surveyed landmarks give the observer's position, with NO compass.

        Wraps `starfix.get_terrestrial_position` (circumscribed-circle method),
        which is the repository's own terrestrial engine.

        p1, p2, p3 : LatLonGeodetic of the identified peaks (p2 is shared).
        angle_12   : observed horizontal angle between p1 and p2 (degrees).
        angle_23   : observed horizontal angle between p2 and p3 (degrees).
        estimate   : LatLonGeodetic DR position, to pick the correct intersection.
    '''
    return get_terrestrial_position(p1, p2, angle_12, p2, p3, angle_23,
                                    estimated_position=estimate)


def horizontal_angle_fix_sigma_km(distance_km: float, sigma_angle_deg: float,
                                  crossing_deg: float = 90.0) -> float:
    ''' Rough position sigma from a horizontal-angle (circle-of-position) fix.

        For an angle subtended by landmarks at range d, a 1-sigma angle error
        displaces the circle of position by ~ d * sigma / sin(subtended), and two
        circles crossing at `crossing_deg` give the usual 1/sin factor. '''
    s = fabs(sin(radians(crossing_deg)))
    if s < 1e-6:
        return float("inf")
    return distance_km * radians(sigma_angle_deg) * sqrt(2.0) / s


# --------------------------------------------------------------------------- #
# Detection: is this landmark usable right now?
# --------------------------------------------------------------------------- #

def visible(distance_km: float, peak_height_m: float, eye_height_m: float = 2.0,
            k: float = K_REFRACTION) -> bool:
    ''' Is the summit above the observer's sea horizon? '''
    return distance_km <= geographic_range_km(eye_height_m, peak_height_m, k)


def apparent_height_above_horizon_deg(distance_km: float, peak_height_m: float,
                                      eye_height_m: float = 2.0,
                                      k: float = K_REFRACTION) -> float:
    ''' How far the summit stands above the visible SEA HORIZON line (degrees).

        At sea this is the natural measurement: the horizon is a sharp, gravity-
        referenced datum in the same frame as the peak, so it beats the IMU. '''
    from starfix import get_dip_of_horizon
    dip_deg = get_dip_of_horizon(eye_height_m) / 60.0        # arcmin -> deg
    return vertical_angle_deg(distance_km, peak_height_m, eye_height_m, k) + dip_deg
