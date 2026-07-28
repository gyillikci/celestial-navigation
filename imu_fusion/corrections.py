''' Apparent-altitude corrections: refraction + topocentric parallax + semidiameter.

    The study, the factor graph and `astro.predicted_altitude` all work in
    GEOMETRIC, GEOCENTRIC altitudes (the clean navigation triangle from a body's
    sub-point).  A phone on the ground, however, measures the APPARENT altitude of
    the disk centre: the light has been bent upward by the atmosphere
    (refraction) and the observer sits on the Earth's surface, not at its centre,
    so a near body (the Moon) appears LOWER by up to ~1 degree (topocentric
    parallax).  This module is the single source of truth for the three transforms
    that bridge measured <-> geometric, and it is mirrored verbatim by the iOS app
    (`ios/CelestialFixMVP/Sources/CelestialMath.swift`).

    MODELS — the same standard closed forms Stellarium / the Nautical Almanac use:
      * Refraction  — Bennett (1982), reused from `starfix.get_refraction` so the
        study and the app share one implementation.  Defined on the APPARENT
        altitude; validated against IAU ERFA `refco` in `validate_ephemeris.py`.
      * Parallax    — exact geometry from the body's distance:
        sin(HP) = R_earth / d,   p(alt) = asin( sin(HP) * cos(alt_geometric) ).
      * Semidiameter — geometry from the distance.  The disk-metrology stage fits
        the FULL circle and returns the centre, so the app needs NO limb
        correction; `semidiameter_deg` is provided only for reducing a
        limb-tangent sight (e.g. a hand sextant) and for reporting.

    DIRECTION.  Refraction raises a body; parallax (for the Moon) lowers it.  So

        measured apparent  --(- refraction)-->  topocentric geometric
                           --(+ parallax)   -->  geocentric geometric  (the triangle)

    and `apparent_from_geometric` is the exact inverse of `geometric_from_apparent`.

    (c) 2026.  MIT License (see LICENSE file).  Part of the celestial-navigation
    project.
'''

from math import sin, cos, asin, tan, radians, degrees

from starfix import get_refraction, EARTH_RADIUS   # EARTH_RADIUS in km

# Standard-atmosphere defaults (Bennett).  The app can override from the device
# barometer/thermometer; at altitudes above ~10 deg the sensitivity is small.
DEFAULT_TEMPERATURE_C = 10.0
DEFAULT_PRESSURE_KPA = 101.0
DEFAULT_HUMIDITY_PCT = 50.0

# Body mean radii (km) for the semidiameter helper.
_BODY_RADIUS_KM = {"moon": 1737.4, "sun": 695700.0}
_MEAN_SUN_DISTANCE_KM = 1.495978707e8


# --------------------------------------------------------------------------- #
# Parallax
# --------------------------------------------------------------------------- #

def horizontal_parallax_deg(distance_km: float) -> float:
    ''' Equatorial horizontal parallax HP (degrees): sin(HP) = R_earth / d.

        The Moon's HP is ~0.95 deg (57'); the Sun's is ~0.0024 deg (8.8").
    '''
    return degrees(asin(EARTH_RADIUS / distance_km))


def parallax_in_altitude_deg(geometric_alt_deg: float, distance_km: float) -> float:
    ''' Parallax in altitude p (degrees) at a given GEOMETRIC altitude:

            p = asin( sin(HP) * cos(alt) )

        Zero at the zenith, maximal (= HP) at the horizon.  A body is seen this
        much LOWER from the surface than from the geocentre.
    '''
    hp = radians(horizontal_parallax_deg(distance_km))
    c = sin(hp) * cos(radians(geometric_alt_deg))
    return degrees(asin(max(-1.0, min(1.0, c))))


def distance_km_from_hp_deg(hp_deg: float) -> float:
    ''' Inverse of `horizontal_parallax_deg` — body distance from its HP. '''
    return EARTH_RADIUS / sin(radians(hp_deg))


# --------------------------------------------------------------------------- #
# Refraction
# --------------------------------------------------------------------------- #

def refraction_deg(apparent_alt_deg: float,
                   temperature_c: float = DEFAULT_TEMPERATURE_C,
                   pressure_kpa: float = DEFAULT_PRESSURE_KPA,
                   humidity_pct: float = DEFAULT_HUMIDITY_PCT) -> float:
    ''' Atmospheric refraction (degrees), Bennett's formula, as a function of the
        APPARENT altitude.  Wraps `starfix.get_refraction` (arc-minutes) so the
        study and the app share one model.  Clamped to zero below the horizon.
    '''
    if apparent_alt_deg <= -1.0:
        return 0.0
    return get_refraction(apparent_alt_deg, temperature_c,
                          pressure_kpa, humidity_pct) / 60.0


# --------------------------------------------------------------------------- #
# Semidiameter (reporting / limb-tangent sights only — the app uses disk centre)
# --------------------------------------------------------------------------- #

def semidiameter_deg(distance_km: float, body: str = "moon") -> float:
    ''' Angular semidiameter (degrees) of a body at a given distance. '''
    r = _BODY_RADIUS_KM.get(body.lower(), _BODY_RADIUS_KM["moon"])
    return degrees(asin(r / distance_km))


# --------------------------------------------------------------------------- #
# The two bridges (exact inverses of one another)
# --------------------------------------------------------------------------- #

def geometric_from_apparent(apparent_alt_deg: float, distance_km: float,
                            temperature_c: float = DEFAULT_TEMPERATURE_C,
                            pressure_kpa: float = DEFAULT_PRESSURE_KPA,
                            humidity_pct: float = DEFAULT_HUMIDITY_PCT) -> float:
    ''' Reduce a MEASURED apparent altitude (disk centre, through the atmosphere,
        from the surface) to the GEOMETRIC GEOCENTRIC altitude that the navigation
        triangle / `astro.predicted_altitude` predicts.  This is the transform the
        device applies to a sight before the position fix.
    '''
    topo = apparent_alt_deg - refraction_deg(apparent_alt_deg, temperature_c,
                                             pressure_kpa, humidity_pct)
    return topo + parallax_in_altitude_deg(topo, distance_km)


def apparent_from_geometric(geocentric_alt_deg: float, distance_km: float,
                            temperature_c: float = DEFAULT_TEMPERATURE_C,
                            pressure_kpa: float = DEFAULT_PRESSURE_KPA,
                            humidity_pct: float = DEFAULT_HUMIDITY_PCT) -> float:
    ''' Forward model: predict the APPARENT altitude a phone would measure, given
        the geometric geocentric altitude from the ephemeris.  Exact inverse of
        `geometric_from_apparent` (used to synthesise sights and golden vectors).
    '''
    # geocentric -> topocentric geometric: solve topo + p(topo) = geocentric.
    topo = geocentric_alt_deg - parallax_in_altitude_deg(geocentric_alt_deg, distance_km)
    for _ in range(3):
        topo = geocentric_alt_deg - parallax_in_altitude_deg(topo, distance_km)
    # topocentric geometric -> apparent: solve app - refraction(app) = topo.
    app = topo
    for _ in range(4):
        app = topo + refraction_deg(app, temperature_c, pressure_kpa, humidity_pct)
    return app


def total_correction_deg(geocentric_alt_deg: float, distance_km: float,
                         **kw) -> float:
    ''' apparent - geocentric (degrees): the net of refraction (up) and parallax
        (down).  For the Moon near the horizon parallax dominates (net negative,
        ~ -0.9 deg); high up, refraction dominates (small positive). '''
    return apparent_from_geometric(geocentric_alt_deg, distance_km, **kw) - geocentric_alt_deg
