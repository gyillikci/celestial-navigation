''' Shared astronomy / geometry helpers for the iPhone IMU-fusion study.

    This module is the single source of truth for
      (a) a celestial body's geographic position (GP / sub-point) at a time, and
      (b) the predicted altitude and azimuth of that body from an observer.

    Both the synthetic-truth generator (`scenario.py`), the factor graph
    (`celestial_factor_graph.py`) and the least-squares baseline use the SAME
    functions here, so a zero-noise run recovers ground truth to machine
    precision.  The hard, accuracy-critical part -- the Sun/Moon ephemeris -- is
    reused verbatim from the repository's `starfix` engine (hourly GHA/Dec from
    the machine-readable nautical almanac, linearly interpolated).

    MODELLING NOTE.  The study works on a spherical Earth (radius
    `starfix.EARTH_RADIUS`).  Observer positions are geocentric lat/lon.  Because
    truth *and* estimate share this convention, the reported quantity -- the
    distance between the estimate and the truth -- is a pure estimation error and
    is unaffected by the (up to ~0.19 deg) geocentric-vs-geodetic offset of the
    real figure of the Earth.

    (c) 2026.  MIT License (see LICENSE file).  Part of the celestial-navigation
    project by August Linnman.
'''

from math import sin, cos, asin, atan2, radians, degrees, sqrt, pi

from starfix import (Sight, LatLonGeodetic, LatLonGeocentric, get_azimuth,
                     EARTH_RADIUS)

# A body's angular semidiameter is already removed by working with the disk
# centre; refraction/parallax are treated as pre-corrected ("observed
# altitude").  The study therefore uses purely geometric altitudes, which keeps
# the measurement model transparent and differentiable.


def body_gp(object_name: str, time_iso: str) -> LatLonGeocentric:
    ''' Return the geographic position (sub-point) of a body at a UTC time.

        Reuses `starfix.Sight`'s ephemeris interpolation.  The GP depends only on
        the interpolated GHA/Dec, not on the (dummy) measured altitude or any
        sextant correction, so those are chosen to be inert.
    '''
    dummy = Sight(object_name=object_name,
                  set_time=time_iso,
                  measured_alt="45:0:0",
                  estimated_position=LatLonGeodetic(0, 0),
                  ho_obs=True,               # skip refraction + dip
                  limb_correction=0,         # disk centre
                  horizontal_parallax=0)     # geometric (no topocentric shift)
    return dummy.get_gp()


def gp_dec_gha(gp: LatLonGeocentric) -> tuple[float, float]:
    ''' Decompose a GP into (declination, Greenwich hour angle) in degrees.

        `starfix` stores a GP as a LatLonGeocentric with lat == declination and
        lon == -(GHA + SHA); hence GHA = -lon.
    '''
    return gp.get_lat(), (-gp.get_lon()) % 360.0


def predicted_altitude(lat: float, lon: float, gp: LatLonGeocentric) -> float:
    ''' Geometric altitude (degrees) of a body (given its GP) from a geocentric
        observer position.  Classic navigation triangle:
            sin(Hc) = sin(L) sin(Dec) + cos(L) cos(Dec) cos(LHA)
    '''
    dec, gha = gp_dec_gha(gp)
    lha = radians(gha + lon)
    latr, decr = radians(lat), radians(dec)
    sin_hc = sin(latr) * sin(decr) + cos(latr) * cos(decr) * cos(lha)
    sin_hc = max(-1.0, min(1.0, sin_hc))
    return degrees(asin(sin_hc))


def predicted_azimuth(lat: float, lon: float, gp: LatLonGeocentric) -> float:
    ''' Geometric azimuth (degrees, 0..360, N=0 E=90) of a body from an
        observer.  Reuses `starfix.get_azimuth` for consistency with the rest of
        the toolkit.
    '''
    observer = LatLonGeocentric(lat, lon)
    return get_azimuth(gp, observer) % 360.0


def altaz(lat: float, lon: float, gp: LatLonGeocentric) -> tuple[float, float]:
    ''' Convenience: (altitude, azimuth) in degrees. '''
    return predicted_altitude(lat, lon, gp), predicted_azimuth(lat, lon, gp)


# --------------------------------------------------------------------------- #
# Local ENU tangent-plane <-> geographic conversions.
# The factor graph works in metres in a local East-North-Up frame anchored at a
# base geographic point; these helpers move between that frame and lat/lon.
# --------------------------------------------------------------------------- #

_R_M = EARTH_RADIUS * 1000.0   # Earth radius in metres


def enu_to_latlon(east_m: float, north_m: float,
                  lat0: float, lon0: float) -> tuple[float, float]:
    ''' Convert a local ENU horizontal offset (metres) to geocentric lat/lon. '''
    lat = lat0 + degrees(north_m / _R_M)
    lon = lon0 + degrees(east_m / (_R_M * cos(radians(lat0))))
    return lat, lon


def latlon_to_enu(lat: float, lon: float,
                  lat0: float, lon0: float) -> tuple[float, float]:
    ''' Convert geocentric lat/lon to a local ENU horizontal offset (metres). '''
    north_m = radians(lat - lat0) * _R_M
    east_m = radians(lon - lon0) * _R_M * cos(radians(lat0))
    return east_m, north_m


def great_circle_km(lat1: float, lon1: float,
                    lat2: float, lon2: float) -> float:
    ''' Great-circle distance (km) on the spherical Earth used by the study. '''
    p1, p2 = radians(lat1), radians(lat2)
    dl = radians(lon2 - lon1)
    a = sin((p2 - p1) / 2) ** 2 + cos(p1) * cos(p2) * sin(dl / 2) ** 2
    return 2 * EARTH_RADIUS * asin(min(1.0, sqrt(a)))
