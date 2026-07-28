''' Mission-planning forecast: what the sky will offer along a route, and what
    position accuracy to expect from it.

    Motivating question: flying A -> B (e.g. Istanbul -> Ankara) in a GPS-denied
    aircraft, does pre-computing the Sun/Moon sky positions at each time along the
    route actually help?

    It does, in a specific and quantifiable way.  The ephemeris itself is
    DETERMINISTIC -- knowing where the Sun will be does not make a sight more
    accurate.  What the forecast buys is OBSERVATION GEOMETRY, known before
    take-off:

      * whether each body is up at all (and, for the Moon, lit enough for a limb),
      * the azimuth separation between the two bodies -- the crossing angle of
        their lines of position, which dominates fix quality (a "celestial DOP"),
      * each body's altitude, which sets the refraction/parallax regime and warns
        about the ill-conditioned near-zenith case,
      * therefore WHEN along the leg the fix will be strong and when it will be
        degenerate -- i.e. when to schedule the shots.

    This module computes that forecast from the project's authoritative ephemeris
    (Stellarium when an export is present -- see `stellarium_source.py` -- else the
    starfix almanac; they agree to well under an arc-minute, so the forecast is
    source-independent), and converts the geometry into a predicted position-error
    scale via the standard two-LOP crossing formula.

    (c) 2026.  MIT License (see LICENSE file).
'''

from datetime import datetime, timedelta, timezone
from math import (sin, cos, tan, asin, acos, atan2, radians, degrees, sqrt,
                  hypot, fabs)

from .astro import body_gp, altaz, gp_dec_gha, body_distance_km, great_circle_km
from . import corrections as C

_R_KM = 6371.0


# --------------------------------------------------------------------------- #
# Route
# --------------------------------------------------------------------------- #

class Leg:
    ''' A great-circle leg flown at constant ground speed. '''

    def __init__(self, name_a, lat_a, lon_a, name_b, lat_b, lon_b,
                 speed_kmh: float = 800.0):
        self.name_a, self.lat_a, self.lon_a = name_a, lat_a, lon_a
        self.name_b, self.lat_b, self.lon_b = name_b, lat_b, lon_b
        self.speed_kmh = speed_kmh
        self.distance_km = great_circle_km(lat_a, lon_a, lat_b, lon_b)
        self.duration_h = self.distance_km / speed_kmh

    def position_at(self, frac: float):
        ''' Great-circle interpolation, frac in [0,1] -> (lat, lon) degrees. '''
        f = max(0.0, min(1.0, frac))
        p1, l1 = radians(self.lat_a), radians(self.lon_a)
        p2, l2 = radians(self.lat_b), radians(self.lon_b)
        d = self.distance_km / _R_KM
        if d < 1e-12:
            return self.lat_a, self.lon_a
        a = sin((1 - f) * d) / sin(d)
        b = sin(f * d) / sin(d)
        x = a * cos(p1) * cos(l1) + b * cos(p2) * cos(l2)
        y = a * cos(p1) * sin(l1) + b * cos(p2) * sin(l2)
        z = a * sin(p1) + b * sin(p2)
        return degrees(atan2(z, hypot(x, y))), degrees(atan2(y, x))


ISTANBUL_ANKARA = Leg("Istanbul", 41.0082, 28.9784,
                      "Ankara", 39.9334, 32.8597, speed_kmh=800.0)


# --------------------------------------------------------------------------- #
# Sky forecast
# --------------------------------------------------------------------------- #

def _iso(dt: datetime) -> str:
    return dt.strftime("%Y-%m-%d %H:%M:%S")


def sky_at(lat: float, lon: float, dt: datetime, bodies=("Sun", "Moon")) -> dict:
    ''' Geometric + apparent altitude and azimuth of each body, from the
        authoritative ephemeris. '''
    out = {}
    iso = _iso(dt)
    for b in bodies:
        gp = body_gp(b, iso)
        alt, az = altaz(lat, lon, gp)
        dist = body_distance_km(b, iso)
        app = C.apparent_from_geometric(alt, dist) if alt > -2 else alt
        out[b] = dict(alt_geom=alt, alt_app=app, az=az, dist_km=dist,
                      refraction_deg=(C.refraction_deg(app) if app > -1 else 0.0),
                      parallax_deg=C.parallax_in_altitude_deg(alt, dist))
    if "Sun" in out and "Moon" in out:
        out["sun_moon_daz"] = _wrap180(out["Moon"]["az"] - out["Sun"]["az"])
    return out


def _wrap180(x: float) -> float:
    return ((x + 180.0) % 360.0) - 180.0


def two_lop_sigma_km(daz_deg: float, sigma_alt_arcmin: float) -> float:
    ''' Predicted 1-sigma position error (km) from crossing TWO altitude LOPs
        whose azimuths differ by `daz_deg`, each sight having `sigma_alt_arcmin`.

        One arc-minute of altitude error = 1 nautical mile (1.852 km) along the
        azimuth.  Two independent LOPs crossing at angle T give an error ellipse
        whose RMS radius is sigma * sqrt(2) / |sin T|; as T -> 0 (bodies in line)
        the fix degenerates.
    '''
    t = radians(fabs(_wrap180(daz_deg)))
    st = fabs(sin(t))
    if st < 1e-6:
        return float("inf")
    return sigma_alt_arcmin * 1.852 * sqrt(2.0) / st


def forecast_leg(leg: Leg, departure: datetime, n_points: int = 13,
                 sigma_alt_arcmin: float = 2.0) -> list:
    ''' Sample the sky along the leg.  Returns a list of dicts, one per waypoint,
        with time, position, both bodies' geometry, the crossing angle and the
        predicted two-LOP fix sigma. '''
    rows = []
    for i in range(n_points):
        f = i / (n_points - 1) if n_points > 1 else 0.0
        dt = departure + timedelta(hours=leg.duration_h * f)
        lat, lon = leg.position_at(f)
        sky = sky_at(lat, lon, dt)
        sun, moon = sky["Sun"], sky["Moon"]
        daz = sky["sun_moon_daz"]
        both_up = sun["alt_app"] > 5.0 and moon["alt_app"] > 5.0
        rows.append(dict(
            t=dt, frac=f, lat=lat, lon=lon,
            sun_alt=sun["alt_app"], sun_az=sun["az"],
            moon_alt=moon["alt_app"], moon_az=moon["az"],
            moon_dist_km=moon["dist_km"],
            sun_refr_arcmin=sun["refraction_deg"] * 60.0,
            moon_parallax_deg=moon["parallax_deg"],
            daz=daz, both_up=both_up,
            fix_sigma_km=(two_lop_sigma_km(daz, sigma_alt_arcmin)
                          if both_up else float("nan")),
        ))
    return rows


def observation_windows(leg: Leg, day: datetime, step_min: int = 15,
                        sigma_alt_arcmin: float = 2.0,
                        min_alt_deg: float = 10.0) -> list:
    ''' Scan a whole day at the leg's mid-point and return the windows in which
        BOTH bodies are usable, with the best crossing angle in each.  This is
        the "when should I schedule the flight / the shots" product. '''
    lat, lon = leg.position_at(0.5)
    start = day.replace(hour=0, minute=0, second=0, microsecond=0)
    rows = []
    for i in range(int(24 * 60 / step_min) + 1):
        dt = start + timedelta(minutes=step_min * i)
        sky = sky_at(lat, lon, dt)
        s, m = sky["Sun"], sky["Moon"]
        ok = s["alt_app"] > min_alt_deg and m["alt_app"] > min_alt_deg
        rows.append(dict(t=dt, sun_alt=s["alt_app"], moon_alt=m["alt_app"],
                         daz=sky["sun_moon_daz"], ok=ok,
                         fix_sigma_km=(two_lop_sigma_km(sky["sun_moon_daz"],
                                                        sigma_alt_arcmin)
                                       if ok else float("nan"))))
    # group consecutive ok rows into windows
    windows, cur = [], None
    for r in rows:
        if r["ok"]:
            cur = cur or []
            cur.append(r)
        elif cur:
            windows.append(cur)
            cur = None
    if cur:
        windows.append(cur)
    return [dict(start=w[0]["t"], end=w[-1]["t"],
                 best=min(w, key=lambda r: r["fix_sigma_km"]),
                 worst=max(w, key=lambda r: r["fix_sigma_km"]))
            for w in windows]
