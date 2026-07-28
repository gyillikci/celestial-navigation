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

# A few named places so the CLI is usable without looking up coordinates.
PLACES = {
    "istanbul": (41.0082, 28.9784),
    "ankara": (39.9334, 32.8597),
    "izmir": (38.4237, 27.1428),
    "antalya": (36.8969, 30.7133),
    "trabzon": (41.0015, 39.7178),
    "london": (51.5072, -0.1276),
    "athens": (37.9838, 23.7275),
    "dubai": (25.2048, 55.2708),
}


def leg_from_names(a: str, b: str, speed_kmh: float = 800.0) -> Leg:
    ''' Build a Leg from two names in `PLACES` (case-insensitive). '''
    ka, kb = a.strip().lower(), b.strip().lower()
    if ka not in PLACES or kb not in PLACES:
        raise KeyError(f"unknown place; known: {', '.join(sorted(PLACES))}")
    return Leg(a.title(), *PLACES[ka], b.title(), *PLACES[kb], speed_kmh=speed_kmh)


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


# --------------------------------------------------------------------------- #
# The mission brief
#
# Design note -- WHY THIS GATES ON AVAILABILITY, NOT ON CROSSING ANGLE.
# A classical two-LOP fix degenerates as the bodies' azimuths converge
# (`two_lop_sigma_km` blows up as 1/sin T).  The study's full stack does NOT:
# ablating the Istanbul->Ankara leg, an ALTITUDE-ONLY fix degraded 76 -> 326 km
# rms between a 89-deg and a 6-deg crossing, while the full stack (IMU +
# differential Sun-Moon dq + parallactic + heading) stayed 6.4 -> 5.7 km, i.e.
# flat.  So the planner's job is NOT to hunt for a 90-deg crossing; it is to
# answer the binary question "is the Moon available at all?" and to hand the crew
# the window.  The crossing angle is reported for information only.
# --------------------------------------------------------------------------- #

# Expected fix accuracy, from the air-regime study runs on this leg (5 seeds,
# 12 shots, 30 km DR prior, full stack).  SIMULATION figures -- the on-device
# tilt floor must still be calibrated before quoting these operationally.
FIX_RMS_SUN_MOON_KM = 5.5      # observed 3.8-6.4 across geometries
FIX_RMS_SUN_ONLY_KM = 5.1      # latitude-strong; longitude leans on the prior
MIN_USABLE_ALT_DEG = 10.0      # below this, refraction model error grows fast
ZENITH_CAUTION_DEG = 85.0      # near-overhead: azimuth (LOP direction) is soft


class MissionBrief:
    ''' The per-leg planning product. '''

    def __init__(self, leg, day, windows, moon_available, fix_mode,
                 expected_rms_km, recommended, warnings, elong_deg, illum):
        self.leg = leg
        self.day = day
        self.windows = windows
        self.moon_available = moon_available
        self.fix_mode = fix_mode                  # "sun+moon" | "sun-only"
        self.expected_rms_km = expected_rms_km
        self.recommended = recommended            # best window dict or None
        self.warnings = warnings
        self.elong_deg = elong_deg
        self.illum = illum


def plan_leg(leg: Leg, day: datetime, step_min: int = 15,
             sigma_alt_arcmin: float = 2.0,
             min_alt_deg: float = MIN_USABLE_ALT_DEG) -> MissionBrief:
    ''' Produce the mission brief for `leg` on `day`.

        The primary output is the binary Moon-availability flag plus the
        both-bodies-up windows; the recommended window is simply the LONGEST one
        (most shot opportunities), not the one with the best crossing angle.
    '''
    from .optical_attitude import moon_elongation_deg, moon_illuminated_fraction

    day = day.replace(hour=0, minute=0, second=0, microsecond=0)
    noon_iso = _iso(day + timedelta(hours=12))
    elong = moon_elongation_deg(noon_iso)
    illum = moon_illuminated_fraction(noon_iso)

    windows = observation_windows(leg, day, step_min=step_min,
                                  sigma_alt_arcmin=sigma_alt_arcmin,
                                  min_alt_deg=min_alt_deg)
    moon_available = bool(windows)
    # Recommend the LONGEST window: more chances to shoot, and the fix quality is
    # geometry-insensitive (see the design note above).
    recommended = (max(windows, key=lambda w: (w["end"] - w["start"]).total_seconds())
                   if windows else None)

    warnings = []
    if not moon_available:
        warnings.append(
            f"Moon NOT usable today (elongation {elong:.0f} deg, illuminated "
            f"{illum:.0%}): it is never above {min_alt_deg:.0f} deg while the Sun "
            f"is. Plan a SUN-ONLY sight: a strong latitude line, with longitude "
            f"resting on the dead-reckoned prior.")
    if illum < 0.03 or illum > 0.98:
        warnings.append(
            f"Moon is near new/full (illuminated {illum:.0%}): the bright-limb "
            f"axis is unusable, so the differential Sun-Moon dq factor -- the "
            f"horizon-free observable -- is degraded.")
    # Sun near the zenith softens the LOP direction (azimuth ill-defined).
    lat, lon = leg.position_at(0.5)
    peak = max((sky_at(lat, lon, day + timedelta(minutes=step_min * i))["Sun"]["alt_app"]
                for i in range(int(24 * 60 / step_min) + 1)))
    if peak > ZENITH_CAUTION_DEG:
        warnings.append(
            f"Sun peaks at {peak:.1f} deg (near zenith): the altitude circle is "
            f"small and its azimuth -- the LOP direction -- is poorly conditioned. "
            f"Prefer shots away from local noon.")

    return MissionBrief(
        leg=leg, day=day, windows=windows, moon_available=moon_available,
        fix_mode="sun+moon" if moon_available else "sun-only",
        expected_rms_km=(FIX_RMS_SUN_MOON_KM if moon_available
                         else FIX_RMS_SUN_ONLY_KM),
        recommended=recommended, warnings=warnings,
        elong_deg=elong, illum=illum)


def find_next_opportunity(leg: Leg, start_day: datetime, max_days: int = 30,
                          min_window_min: int = 60) -> list:
    ''' Scan forward for days offering a usable Sun+Moon window of at least
        `min_window_min`.  Returns a list of (day, brief) for the usable days --
        the "if the mission can slip, fly then" product. '''
    out = []
    for d in range(max_days):
        day = start_day + timedelta(days=d)
        brief = plan_leg(leg, day, step_min=20)
        if brief.recommended is not None:
            w = brief.recommended
            if (w["end"] - w["start"]).total_seconds() / 60.0 >= min_window_min:
                out.append((day, brief))
    return out


def brief_text(brief: MissionBrief) -> str:
    ''' Render the mission brief as a readable block. '''
    leg = brief.leg
    L = []
    L.append(f"MISSION BRIEF  {leg.name_a} -> {leg.name_b}   {brief.day:%Y-%m-%d} (UTC)")
    L.append(f"  Leg: {leg.distance_km:.0f} km, {leg.duration_h*60:.0f} min "
             f"at {leg.speed_kmh:.0f} km/h")
    L.append(f"  Moon: elongation {brief.elong_deg:.0f} deg, illuminated {brief.illum:.0%}")
    L.append("")
    flag = "YES" if brief.moon_available else "NO"
    L.append(f"  MOON AVAILABLE (daytime, with the Sun): {flag}")
    L.append(f"  Fix mode: {brief.fix_mode.upper()}   "
             f"expected ~{brief.expected_rms_km:.1f} km rms (simulation)")
    if brief.fix_mode == "sun-only":
        L.append("    NOTE: one body = one line of position. Latitude is strong; "
                 "longitude comes from the prior, not from the sky.")
    L.append("")
    if brief.windows:
        L.append("  Both-bodies-up windows (UTC):")
        for w in brief.windows:
            mins = (w["end"] - w["start"]).total_seconds() / 60.0
            star = "  <== RECOMMENDED" if w is brief.recommended else ""
            L.append(f"    {w['start']:%H:%M}-{w['end']:%H:%M}  ({mins:>4.0f} min)  "
                     f"crossing {abs(w['best']['daz']):5.1f}-{abs(w['worst']['daz']):5.1f} deg{star}")
        L.append("    (crossing angle is informational: the fused fix is "
                 "geometry-insensitive)")
    else:
        L.append("  Both-bodies-up windows: NONE")
    if brief.warnings:
        L.append("")
        L.append("  WARNINGS:")
        for w in brief.warnings:
            L.append(f"    - {w}")
    return "\n".join(L)


def _main(argv=None):
    import argparse
    p = argparse.ArgumentParser(
        description="Celestial mission planner: is the Moon available on this "
                    "leg, and when?")
    p.add_argument("--from", dest="a", default="Istanbul")
    p.add_argument("--to", dest="b", default="Ankara")
    p.add_argument("--date", default=None, help="YYYY-MM-DD (UTC), default today")
    p.add_argument("--speed-kmh", type=float, default=800.0)
    p.add_argument("--step-min", type=int, default=15)
    p.add_argument("--next", action="store_true",
                   help="if the Moon is unavailable, scan ahead for the next "
                        "usable days")
    p.add_argument("--scan-days", type=int, default=30)
    args = p.parse_args(argv)

    leg = leg_from_names(args.a, args.b, speed_kmh=args.speed_kmh)
    day = (datetime.strptime(args.date, "%Y-%m-%d").replace(tzinfo=timezone.utc)
           if args.date else datetime.now(timezone.utc))
    brief = plan_leg(leg, day, step_min=args.step_min)
    print(brief_text(brief))

    if args.next and not brief.moon_available:
        print("\n  NEXT USABLE DAYS (>= 60 min window):")
        found = find_next_opportunity(leg, day, max_days=args.scan_days)
        if not found:
            print(f"    none within {args.scan_days} days")
        for d, br in found[:8]:
            w = br.recommended
            mins = (w["end"] - w["start"]).total_seconds() / 60.0
            print(f"    {d:%Y-%m-%d}  {w['start']:%H:%M}-{w['end']:%H:%M} UTC "
                  f"({mins:>4.0f} min), illuminated {br.illum:.0%}")


if __name__ == "__main__":
    _main()
