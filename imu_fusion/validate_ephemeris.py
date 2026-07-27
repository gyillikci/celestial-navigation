''' Independent validation of the study's Sun/Moon GROUND TRUTH.

    The whole study treats `astro.body_gp` (-> starfix's SkyAlmanac/Skyfield-
    derived, hourly, 0.1'-quantized, linearly interpolated almanac) as truth.
    A zero-noise run recovers that truth *by construction*, so it is never
    independently checked.  This module recomputes the Sun/Moon geographic
    position (GHA, Dec) and geometric altitude/azimuth from an INDEPENDENT
    ephemeris engine and reports the residuals.

    Engines, auto-selected:
      * "skyfield" -- Skyfield + a JPL DE kernel (de440s.bsp / de421.bsp).  This
        matches the almanac's own source, so it mainly measures the almanac's
        interpolation + quantization ceiling.  Used only if a kernel is available
        locally (the kernel download is often blocked); set SKYFIELD_KERNEL to a
        .bsp path to force it.
      * "astropy" -- astropy + ERFA (IAU SOFA) with the builtin ephemeris.  No
        download, and a genuinely INDEPENDENT implementation (ERFA analytic Sun +
        truncated ELP Moon), so it is the stronger cross-check here.

    A Stellarium export can be compared too -- see stellarium_reference.md.

    (c) 2026.  MIT License (see LICENSE file).
'''

import os
from datetime import datetime, timezone

from .astro import body_gp, gp_dec_gha, predicted_altitude, predicted_azimuth

_ARCSEC = 3600.0


# --------------------------------------------------------------------------- #
# Engine selection
# --------------------------------------------------------------------------- #

def _try_skyfield():
    kernel = os.environ.get("SKYFIELD_KERNEL")
    candidates = [kernel] if kernel else ["de440s.bsp", "de421.bsp"]
    try:
        from skyfield.api import Loader
        here = os.path.dirname(os.path.abspath(__file__))
        load = Loader(here, verbose=False)
        for name in candidates:
            if not name:
                continue
            path = name if os.path.isabs(name) else os.path.join(here, name)
            if os.path.exists(path):
                eph = load(path)
                ts = load.timescale(builtin=True)
                return ("skyfield", (eph, ts))
    except Exception:
        pass
    return None


def _try_astropy():
    try:
        from astropy.coordinates import solar_system_ephemeris
        from astropy.utils import iers
        iers.conf.auto_download = False            # use bundled data, stay offline
        # UT1-UTC beyond the bundled prediction table degrades by ~ms -> sub-
        # arcsecond sidereal time, negligible for an arcsecond-level check.
        iers.conf.iers_degraded_accuracy = "ignore"
        solar_system_ephemeris.set("builtin")
        return ("astropy", None)
    except Exception:
        return None


def select_engine():
    ''' Return (engine_name, handle) or (None, None) if no engine is available. '''
    for probe in (_try_skyfield, _try_astropy):
        got = probe()
        if got:
            return got
    return (None, None)


ENGINE, _HANDLE = select_engine()


# --------------------------------------------------------------------------- #
# Reference GHA / Dec from the independent engine
# --------------------------------------------------------------------------- #

def reference_gha_dec(body: str, dt: datetime):
    ''' Apparent geographic position of `body` at UTC `dt`: (dec_deg, gha_deg)
        from the active independent engine. '''
    if ENGINE == "skyfield":
        eph, ts = _HANDLE
        d = dt.astimezone(timezone.utc)
        t = ts.utc(d.year, d.month, d.day, d.hour, d.minute,
                   d.second + d.microsecond / 1e6)
        earth = eph["earth"]
        target = eph["moon" if body.lower() == "moon" else "sun"]
        astrometric = earth.at(t).observe(target).apparent()
        ra, dec, _ = astrometric.radec(epoch="date")
        gast = t.gast * 15.0                       # hours -> degrees
        gha = (gast - ra._degrees) % 360.0
        return dec.degrees, gha
    if ENGINE == "astropy":
        import warnings
        from astropy.coordinates import get_body, TETE
        from astropy.time import Time
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            t = Time(dt.astimezone(timezone.utc).replace(tzinfo=None),
                     scale="utc")
            app = get_body("moon" if body.lower() == "moon" else "sun",
                           t).transform_to(TETE(obstime=t))
            gast = t.sidereal_time("apparent", "greenwich").deg
            gha = (gast - app.ra.deg) % 360.0
            return app.dec.deg, gha
    raise RuntimeError("no independent ephemeris engine available")


# --------------------------------------------------------------------------- #
# Comparison against the study's ground truth
# --------------------------------------------------------------------------- #

def _iso(dt: datetime) -> str:
    return dt.strftime("%Y-%m-%d %H:%M:%S")


def compare(times, locations=((51.5, 0.0),), bodies=("Sun", "Moon")):
    ''' Residuals (independent engine - starfix almanac) over the given times.

        Returns per-body dicts of GHA/Dec residuals (arc-seconds) and the
        geometric alt/az residual (arc-minutes, computed identically from both
        GPs so refraction/parallax conventions cannot confound it).
    '''
    out = {b: {"gha_as": [], "dec_as": [], "alt_am": [], "az_am": []}
           for b in bodies}
    for dt in times:
        iso = _iso(dt)
        for b in bodies:
            dec_ref, gha_ref = reference_gha_dec(b, dt)
            gp = body_gp(b, iso)
            dec_sf, gha_sf = gp_dec_gha(gp)
            dgha = ((gha_ref - gha_sf + 180) % 360 - 180) * _ARCSEC
            ddec = (dec_ref - dec_sf) * _ARCSEC
            out[b]["gha_as"].append(dgha)
            out[b]["dec_as"].append(ddec)
            # Geometric alt/az from each GP, at each location.
            from starfix import LatLonGeocentric
            gp_ref = LatLonGeocentric(dec_ref, -gha_ref)
            for lat, lon in locations:
                da = (predicted_altitude(lat, lon, gp_ref) -
                      predicted_altitude(lat, lon, gp)) * 60.0
                dz = (((predicted_azimuth(lat, lon, gp_ref) -
                        predicted_azimuth(lat, lon, gp)) + 180) % 360 - 180) * 60.0
                out[b]["alt_am"].append(da)
                out[b]["az_am"].append(dz)
    return out


def summarise(out) -> str:
    def rms(v):
        return (sum(x * x for x in v) / len(v)) ** 0.5 if v else float("nan")

    def mx(v):
        return max(abs(x) for x in v) if v else float("nan")
    lines = [f"Ground-truth validation vs '{ENGINE}' (independent ephemeris):"]
    for b, d in out.items():
        lines.append(
            f"  {b:5}  GHA rms={rms(d['gha_as']):.2f}\" max={mx(d['gha_as']):.2f}\"  "
            f"Dec rms={rms(d['dec_as']):.2f}\" max={mx(d['dec_as']):.2f}\"  "
            f"| geom alt max={mx(d['alt_am']):.3f}' az max={mx(d['az_am']):.3f}'")
    return "\n".join(lines)


def _gast_deg(dt: datetime) -> float:
    ''' Greenwich apparent sidereal time (degrees) from the active engine. '''
    if ENGINE == "skyfield":
        _, ts = _HANDLE
        d = dt.astimezone(timezone.utc)
        t = ts.utc(d.year, d.month, d.day, d.hour, d.minute,
                   d.second + d.microsecond / 1e6)
        return (t.gast * 15.0) % 360.0
    from astropy.time import Time
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t = Time(dt.astimezone(timezone.utc).replace(tzinfo=None), scale="utc")
        return t.sidereal_time("apparent", "greenwich").deg % 360.0


# --------------------------------------------------------------------------- #
# Stellarium reference (a CSV the user exports; see stellarium_reference.md)
# --------------------------------------------------------------------------- #

def load_reference_csv(path: str):
    ''' Read a reference CSV with header `utc,body,ra_deg,dec_deg` (apparent,
        geocentric, of date -- e.g. a Stellarium AstroCalc ephemeris export).
        Returns a list of (datetime, body, ra_deg, dec_deg).  Blank/`#` lines and
        rows with empty values are skipped, so a template with no data rows
        simply yields nothing. '''
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#") or s.lower().startswith("utc"):
                continue
            parts = [p.strip() for p in s.split(",")]
            if len(parts) < 4 or not all(parts[:4]):
                continue
            dt = datetime.fromisoformat(parts[0]).replace(tzinfo=timezone.utc)
            rows.append((dt, parts[1], float(parts[2]), float(parts[3])))
    return rows


def compare_reference_csv(path: str):
    ''' Compare a reference CSV (Stellarium export) against the starfix ground
        truth.  RA(of date) -> GHA via GAST.  Returns a list of per-row residual
        dicts (GHA/Dec in arc-seconds), empty if the file has no data rows. '''
    out = []
    for dt, body, ra_deg, dec_deg in load_reference_csv(path):
        gha_ref = (_gast_deg(dt) - ra_deg) % 360.0
        dec_sf, gha_sf = gp_dec_gha(body_gp(body.capitalize(), _iso(dt)))
        out.append(dict(utc=_iso(dt), body=body,
                        gha_as=((gha_ref - gha_sf + 180) % 360 - 180) * _ARCSEC,
                        dec_as=(dec_deg - dec_sf) * _ARCSEC))
    return out


def default_grid(n_days=60, step_hours=7):
    ''' A grid of UTC times spanning part of the almanac's 2024-2030 range. '''
    from datetime import timedelta
    start = datetime(2026, 3, 24, 12, 0, 0, tzinfo=timezone.utc)
    return [start + timedelta(hours=step_hours * i)
            for i in range(int(n_days * 24 / step_hours))]


if __name__ == "__main__":
    if ENGINE is None:
        print("No independent ephemeris engine available "
              "(install astropy, or place a JPL .bsp for skyfield).")
    else:
        res = compare(default_grid())
        print(summarise(res))
