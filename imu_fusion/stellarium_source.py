''' Stellarium as the authoritative astronomical source.

    When a Stellarium export is present at `AUTHORITATIVE_CSV`, this module makes
    Stellarium the primary ground-truth ephemeris for the study and the device
    table build.  It is the single ingestion point for
    `imu_fusion/tools/stellarium_export.ssc` output (see
    `stellarium_reference.md`).

    Design goals:
      * ENGINE-FREE.  Greenwich apparent sidereal time (to convert Stellarium's
        RA-of-date into a Greenwich hour angle) is computed here with a compact
        IAU formula, so the Stellarium path needs NO astropy/skyfield.
      * TOLERANT.  The loader accepts the rich schema (with distance/alt/az/…) or
        the minimal `utc,body,ra_deg,dec_deg`; blank optional cells, `#` comment
        lines (the exporter's `#SCHEMA`/`#END` markers) and header rows are
        skipped.
      * OPTIONAL.  `get_table()` returns None when no export exists, so the study
        transparently falls back to the starfix almanac until the CSV is dropped
        in.

    (c) 2026.  MIT License (see LICENSE file).
'''

import os
from bisect import bisect_left
from datetime import datetime, timezone
from math import sin, cos, radians

_HERE = os.path.dirname(os.path.abspath(__file__))
AUTHORITATIVE_CSV = os.path.join(_HERE, "sample_data", "stellarium_ephemeris.csv")

_AU_KM = 149_597_870.7


# --------------------------------------------------------------------------- #
# Engine-free sidereal time (so the Stellarium path stands alone)
# --------------------------------------------------------------------------- #

def _julian_day(dt: datetime) -> float:
    ''' Julian Day from a (UTC) datetime. '''
    dt = dt.astimezone(timezone.utc) if dt.tzinfo else dt.replace(tzinfo=timezone.utc)
    return dt.timestamp() / 86400.0 + 2440587.5


def gmst_deg(dt: datetime) -> float:
    ''' Greenwich Mean Sidereal Time (degrees), Meeus 12.4. '''
    jd = _julian_day(dt)
    t = (jd - 2451545.0) / 36525.0
    g = (280.46061837 + 360.98564736629 * (jd - 2451545.0)
         + 0.000387933 * t * t - t * t * t / 38710000.0)
    return g % 360.0


def gast_deg(dt: datetime) -> float:
    ''' Greenwich Apparent Sidereal Time (degrees) = GMST + equation of the
        equinoxes, with the leading nutation terms (arc-second level).  Matches
        Stellarium's apparent RA-of-date convention. '''
    jd = _julian_day(dt)
    t = (jd - 2451545.0) / 36525.0
    omega = radians(125.04452 - 1934.136261 * t)
    lsun = radians(280.4665 + 36000.7698 * t)
    lmoon = radians(218.3165 + 481267.8813 * t)
    d_psi_arcsec = (-17.20 * sin(omega) - 1.32 * sin(2 * lsun)
                    - 0.23 * sin(2 * lmoon) + 0.21 * sin(2 * omega))
    eps = radians(23.439291 - 0.0130042 * t)
    eqeq_deg = (d_psi_arcsec / 3600.0) * cos(eps)
    return (gmst_deg(dt) + eqeq_deg) % 360.0


# --------------------------------------------------------------------------- #
# CSV ingestion
# --------------------------------------------------------------------------- #

def _parse_dt(s: str) -> datetime:
    s = s.strip().replace("T", " ").replace("Z", "")
    return datetime.fromisoformat(s).replace(tzinfo=timezone.utc)


def load_csv(path: str) -> list:
    ''' Parse a Stellarium export into a list of row dicts.  Requires at least
        `utc,body,ra_deg,dec_deg`; `dist_au`/`alt_deg`/`az_deg`/`elong_deg`/
        `phase`/`size_arcsec` are optional.  Returns [] if the file is absent. '''
    if not os.path.exists(path):
        return []
    rows, header = [], None
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = [p.strip() for p in s.split(",")]
            if header is None and parts and parts[0].lower() == "utc":
                header = [c.lower() for c in parts]
                continue
            if header is None or len(parts) < 4:
                continue
            row = dict(zip(header, parts))
            if not row.get("ra_deg") or not row.get("dec_deg"):
                continue
            try:
                rec = dict(dt=_parse_dt(row["utc"]),
                           body=row["body"].strip().lower(),
                           ra_deg=float(row["ra_deg"]) % 360.0,
                           dec_deg=float(row["dec_deg"]))
            except (ValueError, KeyError):
                continue
            for opt in ("dist_au", "alt_deg", "az_deg", "elong_deg",
                        "phase", "size_arcsec"):
                v = row.get(opt, "")
                rec[opt] = float(v) if v not in ("", None) else None
            rows.append(rec)
    return rows


class Table:
    ''' Time-sorted per-body Stellarium ephemeris with linear interpolation. '''

    def __init__(self, rows: list):
        self._by_body: dict = {}
        for r in rows:
            self._by_body.setdefault(r["body"], []).append(r)
        for body in self._by_body:
            self._by_body[body].sort(key=lambda r: r["dt"])
        self._epochs = {b: [r["dt"] for r in rs]
                        for b, rs in self._by_body.items()}

    def has_body(self, body: str) -> bool:
        return body.lower() in self._by_body

    def bodies(self) -> list:
        return list(self._by_body)

    def _interp(self, body: str, dt: datetime):
        ''' Return the bracketing rows and fraction for `dt` (clamped to range). '''
        body = body.lower()
        rs = self._by_body[body]
        ts = self._epochs[body]
        dt = dt.astimezone(timezone.utc) if dt.tzinfo else dt.replace(tzinfo=timezone.utc)
        i = bisect_left(ts, dt)
        if i <= 0:
            return rs[0], rs[0], 0.0
        if i >= len(rs):
            return rs[-1], rs[-1], 0.0
        r0, r1 = rs[i - 1], rs[i]
        span = (r1["dt"] - r0["dt"]).total_seconds()
        frac = 0.0 if span == 0 else (dt - r0["dt"]).total_seconds() / span
        return r0, r1, frac

    def gp_dec_gha(self, body: str, dt: datetime):
        ''' (declination_deg, GHA_deg) at `dt` from Stellarium RA/Dec-of-date. '''
        r0, r1, f = self._interp(body, dt)
        dec = r0["dec_deg"] + f * (r1["dec_deg"] - r0["dec_deg"])
        dra = ((r1["ra_deg"] - r0["ra_deg"] + 180.0) % 360.0) - 180.0
        ra = (r0["ra_deg"] + f * dra) % 360.0
        gha = (gast_deg(dt) - ra) % 360.0
        return dec, gha

    def distance_km(self, body: str, dt: datetime):
        ''' Interpolated geocentric distance (km), or None if not in the export. '''
        r0, r1, f = self._interp(body, dt)
        if r0.get("dist_au") is None or r1.get("dist_au") is None:
            return None
        return (r0["dist_au"] + f * (r1["dist_au"] - r0["dist_au"])) * _AU_KM


# --------------------------------------------------------------------------- #
# Memoised authoritative-table accessor
# --------------------------------------------------------------------------- #

_TABLE_CACHE = {"path": None, "table": None}


def get_table(path: str = AUTHORITATIVE_CSV):
    ''' Return the Stellarium `Table` for `path`, or None if the export is
        missing/empty.  Memoised on the path. '''
    if _TABLE_CACHE["path"] == path and _TABLE_CACHE["table"] is not None:
        return _TABLE_CACHE["table"]
    rows = load_csv(path)
    table = Table(rows) if rows else None
    _TABLE_CACHE["path"] = path
    _TABLE_CACHE["table"] = table
    return table


def clear_cache() -> None:
    _TABLE_CACHE["path"] = None
    _TABLE_CACHE["table"] = None
