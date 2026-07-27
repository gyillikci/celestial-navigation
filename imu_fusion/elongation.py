''' Sun-Moon ELONGATION (angular separation) as a navigation observable, and a
    position error budget that places it against the altitude lines of position.

    In daytime both bodies are up, so the tele can measure each disk centre and
    -- referenced to the shared attitude -- their angular separation E (the
    classic "lunar distance").  Two things make E special:

      * dE/dt is large (~0.5 deg/hr, the Moon's motion against the Sun), so E is
        a strong ABSOLUTE-TIME observable -- the pre-chronometer longitude method.
        With the photo EXIF stripped (no timestamp), E is exactly what recovers
        the time, hence the longitude.
      * dE/d(observer position) is small (only the lunar parallax shifts it), so
        as a DIRECT position line E is weak compared with an altitude sight.

    This module computes E and its Jacobians from the real ephemeris (finite
    differences through astro.body_gp) and prints a position error budget.

    (c) 2026.  MIT License (see LICENSE file).
'''

from datetime import datetime, timedelta, timezone
from math import sin, cos, acos, radians, degrees, sqrt

from .astro import body_gp, gp_dec_gha, predicted_altitude, altaz

_NM_PER_ARCMIN_KM = 1.852            # 1' of altitude = 1 nmi on the ground
_KM_PER_DEG = 111.195                # ground km per degree of great circle


def _sep_from_gp(gp1, gp2):
    ''' Great-circle angle (deg) between two geographic subpoints = the
        GEOCENTRIC angular separation of the two bodies. '''
    d1, g1 = gp_dec_gha(gp1)
    d2, g2 = gp_dec_gha(gp2)
    c = (sin(radians(d1)) * sin(radians(d2))
         + cos(radians(d1)) * cos(radians(d2)) * cos(radians(g1 - g2)))
    return degrees(acos(max(-1.0, min(1.0, c))))


def elongation_deg(time_iso: str) -> float:
    ''' Geocentric Sun-Moon elongation (deg) at UTC `time_iso`. '''
    return _sep_from_gp(body_gp("Sun", time_iso), body_gp("Moon", time_iso))


def topocentric_elongation_deg(time_iso: str, lat: float, lon: float) -> float:
    ''' Elongation as seen from (lat, lon): angle between the observer->Sun and
        observer->Moon directions, using the study's alt/az (which include the
        Moon's parallax through starfix). '''
    a_s, z_s = altaz(lat, lon, body_gp("Sun", time_iso))
    a_m, z_m = altaz(lat, lon, body_gp("Moon", time_iso))
    c = (sin(radians(a_s)) * sin(radians(a_m))
         + cos(radians(a_s)) * cos(radians(a_m)) * cos(radians(z_s - z_m)))
    return degrees(acos(max(-1.0, min(1.0, c))))


def _iso(dt):
    return dt.strftime("%Y-%m-%d %H:%M:%S")


def elongation_time_rate(time_iso: str, dt_min: float = 2.0) -> float:
    ''' dE/dt in deg per hour (central finite difference). '''
    t = datetime.fromisoformat(time_iso).replace(tzinfo=timezone.utc)
    ep = elongation_deg(_iso(t + timedelta(minutes=dt_min)))
    em = elongation_deg(_iso(t - timedelta(minutes=dt_min)))
    return (ep - em) / (2 * dt_min / 60.0)


def elongation_position_rate(time_iso, lat, lon, d_deg=0.5):
    ''' (dE/dlat, dE/dlon) in deg-E per deg of observer move (finite diff). '''
    e_la = (topocentric_elongation_deg(time_iso, lat + d_deg, lon)
            - topocentric_elongation_deg(time_iso, lat - d_deg, lon)) / (2 * d_deg)
    e_lo = (topocentric_elongation_deg(time_iso, lat, lon + d_deg)
            - topocentric_elongation_deg(time_iso, lat, lon - d_deg)) / (2 * d_deg)
    return e_la, e_lo


def altitude_position_km_per_arcmin():
    ''' A body altitude LOP: 1 arcmin of altitude error = 1 nmi on the ground. '''
    return _NM_PER_ARCMIN_KM


def position_budget(time_iso: str, lat: float, lon: float,
                    sigma_alt_arcmin: float = 2.0,
                    sigma_sep_arcmin: float = 2.0,
                    sigma_clock_s: float = 0.0):
    ''' Per-observable POSITION error budget (km, 1-sigma, single line of
        position), plus the elongation-as-time -> longitude channel.

        sigma_alt_arcmin : per-body altitude sight noise.
        sigma_sep_arcmin : Sun-Moon separation measurement noise (two disk
                           centroids referenced to the shared attitude).
        sigma_clock_s    : clock uncertainty (s); if >0, elongation pins time.
    '''
    out = {}
    # altitude lines of position (Sun, Moon) -- the workhorses
    out["alt_sun_km"] = sigma_alt_arcmin * _NM_PER_ARCMIN_KM
    out["alt_moon_km"] = sigma_alt_arcmin * _NM_PER_ARCMIN_KM

    # geometry: the two altitude LOPs cross at the Sun-Moon azimuth difference,
    # so the single-epoch two-body fix is sigma_LOP*sqrt(2)/|sin(dAz)|.
    _, z_s = altaz(lat, lon, body_gp("Sun", time_iso))
    _, z_m = altaz(lat, lon, body_gp("Moon", time_iso))
    d_az = abs(((z_s - z_m + 180) % 360) - 180)
    out["delta_az_deg"] = d_az
    sin_cross = max(abs(sin(radians(d_az))), 1e-3)
    out["two_lop_fix_km"] = out["alt_sun_km"] * sqrt(2.0) / sin_cross

    # elongation as a DIRECT position line (parallax only -> weak)
    e_la, e_lo = elongation_position_rate(time_iso, lat, lon)
    grad_deg_per_deg = sqrt(e_la ** 2 + e_lo ** 2)          # deg-E per deg move
    grad_deg_per_km = grad_deg_per_deg / _KM_PER_DEG
    sig_sep_deg = sigma_sep_arcmin / 60.0
    out["elong_direct_km"] = (sig_sep_deg / grad_deg_per_km
                              if grad_deg_per_km > 1e-6 else float("inf"))
    out["elong_direct_negligible"] = grad_deg_per_km <= 1e-6

    # elongation as an ABSOLUTE-TIME observable -> longitude (lunar distance)
    dEdt = elongation_time_rate(time_iso)                   # deg/hr
    time_sig_hr = sig_sep_deg / abs(dEdt) if dEdt else float("inf")
    out["elong_time_s"] = time_sig_hr * 3600.0
    # Earth turns 15 deg/hr; a time error maps to a longitude error
    lon_sig_deg = time_sig_hr * 15.0
    out["elong_longitude_km"] = lon_sig_deg * _KM_PER_DEG * cos(radians(lat))

    out["dE_dt_deg_per_hr"] = dEdt
    out["dE_dpos_deg_per_deg"] = grad_deg_per_deg
    out["elongation_deg"] = elongation_deg(time_iso)
    return out


def summarise(time_iso, lat=51.5, lon=0.0, **kw):
    b = position_budget(time_iso, lat, lon, **kw)
    direct = ("negligible (parallax-only)" if b["elong_direct_negligible"]
              else f"{b['elong_direct_km']:.0f} km")
    L = [f"Position budget @ {time_iso}  ({lat:.1f},{lon:.1f})  "
         f"elongation={b['elongation_deg']:.1f} deg  dAz={b['delta_az_deg']:.0f} deg",
         f"  altitude LOP per body    : {b['alt_sun_km']:.1f} km  (1 arcmin = 1 nmi)",
         f"  two-body single-epoch fix: {b['two_lop_fix_km']:.1f} km  "
         f"(LOPs cross at dAz={b['delta_az_deg']:.0f} deg)",
         f"  elongation -> TIME       : dE/dt={b['dE_dt_deg_per_hr']:.3f} deg/hr "
         f"-> time to {b['elong_time_s']:.0f} s -> longitude "
         f"{b['elong_longitude_km']:.0f} km  (recovers a stripped clock!)",
         f"  elongation -> DIRECT pos : {direct}"]
    return "\n".join(L)


if __name__ == "__main__":
    print(summarise("2026-03-24 12:00:00"))
