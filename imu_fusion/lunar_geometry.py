''' Where the Moon's face is pointing: libration and the axis position angle.

    `lunar_orientation.render_moon` takes a libration and a pole position angle
    and draws the maria in the right places — but nothing computed them.  They
    were passed in, defaulted to zero.  This module supplies them from the
    ephemeris, which is what makes a rendered template comparable with a real
    photograph.

    THE QUANTITIES

      * **Optical libration (l, b)** — the selenographic longitude and latitude
        of the SUB-OBSERVER POINT, i.e. the point of the Moon's surface at the
        apparent centre of the disk.  The Moon keeps one face toward us, but the
        eccentricity of its orbit and the tilt of its axis rock that face by up
        to about ±8° in longitude and ±7° in latitude, so the disk centre wanders
        over a region a good fraction of the size of a mare.  Identifying a
        feature near the centre without knowing (l, b) is guesswork.

      * **Position angle of the axis (P)** — the direction of the Moon's north
        pole on the sky, measured from celestial north through east.  Combined
        with the parallactic angle it fixes how the whole pattern is rotated in
        a photograph.

      * **Topocentric correction** — the Moon is close enough that an observer on
        the surface sees a slightly different face than a hypothetical observer
        at the Earth's centre.  The shift is of order the parallax, up to ~1°,
        and is applied here to first order.  It is the reason this is computed
        "for Istanbul" rather than for the Earth.

    Reference: Meeus, *Astronomical Algorithms*, 2nd ed., ch. 53 (optical
    libration and position angle of the axis), with the leading nutation terms
    from ch. 22.

    HONEST LIMIT.  Only the OPTICAL libration is computed.  The physical
    libration — a real, forced oscillation of the Moon's body — is a further
    term of order 0.02°, negligible at the precision this study works to, and is
    omitted.  Nothing here has been checked against JPL Horizons, because the
    sandbox this was written in cannot reach it; the magnitudes are checked
    against their known bounds instead.

    (c) 2026.  MIT License (see LICENSE file).
'''

from math import (sin, cos, tan, asin, atan2, radians, degrees, sqrt, hypot)

from .astro import body_gp, gp_dec_gha, body_distance_km
from .stellarium_source import gast_deg, _parse_dt

# Inclination of the lunar equator to the ecliptic (Meeus 53).
I_INC = 1.54242
EARTH_R_KM = 6378.14


def _centuries(dt):
    jd = dt.timestamp() / 86400.0 + 2440587.5
    return (jd - 2451545.0) / 36525.0


def _nutation_obliquity(t):
    ''' Nutation in longitude (deg) and true obliquity (deg), leading terms. '''
    om = radians(125.04452 - 1934.136261 * t)
    ls = radians(280.4665 + 36000.7698 * t)
    lm = radians(218.3165 + 481267.8813 * t)
    dpsi = (-17.20 * sin(om) - 1.32 * sin(2 * ls)
            - 0.23 * sin(2 * lm) + 0.21 * sin(2 * om)) / 3600.0
    deps = (9.20 * cos(om) + 0.57 * cos(2 * ls)
            + 0.10 * cos(2 * lm) - 0.09 * cos(2 * om)) / 3600.0
    eps0 = 23.439291 - 0.0130042 * t
    return dpsi, eps0 + deps


def _moon_radec(time_iso):
    ''' Apparent geocentric right ascension and declination (deg) of the Moon,
        from the project's authoritative ephemeris. '''
    dec, gha = gp_dec_gha(body_gp("Moon", time_iso))
    dt = _parse_dt(time_iso)
    ra = (gast_deg(dt) - gha) % 360.0
    return ra, dec


def _equatorial_to_ecliptic(ra_deg, dec_deg, eps_deg):
    a, d, e = radians(ra_deg), radians(dec_deg), radians(eps_deg)
    lam = atan2(sin(a) * cos(e) + tan(d) * sin(e), cos(a))
    beta = asin(sin(d) * cos(e) - cos(d) * sin(e) * sin(a))
    return degrees(lam) % 360.0, degrees(beta)


def libration_from_ecliptic(lam, beta, ra, t, dpsi, eps):
    ''' The ch.53 geometry alone, given the Moon's apparent position.

        Split out from `geocentric_libration` for one concrete reason: the
        project's almanac only spans 2024-2030, so Meeus's worked example
        (1992 April 12) cannot be run through the ephemeris.  Feeding his own
        lambda/beta/alpha in here checks the geometry against a published
        answer, which is the only independent check available offline.

        lam, beta : apparent geocentric ecliptic longitude/latitude (deg)
        ra        : apparent right ascension (deg)
        t         : Julian centuries from J2000
        dpsi      : nutation in longitude (deg)
        eps       : true obliquity (deg)
    '''
    F = (93.2720950 + 483202.0175233 * t - 0.0036539 * t * t) % 360.0
    Om = (125.0445479 - 1934.1362891 * t + 0.0020754 * t * t) % 360.0

    W = radians(lam - dpsi - Om)
    b_r, i_r = radians(beta), radians(I_INC)
    A = atan2(sin(W) * cos(b_r) * cos(i_r) - sin(b_r) * sin(i_r),
              cos(W) * cos(b_r))
    lib_lon = ((degrees(A) - F + 180.0) % 360.0) - 180.0
    lib_lat = degrees(asin(-sin(W) * cos(b_r) * sin(i_r) - sin(b_r) * cos(i_r)))

    # Position angle of the axis (Meeus 53.5), optical part only.
    V = radians(Om + dpsi)
    e_r = radians(eps)
    X = sin(i_r) * sin(V)
    Y = sin(i_r) * cos(V) * cos(e_r) - cos(i_r) * sin(e_r)
    omega = atan2(X, Y)
    P = degrees(asin(max(-1.0, min(1.0,
        hypot(X, Y) * cos(radians(ra) - omega) / cos(radians(lib_lat))))))
    return lib_lon, lib_lat, P


def geocentric_libration(time_iso: str):
    ''' Optical libration (l, b) in degrees and the axis position angle P.

        l > 0 means the sub-observer point lies EAST on the Moon, so extra
        terrain is visible on the Moon's western (celestial-east) limb.
        Returns dict(lon, lat, pole_pa, ra, dec, lambda, beta).
    '''
    dt = _parse_dt(time_iso)
    t = _centuries(dt)
    dpsi, eps = _nutation_obliquity(t)
    ra, dec = _moon_radec(time_iso)
    lam, beta = _equatorial_to_ecliptic(ra, dec, eps)
    lib_lon, lib_lat, P = libration_from_ecliptic(lam, beta, ra, t, dpsi, eps)
    return dict(lon=lib_lon, lat=lib_lat, pole_pa=P, ra=ra, dec=dec,
                lam=lam, beta=beta, eps=eps)


def topocentric_libration(time_iso: str, obs_lat: float, obs_lon: float,
                          obs_height_m: float = 0.0):
    ''' Libration as seen from a place on the Earth's surface.

        The observer is displaced from the Earth's centre, so the Moon is seen
        from a slightly different direction and a slightly different face is
        presented.  The displacement of the apparent direction is the parallax;
        to first order the sub-observer point moves by the same angle, rotated
        into the lunar disk frame by the axis position angle.

        Returns the geocentric dict plus `lon`/`lat` corrected, `dlon`/`dlat`
        (the correction applied) and `parallax_deg` (its magnitude).
    '''
    g = geocentric_libration(time_iso)
    dt = _parse_dt(time_iso)
    d_km = body_distance_km("Moon", time_iso)

    # Observer's geocentric direction relative to the Moon: hour angle and the
    # standard topocentric parallax in RA/Dec (Meeus 40, small-angle form).
    lst = (gast_deg(dt) + obs_lon) % 360.0
    H = radians(lst - g["ra"])
    phi = radians(obs_lat)
    rho = 1.0 + obs_height_m / 6.378e6
    sin_pi = (EARTH_R_KM / d_km) * rho
    dec_r = radians(g["dec"])
    # displacement of the apparent position (degrees), east and north on the sky
    d_ra = degrees(-sin_pi * cos(phi) * sin(H) / cos(dec_r))
    d_dec = degrees(-sin_pi * (sin(phi) * cos(dec_r)
                               - cos(phi) * sin(H) * 0.0 - cos(phi) * cos(H) * sin(dec_r)))
    east = d_ra * cos(dec_r)
    north = d_dec
    par = hypot(east, north)

    # Rotate the sky-plane displacement into the lunar disk frame.  The Moon's
    # north pole lies at position angle P (from celestial north through east);
    # the sub-observer point moves opposite to the apparent displacement.
    p = radians(g["pole_pa"])
    dlat = -(north * cos(p) + east * sin(p))
    dlon = -(-north * sin(p) + east * cos(p)) / max(cos(radians(g["lat"])), 0.2)

    out = dict(g)
    out.update(lon=g["lon"] + dlon, lat=g["lat"] + dlat,
               dlon=dlon, dlat=dlat, parallax_deg=par,
               geo_lon=g["lon"], geo_lat=g["lat"])
    return out


def disk_position(seleno_lon: float, seleno_lat: float, lib_lon: float,
                  lib_lat: float):
    ''' Where a surface feature falls on the visible disk.

        Returns (u, v, visible): u is toward the Moon's west limb, v toward its
        north pole, both in units of the disk radius, before any rotation into
        image coordinates.  `visible` is False on the far side.
    '''
    lo, la = radians(seleno_lon - lib_lon), radians(seleno_lat)
    b = radians(lib_lat)
    x = cos(la) * sin(lo)
    y = sin(la) * cos(b) - cos(la) * cos(lo) * sin(b)
    z = sin(la) * sin(b) + cos(la) * cos(lo) * cos(b)
    return x, y, z > 0.0


def features_near_centre(lib_lon: float, lib_lat: float, within_frac: float = 0.45,
                         kinds=None):
    ''' Named features whose centres fall within `within_frac` of the disk radius
        of the apparent centre, nearest first — what to look for in the middle of
        a photograph.

        `within_frac` is a fraction of the DISK RADIUS, not an angle: 0.45 of the
        radius is where foreshortening is still mild (the surface is tilted
        27 deg), so features there are recognisable rather than smeared.

        Each entry carries `sep_deg`, the true angular distance along the
        surface from the sub-observer point, and `sep_km` — the honest measure of
        "how far from the centre", which the projected radius `r` understates
        badly near the limb.
    '''
    from .lunar_features import FEATURES, KM_PER_DEG
    out = []
    for name, lon, lat, ang_r, kind in FEATURES:
        if kinds and kind not in kinds:
            continue
        u, v, vis = disk_position(lon, lat, lib_lon, lib_lat)
        r = hypot(u, v)
        if vis and r <= within_frac:
            sep = degrees(asin(max(-1.0, min(1.0, r))))
            out.append(dict(name=name, kind=kind, seleno_lon=lon, seleno_lat=lat,
                            u=u, v=v, r=r, ang_radius=ang_r,
                            sep_deg=sep, sep_km=sep * KM_PER_DEG))
    out.sort(key=lambda f: f["r"])
    return out
