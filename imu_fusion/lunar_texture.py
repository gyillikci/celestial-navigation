''' Render the Moon as it actually looks: real crater texture, real libration,
    real terminator.

    `lunar_orientation.render_moon` draws ten fuzzy ellipses for the maria.  That
    is enough to demonstrate rotational matching, but it is not enough to MATCH A
    PHOTOGRAPH: a photograph resolves crater rays, individual highland patches
    and the exact bite the terminator takes out of the limb, and a ten-blob
    cartoon correlates with all of them about equally badly.

    This module renders from a real albedo map instead.  Stellarium ships one --
    `textures/moon_4k.jpg`, 4096x2048 equirectangular, combined by the Stellarium
    project from USGS Astrogeology and Clementine data, PUBLIC DOMAIN (see
    /usr/share/doc/stellarium-data/copyright, section 4.4).  Using it keeps the
    standing rule that astronomy comes from Stellarium, and it costs nothing in
    licensing.

    WHAT IS MODELLED

      * **Libration** -- the sub-observer selenographic point (l, b) from
        `lunar_geometry`, so the texture is projected from the direction the
        observer actually sees it, not from (0, 0).
      * **Position angle** -- the disk is rotated by the axis position angle P
        (plus whatever camera roll the caller supplies), so the pattern lands at
        the orientation a photograph would show.
      * **The terminator** -- illumination is computed from the SUB-SOLAR
        selenographic point (Meeus 53, "selenographic position of the Sun"), so
        near full Moon the render loses the same crescent of limb the photograph
        does.  This is the part that makes the terminator legible: it is not
        drawn as a phase cut across the disk, it falls out of the surface
        geometry.
      * **Photometry** -- Lommel-Seeliger, I ~ A * mu0 / (mu0 + mu).  The Moon is
        a famously non-Lambertian surface: at full it looks like a flat disk, not
        a shaded ball, and Lommel-Seeliger reproduces that while Lambert does
        not.  Getting this wrong puts a spurious radial brightness gradient into
        the template, which a correlator will happily lock onto.

    WHAT IS NOT

      * No opposition surge, no macroscopic roughness (Hapke), no colour.  All of
        those change brightness, none change WHERE anything is, and position is
        what the matching uses.
      * No topographic shadowing.  Near the terminator real craters throw long
        shadows the smooth-sphere model cannot know about, so the last degree or
        two before the terminator is the least trustworthy part of the render.

    (c) 2026.  MIT License (see LICENSE file).
'''

import os

import numpy as np

from .lunar_geometry import (geocentric_libration, topocentric_libration,
                             libration_from_ecliptic, _centuries,
                             _nutation_obliquity, _equatorial_to_ecliptic)
from .astro import body_gp, gp_dec_gha
from .stellarium_source import gast_deg, _parse_dt

# Stellarium's lunar albedo map.  Ordered best-first; the 4k map is the one worth
# using, `moon.png` is a 1k fallback if only the smaller data package is present.
TEXTURE_CANDIDATES = (
    "/usr/share/stellarium/textures/moon_4k.jpg",
    "/usr/share/stellarium/textures/moon.png",
    os.path.join(os.path.dirname(__file__), "sample_data", "moon_4k.jpg"),
)

_TEX_CACHE = {}


def find_texture(path: str = None):
    ''' Path to the lunar albedo map, or None if Stellarium is not installed. '''
    if path:
        return path if os.path.exists(path) else None
    for p in TEXTURE_CANDIDATES:
        if os.path.exists(p):
            return p
    return None


def load_texture(path: str = None, max_width: int = 4096):
    ''' Equirectangular albedo map as a float array, rows +90 -> -90 latitude,
        columns 0 -> 360 in selenographic longitude with lon 0 at the CENTRE
        column (the usual planetary-texture convention, verified against the
        photograph in `test_imu_fusion`).

        Cached: decoding a 4k JPEG on every render would dominate the runtime.
    '''
    p = find_texture(path)
    if p is None:
        raise FileNotFoundError(
            "no lunar albedo map found; install Stellarium's data package "
            "(textures/moon_4k.jpg) or pass an explicit path")
    key = (p, max_width)
    if key in _TEX_CACHE:
        return _TEX_CACHE[key]
    from PIL import Image
    im = Image.open(p).convert("L")
    if im.width > max_width:
        im = im.resize((max_width, max_width // 2), Image.LANCZOS)
    _TEX_CACHE[key] = np.asarray(im, dtype=float)
    return _TEX_CACHE[key]


# --------------------------------------------------------------------------- #
# Where the Sun is, in selenographic coordinates
# --------------------------------------------------------------------------- #

def subsolar_point(time_iso: str):
    ''' Selenographic longitude and latitude of the Sun (deg).

        The terminator is the great circle 90 deg from this point, so this single
        pair fixes it completely -- there is no separate "phase" parameter.
        Meeus 53: the Sun's selenographic position comes from the same libration
        formulae fed the direction OPPOSITE the Sun as seen from the Moon,

            lam_H = lam_sun + 180 + (d_moon / d_sun) * cos(beta_sun)
                                    * sin(lam_sun - lam_moon)   [radians->deg]
            beta_H = (d_moon / d_sun) * beta_sun

        The correction terms are the parallactic difference between viewing the
        Sun from the Earth and from the Moon: at most about 0.15 deg, which moves
        the terminator by ~1 km on the surface.  Small, but free.

        Also returns `colongitude` = (90 - l0) mod 360, the quantity lunar
        observers actually quote, and the sub-solar unit vector in the
        Moon-fixed frame.
    '''
    dt = _parse_dt(time_iso)
    t = _centuries(dt)
    dpsi, eps = _nutation_obliquity(t)

    dec_m, gha_m = gp_dec_gha(body_gp("Moon", time_iso))
    ra_m = (gast_deg(dt) - gha_m) % 360.0
    lam_m, _ = _equatorial_to_ecliptic(ra_m, dec_m, eps)

    dec_s, gha_s = gp_dec_gha(body_gp("Sun", time_iso))
    ra_s = (gast_deg(dt) - gha_s) % 360.0
    lam_s, beta_s = _equatorial_to_ecliptic(ra_s, dec_s, eps)

    from .astro import body_distance_km
    d_moon = body_distance_km("Moon", time_iso)
    d_sun = _sun_distance_km(t)
    ratio = d_moon / d_sun

    lam_h = (lam_s + 180.0
             + np.degrees(ratio * np.cos(np.radians(beta_s))
                          * np.sin(np.radians(lam_s - lam_m)))) % 360.0
    beta_h = ratio * beta_s

    l0, b0, _ = libration_from_ecliptic(lam_h, beta_h, ra_m, t, dpsi, eps)
    return dict(lon=l0, lat=b0, colongitude=(90.0 - l0) % 360.0,
                vec=_seleno_vec(l0, b0), lam_sun=lam_s, lam_moon=lam_m)


def _sun_distance_km(t):
    ''' Earth-Sun distance (km) from the low-precision solar theory (Meeus 25).
        Only the 1.7% eccentricity variation matters here. '''
    m = np.radians((357.52911 + 35999.05029 * t) % 360.0)
    e = 0.016708634 - 0.000042037 * t
    v = m + np.radians((1.914602 - 0.004817 * t) * np.sin(m)
                       + 0.019993 * np.sin(2 * m) + 0.000289 * np.sin(3 * m))
    r_au = 1.000001018 * (1 - e * e) / (1 + e * np.cos(v))
    return r_au * 1.495978707e8


def _seleno_vec(lon_deg, lat_deg):
    lo, la = np.radians(lon_deg), np.radians(lat_deg)
    return np.array([np.cos(la) * np.sin(lo), np.sin(la), np.cos(la) * np.cos(lo)])


# --------------------------------------------------------------------------- #
# The render
# --------------------------------------------------------------------------- #

def _libration_R(lib_lon, lib_lat):
    ''' Moon-fixed -> view frame, where view +z points at the observer, +y is the
        Moon's north pole and +x is selenographic EAST (which appears on the
        RIGHT of a north-up view, i.e. toward celestial west). '''
    a, b = np.radians(lib_lon), np.radians(lib_lat)
    Ry = np.array([[np.cos(a), 0, -np.sin(a)], [0, 1, 0], [np.sin(a), 0, np.cos(a)]])
    Rx = np.array([[1, 0, 0], [0, np.cos(b), -np.sin(b)], [0, np.sin(b), np.cos(b)]])
    return Rx @ Ry


def render(size=512, radius_frac=1.0, libration=(0.0, 0.0), subsolar=None,
           rotation_deg=0.0, texture=None, flip_east=False, sky=0.0,
           gain=1.0, gamma=1.0):
    ''' Orthographic render of the illuminated Moon.

        libration    : (lon, lat) of the sub-observer point, degrees.
        subsolar     : (lon, lat) of the sub-solar point, degrees.  None -> the
                       whole disk lit (useful for pure geometry checks).
        rotation_deg : in-image rotation applied to the disk, counter-clockwise,
                       i.e. the position angle of the Moon's north pole in the
                       rendered frame.  Pass P for a north-up sky view, or
                       P + q + roll to imitate a hand-held camera.
        flip_east    : mirror the disk left-right.  A view through a telescope
                       diagonal is mirrored; a camera is not.  Exposed so the
                       matcher can test the hypothesis rather than assume it.

        Returns (image, (cx, cy, r_px)).  Pixels outside the disk are `sky`;
        pixels on the unlit side go to `sky` as well, which is what makes the
        terminator a real edge in the array.
    '''
    c = (size - 1) / 2.0
    r_px = radius_frac * (size / 2.0)
    yy, xx = np.mgrid[0:size, 0:size].astype(float)
    img = render_on_grid(xx, yy, c, c, r_px, libration=libration,
                         subsolar=subsolar, rotation_deg=rotation_deg,
                         texture=texture, flip_east=flip_east, sky=sky,
                         gain=gain, gamma=gamma)
    return img, (c, c, r_px)


def render_on_grid(xx, yy, cx, cy, r_px, libration=(0.0, 0.0), subsolar=None,
                   rotation_deg=0.0, texture=None, flip_east=False, sky=0.0,
                   gain=1.0, gamma=1.0):
    ''' The same render, evaluated at arbitrary pixel coordinates.

        `render` centres the disk in its own square array, which forces whoever
        wants it aligned with a PHOTOGRAPH to paste it at an integer offset --
        and half a pixel of paste error is twice the tie-point residual this
        study is trying to measure.  Evaluating directly on the photograph's grid
        removes the resampling and the rounding together: `cx`, `cy` and `r_px`
        may be fractional and are honoured exactly.
    '''
    tex = load_texture() if texture is None else texture
    nrow, ncol = tex.shape

    X = (xx - cx) / r_px
    Y = -(yy - cy) / r_px                      # image row 0 is the top
    if flip_east:
        X = -X

    # undo the in-image rotation to get disk coordinates with Moon north up
    th = np.radians(rotation_deg)
    ct, st = np.cos(th), np.sin(th)
    xd = X * ct + Y * st
    yd = -X * st + Y * ct

    rr = xd * xd + yd * yd
    disk = rr <= 1.0
    zd = np.sqrt(np.clip(1.0 - rr, 0.0, None))

    # view frame -> Moon-fixed frame
    R = _libration_R(*libration)
    v = np.stack([xd, yd, zd], axis=-1)
    f = v @ R                                   # == (R.T @ v) per pixel

    lat = np.degrees(np.arcsin(np.clip(f[..., 1], -1, 1)))
    lon = np.degrees(np.arctan2(f[..., 0], f[..., 2]))

    # equirectangular lookup: lon 0 at the centre column, lat +90 at row 0
    col = (lon + 180.0) / 360.0 * (ncol - 1)
    row = (90.0 - lat) / 180.0 * (nrow - 1)
    albedo = _bilinear(tex, np.clip(col, 0, ncol - 1), np.clip(row, 0, nrow - 1))

    if subsolar is None:
        img = np.where(disk, albedo, sky)
    else:
        s = _seleno_vec(*subsolar)
        mu0 = f @ s                             # cos(incidence)
        mu = zd                                 # cos(emission)
        lit = disk & (mu0 > 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            shade = np.where(lit, mu0 / np.maximum(mu0 + mu, 1e-6), 0.0)
        img = np.where(lit, albedo * 2.0 * shade, sky)

    if gain != 1.0:
        img = img * gain
    if gamma != 1.0:
        img = np.clip(img, 0, None) ** gamma
    return img


def _bilinear(g, x, y):
    x0 = np.clip(np.floor(x).astype(int), 0, g.shape[1] - 2)
    y0 = np.clip(np.floor(y).astype(int), 0, g.shape[0] - 2)
    fx, fy = x - x0, y - y0
    return (g[y0, x0] * (1 - fx) * (1 - fy) + g[y0, x0 + 1] * fx * (1 - fy)
            + g[y0 + 1, x0] * (1 - fx) * fy + g[y0 + 1, x0 + 1] * fx * fy)


def render_for(time_iso: str, obs_lat: float, obs_lon: float, size=512,
               radius_frac=1.0, roll_deg=0.0, include_parallactic=True,
               obs_height_m=0.0, **kw):
    ''' Render the Moon as it appears from a place at a time.

        Everything the render needs -- libration, axis position angle, sub-solar
        point -- is computed here from the ephemeris; the caller supplies only
        the observer and the camera roll.  With `include_parallactic` the disk is
        rotated by (P - q), which is the orientation in a frame whose "up" is the
        observer's zenith -- i.e. what a levelled camera records.  See
        `sky_rotation` for why the sign is a MINUS.

        Returns (image, geometry_dict).
    '''
    lib = topocentric_libration(time_iso, obs_lat, obs_lon, obs_height_m)
    sun = subsolar_point(time_iso)
    q = parallactic_angle(time_iso, obs_lat, obs_lon) if include_parallactic else 0.0
    rot = lib["pole_pa"] - q + roll_deg
    img, geom = render(size=size, radius_frac=radius_frac,
                       libration=(lib["lon"], lib["lat"]),
                       subsolar=(sun["lon"], sun["lat"]),
                       rotation_deg=rot, **kw)
    return img, dict(libration=lib, subsolar=sun, parallactic=q, rotation=rot,
                     geometry=geom)


def parallactic_angle(time_iso: str, obs_lat: float, obs_lon: float,
                      body: str = "Moon"):
    ''' Parallactic angle q (deg): the POSITION ANGLE OF THE ZENITH at the body,
        measured from celestial north through east. '''
    dec, gha = gp_dec_gha(body_gp(body, time_iso))
    H = np.radians((gha + obs_lon) % 360.0)
    d, p = np.radians(dec), np.radians(obs_lat)
    return float(np.degrees(np.arctan2(
        np.sin(H), np.tan(p) * np.cos(d) - np.sin(d) * np.cos(H))))


def sky_rotation(time_iso: str, obs_lat: float, obs_lon: float,
                 roll_deg: float = 0.0):
    ''' In-image counter-clockwise angle of the Moon's north pole, for a camera
        whose "up" is the zenith.  Equals P - q + roll.

        WHY MINUS.  Both angles are position angles measured from celestial north
        through EAST, and in an un-mirrored photograph of the sky east lies
        counter-clockwise from north.  The zenith sits at PA q, so with the
        zenith placed at image-up, celestial north falls at image angle -q; the
        Moon's pole, at PA P from north, therefore lands at P - q.

        The plus sign is the natural-looking mistake and it is wrong: it makes
        the recovered camera roll drift by 2q, up to ~90 deg across a night,
        which then gets absorbed as a bogus "roll" and hides the real one.  This
        was settled here by computing the image angle of the pole directly from
        3-vectors rather than by trusting the algebra -- the same check lives in
        `test_imu_fusion.TestLunarTexture`.
    '''
    from .lunar_geometry import geocentric_libration
    P = geocentric_libration(time_iso)["pole_pa"]
    return P - parallactic_angle(time_iso, obs_lat, obs_lon) + roll_deg
