''' What the camera can actually SEE, as opposed to what the DEM says is there.

    A skyline render marches until the earth curves away, so it answers *what is
    GEOMETRICALLY visible*.  A photograph answers *what is ATMOSPHERICALLY
    visible*, and on a hazy day the two differ by tens of kilometres.

    THE FAILURE THIS MODULE EXISTS FOR.  On the Istanbul ultrawide the render put
    the Samanli mountains on the far shore of the Marmara -- 680 m at 39.5 km --
    into the skyline at +50 arcmin.  That is correct geometry: a 680 m peak
    clears a 5 m eye's horizon out to 100 km.  The photograph shows nothing there
    but sea.  The extractor therefore traced the WATER horizon, and every one of
    those samples entered the fit as a 50 arcmin error.  The residual was 32.8
    arcmin and no position could lower it, because the error was in the forward
    model, not in the position.

    THE PHYSICS.  Koschmieder: a dark object seen against the horizon sky has
    apparent contrast

        C(d) = C0 * exp(-tau)

    and is detectable while C exceeds a threshold eps (0.02 by the meteorological
    convention, which is what makes the "meteorological visual range" V the
    distance at which a black object against the sky disappears:
    tau = sigma*V = ln(1/0.02) = 3.912).

    WHY A FLAT `d_max` CAP IS NOT ENOUGH.  Aerosol is concentrated near the
    ground -- scale height of order a kilometre -- so a sight line that CLIMBS to
    a summit accumulates far less optical depth than one that runs horizontally
    through the murk.  A 680 m peak at 40 km and a 5 m islet at 40 km are not
    equally visible, and a single `d_max` cannot tell them apart.  This module
    integrates the extinction along the actual slant path, so the summit is
    correctly favoured over the islet.

    HONEST LIMITS.  This is a grey-atmosphere, single-scattering, horizontally
    homogeneous model.  It has nothing to say about a bank of sea fog on one
    bearing and clear air on the next, about the strong wavelength dependence of
    aerosol scattering (a red channel sees notably further than a blue one), or
    about terrain that is bright rather than dark against the sky.  Its job is to
    keep a forward model from asserting mountains the camera never recorded, and
    for that a crude model beats none at all.

    (c) 2026.  MIT License (see LICENSE file).
'''

import math

import numpy as np

# Koschmieder's constant: ln(1/0.02).  The optical depth at which a black object
# against the horizon sky reaches the conventional 2% contrast threshold, which
# is what DEFINES the meteorological visual range.
KOSCHMIEDER = 3.912023005428146

# Contrast threshold.  0.02 is the meteorological convention, set by the human
# eye.  A camera with a clean sensor and post-processing does better, so pass a
# smaller value if the extractor is finding faint edges; the DP extractor here
# locks onto steps of a few units in 255, which is nearer 0.01.
CONTRAST_THRESHOLD = 0.02

# Aerosol scale height, metres.  Boundary-layer haze is far more concentrated
# than the Rayleigh atmosphere: 1-1.5 km is the usual range, against ~8 km for
# molecular scattering.  It is the small number that makes summits punch through
# haze that hides the sea surface at the same distance.
AEROSOL_SCALE_HEIGHT_M = 1200.0


def extinction_per_km(visibility_km):
    ''' Extinction coefficient sigma, per km, from meteorological visual range. '''
    v = float(visibility_km)
    if v <= 0.0:
        raise ValueError("visibility_km must be positive")
    return KOSCHMIEDER / v


def slant_optical_depth(distance_km, eye_m, target_m, visibility_km,
                        scale_height_m=AEROSOL_SCALE_HEIGHT_M):
    ''' Optical depth along a sight line that climbs from `eye_m` to `target_m`.

        The visual range `visibility_km` is quoted for a HORIZONTAL path at the
        surface, so sigma is the sea-level value and the density falls as
        exp(-h/H) along the way.  Treating height as linear in path length (true
        to well under a percent for these geometries),

            tau = sigma * d * (H / (h1 - h0)) * (exp(-h0/H) - exp(-h1/H))

        which collapses to the horizontal result sigma*d*exp(-h0/H) as h1 -> h0.

        Vectorised over any broadcastable combination of the arguments.
    '''
    d = np.asarray(distance_km, dtype=float)
    h0 = np.asarray(eye_m, dtype=float)
    h1 = np.asarray(target_m, dtype=float)
    hs = float(scale_height_m)
    sigma = extinction_per_km(visibility_km)
    dh = h1 - h0
    # the h1 -> h0 limit is exp(-h0/H); use it wherever the height change is
    # small enough that the closed form would lose precision to cancellation
    small = np.abs(dh) < 1e-6
    dh_safe = np.where(small, 1.0, dh)
    mean_density = np.where(
        small,
        np.exp(-h0 / hs),
        (hs / dh_safe) * (np.exp(-h0 / hs) - np.exp(-h1 / hs)))
    return sigma * d * mean_density


def contrast(distance_km, eye_m, target_m, visibility_km,
             scale_height_m=AEROSOL_SCALE_HEIGHT_M, contrast_0=1.0):
    ''' Apparent contrast of dark terrain against the horizon sky. '''
    tau = slant_optical_depth(distance_km, eye_m, target_m, visibility_km,
                              scale_height_m)
    return float(contrast_0) * np.exp(-tau)


def is_detectable(distance_km, eye_m, target_m, visibility_km,
                  scale_height_m=AEROSOL_SCALE_HEIGHT_M, contrast_0=1.0,
                  threshold=CONTRAST_THRESHOLD):
    ''' Boolean mask: would this terrain register in the photograph at all? '''
    return contrast(distance_km, eye_m, target_m, visibility_km,
                    scale_height_m, contrast_0) >= float(threshold)


def detection_range_km(visibility_km, eye_m=0.0, target_m=0.0,
                       scale_height_m=AEROSOL_SCALE_HEIGHT_M, contrast_0=1.0,
                       threshold=CONTRAST_THRESHOLD):
    ''' Range at which terrain of this height stops being detectable.

        For a sea-level target this returns the visual range itself (that is the
        definition).  For a summit it returns considerably more, which is the
        whole point of carrying the slant path.
    '''
    lo, hi = 1e-4, 2000.0
    want = math.log(float(contrast_0) / float(threshold))
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        tau = float(slant_optical_depth(mid, eye_m, target_m, visibility_km,
                                        scale_height_m))
        if tau < want:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def visibility_from_bracket(seen_km, seen_height_m, unseen_km, unseen_height_m,
                            eye_m=0.0, scale_height_m=AEROSOL_SCALE_HEIGHT_M,
                            threshold=CONTRAST_THRESHOLD):
    ''' Bracket the day's visibility from what the photograph does and does not show.

        The only visibility measurement available after the fact is the picture
        itself: something at `seen_km` IS in the frame, something at
        `unseen_km` is NOT.  Returns (v_lo, v_hi) -- the visual range must be at
        least enough to show the first and little enough to hide the second.

        On the Istanbul frame: Buyukada, 169 m at 6.6 km, sharply visible; the
        Samanli mountains, 680 m at 39.5 km, absent.  Returns a wide but real
        bracket, and the solve was insensitive across it -- which is the useful
        thing to know, because it says the exact number does not matter.

        Returns v_hi = inf if the "unseen" object would be invisible at any
        visibility (it never constrains from above), and v_lo = 0 likewise.
    '''
    def v_needed(d, h):
        lo, hi = 0.05, 5000.0
        for _ in range(90):
            mid = 0.5 * (lo + hi)
            c = float(contrast(d, eye_m, h, mid, scale_height_m))
            if c < threshold:
                lo = mid          # too hazy, need more visibility
            else:
                hi = mid
        return 0.5 * (lo + hi)
    v_lo = v_needed(seen_km, seen_height_m)          # at least this clear
    v_hi = v_needed(unseen_km, unseen_height_m)      # less clear than this
    return (float(v_lo), float(v_hi))
