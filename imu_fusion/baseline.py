''' Incumbent-tool baselines, for context against the factor graph.

    Reuses the repository's own `starfix` engine on the SAME synthetic
    measurements the factor graph consumes:

      * `starfix_single_fix`   -- one epoch's Sun+Moon fix via
        SightCollection.get_intersections (the classic two-LOP intersection,
        no IMU, no multi-epoch smoothing).
      * `starfix_mc_sigma`     -- the engine's Monte-Carlo instability metric
        (SightCollection.get_intersections_mc), the incumbent accuracy proxy.

    These have no IMU and no cross-epoch fusion, so they are the "before"
    picture the factor graph improves on.

    (c) 2026.  MIT License (see LICENSE file).
'''

from starfix import (Sight, SightCollection, LatLonGeodetic, IntersectError)

from .astro import great_circle_km


def _dms(angle_deg: float) -> str:
    ''' Format a decimal degree angle as a "D:M:S" string for starfix. '''
    sign = "-" if angle_deg < 0 else ""
    a = abs(angle_deg)
    d = int(a)
    m = int((a - d) * 60)
    s = (a - d - m / 60) * 3600
    return f"{sign}{d}:{m}:{s:.3f}"


def _sight(obs, drp):
    ''' Build a starfix Sight from a synthetic Observation (geometric altitude,
        so refraction/dip/parallax/semidiameter are all disabled). '''
    return Sight(object_name=obs.body,
                 set_time=obs.time_iso,
                 measured_alt=_dms(obs.meas_alt),
                 estimated_position=drp,
                 ho_obs=True, limb_correction=0, horizontal_parallax=0)


def starfix_single_fix(keyframe, anchor_lat, anchor_lon):
    ''' Deterministic Sun+Moon fix for one keyframe via starfix.  Returns
        (error_km, lat, lon) or (None, None, None) if it cannot intersect. '''
    drp = LatLonGeodetic(anchor_lat, anchor_lon)
    try:
        sights = [_sight(o, drp) for o in keyframe.observations]
        coll = SightCollection(sights)
        res = coll.get_intersections(return_geodetic=True,
                                     estimated_position=drp)
    except (IntersectError, ValueError, AssertionError):
        return None, None, None
    coord = res[0] if isinstance(res, tuple) else res
    # coord may itself be a pair of candidate points; pick the first.
    if isinstance(coord, tuple):
        coord = coord[0]
    lat, lon = coord.get_lat(), coord.get_lon()
    return great_circle_km(lat, lon, keyframe.true_lat, keyframe.true_lon), \
        lat, lon
