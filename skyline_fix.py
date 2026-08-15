''' Bridge from horizonator skyline fixes to the celestial-navigation toolkit.

    A "skyline fix" is a position estimated by matching a photographed
    island/coast skyline against a digital elevation model (see the
    horizonator repo, experiments/skyfix.py). skyfix writes a JSON result
    with the estimated position and a covariance from the local curvature
    of the match cost. This module converts that JSON into toolkit objects:

    - as_position():  a LatLonGeodetic, usable as the estimated_position
      argument of get_intersections()/SightCollection reductions -- a
      skyline fix is typically far tighter than dead reckoning, so it makes
      an excellent intersection disambiguator and linearization point.
    - as_circle():    a Circle centered on the fix whose radius is the
      n-sigma horizontal uncertainty, for plotting alongside celestial
      circles of position or for approximate fusion by intersection.

    Example:
        from skyline_fix import SkylineFix
        fix = SkylineFix.from_json ("bafa_fix.json")
        intersection, fitness, diag = get_intersections (
            circle1, circle2, estimated_position=fix.as_position ())

    (c) 2026, part of the horizonator skyline-matching study, MIT License
'''

import json
from math import hypot, pi

from starfix import LatLonGeodetic, Circle, EARTH_CIRCUMFERENCE


class SkylineFix:
    ''' A position fix from DEM skyline matching, with uncertainty '''

    def __init__ (self, lat : float, lon : float,
                  sigma_north_m : float, sigma_east_m : float,
                  rms_mrad : float = 0.0):
        self.lat            = lat
        self.lon            = lon
        self.sigma_north_m  = sigma_north_m
        self.sigma_east_m   = sigma_east_m
        self.rms_mrad       = rms_mrad

    @staticmethod
    def from_json (path : str) -> "SkylineFix":
        ''' Load a fix written by horizonator experiments/skyfix.py '''
        with open (path, encoding="utf-8") as f:
            d = json.load (f)
        return SkylineFix (d["lat"], d["lon"],
                           d.get("sigma_n_m", 100.0),
                           d.get("sigma_e_m", 100.0),
                           d.get("rms_mrad", 0.0))

    def as_position (self) -> LatLonGeodetic:
        ''' The fix as a geodetic coordinate '''
        return LatLonGeodetic (lat=self.lat, lon=self.lon)

    def sigma_m (self) -> float:
        ''' Scalar (RMS-combined) horizontal uncertainty in meters '''
        return hypot (self.sigma_north_m, self.sigma_east_m)

    def as_circle (self, n_sigma : float = 2.0) -> Circle:
        ''' The fix as a Circle of position with an n-sigma radius '''
        radius_m  = max (n_sigma * self.sigma_m (), 1.0)
        angle_deg = (radius_m / (EARTH_CIRCUMFERENCE * 1000.0)) * 360.0
        return Circle (self.as_position (), angle_deg)
