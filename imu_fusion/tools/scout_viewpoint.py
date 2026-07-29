''' Before you drive there: can this viewpoint fix your position?

    Every terrain-resection failure in this study was predictable from the
    geometry alone, and both took hours of grid search to discover instead.  This
    asks the DEM the question in about a second:

        python -m imu_fusion.tools.scout_viewpoint --lat 40.9046 --lon 29.2094 \\
            --sector 126 144 --sigma 0.14

    It ray-marches the sector, records how far away each skyline feature actually
    is, and rules on the spread -- because after a compass bias and a pitch bias
    have absorbed the two first-order signals, RANGE SPREAD is the only thing
    left carrying position.  A tack-sharp view of one distant ridge scores worse
    than a soft view that also contains something close.

    Use it to choose where to stand and which way to look, and to decide whether
    a second sector is needed.  It reports the error ellipse, so a 0.3 x 16 km
    answer tells you plainly that you have a line of position and must cross it.

    (c) 2026.  MIT License (see LICENSE file).
'''

import argparse

from ..terrain_resection import DemTiles
from ..resection_geometry import scene_report


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--lat", type=float, required=True)
    p.add_argument("--lon", type=float, required=True)
    p.add_argument("--sector", nargs=2, type=float, required=True,
                   metavar=("AZ_START", "AZ_END"))
    p.add_argument("--eye", type=float, default=None,
                   help="camera height in metres AMSL (default: DEM + 1.6)")
    p.add_argument("--sigma", type=float, default=0.05,
                   help="skyline match accuracy in degrees; use what you "
                        "actually achieve, not the pixel scale (real "
                        "photographs here gave 0.019 to 0.14)")
    p.add_argument("--max-km", type=float, default=150.0)
    p.add_argument("--near", nargs=3, type=float, action="append", default=[],
                   metavar=("LAT", "LON", "HEIGHT_M"),
                   help="a landmark the DEM does not contain (tower, mast, "
                        "lighthouse). Repeatable. Including something CLOSE is "
                        "the single highest-leverage change you can make.")
    p.add_argument("--dem", default="imu_fusion/dem")
    a = p.parse_args(argv)

    dem = DemTiles(a.dem)
    eye = a.eye
    if eye is None:
        import numpy as np
        eye = float(dem.elevation(np.array([a.lat]), np.array([a.lon]))[0]) + 1.6
    r = scene_report(dem, a.lat, a.lon, eye, a.sector[0], a.sector[1],
                     sigma_deg=a.sigma, d_max_km=a.max_km,
                     extra_landmarks=[tuple(v) for v in a.near])

    print(f"viewpoint {a.lat:.5f} {a.lon:.5f}, eye {eye:.0f} m, "
          f"sector {a.sector[0]:.0f}-{a.sector[1]:.0f} deg")
    print(f"relief across the sector: {r['relief_deg']:.2f} deg")
    if not r["features"]:
        print("  " + r["message"])
        return
    print(f"\n{'azimuth':>9}{'range km':>10}{'height m':>10}{'elev deg':>10}")
    for f in sorted(r["features"], key=lambda q: q["azimuth"]):
        tag = "  (given)" if f.get("man_made") else ""
        print(f"{f['azimuth']:9.1f}{f['distance_km']:10.1f}"
              f"{f['height_m']:10.0f}{f['elevation_deg']:10.2f}{tag}")
    s = r["sensitivity"]
    print(f"\nrange spread {min(s['ranges_km']):.0f}-{max(s['ranges_km']):.0f} km")
    print(f"  lateral {s['lateral_raw']:.3f} -> {s['lateral_absorbed']:.3f} deg/km "
          f"once a compass bias is free")
    print(f"  radial  {s['radial_raw']:.4f} -> {s['radial_absorbed']:.4f} deg/km "
          f"once a pitch bias is free")
    print(f"\nat {a.sigma:.3f} deg match accuracy: "
          f"{r['across_km']:.2f} km across x {r['along_km']:.1f} km along")
    print(f"VERDICT: {r['message']}")


if __name__ == "__main__":
    main()
