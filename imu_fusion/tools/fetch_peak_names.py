''' Attach OpenStreetMap peak names to DEM-derived skyline summits.

    `terrain_resection.skyline_peaks` finds summits geometrically — it knows a
    summit's bearing, elevation angle, distance and height, but not what it is
    CALLED.  Names come from a gazetteer, and this fetches them from
    OpenStreetMap via the Overpass API.

    NETWORK REQUIRED.  Overpass and the GeoNames dumps are unreachable from the
    sandbox this study was written in, so the naming step has never been run
    here; the panorama viewer ships with names only for summits that could be
    identified independently.  Run this yourself where the network allows.

        python -m imu_fusion.tools.fetch_peak_names --bbox 36.6 27.0 37.2 27.9 \\
            --out imu_fusion/results/peak_names.json

    Then match them to the summits of a viewpoint:

        from imu_fusion.tools.fetch_peak_names import load, name_summits
        named = name_summits(summits, load("imu_fusion/results/peak_names.json"))

    OSM data is © OpenStreetMap contributors, ODbL.  Attribute it if you ship it.

    (c) 2026.  MIT License (see LICENSE file).
'''

import argparse
import json
import math
import os

OVERPASS = "https://overpass-api.de/api/interpreter"

# natural=peak covers summits; natural=volcano and natural=saddle are useful too.
_QUERY = """[out:json][timeout:180];
(
  node["natural"="peak"]({s},{w},{n},{e});
  node["natural"="volcano"]({s},{w},{n},{e});
);
out body;"""


def fetch(south, west, north, east, url: str = OVERPASS, timeout: int = 180):
    ''' Query Overpass for named peaks in a bounding box.

        Returns [{name, lat, lon, ele_m}], skipping unnamed nodes.  Raises on a
        network or HTTP failure — the caller decides whether that is fatal.
    '''
    import urllib.request
    import urllib.parse
    q = _QUERY.format(s=south, w=west, n=north, e=east)
    data = urllib.parse.urlencode({"data": q}).encode()
    req = urllib.request.Request(url, data=data,
                                 headers={"User-Agent": "celestial-navigation/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        payload = json.loads(r.read().decode())
    out = []
    for el in payload.get("elements", []):
        tags = el.get("tags", {})
        name = tags.get("name") or tags.get("name:en")
        if not name:
            continue
        ele = tags.get("ele")
        try:
            ele = float(str(ele).split()[0]) if ele else None
        except ValueError:
            ele = None
        out.append(dict(name=name, lat=el["lat"], lon=el["lon"], ele_m=ele))
    return out


def load(path: str):
    ''' Read a previously fetched peak list; [] if the file is absent. '''
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return json.load(f)


# --------------------------------------------------------------------------- #
# Matching gazetteer entries to DEM summits
# --------------------------------------------------------------------------- #

def summit_position(view_lat, view_lon, bearing_deg, distance_km):
    ''' Where a summit seen at this bearing/distance actually is (great circle). '''
    R = 6371.0
    br, d = math.radians(bearing_deg), distance_km / R
    p1, l1 = math.radians(view_lat), math.radians(view_lon)
    p2 = math.asin(math.sin(p1) * math.cos(d) + math.cos(p1) * math.sin(d) * math.cos(br))
    l2 = l1 + math.atan2(math.sin(br) * math.sin(d) * math.cos(p1),
                         math.cos(d) - math.sin(p1) * math.sin(p2))
    return math.degrees(p2), (math.degrees(l2) + 540) % 360 - 180


def _sep_km(a_lat, a_lon, b_lat, b_lon):
    R = 6371.0
    p1, p2 = math.radians(a_lat), math.radians(b_lat)
    dp, dl = p2 - p1, math.radians(b_lon - a_lon)
    h = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * R * math.asin(min(1.0, math.sqrt(h)))


def name_summits(summits, gazetteer, view_lat, view_lon,
                 tolerance_km: float = 1.0, height_tolerance_m: float = 250.0):
    ''' Attach a `name` to each summit dict that has one nearby in the gazetteer.

        summits : [{az, el, dist, h}, ...] as produced for the panorama viewer.
        Matching is by GROUND POSITION, not bearing alone: the summit's lat/lon
        is reconstructed from its bearing and distance, then the nearest named
        peak within `tolerance_km` (and, if the gazetteer gives an elevation,
        within `height_tolerance_m`) wins.

        The tolerance matters: a DEM summit is the highest CELL of a massif,
        which can sit a few hundred metres from the surveyed spot height, and
        1 arc-second is ~30 m on the ground.  Too tight and nothing matches; too
        loose and a neighbouring peak steals the name.
    '''
    out = []
    for s in summits:
        s = dict(s)
        lat, lon = summit_position(view_lat, view_lon, s["az"], s["dist"])
        best, best_d = None, tolerance_km
        for g in gazetteer:
            d = _sep_km(lat, lon, g["lat"], g["lon"])
            if d > best_d:
                continue
            if (g.get("ele_m") is not None and s.get("h") is not None
                    and abs(g["ele_m"] - s["h"]) > height_tolerance_m):
                continue
            best, best_d = g, d
        if best:
            s["name"] = best["name"]
            s["match_km"] = round(best_d, 3)
        out.append(s)
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description="Fetch OSM peak names for a bounding box")
    p.add_argument("--bbox", nargs=4, type=float, required=True,
                   metavar=("SOUTH", "WEST", "NORTH", "EAST"))
    p.add_argument("--out", required=True)
    p.add_argument("--url", default=OVERPASS)
    args = p.parse_args(argv)
    peaks = fetch(*args.bbox, url=args.url)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(peaks, f, indent=1)
    print(f"{len(peaks)} named peaks -> {args.out}")


if __name__ == "__main__":
    main()
