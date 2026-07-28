''' Put selenographic coordinates on a photograph of the Moon.

    Given a photograph, an observing site and an approximate time, this fits the
    limb, matches the disk against a Stellarium-textured render, solves for the
    disk geometry from crater tie-points (`lunar_match`), and then draws:

      * the SUB-EARTH POINT -- the selenographic longitude and latitude the
        camera was looking straight down at;
      * a selenographic graticule, so any pixel can be read as a coordinate;
      * named craters and maria from `lunar_features`;
      * the terminator, computed from the sub-solar point rather than from a
        phase angle, so it is a real great circle on a real sphere;
      * a magnified comparison of the disk CENTRE and of the TERMINATOR, side by
        side with the render -- the two places worth looking at closely.

    Run it:

        python -m imu_fusion.tools.annotate_moon PHOTO.JPG \\
            --time 2026-07-29T01:20:00 --lat 41.0082 --lon 28.9784 \\
            --out imu_fusion/results/fig_moon_match.png

    The rotation and the libration are MEASURED, not assumed, so the overlay
    landing on the right craters is a real check of the ephemeris chain and not
    a tautology.  The printed summary reports the ephemeris-minus-measured
    difference, which is the number that says whether it worked.

    (c) 2026.  MIT License (see LICENSE file).
'''

import argparse
import os

import numpy as np

from ..lunar_features import FEATURES
from ..lunar_geometry import features_near_centre
from ..lunar_match import project, solve_photo
from ..lunar_texture import render, parallactic_angle, _seleno_vec

GRID = "#39d0ff"
TERM = "#ffb300"
HOT = "#ff3b6b"
COOL = "#8ef58e"


def load_gray(path):
    from PIL import Image
    return np.asarray(Image.open(path).convert("L"), dtype=float)


def terminator_points(sun_lon, sun_lat, n=721):
    ''' The great circle 90 deg from the sub-solar point. '''
    s = _seleno_vec(sun_lon, sun_lat)
    a = np.cross(s, [0.0, 1.0, 0.0])
    if np.linalg.norm(a) < 1e-6:
        a = np.cross(s, [1.0, 0.0, 0.0])
    a = a / np.linalg.norm(a)
    b = np.cross(s, a)
    t = np.linspace(0, 2 * np.pi, n)
    p = np.outer(np.cos(t), a) + np.outer(np.sin(t), b)
    return (np.degrees(np.arctan2(p[:, 0], p[:, 2])),
            np.degrees(np.arcsin(np.clip(p[:, 1], -1, 1))))


def _polyline(ax, lons, lats, lib, rot, cx, cy, R, **kw):
    xs, ys, vs = zip(*[project(a, b, lib, rot, cx, cy, R)
                       for a, b in zip(lons, lats)])
    xs = np.array(xs, float)
    xs[np.array(vs) <= 0] = np.nan
    ax.plot(xs, np.array(ys), **kw)


def _draw_overlay(ax, lib, rot, cx, cy, R, sun, labels, label_frac=0.97,
                  min_sep_px=34.0, fontsize=6.4):
    for lon0 in range(-180, 180, 15):
        la = np.linspace(-88, 88, 179)
        _polyline(ax, np.full_like(la, lon0), la, lib, rot, cx, cy, R,
                  color=GRID, lw=0.4, alpha=0.5)
    for lat0 in range(-75, 90, 15):
        lo = np.linspace(-180, 180, 361)
        _polyline(ax, lo, np.full_like(lo, lat0), lib, rot, cx, cy, R,
                  color=GRID, lw=0.4, alpha=0.5)
    tl, tb = terminator_points(*sun)
    _polyline(ax, tl, tb, lib, rot, cx, cy, R, color=TERM, lw=2.2)

    x0, y0, _ = project(lib[0], lib[1], lib, rot, cx, cy, R)
    ax.plot([x0], [y0], "+", color=HOT, ms=20, mew=2.4)
    ax.plot([x0], [y0], "o", mfc="none", mec=HOT, ms=30, mew=1.5)

    if not labels:
        return
    near = {f["name"] for f in features_near_centre(*lib, 0.30)}
    placed = []
    for name, lo, la, _, _ in FEATURES:
        x, y, vis = project(lo, la, lib, rot, cx, cy, R)
        if vis <= 0 or (x - cx) ** 2 + (y - cy) ** 2 > (label_frac * R) ** 2:
            continue
        hot = name in near
        # Crowding is worst exactly where the interesting features are, so the
        # separation test applies to the near-centre labels too -- letting them
        # through unconditionally turns the middle of the disk into a smear.
        sep = min_sep_px * (0.7 if hot else 1.0)
        if any((x - px) ** 2 + (y - py) ** 2 < sep ** 2 for px, py in placed):
            continue
        placed.append((x, y))
        c = HOT if hot else COOL
        ax.plot([x], [y], ".", color=c, ms=5 if hot else 3.2)
        ax.annotate(name, (x, y), xytext=(5, 4), textcoords="offset points",
                    color=c, fontsize=fontsize + (1.2 if hot else 0.0),
                    fontweight="bold" if hot else "normal")


def annotate(photo, time_iso, obs_lat, obs_lon, out, title=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    gray = load_gray(photo)
    sol = solve_photo(gray, time_iso, obs_lat, obs_lon)
    cx, cy, R, rot, ll, lb = sol["fit"]["params"]
    lib = (ll, lb)
    sund = sol["subsolar"]
    sun = (sund["lon"], sund["lat"])
    eph = sol["ephem_libration"]
    q = parallactic_angle(time_iso, obs_lat, obs_lon)

    ref, _ = render(size=int(2 * R * 1.12), radius_frac=1.0 / 1.12, libration=lib,
                    subsolar=sun, rotation_deg=rot)
    pad = 1.12 * R
    ext = [cx - pad, cx + pad, cy + pad, cy - pad]

    fig = plt.figure(figsize=(16.0, 11.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.55, 1.0], hspace=0.10, wspace=0.05)
    ax_p = fig.add_subplot(gs[0, 0])
    ax_r = fig.add_subplot(gs[0, 1])
    ax_q = fig.add_subplot(gs[0, 2])
    zooms = [fig.add_subplot(gs[1, i]) for i in range(3)]

    ax_p.imshow(gray, cmap="gray", origin="upper")
    _draw_overlay(ax_p, lib, rot, cx, cy, R, sun, True)
    ax_p.set_title("photograph + selenographic overlay", fontsize=10)

    ax_r.imshow(ref, cmap="gray", origin="upper", extent=ext)
    _draw_overlay(ax_r, lib, rot, cx, cy, R, sun, False)
    ax_r.set_title("Stellarium albedo map rendered for this geometry", fontsize=10)

    # tie-point residuals, magnified
    ax_q.imshow(gray, cmap="gray", origin="upper", alpha=0.45)
    res = sol["fit"]["residuals"]
    for k, t in enumerate(sol["ties"]):
        rx, ry = res[2 * k], res[2 * k + 1]
        ax_q.arrow(t["x"], t["y"], -rx * 120, -ry * 120, color=HOT,
                   width=1.2, head_width=6, length_includes_head=True)
        ax_q.plot([t["x"]], [t["y"]], ".", color=COOL, ms=3)
    ax_q.plot([cx + 0.55 * R], [cy + 0.86 * R], ".", color=COOL, ms=3)
    ax_q.arrow(cx + 0.55 * R, cy + 0.86 * R, 120, 0, color=HOT, width=1.2,
               head_width=6, length_includes_head=True)
    ax_q.text(cx + 0.55 * R, cy + 0.80 * R, "1 px residual, x120",
              color="w", fontsize=8)
    ax_q.set_title(f"tie-point residuals — {sol['fit']['rms_px']:.3f} px rms",
                   fontsize=10)

    for ax in (ax_p, ax_r, ax_q):
        ax.set_xlim(cx - pad, cx + pad)
        ax.set_ylim(cy + pad, cy - pad)
        ax.set_xticks([]); ax.set_yticks([])

    # --- magnified: the disk centre (photo | render) and the terminator ------
    zc = 0.22 * R
    zx, zy, _ = project(lib[0], lib[1], lib, rot, cx, cy, R)
    # terminator crop: the point of the terminator nearest the disk centre
    tl, tb = terminator_points(*sun)
    pts = [project(a, b, lib, rot, cx, cy, R) for a, b in zip(tl, tb)]
    pts = [(x, y) for x, y, v in pts if v > 0]
    tx, ty = min(pts, key=lambda p: (p[0] - cx) ** 2 + (p[1] - cy) ** 2)

    specs = [(gray, None, zx, zy, zc, "photograph — disk centre"),
             (ref, ext, zx, zy, zc, "render — disk centre"),
             (gray, None, tx, ty, 0.30 * R, "photograph — terminator")]
    for ax, (im, extent, ax_, ay_, half, ttl) in zip(zooms, specs):
        if extent is None:
            ax.imshow(im, cmap="gray", origin="upper")
        else:
            ax.imshow(im, cmap="gray", origin="upper", extent=extent)
        _draw_overlay(ax, lib, rot, cx, cy, R, sun, "centre" in ttl,
                      label_frac=1.30, min_sep_px=16.0, fontsize=8.0)
        ax.set_xlim(ax_ - half, ax_ + half)
        ax.set_ylim(ay_ + half, ay_ - half)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(ttl, fontsize=10)

    dlon, dlat = ll - eph["lon"], lb - eph["lat"]
    head = title or os.path.basename(photo)
    fig.suptitle(
        f"{head}    site {obs_lat:.4f}N {obs_lon:.4f}E    ephemeris epoch {time_iso}Z\n"
        f"sub-Earth point MEASURED {ll:+.2f}°E {lb:+.2f}°N   "
        f"(ephemeris {eph['lon']:+.2f}, {eph['lat']:+.2f}; "
        f"difference {dlon:+.2f}°, {dlat:+.2f}°)      "
        f"sub-solar {sun[0]:+.2f}°E {sun[1]:+.2f}°N, colongitude {sund['colongitude']:.1f}°\n"
        f"disk radius {R:.1f} px    in-image rotation {rot:+.2f}°    "
        f"axis PA P={eph['pole_pa']:+.2f}°  parallactic q={q:+.2f}°  →  "
        f"implied camera roll {((rot - (eph['pole_pa'] - q) + 180) % 360) - 180:+.1f}°     "
        f"disk NCC {sol['disk_ncc']:.4f},  {len(sol['ties'])} tie-points, "
        f"{sol['fit']['rms_px']:.3f} px rms",
        fontsize=9.5)
    fig.tight_layout(rect=[0, 0, 1, 0.935])
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, dpi=125)
    plt.close(fig)
    sol["out"] = out
    return sol


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("photo")
    p.add_argument("--time", required=True, help="UTC, ISO 8601")
    p.add_argument("--lat", type=float, required=True)
    p.add_argument("--lon", type=float, required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--title")
    a = p.parse_args(argv)
    r = annotate(a.photo, a.time, a.lat, a.lon, a.out, title=a.title)
    cx, cy, R, rot, ll, lb = r["fit"]["params"]
    e = r["ephem_libration"]
    print(f"tie-points {len(r['ties'])}   rms {r['fit']['rms_px']:.3f} px"
          f"   disk NCC {r['disk_ncc']:.4f}")
    print(f"sub-Earth point measured {ll:+.3f}E {lb:+.3f}N"
          f"   ephemeris {e['lon']:+.3f}E {e['lat']:+.3f}N"
          f"   d=({ll - e['lon']:+.3f}, {lb - e['lat']:+.3f})")
    print(f"rotation {rot:+.3f} deg    radius {R:.2f} px")
    print(f"wrote {r['out']}")


if __name__ == "__main__":
    main()
