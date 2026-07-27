<!--- © 2026.  MIT License (see LICENSE file). -->

# Stellarium as a ground-truth reference

The study's Sun/Moon ground truth comes from `starfix`'s machine-readable
almanac (Skyfield/JPL under the hood). `validate_ephemeris.py` already
cross-checks it against an **independent** engine (astropy/ERFA — a separate
implementation — or Skyfield/JPL if a kernel is available). You can add
**Stellarium** (VSOP87 + ELP2000/MPP02) as a third, independent witness.

Stellarium is a desktop GUI app, so this can't run inside the headless sandbox —
you export a small CSV from Stellarium once, drop it in, and the validator
compares it.

## Export from Stellarium

1. Open **Stellarium**. Set the location/date if you like (the comparison below
   uses **geocentric** apparent RA/Dec of date, so the observer location does not
   matter for it).
2. Open **Astronomical Calculations** (`F10`) → **Ephemeris** tab.
3. **Celestial body**: choose *Sun* (run once), then *Moon* (run again) — or use
   the two-object mode if your version supports it.
4. Set the **From/To** dates and a **time step** (e.g. 6 hours) covering the
   times you want to check.
5. Make sure coordinates are **RA/Dec (of date)**, and set positions to
   **geocentric** (Configuration → *Tools*/*Astronomical Calculations*: enable
   "topocentric coordinates" **off** — the Moon's geocentric position is what the
   almanac stores).
6. Click **Calculate**, then **Save** / export to CSV.

## CSV format the validator expects

Reformat the export (or hand-enter a few rows) into this simple schema and save
it as `imu_fusion/results/stellarium_reference.csv`:

```
utc,body,ra_deg,dec_deg
2026-03-24 12:00:00,Sun,0.4523,1.5333
2026-03-24 12:00:00,Moon,72.19,27.9197
```

- `utc` — ISO time, UTC (`YYYY-MM-DD HH:MM:SS`).
- `body` — `Sun` or `Moon`.
- `ra_deg`, `dec_deg` — **apparent, geocentric, of date**, in **degrees**
  (convert Stellarium's `HH MM SS` RA to degrees ×15 if needed).

Lines that are blank, start with `#`, or have empty fields are ignored — so the
shipped `stellarium_reference_template.csv` (header only) is skipped until you
add real rows.

## Run the comparison

```python
from imu_fusion.validate_ephemeris import compare_reference_csv
for r in compare_reference_csv("imu_fusion/results/stellarium_reference.csv"):
    print(f"{r['utc']} {r['body']:5} GHA {r['gha_as']:+.1f}\"  Dec {r['dec_as']:+.1f}\"")
```

Residuals should be a handful of arc-seconds (Stellarium's VSOP87/ELP vs the
almanac's Skyfield/JPL), well under the study's arc-minute-level measurement
noise. The validator converts Stellarium's RA-of-date to GHA with the same
Greenwich apparent sidereal time it uses everywhere, so only the ephemeris
differs.
