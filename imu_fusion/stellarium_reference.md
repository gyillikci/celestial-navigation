<!--- © 2026.  MIT License (see LICENSE file). -->

# Stellarium as the authoritative astronomical source

The project's Sun/Moon astronomy is being re-anchored so that it **comes from
Stellarium** (VSOP87 + ELP2000/MPP02). Stellarium's own Sun/Moon numbers are, in
turn, validated against JPL DE — so this makes JPL-grade positions traceable to
Stellarium the reference the whole team trusts.

Stellarium is a desktop GUI app and **cannot render in the headless CI sandbox**
(it segfaults in the software-GL init there — verified). So the flow is: **you
run one export in Stellarium on your Mac, drop the CSV in, and the Python side
ingests it** as the primary ephemeris (device tables are then baked from it).

## Recommended: run the export script (full fields, guaranteed schema)

The script `imu_fusion/tools/stellarium_export.ssc` writes exactly the CSV the
loader expects, including **distance** (needed for parallax/semidiameter) and
apparent **altitude** (to validate the refraction+parallax chain against
Stellarium). It is **self-describing** — it first prints a `#SCHEMA` block with
the real field names/units of your Stellarium version.

1. In Stellarium: **Configuration (F2) → Plugins → Script Console** (enable it,
   restart if prompted). Or use the menu **Scripts**.
2. Open the **Script Console**, paste the contents of
   `imu_fusion/tools/stellarium_export.ssc`, and press **Run** (▶).
   - Alternatively from a terminal:
     `/Applications/Stellarium.app/Contents/MacOS/stellarium --startup-script /path/to/stellarium_export.ssc`
3. **Smoke-test first.** Before the full run, edit the top of the script to a tiny
   window (e.g. `START_UTC="2026-03-24T00:00:00"`, `END_UTC="2026-03-25T00:00:00"`,
   `STEP_HOURS=6`) and run once. Open the output file and **paste me the
   `#SCHEMA` block** so I can confirm the field names/units for your version.
4. Then set the window you want and run for real. **Keep the step HOURLY**
   (the default): the Moon moves ~0.5°/hr, so linear interpolation needs an
   hourly grid to stay arc-second-accurate (this matches the nautical almanac's
   hourly cadence). Narrow the date window rather than coarsen the step — a full
   year hourly is ~8760 rows/body, which is fine; a voyage of a few days is tiny.
5. The output file path is printed in the Script Console log — typically
   **`~/Library/Application Support/Stellarium/output.txt`** on macOS. Copy it to:

   ```
   imu_fusion/sample_data/stellarium_ephemeris.csv
   ```

### Schema the loader reads

```
utc,body,ra_deg,dec_deg,dist_au,alt_deg,az_deg,elong_deg,phase,size_arcsec
2026-03-24T12:00:00,Sun,358.44,1.5333,0.9959,51.07,169.5,,,1918.6
2026-03-24T12:00:00,Moon,72.19,27.9197,0.00246,26.29,74.5,44.7,0.98,1866.0
```

- `utc` — ISO UTC.  `body` — `Sun`/`Moon`.
- `ra_deg`,`dec_deg` — **apparent, of date** (Stellarium script returns degrees).
- `dist_au` — geocentric distance in AU (Moon → HP/parallax & angular size).
- `alt_deg`,`az_deg` — topocentric apparent (atmosphere ON) for the script's
  observer; used only to witness the correction chain, not for GHA.
- `elong_deg`,`phase`,`size_arcsec` — optional extras (elongation, illuminated
  fraction, angular diameter).  Blank cells are fine.

Lines that are blank, start with `#`, or lack `ra_deg`/`dec_deg` are ignored, so
the `#SCHEMA`/`#END` markers and comments pass through harmlessly.

## Fallback: the AstroCalc Ephemeris GUI (RA/Dec only)

If the script won't run, use **Astronomical Calculations (F10) → Ephemeris**:
choose Sun then Moon, set From/To + step, coordinates **RA/Dec (of date)**,
topocentric **off**, **Calculate → Save CSV**. Reformat to at least
`utc,body,ra_deg,dec_deg` (RA `HH MM SS` → degrees ×15). Distance/altitude will
be absent; the Moon's distance then falls back to the almanac HP until a
full-field export is provided.

## What happens on the Python side

`validate_ephemeris.py` loads the CSV, converts RA-of-date → GHA with an
**engine-free** Greenwich apparent sidereal time (no astropy needed), and — when
`imu_fusion/sample_data/stellarium_ephemeris.csv` is present — makes Stellarium
the **primary** ground-truth engine (`ENGINE == "stellarium"`). `astro.body_gp` /
`body_distance_km` then read Stellarium; `export_golden.py` bakes the device
tables from it. Without the file, the pipeline falls back to the starfix almanac,
so the study stays green until you drop the export in.
