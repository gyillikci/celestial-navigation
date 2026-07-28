''' A small catalogue of named lunar features, for labelling a matched disk.

    `lunar_orientation.MARIA` has ten maria, which is enough to draw a schematic
    and not enough to say WHAT IS AT THE CENTRE of a photograph -- the disk
    centre lands in the highlands between the maria, so a maria-only catalogue
    answers "nothing" to the one question worth asking.  This adds the named
    craters, bays and smaller seas that a 300-pixel-radius image resolves.

    Coordinates are selenographic longitude (+E) and latitude (+N) of the
    feature centre, and an approximate angular radius in units of the Moon's
    radius (so 0.01 is about 17 km, roughly Triesnecker).

    ACCURACY.  These are catalogue centres rounded to about a tenth of a degree,
    and a few of the mare "centres" are eyeballed -- a mare has no well-defined
    centre.  Good enough to put a label within a crater diameter of the crater;
    NOT an astrometric reference.  If you need better, the USGS Gazetteer of
    Planetary Nomenclature is the authority.

    1 degree of selenographic arc at the sub-observer point is 30.3 km, and about
    1.7% of the apparent disk radius -- so on a photograph with a 322 px radius,
    a degree is 5.6 px.  That sets the standard for whether a label has landed.

    (c) 2026.  MIT License (see LICENSE file).
'''

# (name, lon_E, lat_N, angular radius in Moon radii, kind)
FEATURES = [
    # --- seas and bays -----------------------------------------------------
    ("Oceanus Procellarum", -57.0, 18.0, 0.34, "mare"),
    ("Mare Imbrium",        -16.0, 33.0, 0.26, "mare"),
    ("Mare Serenitatis",     18.0, 28.0, 0.18, "mare"),
    ("Mare Tranquillitatis", 28.0,  8.0, 0.22, "mare"),
    ("Mare Crisium",         59.0, 17.0, 0.16, "mare"),
    ("Mare Fecunditatis",    52.0, -8.0, 0.17, "mare"),
    ("Mare Nectaris",        35.0,-15.0, 0.12, "mare"),
    ("Mare Humorum",        -39.0,-24.0, 0.13, "mare"),
    ("Mare Nubium",         -17.0,-21.0, 0.17, "mare"),
    ("Mare Vaporum",          3.0, 13.0, 0.10, "mare"),
    ("Mare Frigoris",         1.0, 56.0, 0.30, "mare"),
    ("Mare Smythii",         87.0,  2.0, 0.10, "mare"),
    ("Mare Humboldtianum",   81.0, 57.0, 0.08, "mare"),
    ("Sinus Iridum",        -32.0, 45.0, 0.07, "mare"),
    ("Sinus Aestuum",        -9.0, 11.0, 0.06, "mare"),
    ("Sinus Medii",           0.5,  0.5, 0.04, "mare"),
    ("Palus Putredinis",      0.4, 26.5, 0.05, "mare"),
    # --- craters near the disk centre (the ones this study cares about) ----
    ("Triesnecker",           3.6,  4.2, 0.007, "crater"),
    ("Rhaeticus",             4.9,  0.0, 0.013, "crater"),
    ("Reaumur",               0.7, -2.4, 0.014, "crater"),
    ("Herschel",             -2.1, -5.7, 0.011, "crater"),
    ("Ptolemaeus",           -1.8, -9.3, 0.045, "crater"),
    ("Alphonsus",            -2.8,-13.7, 0.033, "crater"),
    ("Arzachel",             -1.9,-18.3, 0.028, "crater"),
    ("Hipparchus",            4.8, -5.5, 0.043, "crater"),
    ("Albategnius",           4.1,-11.2, 0.038, "crater"),
    ("Agrippa",              10.5,  4.1, 0.013, "crater"),
    ("Godin",                10.2,  1.8, 0.010, "crater"),
    ("Manilius",              9.1, 14.5, 0.011, "crater"),
    ("Julius Caesar",        15.4,  9.0, 0.025, "crater"),
    # --- prominent craters elsewhere on the disk ---------------------------
    ("Copernicus",          -20.1,  9.6, 0.027, "crater"),
    ("Eratosthenes",        -11.3, 14.5, 0.017, "crater"),
    ("Archimedes",           -4.0, 29.7, 0.023, "crater"),
    ("Autolycus",             1.5, 30.7, 0.011, "crater"),
    ("Plato",                -9.4, 51.6, 0.029, "crater"),
    ("Aristarchus",         -47.5, 23.7, 0.012, "crater"),
    ("Kepler",              -38.0,  8.1, 0.009, "crater"),
    ("Grimaldi",            -68.3, -5.4, 0.049, "crater"),
    ("Gassendi",            -40.1,-17.6, 0.032, "crater"),
    ("Bullialdus",          -22.2,-20.7, 0.017, "crater"),
    ("Tycho",               -11.4,-43.3, 0.024, "crater"),
    ("Clavius",             -14.4,-58.6, 0.065, "crater"),
    ("Schickard",           -55.3,-44.3, 0.064, "crater"),
    ("Pythagoras",          -63.0, 63.5, 0.020, "crater"),
    ("Theophilus",           26.4,-11.4, 0.028, "crater"),
    ("Langrenus",            61.1, -8.9, 0.038, "crater"),
    ("Petavius",             60.4,-25.3, 0.050, "crater"),
    ("Posidonius",           29.9, 31.8, 0.027, "crater"),
    ("Aristoteles",          17.4, 50.2, 0.025, "crater"),
    ("Eudoxus",              16.3, 44.3, 0.019, "crater"),
    ("Plinius",              23.7, 15.4, 0.012, "crater"),
    ("Menelaus",             16.0, 16.3, 0.008, "crater"),
    ("Proclus",              46.8, 16.1, 0.008, "crater"),
]

# Kilometres per degree of selenographic arc at the sub-observer point.
KM_PER_DEG = 30.31
