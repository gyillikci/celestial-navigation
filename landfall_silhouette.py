''' Sample for landfall via silhouette reading (terrestrial navigation from a photograph)
    © August Linnman, 2025, email: august@linnman.net
    MIT License (see LICENSE file)

    Classic landfall navigation (see terrestrial.py) requires a sextant held horizontally
    to measure the angle between pairs of identified landmarks. This sample shows how the
    same angles can instead be "read" directly out of a photographed skyline/silhouette
    (a distant coastline or mountain range with a handful of identifiable peaks), using
    only the pixel positions of the landmarks, the image width and the camera's horizontal
    field of view -- no compass heading is required, since only the *relative* angles
    between landmarks matter, exactly as with a sextant.

    This is a lightweight, purely geometric take on the idea explored by the CrossLocate
    project (matching a photographed natural skyline against rendered terrain silhouettes
    to geo-localize the photographer): instead of a trained image-retrieval model and a
    large database of rendered silhouettes, a handful of manually (or automatically)
    identified charted landmarks are used to resect the position directly.

    It also demonstrates TerrestrialFixCollection, which generalizes
    get_terrestrial_position() from exactly three landmarks to three-or-more, with
    redundant landmarks improving the accuracy of the fix.
'''
from time import time
from math import tan
from starfix import LatLonGeodetic, get_terrestrial_position, get_google_map_string,\
      spherical_distance, angle_between_points, get_angles_from_silhouette,\
      TerrestrialFixCollection, km_to_nm, deg_to_rad

def main ():
    ''' Main body of script '''

    starttime = time ()

    # The same three lighthouses used in terrestrial.py, ordered left to right
    # as they would appear in a photograph of the coastline.
    p1 = LatLonGeodetic (58.7396,   17.8656)   # Landsort
    p2 = LatLonGeodetic (58.594091, 17.467489) # Gustaf Dalén
    p3 = LatLonGeodetic (58.60355,  17.316041) # Hävringe

    # A fourth, synthetic reference point, added purely to demonstrate that
    # TerrestrialFixCollection supports more than three landmarks (redundant
    # observations). It is NOT a real charted landmark.
    p4 = LatLonGeodetic (58.55, 17.05)

    # Establish a "true" observer position using the same real, previously measured
    # sextant angles as in terrestrial.py (20 degrees between p1/p2, 45 degrees
    # between p2/p3).
    angle_1 = 20
    angle_2 = 45
    candidates, _, _, _, _ = \
        get_terrestrial_position (p3, p2, angle_1, p2, p1, angle_2)
    assert isinstance (candidates, tuple)
    # Eliminate the false intersection (located at one of the lighthouses)
    true_position_gc = None
    for candidate in candidates:
        if all (spherical_distance (candidate, lh) > 0.001 for lh in (p1, p2, p3)):
            true_position_gc = candidate
            break
    assert true_position_gc is not None

    # Simulate a photograph taken from the true position. Only the relative bearing
    # of each landmark (from the true angle between it and its neighbour) is needed
    # to place it on a synthetic image -- exactly what a camera would record.
    image_width_px = 4000
    horizontal_fov_deg = 60

    angle_p1_p2 = angle_between_points (true_position_gc, p1, p2)
    angle_p2_p3 = angle_between_points (true_position_gc, p2, p3)
    angle_p3_p4 = angle_between_points (true_position_gc, p3, p4)

    focal_length_px = (image_width_px / 2) / tan (deg_to_rad (horizontal_fov_deg / 2))

    def pixel_for_offset (offset_deg : float) -> float:
        return image_width_px / 2 + focal_length_px * tan (deg_to_rad (offset_deg))

    # p2 is placed at the center of the image (offset 0); p1 is to its left,
    # p3 and p4 are to its right, matching their real bearing order.
    pixel_p1 = pixel_for_offset (-angle_p1_p2)
    pixel_p2 = pixel_for_offset (0)
    pixel_p3 = pixel_for_offset (angle_p2_p3)
    pixel_p4 = pixel_for_offset (angle_p2_p3 + angle_p3_p4)

    # Simulate small measurement errors when identifying the landmarks in the photo
    # (e.g. sub-pixel localization error), to show how a redundant landmark improves
    # the fix. The errors below are fixed (not random) so the sample is reproducible.
    pixel_p1 += 3
    pixel_p2 -= 2
    pixel_p3 += 4
    pixel_p4 -= 3

    print ("--- Fix using 3 landmarks read from a silhouette photograph ---")
    angles_3 = get_angles_from_silhouette \
        ([pixel_p1, pixel_p2, pixel_p3], image_width_px, horizontal_fov_deg)
    fix_3 = TerrestrialFixCollection ([p1, p2, p3], angles_3)
    position_3, fitness_3 = fix_3.get_position (estimated_position=true_position_gc)
    error_3_m = round (spherical_distance (position_3, true_position_gc) * 1000, 1)
    print ("Your location = " + get_google_map_string (position_3, 4))
    print ("Fitness = " + str (round (fitness_3, 6)))
    print ("Error vs. true position = " + str (error_3_m) + " m")

    print ("\n--- Fix using 4 landmarks (one redundant) read from the same photograph ---")
    angles_4 = get_angles_from_silhouette \
        ([pixel_p1, pixel_p2, pixel_p3, pixel_p4], image_width_px, horizontal_fov_deg)
    fix_4 = TerrestrialFixCollection ([p1, p2, p3, p4], angles_4)
    position_4, fitness_4 = fix_4.get_position (estimated_position=true_position_gc)
    error_4_m = round (spherical_distance (position_4, true_position_gc) * 1000, 1)
    print ("Your location = " + get_google_map_string (position_4, 4))
    print ("Fitness = " + str (round (fitness_4, 6)))
    print ("Error vs. true position = " + str (error_4_m) + " m")

    print ("\nTrue position (from sextant angles in terrestrial.py) = " +\
           get_google_map_string (true_position_gc, 4))
    print ("Distance in nautical miles between the two fixes = " +\
           str (round (km_to_nm (spherical_distance (position_3, position_4)), 4)))

    endtime = time ()
    taken_ms = round ((endtime-starttime)*1000, 2)
    print ("\nTime taken = " +str (taken_ms)+" ms")

if __name__ == '__main__':
    main()
