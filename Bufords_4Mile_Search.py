"""
Buford's Tavern 4-mile radius search.

Tests the literal reading of Cipher 2's "deposited in Bedford County
about four miles from Buford's." Assumes Buford's = Buford's Tavern at
Locust Level (37°23'11" N, 79°44'14" W) — the most common interpretation
in Beale literature.

Procedure:
  1. Enumerate every DMS coordinate at 1-arcsecond resolution inside a
     bounding box around Buford's that fully contains a 4-mile circle.
  2. Compute great-circle distance to Buford's via haversine; filter to
     those actually within 4 miles.
  3. For each remaining coordinate, check whether BOTH its latitude and
     longitude digit triplets appear in Beale Cipher 1's 2-digit-pair
     sequence.
  4. Report: total candidates, candidates within 4 miles, candidates
     within 4 miles with both triplets in the cipher.

This test specifically corrects for the geographic clip in the previous
Geographic_Null_Test.py annulus: the 79°45' W western limit there
excludes part of the 4-mile circle around Buford's Tavern.

Run:
    python3 Bufords_4Mile_Search.py
"""

import math
from collections import defaultdict


CIPHER_NUMBERS = [
    71, 194, 38, 1701, 89, 76, 11, 83, 1629, 48, 94, 63, 132, 16, 111, 95, 84, 341, 975, 14, 40, 64,
    27, 81, 139, 213, 63, 90, 1120, 8, 15, 3, 126, 2018, 40, 74, 758, 485, 604, 230, 436, 664, 582,
    150, 251, 284, 308, 231, 124, 211, 486, 225, 401, 370, 11, 101, 305, 139, 189, 17, 33, 88, 208,
    193, 145, 1, 94, 73, 416, 918, 263, 28, 500, 538, 356, 117, 136, 219, 27, 176, 130, 10, 460, 25,
    485, 18, 436, 65, 84, 200, 283, 118, 320, 138, 36, 416, 280, 15, 71, 224, 961, 44, 16, 401, 39,
    88, 61, 304, 12, 21, 24, 283, 134, 92, 63, 246, 486, 682, 7, 219, 184, 360, 780, 18, 64, 463,
    474, 131, 160, 79, 73, 440, 95, 18, 64, 581, 34, 69, 128, 367, 460, 17, 81, 12, 103, 820, 62,
    116, 97, 103, 862, 70, 60, 1317, 471, 540, 208, 121, 890, 346, 36, 150, 59, 568, 614, 13, 120,
    63, 219, 812, 2160, 1780, 99, 35, 18, 21, 136, 872, 15, 28, 170, 88, 4, 30, 44, 112, 18, 147,
    436, 195, 320, 37, 122, 113, 6, 140, 8, 120, 305, 42, 58, 461, 44, 106, 301, 13, 408, 680, 93,
    86, 116, 530, 82, 568, 9, 102, 38, 416, 89, 71, 216, 728, 965, 818, 2, 38, 121, 195, 14, 326,
    148, 234, 18, 55, 131, 234, 361, 824, 5, 81, 623, 48, 961, 19, 26, 33, 10, 1101, 365, 92, 88,
    181, 275, 346, 201, 206, 86, 36, 219, 324, 829, 840, 64, 326, 19, 48, 122, 85, 216, 284, 919,
    861, 326, 985, 233, 64, 68, 232, 431, 960, 50, 29, 81, 216, 321, 603, 14, 612, 81, 360, 36, 51,
    62, 194, 78, 60, 200, 314, 676, 112, 4, 28, 18, 61, 136, 247, 819, 921, 1060, 464, 895, 10, 6,
    66, 119, 38, 41, 49, 602, 423, 962, 302, 294, 875, 78, 14, 23, 111, 109, 62, 31, 501, 823, 216,
    280, 34, 24, 150, 1000, 162, 286, 19, 21, 17, 340, 19, 242, 31, 86, 234, 140, 607, 115, 33, 191,
    67, 104, 86, 52, 88, 16, 80, 121, 67, 95, 122, 216, 548, 96, 11, 201, 77, 364, 218, 65, 667,
    890, 236, 154, 211, 10, 98, 34, 119, 56, 216, 119, 71, 218, 1164, 1496, 1817, 51, 39, 210, 36,
    3, 19, 540, 232, 22, 141, 617, 84, 290, 80, 46, 207, 411, 150, 29, 38, 46, 172, 85, 194, 39,
    261, 543, 897, 624, 18, 212, 416, 127, 931, 19, 4, 63, 96, 12, 101, 418, 16, 140, 230, 460,
    538, 19, 27, 88, 612, 1431, 90, 716, 275, 74, 83, 11, 426, 89, 72, 84, 1300, 1706, 814, 221,
    132, 40, 102, 34, 868, 975, 1101, 84, 16, 79, 23, 16, 81, 122, 324, 403, 912, 227, 936, 447,
    55, 86, 34, 43, 212, 107, 96, 314, 264, 1065, 323, 428, 601, 203, 124, 95, 216, 814, 2906, 654,
    820, 2, 301, 112, 176, 213, 71, 87, 96, 202, 35, 10, 2, 41, 17, 84, 221, 736, 820, 214, 11,
    60, 760,
]

# Buford's Tavern (Locust Level) coordinates
BUFORDS_LAT_DMS = (37, 23, 11)   # 37°23'11" N
BUFORDS_LON_DMS = (79, 44, 14)   # 79°44'14" W
BUFORDS_LAT_DEG = BUFORDS_LAT_DMS[0] + BUFORDS_LAT_DMS[1] / 60 + BUFORDS_LAT_DMS[2] / 3600
BUFORDS_LON_DEG = BUFORDS_LON_DMS[0] + BUFORDS_LON_DMS[1] / 60 + BUFORDS_LON_DMS[2] / 3600

# 4-mile bounding box (with extra arcmin margin so we don't clip the circle)
RADIUS_MILES = 4.0
LAT_MARGIN_ARCMIN = 6     # ~5.7 mi at 1 arcmin/0.92 mi (safe margin)
LON_MARGIN_ARCMIN = 8     # ~7.3 mi at this latitude


# ---------------------------------------------------------------------------
# Cipher pair building + triplet lookup
# ---------------------------------------------------------------------------

def build_pairs():
    s = "".join(str(n) for n in CIPHER_NUMBERS)
    return [s[i:i + 2] for i in range(0, len(s) - 1, 2)]


def build_triplet_positions(pairs):
    positions = defaultdict(list)
    for i in range(len(pairs) - 2):
        t = (pairs[i], pairs[i + 1], pairs[i + 2])
        positions[t].append(i)
    return positions


# ---------------------------------------------------------------------------
# Coordinate utilities
# ---------------------------------------------------------------------------

def haversine_miles(lat1_deg, lon1_deg, lat2_deg, lon2_deg):
    """Great-circle distance in statute miles. lon as positive west."""
    R_miles = 3958.7613
    lat1, lat2 = math.radians(lat1_deg), math.radians(lat2_deg)
    # Treat both as west longitude so sign is consistent
    dlat = lat2 - lat1
    dlon = math.radians(lon2_deg - lon1_deg)
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return 2 * R_miles * math.asin(math.sqrt(a))


def dms_to_deg(deg, minutes, seconds):
    return deg + minutes / 60 + seconds / 3600


def deg_to_dms(deg_total):
    deg = int(deg_total)
    rem = (deg_total - deg) * 60
    minutes = int(rem)
    seconds = round((rem - minutes) * 60)
    if seconds == 60:
        seconds = 0
        minutes += 1
    if minutes == 60:
        minutes = 0
        deg += 1
    return deg, minutes, seconds


def dms_to_triplet(deg, minutes, seconds):
    """Convert DMS to the 2-character pair triplet used in cipher matching."""
    return f"{deg:02d}", f"{minutes:02d}", f"{seconds:02d}"


# ---------------------------------------------------------------------------
# Main enumeration
# ---------------------------------------------------------------------------

def main():
    pairs = build_pairs()
    triplet_positions = build_triplet_positions(pairs)
    triplet_set = set(triplet_positions.keys())

    print("=" * 80)
    print("BUFORD'S TAVERN 4-MILE RADIUS SEARCH")
    print("=" * 80)
    print(f"  Buford's Tavern reference: 37°23'11\" N, 79°44'14\" W")
    print(f"                              ({BUFORDS_LAT_DEG:.6f}, {BUFORDS_LON_DEG:.6f})")
    print(f"  Search radius: {RADIUS_MILES:.1f} miles")
    print()

    # Build the bounding-box search region in arcseconds, with margin
    lat_center_arcsec = BUFORDS_LAT_DMS[0] * 3600 + BUFORDS_LAT_DMS[1] * 60 + BUFORDS_LAT_DMS[2]
    lon_center_arcsec = BUFORDS_LON_DMS[0] * 3600 + BUFORDS_LON_DMS[1] * 60 + BUFORDS_LON_DMS[2]
    lat_arcsec_min = lat_center_arcsec - LAT_MARGIN_ARCMIN * 60
    lat_arcsec_max = lat_center_arcsec + LAT_MARGIN_ARCMIN * 60
    lon_arcsec_min = lon_center_arcsec - LON_MARGIN_ARCMIN * 60
    lon_arcsec_max = lon_center_arcsec + LON_MARGIN_ARCMIN * 60

    n_lat = lat_arcsec_max - lat_arcsec_min + 1
    n_lon = lon_arcsec_max - lon_arcsec_min + 1
    n_total_box = n_lat * n_lon

    # Convert arcsec range to DMS for display
    bbox_lat_min = deg_to_dms(lat_arcsec_min / 3600.0)
    bbox_lat_max = deg_to_dms(lat_arcsec_max / 3600.0)
    bbox_lon_min = deg_to_dms(lon_arcsec_min / 3600.0)
    bbox_lon_max = deg_to_dms(lon_arcsec_max / 3600.0)
    print(f"  Bounding box (full circle coverage):")
    print(f"    Latitude  range: {bbox_lat_min[0]:02d}°{bbox_lat_min[1]:02d}'{bbox_lat_min[2]:02d}\" N "
          f"to {bbox_lat_max[0]:02d}°{bbox_lat_max[1]:02d}'{bbox_lat_max[2]:02d}\" N "
          f"({n_lat} values)")
    print(f"    Longitude range: {bbox_lon_min[0]:02d}°{bbox_lon_min[1]:02d}'{bbox_lon_min[2]:02d}\" W "
          f"to {bbox_lon_max[0]:02d}°{bbox_lon_max[1]:02d}'{bbox_lon_max[2]:02d}\" W "
          f"({n_lon} values)")
    print(f"    Total bounding-box candidates: {n_total_box:,}")
    print()

    # First pass: identify which lat/lon triplets in the bounding box exist in the cipher
    lat_in_cipher = []
    for arcsec in range(lat_arcsec_min, lat_arcsec_max + 1):
        deg, m, s = deg_to_dms(arcsec / 3600.0)
        t = dms_to_triplet(deg, m, s)
        if t in triplet_set:
            lat_in_cipher.append((deg, m, s, t))

    lon_in_cipher = []
    for arcsec in range(lon_arcsec_min, lon_arcsec_max + 1):
        deg, m, s = deg_to_dms(arcsec / 3600.0)
        t = dms_to_triplet(deg, m, s)
        if t in triplet_set:
            lon_in_cipher.append((deg, m, s, t))

    print(f"  Latitudes in bounding box with cipher-resident triplets:  "
          f"{len(lat_in_cipher)} of {n_lat}")
    print(f"  Longitudes in bounding box with cipher-resident triplets: "
          f"{len(lon_in_cipher)} of {n_lon}")
    print()

    # Second pass: compute joint candidates and filter by 4-mile haversine
    inside_circle_total = 0
    inside_circle_with_both = []

    # To compute total inside the circle (regardless of cipher match), we'd
    # have to iterate every box cell — for a ~200k cell box this is fast
    # enough but we can sample. We do compute it for context.
    for lat_arcsec in range(lat_arcsec_min, lat_arcsec_max + 1):
        lat_deg = lat_arcsec / 3600.0
        for lon_arcsec in range(lon_arcsec_min, lon_arcsec_max + 1):
            lon_deg = lon_arcsec / 3600.0
            d = haversine_miles(lat_deg, lon_deg, BUFORDS_LAT_DEG, BUFORDS_LON_DEG)
            if d <= RADIUS_MILES:
                inside_circle_total += 1

    print(f"  Total DMS candidates inside the 4-mile circle: "
          f"{inside_circle_total:,}")
    print()

    # Now compute joint candidates with cipher triplet matches inside the circle
    for lat_d, lat_m, lat_s, lat_t in lat_in_cipher:
        lat_deg = dms_to_deg(lat_d, lat_m, lat_s)
        for lon_d, lon_m, lon_s, lon_t in lon_in_cipher:
            lon_deg = dms_to_deg(lon_d, lon_m, lon_s)
            d = haversine_miles(lat_deg, lon_deg, BUFORDS_LAT_DEG, BUFORDS_LON_DEG)
            if d <= RADIUS_MILES:
                inside_circle_with_both.append({
                    "lat_dms": (lat_d, lat_m, lat_s),
                    "lon_dms": (lon_d, lon_m, lon_s),
                    "lat_triplet": lat_t,
                    "lon_triplet": lon_t,
                    "lat_positions": triplet_positions[lat_t],
                    "lon_positions": triplet_positions[lon_t],
                    "distance_mi": d,
                })

    print("=" * 80)
    print("RESULT")
    print("=" * 80)
    print(f"  Candidates inside the 4-mile circle:                    "
          f"{inside_circle_total:,}")
    print(f"  Candidates inside the 4-mile circle with BOTH triplets: "
          f"{len(inside_circle_with_both)}")
    print()

    if not inside_circle_with_both:
        print("  NO coordinate inside the 4-mile circle around Buford's Tavern has")
        print("  both its latitude and longitude triplets present in Beale Cipher 1.")
        print()
        print("  Implication: under the Locust Level / Buford's Tavern interpretation,")
        print("  no candidate satisfies BOTH the literal Cipher 2 constraint")
        print("  ('about four miles from Buford's') AND the digit-pattern test on")
        print("  Cipher 1. The published target (37°12'21\" N, 79°23'16\" W) is")
        print("  outside this 4-mile radius.")
        print()
    else:
        print("  Candidates within 4 miles of Buford's Tavern that have BOTH")
        print("  digit triplets present in Cipher 1:")
        print()
        for c in inside_circle_with_both:
            ld, lm, ls = c["lat_dms"]
            od, om, os = c["lon_dms"]
            print(f"    {ld:02d}°{lm:02d}'{ls:02d}\" N, {od:02d}°{om:02d}'{os:02d}\" W")
            print(f"      Lat triplet {c['lat_triplet']} at pair index(es) {c['lat_positions']}")
            print(f"      Lon triplet {c['lon_triplet']} at pair index(es) {c['lon_positions']}")
            print(f"      Great-circle distance from Buford's: {c['distance_mi']:.2f} mi")
            print()


if __name__ == "__main__":
    main()
