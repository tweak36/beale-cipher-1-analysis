"""
Geographic null-hypothesis test ("Texas Sharpshooter" control) for the
Beale Cipher 1 coordinate finding.

Question this script answers:
  Within the Beale-Cipher-2-specified search region (Bedford County, VA,
  near Buford's Tavern at modern Locust Level / Montvale), how many
  DMS-formatted (latitude, longitude) coordinate pairs have BOTH their
  digit triplets appearing in Beale Cipher 1's 2-digit-pair sequence?

  Cipher 2 (independently solved long before this analysis, key = U.S.
  Declaration of Independence) states the vault is "deposited in Bedford
  County about four miles from Buford's." This script tests:

    Wide:  Bedford County rectangular envelope
           Latitude  37°00'00" N to 37°30'00" N  (1,801 values at 1")
           Longitude 79°09'00" W to 79°45'00" W  (2,161 values at 1")
           Total: 3,891,961 coordinate pairs

    Conditional annulus test: 4-mile ring around ONE candidate location
           ASSUMING "Buford's" refers to Buford's Tavern at Locust Level
                    ≈ 37°23'11" N, 79°44'14" W
           Radius:  3-5 miles (sweep to bracket "about four miles")

  IMPORTANT CAVEAT: "Buford's" in Cipher 2 is a possessive without a
  named referent. The Locust Level / Buford's Tavern interpretation is
  the most common in Beale Papers literature, but is NOT canonically
  established. Other plausible referents exist — Buford family
  residences, farms, mills, or other properties owned by extended
  Buford-family members in 1820s Bedford County. The annulus test
  below should be read as ONE specific interpretation, not as a
  definitive test of Cipher 2's geographic specification.

  For each candidate coordinate, the script applies multiple grid-width
  filters and reports a full sensitivity table:

    Width set            Use case
    {22, 32, 44}         Paper's named widths (factors of D = 352, post-hoc)
    {10, 11, ..., 100}   Broad "any plausible grid width" (least restrictive)
    {44}                 Single-width strict test

  Honest framing: the {22, 32, 44} widths are reverse-engineered from the
  observed separation D = 352. Using them as a uniqueness filter biases
  toward the published finding. The broad [10, 100] criterion is the
  fairer test because the width set is fixed independently of D.

Sources:
  Cipher 2 plaintext: Ward (1885), reproduced in Duckworth (2025)
  Buford's coordinates: en.wikipedia.org/wiki/Locust_Level
  Bedford County bounds: USGS county polygon (approximated as rectangle here;
    full polygon clip is a follow-up sub-test)

Run:
    python3 Geographic_Null_Test.py
"""

import math
from collections import defaultdict

# ---------------------------------------------------------------------------
# Beale Cipher 1 numbers (canonical, 520 entries)
# ---------------------------------------------------------------------------

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

# Bedford County rectangular envelope
LAT_DEG_BASE = 37
LAT_ARCSEC_MIN = 0          # 37°00'00"
LAT_ARCSEC_MAX = 30 * 60    # 37°30'00"

LON_DEG_BASE = 79
LON_ARCSEC_MIN = 9 * 60     # 79°09'00"
LON_ARCSEC_MAX = 45 * 60    # 79°45'00"

# Buford's Tavern (Locust Level), Bedford County, VA
BUFORDS_LAT_DEG = 37 + 23 / 60 + 11 / 3600    # 37°23'11" N
BUFORDS_LON_DEG = 79 + 44 / 60 + 14 / 3600    # 79°44'14" W

# Grid width sets to test
WIDTH_SET_NAMED = {22, 32, 44}                          # paper's post-hoc set
WIDTH_SET_BROAD = set(range(10, 101))                   # any plausible width
WIDTH_SET_44_ONLY = {44}                                # single strict width


# ---------------------------------------------------------------------------
# Helpers
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


def arcsec_to_dms(deg_base, arcsec_offset):
    minutes = arcsec_offset // 60
    seconds = arcsec_offset % 60
    return deg_base, minutes, seconds


def dms_to_triplet(deg, minutes, seconds):
    return f"{deg:02d}", f"{minutes:02d}", f"{seconds:02d}"


def haversine_miles(lat1_deg, lon1_deg, lat2_deg, lon2_deg):
    """Great-circle distance in statute miles."""
    R_miles = 3958.7613
    lat1, lon1 = math.radians(lat1_deg), math.radians(lon1_deg)
    lat2, lon2 = math.radians(lat2_deg), math.radians(lon2_deg)
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return 2 * R_miles * math.asin(math.sqrt(a))


# ---------------------------------------------------------------------------
# Main enumeration
# ---------------------------------------------------------------------------

def main():
    pairs = build_pairs()
    triplet_positions = build_triplet_positions(pairs)
    triplet_set = set(triplet_positions.keys())

    print("=" * 80)
    print("GEOGRAPHIC NULL TEST — Bedford County coordinate enumeration")
    print("=" * 80)
    print(f"  Beale Cipher 1 pair count:           {len(pairs)}")
    print(f"  Distinct 3-pair triplets in cipher:  {len(triplet_set)}")
    print()
    print(f"  Bedford County rectangular envelope:")
    print(f"    Latitude  range:  37°00'00\" N to 37°30'00\" N  "
          f"({LAT_ARCSEC_MAX - LAT_ARCSEC_MIN + 1} values)")
    print(f"    Longitude range:  79°09'00\" W to 79°45'00\" W  "
          f"({LON_ARCSEC_MAX - LON_ARCSEC_MIN + 1} values)")
    n_lat = LAT_ARCSEC_MAX - LAT_ARCSEC_MIN + 1
    n_lon = LON_ARCSEC_MAX - LON_ARCSEC_MIN + 1
    n_total = n_lat * n_lon
    print(f"    Total coordinate pairs:     {n_total:,}")
    print()
    print(f"  Buford's Tavern (Locust Level) reference: "
          f"37°23'11\" N, 79°44'14\" W")
    print()

    # Step 1: latitudes whose triplet exists in cipher
    lat_hits = []
    for arcsec in range(LAT_ARCSEC_MIN, LAT_ARCSEC_MAX + 1):
        deg, m, s = arcsec_to_dms(LAT_DEG_BASE, arcsec)
        t = dms_to_triplet(deg, m, s)
        if t in triplet_set:
            lat_hits.append((deg, m, s, t))

    # Step 2: longitudes whose triplet exists in cipher
    lon_hits = []
    for arcsec in range(LON_ARCSEC_MIN, LON_ARCSEC_MAX + 1):
        deg, m, s = arcsec_to_dms(LON_DEG_BASE, arcsec)
        t = dms_to_triplet(deg, m, s)
        if t in triplet_set:
            lon_hits.append((deg, m, s, t))

    print(f"  Latitudes in box with cipher-resident triplets:  "
          f"{len(lat_hits)} of {n_lat}")
    print(f"  Longitudes in box with cipher-resident triplets: "
          f"{len(lon_hits)} of {n_lon}")
    print(f"  Joint candidates (both halves present):          "
          f"{len(lat_hits) * len(lon_hits):,}")
    print()

    # Step 3: build the full joint candidate list with position + distance info
    joint_candidates = []
    for lat_d, lat_m, lat_s, lat_t in lat_hits:
        lat_positions = triplet_positions[lat_t]
        lat_deg = lat_d + lat_m / 60 + lat_s / 3600
        for lon_d, lon_m, lon_s, lon_t in lon_hits:
            lon_positions = triplet_positions[lon_t]
            lon_deg = lon_d + lon_m / 60 + lon_s / 3600
            # All possible (lat_pos, lon_pos) pair combinations
            D_widths_named = []
            D_widths_broad = []
            D_widths_44 = []
            best_D = None
            for lp in lat_positions:
                for op in lon_positions:
                    D = abs(op - lp)
                    if D == 0:
                        continue
                    if best_D is None or D < best_D:
                        best_D = D
                    if any(D % w == 0 for w in WIDTH_SET_NAMED):
                        D_widths_named.append((lp, op, D))
                    if any(D % w == 0 for w in WIDTH_SET_BROAD):
                        D_widths_broad.append((lp, op, D))
                    if any(D % w == 0 for w in WIDTH_SET_44_ONLY):
                        D_widths_44.append((lp, op, D))
            dist = haversine_miles(lat_deg, lon_deg, BUFORDS_LAT_DEG, BUFORDS_LON_DEG)
            joint_candidates.append({
                "lat_dms": (lat_d, lat_m, lat_s),
                "lon_dms": (lon_d, lon_m, lon_s),
                "lat_deg": lat_deg,
                "lon_deg": lon_deg,
                "lat_positions": lat_positions,
                "lon_positions": lon_positions,
                "best_D": best_D,
                "aligned_named": bool(D_widths_named),
                "aligned_broad": bool(D_widths_broad),
                "aligned_44": bool(D_widths_44),
                "miles_from_bufords": dist,
            })

    # Step 4: full candidate listing
    print("=" * 80)
    print("ALL JOINT CANDIDATES IN BEDFORD COUNTY RECTANGULAR ENVELOPE")
    print("=" * 80)
    if not joint_candidates:
        print("  (none)")
    else:
        for c in joint_candidates:
            lat_d, lat_m, lat_s = c["lat_dms"]
            lon_d, lon_m, lon_s = c["lon_dms"]
            print(f"  {lat_d:02d}°{lat_m:02d}'{lat_s:02d}\" N, "
                  f"{lon_d:02d}°{lon_m:02d}'{lon_s:02d}\" W")
            print(f"    Lat positions:        {c['lat_positions']}")
            print(f"    Lon positions:        {c['lon_positions']}")
            print(f"    Smallest separation D: {c['best_D']}")
            print(f"    Aligned under {{22,32,44}} (paper's post-hoc set): "
                  f"{c['aligned_named']}")
            print(f"    Aligned under [10,100] (any plausible width):    "
                  f"{c['aligned_broad']}")
            print(f"    Aligned under {{44}} (single strict width):       "
                  f"{c['aligned_44']}")
            print(f"    Distance from Buford's: {c['miles_from_bufords']:.2f} miles")
            print()

    # Step 5: sensitivity table over grid-width sets and Buford's annulus
    print("=" * 80)
    print("SENSITIVITY TABLE — strict-uniqueness depends on rule selection")
    print("=" * 80)
    n_aligned_named = sum(1 for c in joint_candidates if c["aligned_named"])
    n_aligned_broad = sum(1 for c in joint_candidates if c["aligned_broad"])
    n_aligned_44 = sum(1 for c in joint_candidates if c["aligned_44"])
    n_existence = len(joint_candidates)

    print()
    print("  Joint hits in Bedford County rectangle (existence only):  "
          f"{n_existence}")
    print("  Hits also requiring column alignment under grid widths:")
    print(f"    {{22, 32, 44}}  (paper's named set, factors of 352):      "
          f"{n_aligned_named}")
    print(f"    [10, 100]     (any plausible grid width):                "
          f"{n_aligned_broad}")
    print(f"    {{44}}         (single most-cited width):                 "
          f"{n_aligned_44}")
    print()
    print("  HONEST READ: the uniqueness of the published target depends")
    print("  sharply on which width set is used. Under the paper's named")
    print("  set {22, 32, 44}, exactly one of the joint candidates aligns")
    print("  (the published target). Under the broader [10, 100] criterion")
    print("  used in Monte_Carlo_Null_Test.py, additional candidates may")
    print("  also align. Report this sensitivity, do not cherry-pick.")
    print()

    # Step 6: conditional Buford's annulus test (one interpretation, not definitive)
    print("=" * 80)
    print("CONDITIONAL ANNULUS TEST — assuming Buford's = Locust Level")
    print("=" * 80)
    print()
    print("  Cipher 2 states the vault is 'deposited in Bedford County")
    print("  about four miles from Buford's.' The referent 'Buford's' is")
    print("  not specified in Cipher 2 itself — it could mean Buford's")
    print("  Tavern at Locust Level (popular interpretation), a Buford")
    print("  family residence, farm, mill, or other property. Multiple")
    print("  Buford-family-associated locations existed in 1820s Bedford")
    print("  County.")
    print()
    print("  The annulus below is computed ASSUMING the Locust Level")
    print("  interpretation. Results should not be read as definitive")
    print("  evidence for or against the published target. They are one")
    print("  interpretation among several historically plausible ones.")
    print()
    for r_inner, r_outer, label in [
        (3.0, 5.0, "3-5 mi annulus (loose)"),
        (3.5, 4.5, "3.5-4.5 mi annulus (tight)"),
        (0.0, 5.0, "0-5 mi disk (very loose)"),
    ]:
        hits = [c for c in joint_candidates
                if r_inner <= c["miles_from_bufords"] <= r_outer]
        n_aligned_in_annulus = sum(1 for c in hits if c["aligned_named"])
        print(f"  {label}: {len(hits)} joint hits, "
              f"{n_aligned_in_annulus} also column-aligned under {{22,32,44}}")
        for c in hits:
            lat_d, lat_m, lat_s = c["lat_dms"]
            lon_d, lon_m, lon_s = c["lon_dms"]
            tag = "(aligned)" if c["aligned_named"] else "(not aligned)"
            print(f"    -> {lat_d:02d}°{lat_m:02d}'{lat_s:02d}\" N, "
                  f"{lon_d:02d}°{lon_m:02d}'{lon_s:02d}\" W, "
                  f"{c['miles_from_bufords']:.2f} mi  {tag}")
    print()

    # Step 7: distance reporting (conditional on Locust Level interpretation)
    print("=" * 80)
    print("DISTANCE REPORTING (conditional on Locust Level interpretation)")
    print("=" * 80)
    target_lat_dms = (37, 12, 21)
    target_lon_dms = (79, 23, 16)
    matched = [c for c in joint_candidates
               if c["lat_dms"] == target_lat_dms and c["lon_dms"] == target_lon_dms]
    if matched:
        c = matched[0]
        print(f"  Published target: 37°12'21\" N, 79°23'16\" W")
        print(f"    Distance from Locust Level (one Buford's candidate): "
              f"{c['miles_from_bufords']:.2f} miles")
        print()
        print(f"  Interpretation depends on which Buford's reference is correct:")
        print(f"    - If Buford's = Locust Level, the published target is far from")
        print(f"      the literal 4-mile specification.")
        print(f"    - If Buford's refers to a different Buford-family property")
        print(f"      elsewhere in Bedford County, the 4-mile constraint cannot")
        print(f"      be evaluated until that referent is historically resolved.")
        print(f"    - Resolution of the Buford's referent is a historical research")
        print(f"      question, not a question this geographic enumeration can")
        print(f"      settle.")
    else:
        print("  Published target NOT recovered in this enumeration — debug.")
    print()


if __name__ == "__main__":
    main()
