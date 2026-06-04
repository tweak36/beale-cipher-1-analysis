"""
Segmentation sensitivity test for the Beale Cipher 1 coordinate finding.

The published analysis pairs the cipher's digit string at offset 0 with
2-character segments. This script tests how the finding behaves under
alternative segmentations that a hostile reviewer would (correctly)
flag as researcher degrees of freedom:

  - Segment length:  1, 2, 3, 4 digits per group
  - Starting offset: 0 and 1 (phase shift to confirm offset choice
                              isn't load-bearing)

For each (length, offset) configuration the script checks whether the
published targets are still recoverable as 6-digit DMS-formatted
triplets matching:

    Latitude  37° 12' 21" N  -> digit string "371221"
    Longitude 79° 23' 16" W  -> digit string "792316"

The target SUBSTRINGS are searched in the resulting segment stream as
contiguous segment sequences whose concatenation matches the target.
For example, with length=3 we look for three 2-digit-equivalent
sequences that concatenate to "371221" (i.e., segments ["371", "221"]
under length=3 produce the latitude). With length=1 each segment is a
single digit, etc.

Honest scope:
  - This test addresses one specific axis of researcher freedom
    (segment length + offset). It does NOT exhaust all transformations
    a hostile reviewer might propose (sliding window, base conversion,
    nth-letter book-cipher decoding against alternate keys, etc.).
  - A result of "target recoverable under only one segmentation" is
    moderate evidence against segmentation cherry-picking; a result of
    "target recoverable under many segmentations" weakens uniqueness.

Run:
    python3 Segmentation_Sensitivity_Test.py
"""

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

LAT_DIGITS = "371221"
LON_DIGITS = "792316"

SEGMENT_LENGTHS = [1, 2, 3, 4]
OFFSETS = [0, 1]


def build_digit_string(numbers):
    return "".join(str(n) for n in numbers)


def find_substring_positions(text, target):
    """Return list of every starting index where target occurs in text."""
    positions = []
    start = 0
    while True:
        idx = text.find(target, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions


def find_in_segmented(digit_string, target, seg_len, offset):
    """
    Build the segment stream starting at `offset` with segments of length
    `seg_len`. Return list of segment indices where the target string
    matches a contiguous concatenation of segments (i.e., the target must
    start exactly at a segment boundary and span exactly target_len/seg_len
    segments — only possible if seg_len evenly divides target length).
    """
    if len(target) % seg_len != 0:
        return None      # target does not fit cleanly into segments
    target_segments = [target[i:i + seg_len] for i in range(0, len(target), seg_len)]

    segments = []
    seg_origins = []   # starting digit-position of each segment
    pos = offset
    while pos + seg_len <= len(digit_string):
        segments.append(digit_string[pos:pos + seg_len])
        seg_origins.append(pos)
        pos += seg_len

    hits = []
    for i in range(len(segments) - len(target_segments) + 1):
        if all(segments[i + j] == target_segments[j] for j in range(len(target_segments))):
            hits.append((i, seg_origins[i]))
    return hits


def main():
    text = build_digit_string(CIPHER_NUMBERS)
    n_digits = len(text)

    print("=" * 80)
    print("SEGMENTATION SENSITIVITY TEST")
    print("=" * 80)
    print(f"  Total digit string length: {n_digits} characters")
    print(f"  Targets:  latitude  '{LAT_DIGITS}'  (DMS 37°12'21\" N)")
    print(f"            longitude '{LON_DIGITS}'  (DMS 79°23'16\" W)")
    print()
    print(f"  Substrings searched IN PLACE (no segmentation):")
    lat_subs = find_substring_positions(text, LAT_DIGITS)
    lon_subs = find_substring_positions(text, LON_DIGITS)
    print(f"    Latitude  substring '{LAT_DIGITS}' positions in digit string: {lat_subs}")
    print(f"    Longitude substring '{LON_DIGITS}' positions in digit string: {lon_subs}")
    print()
    print("  Note: substring search ignores segmentation entirely. If a target")
    print("  appears at digit position p that is not on a segment boundary for")
    print("  some chosen (length, offset), it would be invisible to that")
    print("  segmentation but still 'in the string.'")
    print()

    # Tabulate by (segment_length, offset)
    print("=" * 80)
    print("TARGET RECOVERABILITY BY (SEGMENT LENGTH, OFFSET)")
    print("=" * 80)
    print()
    header = f"  {'len':>3}  {'offset':>6}  {'lat hits':>10}  {'lon hits':>10}  {'both found':>10}"
    print(header)
    print("  " + "-" * (len(header) - 2))

    summary = []
    for seg_len in SEGMENT_LENGTHS:
        for offset in OFFSETS:
            lat_hits = find_in_segmented(text, LAT_DIGITS, seg_len, offset)
            lon_hits = find_in_segmented(text, LON_DIGITS, seg_len, offset)
            if lat_hits is None or lon_hits is None:
                continue
            both = bool(lat_hits) and bool(lon_hits)
            summary.append((seg_len, offset, lat_hits, lon_hits, both))
            print(f"  {seg_len:>3}  {offset:>6}  {len(lat_hits):>10}  "
                  f"{len(lon_hits):>10}  {str(both):>10}")
    print()

    # Detail any configuration where both targets are found
    print("=" * 80)
    print("DETAIL OF CONFIGURATIONS WHERE BOTH TARGETS RECOVERED")
    print("=" * 80)
    print()
    any_both = False
    for seg_len, offset, lat_hits, lon_hits, both in summary:
        if not both:
            continue
        any_both = True
        print(f"  segment length = {seg_len}, offset = {offset}:")
        print(f"    Latitude  hits (segment_idx, digit_pos): {lat_hits}")
        print(f"    Longitude hits (segment_idx, digit_pos): {lon_hits}")
        # Compute segment-index separation D for the closest lat/lon pair
        Ds = [abs(o[0] - l[0]) for l in lat_hits for o in lon_hits if abs(o[0] - l[0]) > 0]
        if Ds:
            print(f"    Smallest segment-index separation D: {min(Ds)}")
            print(f"    All non-zero separations: {sorted(set(Ds))}")
        print()

    if not any_both:
        print("  (no configuration recovered both targets)")
        print()

    # Verdict
    print("=" * 80)
    print("VERDICT")
    print("=" * 80)
    print()
    n_both = sum(1 for *_, both in summary if both)
    print(f"  Configurations tested: {len(summary)}")
    print(f"  Configurations where BOTH targets recovered: {n_both}")
    print()
    if n_both == 1:
        print("  Only the published (length=2, offset=0) segmentation recovers both")
        print("  targets. This is moderate evidence against segmentation cherry-picking:")
        print("  among the tested variants, only one works.")
    elif n_both > 1:
        print(f"  {n_both} out of {len(summary)} segmentation variants recover both targets.")
        print("  The published (length=2, offset=0) is not unique in this dimension.")
        print("  Researcher-selection concern on segmentation is NOT fully defeated.")
    else:
        print("  No segmentation recovered both targets — sanity check failed.")
    print()


if __name__ == "__main__":
    main()
