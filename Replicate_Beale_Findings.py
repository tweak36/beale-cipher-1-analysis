"""
Replicate the digit-segment search from "Deciphering Beale Cipher 1"
(Duckworth 2025).

What this script does:
  - Concatenates the canonical 520 numbers of Beale Cipher 1 into one
    digit string (1,336 characters total).
  - Segments into sequential 2-digit pairs from offset 0 (668 pairs).
  - Reports whether the target sequences "37 12 21" (latitude) and
    "79 23 16" (longitude) appear at any position, and which source
    cipher numbers they were drawn from.
  - Reports the pair-index separation D between the two triplets and
    the divisor structure of D relevant to Cardan-grille column-wrap
    analysis.

What this script does NOT do:
  - Establish statistical significance of finding these sequences
    (see Monte_Carlo_Null_Test.py and Geographic_Null_Test.py for that)
  - Defeat the p-hacking / multiple-comparison concerns the paper
    itself acknowledges
  - Make any claim about the historical authenticity of the Beale Papers

Source: Canonical Beale Cipher 1 is the 520-number sequence used in
this script (1,336 total digit characters). Reproduced from the Beale
Papers (Ward 1885) and matches independent transcriptions.

Run: python3 Replicate_Beale_Findings.py
"""


cipher_numbers = [
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
    60, 760
]


def divisors_of(n):
    """Return sorted list of all positive divisors of n."""
    return [k for k in range(1, n + 1) if n % k == 0]


def run_replication():
    # Step 1: Concatenate numbers into a single digit string and map character indices
    full_string = ""
    char_to_num_map = []  # Maps each character index back to the index of the original source number

    for num_idx, num in enumerate(cipher_numbers):
        num_str = str(num)
        for char in num_str:
            full_string += char
            char_to_num_map.append(num_idx)

    print(f"Total cipher numbers:                {len(cipher_numbers)}")
    print(f"Total concatenated digit length:     {len(full_string)} characters")

    # Step 2: Segment into two-digit pairs
    pairs = []
    pair_origins = []  # Stores the original numbers that contributed to each pair

    for i in range(0, len(full_string) - 1, 2):
        pair = full_string[i:i+2]
        pairs.append(pair)
        pair_origins.append(char_to_num_map[i])

    print(f"Total 2-digit pairs at offset 0:     {len(pairs)} pairs")

    # Step 3: Scan for all occurrences of the latitude triplet [37, 12, 21]
    lat_positions = []
    for idx in range(len(pairs) - 2):
        if pairs[idx] == "37" and pairs[idx+1] == "12" and pairs[idx+2] == "21":
            lat_positions.append(idx)

    # Step 4: Scan for all occurrences of the longitude triplet [79, 23, 16]
    lon_positions = []
    for idx in range(len(pairs) - 2):
        if pairs[idx] == "79" and pairs[idx+1] == "23" and pairs[idx+2] == "16":
            lon_positions.append(idx)

    # Step 5: Report findings
    print()
    print("=" * 68)
    print("           BEALE CIPHER 1 — TARGET TRIPLET SEARCH")
    print("=" * 68)

    def contributing_source_numbers(pair_pos):
        """Unique cipher numbers (and their list indices) whose digits contribute
        to the triplet of three pairs starting at pair_pos."""
        char_start = pair_pos * 2
        char_end = char_start + 6  # three 2-digit pairs = six digits
        unique_indices = sorted(set(char_to_num_map[char_start:char_end]))
        return unique_indices, [cipher_numbers[i] for i in unique_indices]

    if lat_positions:
        for pos in lat_positions:
            indices, source_nums = contributing_source_numbers(pos)
            print(f"  Latitude  [37,12,21] FOUND at pair index {pos}")
            print(f"    drawn from source numbers {source_nums} "
                  f"at list indices {indices}")
    else:
        print("  Latitude  [37,12,21] NOT FOUND")

    if lon_positions:
        for pos in lon_positions:
            indices, source_nums = contributing_source_numbers(pos)
            print(f"  Longitude [79,23,16] FOUND at pair index {pos}")
            print(f"    drawn from source numbers {source_nums} "
                  f"at list indices {indices}")
            # Plain observation about the longitude origin: per-token vs cross-token
            if source_nums == [79, 23, 16]:
                print(f"    NOTE: This triplet appears as three CONSECUTIVE cipher numbers")
                print(f"          at list indices {indices} (not an emergent")
                print(f"          cross-token artifact of 2-digit pairing).")
    else:
        print("  Longitude [79,23,16] NOT FOUND")

    if not lat_positions or not lon_positions:
        return

    # Step 6: Separation and divisibility analysis
    lat_pos = lat_positions[0]
    lon_pos = lon_positions[0]
    D = abs(lon_pos - lat_pos)
    print()
    print("=" * 68)
    print("                  SEPARATION ANALYSIS")
    print("=" * 68)
    print(f"  Latitude  pair index:   {lat_pos}")
    print(f"  Longitude pair index:   {lon_pos}")
    print(f"  Separation D = |lon_pos - lat_pos| = {D}")
    print()
    divs = divisors_of(D)
    print(f"  All divisors of D = {D}:")
    print(f"    {divs}")
    print(f"  Number of divisors: {len(divs)}")
    print()
    print("  Interpretation: if the cipher is laid out as a grid of width W and")
    print("  W divides D, both triplets fall in the same column of the grid.")
    print("  The published paper highlights this property for W in {22, 32, 44},")
    print("  which are all divisors of 352. Note that these specific widths were")
    print("  selected after observing D = 352, so this is a structural")
    print("  observation, not an independent statistical test.")
    print()

    # Step 7: Raw triplet display
    print("=" * 68)
    print("                    TRIPLET DISPLAY")
    print("=" * 68)
    print(f"  Latitude  triplet (pair index {lat_pos}):  "
          f"[ {pairs[lat_pos]}  {pairs[lat_pos+1]}  {pairs[lat_pos+2]} ]   "
          f"-> 37° 12' 21\" N")
    print(f"  Longitude triplet (pair index {lon_pos}):  "
          f"[ {pairs[lon_pos]}  {pairs[lon_pos+1]}  {pairs[lon_pos+2]} ]   "
          f"-> 79° 23' 16\" W")
    print()


if __name__ == "__main__":
    run_replication()
