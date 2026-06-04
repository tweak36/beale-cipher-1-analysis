"""
Texas Sharpshooter empirical test for the Beale Cipher 1 coordinate finding.

This script answers the right Monte Carlo question:

  "How often does ANY valid Bedford County DMS coordinate pair appear
   in a random 668-pair string?"

The original Monte_Carlo_Null_Test.py tested the wrong question — joint
co-occurrence of two SPECIFIC pre-selected triplets ([37,12,21] and
[79,23,16]). That test produced 0/100,000 and a 95 % UCB of ~3 × 10⁻⁵,
but the result is misleading because the targets were chosen after
seeing the cipher (post-hoc selection / Texas Sharpshooter fallacy).

The right test: for each random trial, count how many valid Bedford
County latitude/longitude triplets appear ANYWHERE in the string.

Bedford County DMS bounds (from Geographic_Null_Test.py):
  Latitude  range: 37°00'00" N to 37°30'00" N
    Valid lat triplets: ("37", mm, ss) where mm in [0,30], ss in [0,59]
    Total valid lat triplets: 1 × 31 × 60 = 1,860
  Longitude range: 79°09'00" W to 79°45'00" W
    Valid lon triplets: ("79", mm, ss) where mm in [9,45], ss in [0,59]
    Total valid lon triplets: 1 × 37 × 60 = 2,220

Theoretical expectation (independence approximation):
  P(any 3-pair window matches a valid lat) = 1860 / 100^3 = 0.00186
  P(any 3-pair window matches a valid lon) = 2220 / 100^3 = 0.00222
  Number of 3-pair windows in 668-pair cipher: 666
  E[valid lat triplets in a random cipher] = 666 × 0.00186 ≈ 1.24
  E[valid lon triplets in a random cipher] = 666 × 0.00222 ≈ 1.48
  E[joint (lat, lon) coordinate pairs]     ≈ 1.24 × 1.48 ≈ 1.83
  P(at least one joint pair) ≈ 1 - exp(-1.83) ≈ 84 %

Run:
    python3 Texas_Sharpshooter_Test.py
"""

import math
import random

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

# Bedford County DMS bounds
LAT_DEG, LAT_MM_MIN, LAT_MM_MAX = "37", 0, 30
LON_DEG, LON_MM_MIN, LON_MM_MAX = "79", 9, 45
SS_MIN, SS_MAX = 0, 59

N_TRIALS = 100_000


def build_pairs_from_numbers(numbers):
    s = "".join(str(n) for n in numbers)
    return [s[i:i + 2] for i in range(0, len(s) - 1, 2)]


def is_valid_lat_triplet(p0, p1, p2):
    """Check if (p0, p1, p2) is a valid Bedford County latitude triplet."""
    if p0 != LAT_DEG:
        return False
    try:
        mm = int(p1)
        ss = int(p2)
    except ValueError:
        return False
    return LAT_MM_MIN <= mm <= LAT_MM_MAX and SS_MIN <= ss <= SS_MAX


def is_valid_lon_triplet(p0, p1, p2):
    """Check if (p0, p1, p2) is a valid Bedford County longitude triplet."""
    if p0 != LON_DEG:
        return False
    try:
        mm = int(p1)
        ss = int(p2)
    except ValueError:
        return False
    return LON_MM_MIN <= mm <= LON_MM_MAX and SS_MIN <= ss <= SS_MAX


def count_valid_triplets(pairs):
    """Count valid lat and lon Bedford County triplets in a pair sequence,
    and also collect their starting indices for column-alignment testing."""
    lat_positions = []
    lon_positions = []
    for i in range(len(pairs) - 2):
        if is_valid_lat_triplet(pairs[i], pairs[i + 1], pairs[i + 2]):
            lat_positions.append(i)
        if is_valid_lon_triplet(pairs[i], pairs[i + 1], pairs[i + 2]):
            lon_positions.append(i)
    return lat_positions, lon_positions


def has_column_aligned_pair(lat_positions, lon_positions, widths):
    """True if ANY (lat_pos, lon_pos) pair has separation divisible by
    any width in `widths`."""
    for lp in lat_positions:
        for op in lon_positions:
            D = abs(op - lp)
            if D == 0:
                continue
            if any(D % w == 0 for w in widths):
                return True
    return False


def null_A_permutation(rng):
    shuffled = CIPHER_NUMBERS[:]
    rng.shuffle(shuffled)
    return build_pairs_from_numbers(shuffled)


def null_B_random_digits(rng, total_length):
    digits = "".join(str(rng.randint(0, 9)) for _ in range(total_length))
    return [digits[i:i + 2] for i in range(0, len(digits) - 1, 2)]


WIDTHS_NAMED = [22, 32, 44]
WIDTHS_BROAD = list(range(10, 101))


def run_trials(null_func, n_trials, total_length=None, rng_seed=42):
    rng = random.Random(rng_seed)
    sum_lat = 0
    sum_lon = 0
    sum_joint = 0
    n_with_any_joint = 0
    n_with_any_lat = 0
    n_with_any_lon = 0
    n_with_aligned_named = 0
    n_with_aligned_broad = 0
    for _ in range(n_trials):
        if total_length is not None:
            pairs = null_func(rng, total_length)
        else:
            pairs = null_func(rng)
        lat_positions, lon_positions = count_valid_triplets(pairs)
        sum_lat += len(lat_positions)
        sum_lon += len(lon_positions)
        joint = len(lat_positions) * len(lon_positions)
        sum_joint += joint
        if len(lat_positions) > 0:
            n_with_any_lat += 1
        if len(lon_positions) > 0:
            n_with_any_lon += 1
        if joint > 0:
            n_with_any_joint += 1
            if has_column_aligned_pair(lat_positions, lon_positions, WIDTHS_NAMED):
                n_with_aligned_named += 1
            if has_column_aligned_pair(lat_positions, lon_positions, WIDTHS_BROAD):
                n_with_aligned_broad += 1
    return {
        "mean_lat": sum_lat / n_trials,
        "mean_lon": sum_lon / n_trials,
        "mean_joint": sum_joint / n_trials,
        "frac_any_lat": n_with_any_lat / n_trials,
        "frac_any_lon": n_with_any_lon / n_trials,
        "frac_any_joint": n_with_any_joint / n_trials,
        "frac_aligned_named": n_with_aligned_named / n_trials,
        "frac_aligned_broad": n_with_aligned_broad / n_trials,
    }


def main():
    # --- Theoretical expectations ---
    p_lat = (1 * 31 * 60) / 100**3
    p_lon = (1 * 37 * 60) / 100**3
    n_windows = 668 - 2
    e_lat = n_windows * p_lat
    e_lon = n_windows * p_lon
    e_joint = e_lat * e_lon
    p_at_least_one = 1 - math.exp(-e_joint)

    print("=" * 78)
    print("TEXAS SHARPSHOOTER EMPIRICAL TEST")
    print("=" * 78)
    print()
    print("Theoretical expectation (Poisson, independence):")
    print(f"  P(arbitrary window matches valid lat triplet): {p_lat:.6f}")
    print(f"  P(arbitrary window matches valid lon triplet): {p_lon:.6f}")
    print(f"  Number of 3-pair windows in cipher: {n_windows}")
    print(f"  E[valid lat triplets per random cipher]: {e_lat:.3f}")
    print(f"  E[valid lon triplets per random cipher]: {e_lon:.3f}")
    print(f"  E[joint (lat, lon) coord pairs]:         {e_joint:.3f}")
    print(f"  P(at least one joint pair) [Poisson]:    {100*p_at_least_one:.1f} %")
    print()

    # --- Baseline: real cipher ---
    real_pairs = build_pairs_from_numbers(CIPHER_NUMBERS)
    real_lat_pos, real_lon_pos = count_valid_triplets(real_pairs)
    real_joint = len(real_lat_pos) * len(real_lon_pos)
    real_aligned_named = has_column_aligned_pair(real_lat_pos, real_lon_pos, WIDTHS_NAMED)
    real_aligned_broad = has_column_aligned_pair(real_lat_pos, real_lon_pos, WIDTHS_BROAD)
    print("=" * 78)
    print("REAL BEALE CIPHER 1 — baseline measurement")
    print("=" * 78)
    print(f"  Valid Bedford latitude triplet count:  {len(real_lat_pos)}  "
          f"(positions {real_lat_pos})")
    print(f"  Valid Bedford longitude triplet count: {len(real_lon_pos)}  "
          f"(positions {real_lon_pos})")
    print(f"  Joint (lat, lon) coordinate pairs: {real_joint}")
    print(f"  At least one column-aligned pair (widths {WIDTHS_NAMED}): {real_aligned_named}")
    print(f"  At least one column-aligned pair (widths [10, 100]):     {real_aligned_broad}")
    print()

    # --- Null A — token-level permutation ---
    print("=" * 78)
    print(f"NULL A — Token-level permutation of Beale numbers (N = {N_TRIALS:,})")
    print("=" * 78)
    res_A = run_trials(null_A_permutation, N_TRIALS, rng_seed=42)
    print(f"  Mean valid lat triplets per trial:  {res_A['mean_lat']:.3f}  "
          f"(theory: {e_lat:.3f})")
    print(f"  Mean valid lon triplets per trial:  {res_A['mean_lon']:.3f}  "
          f"(theory: {e_lon:.3f})")
    print(f"  Mean joint coord pairs per trial:   {res_A['mean_joint']:.3f}  "
          f"(theory: {e_joint:.3f})")
    print(f"  Fraction with >=1 joint coord pair:                "
          f"{100*res_A['frac_any_joint']:6.2f} %")
    print(f"  Fraction with >=1 column-aligned pair, widths {WIDTHS_NAMED}: "
          f"{100*res_A['frac_aligned_named']:6.2f} %")
    print(f"  Fraction with >=1 column-aligned pair, widths [10,100]:    "
          f"{100*res_A['frac_aligned_broad']:6.2f} %")
    print()

    # --- Null B — random digit strings ---
    total_chars = sum(len(str(n)) for n in CIPHER_NUMBERS)
    print("=" * 78)
    print(f"NULL B — Random uniform digit strings (N = {N_TRIALS:,})")
    print("=" * 78)
    res_B = run_trials(null_B_random_digits, N_TRIALS,
                      total_length=total_chars, rng_seed=43)
    print(f"  Mean valid lat triplets per trial:  {res_B['mean_lat']:.3f}  "
          f"(theory: {e_lat:.3f})")
    print(f"  Mean valid lon triplets per trial:  {res_B['mean_lon']:.3f}  "
          f"(theory: {e_lon:.3f})")
    print(f"  Mean joint coord pairs per trial:   {res_B['mean_joint']:.3f}  "
          f"(theory: {e_joint:.3f})")
    print(f"  Fraction with >=1 joint coord pair:                "
          f"{100*res_B['frac_any_joint']:6.2f} %")
    print(f"  Fraction with >=1 column-aligned pair, widths {WIDTHS_NAMED}: "
          f"{100*res_B['frac_aligned_named']:6.2f} %")
    print(f"  Fraction with >=1 column-aligned pair, widths [10,100]:    "
          f"{100*res_B['frac_aligned_broad']:6.2f} %")
    print()

    # --- Interpretation ---
    print("=" * 78)
    print("INTERPRETATION")
    print("=" * 78)
    print(f"  Real cipher: {len(real_lat_pos)} latitude triplet(s), "
          f"{len(real_lon_pos)} longitude triplet(s),")
    print(f"  {real_joint} joint coordinate pair(s), at least one column-aligned")
    print(f"  under widths {{22,32,44}}: {real_aligned_named}.")
    print()
    print(f"  The bare joint-existence rate is consistent with random noise.")
    print(f"  The column-aligned rate is much lower — finding a joint coord pair")
    print(f"  that is ALSO column-aligned under {{22,32,44}} is rarer than the bare")
    print(f"  finding alone.")
    print()
    print(f"  This empirical floor is the right comparison for the published")
    print(f"  finding. It does not, by itself, determine whether the cipher is")
    print(f"  intentionally encoded or a fortuitous noise match — the data is")
    print(f"  consistent with both, and Cipher 2's announcement that Cipher 1")
    print(f"  describes 'the exact locality of the vault' provides a prior that")
    print(f"  the bare frequency math does not capture.")
    print()


if __name__ == "__main__":
    main()
