# beale_coordinate_grid_view.py

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


def make_pairs(numbers):
    full_string = "".join(str(n) for n in numbers)
    pairs = [full_string[i:i + 2] for i in range(0, len(full_string) - 1, 2)]
    return full_string, pairs


def find_triplet(pairs, target):
    for i in range(len(pairs) - 2):
        if pairs[i:i + 3] == target:
            return i
    return -1


def print_grid_relationship(pairs, lat_idx, lon_idx, width, context=4):
    lat_row = lat_idx // width
    lat_col = lat_idx % width

    lon_row = lon_idx // width
    lon_col = lon_idx % width

    print("\n" + "=" * 80)
    print(f"GRID WIDTH: {width}")
    print("=" * 80)

    print(f"Latitude starts at pair index : {lat_idx}")
    print(f"Longitude starts at pair index: {lon_idx}")
    print(f"Distance between starts       : {lon_idx - lat_idx}")
    print()
    print(f"Latitude grid position        : row {lat_row}, column {lat_col}")
    print(f"Longitude grid position       : row {lon_row}, column {lon_col}")
    print(f"Rows apart                    : {lon_row - lat_row}")
    print(f"Same starting column?         : {lat_col == lon_col}")

    if lat_col == lon_col:
        print("\nRESULT: The DMS triplets are vertically aligned in this grid.")
    else:
        print("\nRESULT: The DMS triplets are NOT vertically aligned in this grid.")

    # Print only the useful slice of each row around the target columns
    start_col = max(0, lat_col - context)
    end_col = min(width, lat_col + 3 + context)

    lat_row_values = pairs[lat_row * width:(lat_row + 1) * width]
    lon_row_values = pairs[lon_row * width:(lon_row + 1) * width]

    print("\nColumn view around the coordinate triplets:")
    print()

    # Column numbers
    col_line = "COL:      "
    for c in range(start_col, end_col):
        col_line += f"{c:>4}"
    print(col_line)

    # Marker row
    marker_line = "          "
    for c in range(start_col, end_col):
        if lat_col <= c <= lat_col + 2:
            marker_line += "   ▲"
        else:
            marker_line += "    "
    print(marker_line)

    # Latitude row
    lat_line = f"ROW {lat_row:<3} : "
    for c in range(start_col, end_col):
        value = lat_row_values[c] if c < len(lat_row_values) else "--"
        if lat_col <= c <= lat_col + 2:
            lat_line += f"[{value:>2}]"
        else:
            lat_line += f" {value:>2} "
    print(lat_line + "   <-- Latitude")

    # Longitude row
    lon_line = f"ROW {lon_row:<3} : "
    for c in range(start_col, end_col):
        value = lon_row_values[c] if c < len(lon_row_values) else "--"
        if lon_col <= c <= lon_col + 2:
            lon_line += f"[{value:>2}]"
        else:
            lon_line += f" {value:>2} "
    print(lon_line + "   <-- Longitude")

    print(marker_line.replace("▲", "▼"))
    print()
    print("Stacked DMS alignment:")
    print(f"          Latitude : [{pairs[lat_idx]}] [{pairs[lat_idx + 1]}] [{pairs[lat_idx + 2]}]")
    print(f"          Longitude: [{pairs[lon_idx]}] [{pairs[lon_idx + 1]}] [{pairs[lon_idx + 2]}]")


def main():
    latitude_target = ["37", "12", "21"]
    longitude_target = ["79", "23", "16"]

    full_string, pairs = make_pairs(cipher_numbers)

    lat_idx = find_triplet(pairs, latitude_target)
    lon_idx = find_triplet(pairs, longitude_target)

    print("=" * 80)
    print("BEALE CIPHER 1 — DMS COORDINATE GRID RELATIONSHIP")
    print("=" * 80)

    print(f"Total concatenated digit length : {len(full_string)}")
    print(f"Total two-digit pairs           : {len(pairs)}")
    print()

    if lat_idx == -1:
        print("Latitude triplet 37 12 21 was NOT found.")
        return

    if lon_idx == -1:
        print("Longitude triplet 79 23 16 was NOT found.")
        return

    distance = lon_idx - lat_idx

    print(f"Latitude found  : {pairs[lat_idx]}° {pairs[lat_idx + 1]}' {pairs[lat_idx + 2]}\" N")
    print(f"  Pair index    : {lat_idx}")
    print()
    print(f"Longitude found : {pairs[lon_idx]}° {pairs[lon_idx + 1]}' {pairs[lon_idx + 2]}\" W")
    print(f"  Pair index    : {lon_idx}")
    print()
    print(f"Distance between triplet starts: {distance} two-digit cells")

    print("\nGrid widths from 10 to 100 that vertically align the two triplets:")
    valid_widths = []
    for width in range(10, 101):
        if distance % width == 0:
            valid_widths.append(width)

    print(valid_widths)

    # Print detailed views for the most important widths
    for width in [22, 32, 44]:
        if distance % width == 0:
            print_grid_relationship(pairs, lat_idx, lon_idx, width)


if __name__ == "__main__":
    main()
