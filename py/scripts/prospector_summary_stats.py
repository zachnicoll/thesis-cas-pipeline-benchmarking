from typing import List, Tuple, Dict

# x = chunk size, y = threshold score
# 5
# 10
# 15
# 20
# 25
# 50
#     100   150	200	250	300	1000

thresh_map: Dict[int, int] = {
    5: 0,
    10: 1,
    15: 2,
    20: 3,
    25: 4,
    50: 5
}


def parse_file(filename: str):
    data_map: Dict[int, Dict[int, Tuple[float, float, float]]] = {}
    current_chunk_thresh = [0, 0]
    rows: List[Tuple[float, float, float]] = []

    for line in open(filename):
        if line.startswith("==>"):
            current_chunk_thresh = list(map(lambda x: int(x.replace(".csv <==", "")), line.split('_')[-2:]))
            if current_chunk_thresh[0] not in data_map:
                data_map[current_chunk_thresh[0]] = {}
        elif not line.startswith("Family") and len(line) > 1:
            p, r, a = map(lambda x: float(x), line.split(',')[1:])
            rows.append((p, r, a))
        elif len(line) == 1 and len(rows) > 0:
            p_sum, r_sum, a_sum = 0.0, 0.0, 0.0
            for row in rows:
                p_sum += row[0]
                r_sum += row[1]
                a_sum += row[2]

            num_rows = len(rows)
            av_p = p_sum / num_rows
            av_r = r_sum / num_rows
            av_a = a_sum / num_rows

            data_map[current_chunk_thresh[0]][current_chunk_thresh[1]] = (av_p, av_r, av_a)
            rows = []

    return data_map


def to_csv(file_id: str, data_map: Dict[int, Dict[int, Tuple[float, float, float]]], metric: int):
    f = open(f"{file_id}_summarized.csv", "w+")

    f.writelines(",5,10,15,20,25,50\n")

    for chunk in data_map:
        lines = [""] * 6

        for threshold in data_map[chunk]:
            m = data_map[chunk][threshold][metric]
            lines[thresh_map[threshold]] = f"{m}"

        lines = f"{chunk}," + ",".join(lines) + "\n"
        f.writelines(lines)


# prosp_50 = parse_file(f"prospector_output/summary/MERGED_STATS_NEW_PROSP_50.csv")
#
# for i in range(3):
#     to_csv(f"prosp_50_{i}", prosp_50, i)

prosp_5 = parse_file(f"prospector_output/summary/MERGED_STATS_NEW_PROSP_5.csv")

for i in range(3):
    to_csv(f"prosp_5_{i}", prosp_5, i)
