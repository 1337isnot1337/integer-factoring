#!/usr/bin/env python3
"""Complete F4 coefficient-space tests for selected cubic complementary maps."""

import argparse
import hashlib
from itertools import product
import json
from math import isqrt, sqrt
from pathlib import Path
import random

from explore_complement_edge_pencils import minimum_decomposition
from explore_cubic_complement_outputs import ComplementOutputs


MULTIPLY = ((0, 0, 0, 0), (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2))
INVERSE = (0, 1, 3, 2)


def scale(value, scalar, mask):
    if scalar == 0:
        return 0
    if scalar == 1:
        return value
    low, high = value & mask, (value >> 1) & mask
    if scalar == 2:
        return high | ((low ^ high) << 1)
    return (low ^ high) | (low << 1)


def rank4(rows, width, low_pivot=False):
    mask = ((1 << (2 * width)) - 1) // 3
    pivots = {}
    for row in rows:
        while row:
            pivot = ((row & -row).bit_length() - 1) // 2 if low_pivot else (row.bit_length() - 1) // 2
            scalar = (row >> (2 * pivot)) & 3
            if pivot not in pivots:
                pivots[pivot] = scale(row, INVERSE[scalar], mask)
                break
            row ^= scale(pivots[pivot], scalar, mask)
    return len(pivots)


def matrices(family, complement):
    d = family.pairs
    output = [[(i, 1)] + [(d + j, value) for j, value in enumerate(complement[i]) if value]
              for i in range(d)]
    planes = [(sum((value & 1) << j for j, value in row),
               sum(((value >> 1) & 1) << j for j, value in row)) for row in output]
    rows = [[0] * (d * d) for _ in range(5)]
    columns = [[0] * (d * d) for _ in range(5)]
    for k, y, v in product(range(5), range(d), range(d)):
        parent = family.parent[k][y][v]
        for x in range(d):
            low = high = 0
            for index, scalar in output[x]:
                if scalar & 1:
                    low ^= parent[index]
                if scalar & 2:
                    high ^= parent[index]
            for u, (other_low, other_high) in enumerate(planes):
                hh = (high & other_high).bit_count() & 1
                first = ((low & other_low).bit_count() & 1) ^ hh
                second = ((low & other_high).bit_count() & 1) ^ ((high & other_low).bit_count() & 1) ^ hh
                value = first | (second << 1)
                rows[k][x * d + y] ^= value << (2 * (v * d + u))
                columns[k][x * d + u] ^= value << (2 * (v * d + y))
    assert rows[0] == columns[0] == [1 << (2 * i) for i in range(d * d)]
    for x, y, v, u in product(range(d), repeat=4):
        h = MULTIPLY[complement[x][y]][complement[u][v]]
        gram = int(x == u)
        for j in range(d):
            gram ^= MULTIPLY[complement[x][j]][complement[u][j]]
        assert ((rows[4][x * d + y] >> (2 * (v * d + u))) & 3) == h
        assert (((rows[1][x * d + y] ^ rows[3][x * d + y]) >> (2 * (v * d + u))) & 3) == (h ^ (gram if y == v else 0))
    return rows, columns


def profile(family, complement):
    rows, columns = matrices(family, complement)
    d, width = family.pairs, family.pairs ** 2
    tensor_width = width * width
    full_mask = ((1 << (2 * tensor_width)) - 1) // 3
    pivots, selected = {}, []
    for index, matrix in enumerate(rows):
        value = sum(row << (2 * width * i) for i, row in enumerate(matrix))
        while value:
            pivot = (value.bit_length() - 1) // 2
            scalar = (value >> (2 * pivot)) & 3
            if pivot not in pivots:
                pivots[pivot] = scale(value, INVERSE[scalar], full_mask)
                selected.append(index)
                break
            value ^= scale(pivots[pivot], scalar, full_mask)
    assert selected[0] == 0
    rows, columns = [[cut[i] for i in selected] for cut in (rows, columns)]
    dimension = len(selected)
    ranks = [[0, 0] for _ in range(4 ** dimension)]
    mask = ((1 << (2 * width)) - 1) // 3
    coefficient_mask = ((1 << (2 * dimension)) - 1) // 3
    checked = 0
    for pivot in range(dimension):
        for tail in range(4 ** pivot):
            code = (1 << (2 * pivot)) | tail
            pair = []
            for cut in (rows, columns):
                matrix = [0] * width
                for i, source in enumerate(cut):
                    scalar = (code >> (2 * i)) & 3
                    if scalar:
                        matrix = [a ^ scale(b, scalar, mask) for a, b in zip(matrix, source)]
                value = rank4(matrix, width)
                # Both pivot orders are cheap at 49 rows and are checked for every point.
                assert value == rank4(matrix, width, low_pivot=True)
                pair.append(value)
            for scalar in (1, 2, 3):
                ranks[scale(code, scalar, coefficient_mask)] = pair
            checked += 1
    upper, witness = minimum_decomposition([sqrt(a * b) for a, b in ranks], target=1)
    lower, _ = minimum_decomposition([isqrt(a * b) for a, b in ranks], target=1)
    return {"complement_matrix": complement, "dimension": dimension,
            "basis_degrees": [family.degrees[i] for i in selected],
            "projective_points_checked": checked,
            "found_decomposition_cost": upper, "integer_shortest_lower_bound": lower,
            "witness": witness, "rank_table": ranks,
            "field4_seed_criterion_excluded": lower ** 70 > d ** 99}


def candidates(d, seed):
    for scalar in (0, 1, 2, 3):
        yield "scalar_identity", [[scalar * int(i == j) for j in range(d)] for i in range(d)]
    for scalar in (2, 3):
        for rank in range(1, d):
            yield "scalar_projection", [[scalar * int(i == j and i < rank) for j in range(d)] for i in range(d)]
        for shift in (1, 2, 3):
            yield "scalar_nilpotent", [[scalar * int(j == i + shift) for j in range(d)] for i in range(d)]
        for vector in (1, 3, 7, 15, (1 << d) - 1):
            yield "scalar_rank_one", [[scalar * int((vector >> i & 1) and (vector >> j & 1))
                                       for j in range(d)] for i in range(d)]
    rng = random.Random(seed)
    for case in range(24):
        if case % 3 == 0:
            diagonal = [rng.randrange(4) for _ in range(d)]
            yield "random_diagonal", [[diagonal[i] * int(i == j) for j in range(d)] for i in range(d)]
        elif case % 3 == 1:
            first, second = [[rng.randrange(4) for _ in range(d)] for _ in range(2)]
            yield "random_rank_one", [[MULTIPLY[first[i]][second[j]] for j in range(d)] for i in range(d)]
        else:
            yield "random_dense", [[rng.randrange(4) for j in range(d)] for i in range(d)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=int, default=7)
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    family = ComplementOutputs(args.pairs)
    result = {"pairs": args.pairs, "field": "F2[w]/(w^2+w+1)", "seed": args.seed,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "cases": [], "complete": False}
    for name, matrix in candidates(args.pairs, args.seed):
        record = dict(profile(family, matrix), family=name)
        result["cases"].append(record)
        print(json.dumps({"cases": len(result["cases"]), "family": name,
                          "cost": record["found_decomposition_cost"],
                          "lower": record["integer_shortest_lower_bound"],
                          "dimension": record["dimension"]}), flush=True)
        args.output.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")
    result["complete"] = True
    result["all_recorded_spaces_excluded_over_F4"] = all(r["field4_seed_criterion_excluded"] for r in result["cases"])
    args.output.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
