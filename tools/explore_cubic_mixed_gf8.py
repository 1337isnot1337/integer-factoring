#!/usr/bin/env python3
"""Complete F8 searches in five structured unions of cubic coefficient spaces."""

import argparse
import hashlib
import json
from math import isqrt
from pathlib import Path
import random
import subprocess

from explore_cubic_complement_outputs import ComplementOutputs
from explore_cubic_complement_gf8 import matrices, scale, INVERSE, TABLE, rank8


def cases():
    def scaled(scalar, matrix):
        return [[TABLE[scalar][value] for value in row] for row in matrix]
    shift = [[int(j == (i + 1) % 7) for j in range(7)] for i in range(7)]
    jordan = [[int(j == i + 1) for j in range(7)] for i in range(7)]
    first, second = [1, 2, 4, 3, 6, 7, 5], [1, 1, 2, 2, 4, 4, 3]
    rank_one = [[TABLE[a][b] for b in second] for a in first]
    one_row = [[2 * int(i == 0 and j == 1) for j in range(7)] for i in range(7)]
    other_row = [[4 * int(i == 2 and j == 3) for j in range(7)] for i in range(7)]
    return [
        ("cyclic_shift_1_and_z", shift, scaled(2, shift)),
        ("cyclic_shift_z_and_z_squared", scaled(2, shift), scaled(4, shift)),
        ("jordan_shift_1_and_z", jordan, scaled(2, jordan)),
        ("dense_rank_one_1_and_z", rank_one, scaled(2, rank_one)),
        ("two_distinct_changed_output_rows", one_row, other_row),
    ]


def mixed_basis(family, first, second):
    rows, columns = [], []
    for complement in (first, second):
        r, c = matrices(family, complement)
        rows.extend(r)
        columns.extend(c)
    selected, pivots = [], {}
    for index, matrix in enumerate(rows):
        value = sum(row << (147 * i) for i, row in enumerate(matrix))
        while value:
            pivot = (value.bit_length() - 1) // 3
            scalar = (value >> (3 * pivot)) & 7
            if pivot not in pivots:
                pivots[pivot] = scale(value, INVERSE[scalar], 2401)
                selected.append(index)
                break
            value ^= scale(pivots[pivot], scalar, 2401)
    assert selected[0] == 0
    return [rows[i] for i in selected], [columns[i] for i in selected], selected


def native_input(rows, columns):
    lines = [f"7 {len(rows)} 16"]
    for cut in (rows, columns):
        for matrix in cut:
            for row in matrix:
                planes = [sum(((row >> (3 * j + bit)) & 1) << j for j in range(49))
                          for bit in range(3)]
                lines.append(" ".join(format(plane, "x") for plane in planes))
    return "\n".join(lines) + "\n"


def decision(points, dimension):
    pivots = {}
    for point in points:
        value = point["code"]
        while value:
            pivot = (value.bit_length() - 1) // 3
            scalar = (value >> (3 * pivot)) & 7
            if pivot not in pivots:
                pivots[pivot] = scale(value, INVERSE[scalar], dimension)
                break
            value ^= scale(pivots[pivot], scalar, dimension)
    residual = 1
    while residual:
        pivot = (residual.bit_length() - 1) // 3
        if pivot not in pivots:
            break
        residual ^= scale(pivots[pivot], (residual >> (3 * pivot)) & 7, dimension)
    result = {"affordable_span_dimension": len(pivots), "endpoint_remainder": residual}
    if residual:
        return dict(result, status="excluded: affordable span misses the crossing")
    products = {}
    for point in points:
        for scalar in range(1, 8):
            products[scale(point["code"], scalar, dimension)] = point["product"]
    best, pair = None, None
    denominator = 10 ** 12
    for code, first in products.items():
        if (1 ^ code) not in products:
            continue
        second = products[1 ^ code]
        roots = [isqrt(value * denominator * denominator) for value in (first, second)]
        upper = sum(root + int(root * root != value * denominator * denominator)
                    for root, value in zip(roots, (first, second)))
        if best is None or upper < best:
            best, pair = upper, [code, 1 ^ code]
    if best is not None and best ** 985 < 7 ** 1393 * denominator ** 985:
        return dict(result, status="strict two-piece candidate: independent full audit required",
                    witness=pair, witness_upper_numerator=best, cost_denominator=denominator)
    return dict(result, status="unresolved: affordable span contains crossing; no certified two-piece seed",
                best_two_piece_upper_numerator=best, best_two_piece_witness=pair,
                cost_denominator=denominator)


def verify(rows, columns, screen, seed):
    dimension = len(rows)
    events = [json.loads(line) for line in screen.splitlines()]
    assert events[0] == {"event": "start", "labels": 7, "dimension": dimension,
                         "field_size": 8, "field_polynomial": 11,
                         "integer_cost_gate": 16, "product_cap": 255}
    total = (8 ** dimension - 1) // 7
    assert events[-1]["event"] == "complete" and events[-1]["projective_tested"] == total
    points = [event for event in events if event["event"] == "eligible"]
    by_code = {point["code"]: point for point in points}
    assert len(points) == len(by_code) == events[-1]["eligible"]
    selected = set(by_code)
    rng = random.Random(seed)
    for _ in range(12):
        code = rng.randrange(1, 8 ** dimension)
        pivot = (code.bit_length() - 1) // 3
        selected.add(scale(code, INVERSE[(code >> (3 * pivot)) & 7], dimension))
    for code in selected:
        pair = []
        for cut in (rows, columns):
            matrix = [0] * 49
            for i, source in enumerate(cut):
                scalar = (code >> (3 * i)) & 7
                if scalar:
                    matrix = [a ^ scale(b, scalar, 49) for a, b in zip(matrix, source)]
            pair.append(rank8(matrix, 49))
        assert (code in by_code) == (pair[0] * pair[1] < 256)
        if code in by_code:
            point = by_code[code]
            assert pair == [point["row_rank"], point["column_rank"]]
            assert point["product"] == pair[0] * pair[1]
    return dict(decision(points, dimension), projective_vectors_tested=total,
                affordable_count=len(points), independent_rank_vectors=len(selected),
                native_screen_sha256=hashlib.sha256(screen.encode()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path)
    parser.add_argument("--replay", type=Path)
    parser.add_argument("--case", type=str)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.native and not args.replay:
        parser.error("a native evaluator is required for new enumeration")
    family = ComplementOutputs(7)
    hashes = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
              for name in (Path(__file__).name, "screen_cubic_gf8.cpp",
                           "explore_cubic_complement_gf8.py", "explore_cubic_complement_outputs.py")}
    result = {"field": "F2[z]/(z^3+z+1)", "source_sha256": hashes,
              "cases": [], "complete": False}
    saved = None
    if args.replay:
        saved = json.loads(args.replay.read_text())
        assert saved["complete"] and saved["source_sha256"] == hashes
        existing = {record["name"]: record for record in saved["cases"]}
    for index, (name, first, second) in enumerate(cases()):
        if args.case and name != args.case:
            continue
        rows, columns, selected = mixed_basis(family, first, second)
        assert len(rows) <= 8
        payload = native_input(rows, columns)
        if saved:
            screen = existing[name]["native_screen"]
        else:
            screen = subprocess.run([str(args.native)], input=payload, text=True,
                                    stdout=subprocess.PIPE, check=True).stdout
        record = dict(name=name, complement_matrices=[first, second],
                      full_union_dimension=len(rows), selected_raw_coefficient_indices=selected,
                      native_input_sha256=hashlib.sha256(payload.encode()).hexdigest(),
                      native_screen=screen, **verify(rows, columns, screen, 20261030 + index))
        if saved:
            assert record == existing[name]
        result["cases"].append(record)
        args.output.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")
        print(json.dumps({key: value for key, value in record.items()
                          if key not in ("native_screen", "complement_matrices")}), flush=True)
    result["complete"] = len(result["cases"]) == (1 if args.case else 5)
    result["all_recorded_mixed_spaces_excluded_over_F8"] = all(
        record["status"].startswith("excluded:") for record in result["cases"])
    result["projective_vectors_tested"] = sum(record["projective_vectors_tested"] for record in result["cases"])
    result["scope"] = "Each recorded space is the full F8 span of all five raw coefficient directions from both stated output maps."
    args.output.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
