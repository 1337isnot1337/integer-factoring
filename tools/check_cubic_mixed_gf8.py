#!/usr/bin/env python3
"""Audit complete mixed F8 spaces, separating functionals, and scalar closure."""

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import subprocess

from explore_cubic_complement_outputs import ComplementOutputs
from explore_cubic_complement_gf8 import INVERSE, TABLE, matrices, rank8, scale
from explore_cubic_mixed_gf8 import cases, mixed_basis, native_input


def flatten(matrix):
    return sum(row << (147 * i) for i, row in enumerate(matrix))


class Coordinates:
    """Low-pivot elimination, opposite to the mixed explorer's basis selection."""

    def __init__(self, basis):
        self.dimension = len(basis)
        self.pivots = {}
        for i, matrix in enumerate(basis):
            value, tag = flatten(matrix), 1 << (3 * i)
            while value:
                pivot = ((value & -value).bit_length() - 1) // 3
                scalar = (value >> (3 * pivot)) & 7
                if pivot not in self.pivots:
                    self.pivots[pivot] = (scale(value, INVERSE[scalar], 2401),
                                          scale(tag, INVERSE[scalar], self.dimension))
                    break
                old, old_tag = self.pivots[pivot]
                value ^= scale(old, scalar, 2401)
                tag ^= scale(old_tag, scalar, self.dimension)
            assert value, "The selected coefficient matrices must be independent."

    def express(self, matrix):
        value, tag = flatten(matrix), 0
        while value:
            pivot = ((value & -value).bit_length() - 1) // 3
            if pivot not in self.pivots:
                return None
            scalar = (value >> (3 * pivot)) & 7
            old, old_tag = self.pivots[pivot]
            value ^= scale(old, scalar, 2401)
            tag ^= scale(old_tag, scalar, self.dimension)
        return tag


def reconstruct(basis, coordinates):
    result = 0
    for i, matrix in enumerate(basis):
        result ^= scale(flatten(matrix), (coordinates >> (3 * i)) & 7, 2401)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("evidence/cubic-mixed-F8.json"))
    parser.add_argument("--native", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    assert data["complete"] and data["all_recorded_mixed_spaces_excluded_over_F8"]
    assert len(data["cases"]) == 5 and data["projective_vectors_tested"] == 5692269
    for name, digest in data["source_sha256"].items():
        assert hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == digest
    assert 16 ** 70 > 7 ** 99 and 99 ** 2 > 2 * 70 ** 2

    family = ComplementOutputs(7)
    zero = [[0] * 7 for _ in range(7)]
    zero_rows, zero_columns = matrices(family, zero)
    total = affordable = independent = native_replays = 0
    checks, bases, raw_cache = [], [], {}
    for index, ((name, first, second), record) in enumerate(zip(cases(), data["cases"])):
        assert record["name"] == name and record["complement_matrices"] == [first, second]
        rows, columns, selected = mixed_basis(family, first, second)
        dimension = len(rows)
        assert dimension == record["full_union_dimension"]
        assert selected == record["selected_raw_coefficient_indices"]
        payload = native_input(rows, columns)
        assert hashlib.sha256(payload.encode()).hexdigest() == record["native_input_sha256"]
        screen = record["native_screen"]
        assert hashlib.sha256(screen.encode()).hexdigest() == record["native_screen_sha256"]
        events = [json.loads(line) for line in screen.splitlines()]
        assert events[0] == {"event": "start", "labels": 7, "dimension": dimension,
                             "field_size": 8, "field_polynomial": 11,
                             "integer_cost_gate": 16, "product_cap": 255}
        points = [event for event in events if event["event"] == "eligible"]
        count = (8 ** dimension - 1) // 7
        assert events[-1]["event"] == "complete"
        assert events[-1]["projective_tested"] == record["projective_vectors_tested"] == count
        assert events[-1]["eligible"] == record["affordable_count"] == len(points)
        codes = [point["code"] for point in points]
        assert len(codes) == len(set(codes))
        for point in points:
            code = point["code"]
            assert code >> (3 * ((code.bit_length() - 1) // 3)) == 1
            assert point["product"] == point["row_rank"] * point["column_rank"] < 256
            assert code & 7 == 0
        assert rank8(codes, dimension) == dimension - 1 == record["affordable_span_dimension"]
        assert rank8(codes + [1], dimension) == dimension
        assert record["endpoint_remainder"] == 1
        coordinate_map = Coordinates(rows)
        raw_coordinates = []
        for complement in (first, second):
            key = tuple(map(tuple, complement))
            if key not in raw_cache:
                raw_cache[key] = matrices(family, complement)
            raw_rows, raw_columns = raw_cache[key]
            for row_matrix, column_matrix in zip(raw_rows, raw_columns):
                coordinate = coordinate_map.express(row_matrix)
                assert coordinate is not None
                assert reconstruct(rows, coordinate) == flatten(row_matrix)
                assert reconstruct(columns, coordinate) == flatten(column_matrix)
                raw_coordinates.append(coordinate)
        zero_coordinates = [coordinate_map.express(matrix) for matrix in zero_rows]
        contains_zero_section = all(value is not None for value in zero_coordinates)
        if contains_zero_section:
            for matrix, coordinate in zip(zero_columns, zero_coordinates):
                assert reconstruct(columns, coordinate) == flatten(matrix)
        if index < 4:
            assert contains_zero_section
            assert first != zero
            assert second == [[TABLE[2][value] for value in row] for row in first]
        if args.native and index == 0:
            replay = subprocess.run([str(args.native)], input=payload, text=True,
                                    stdout=subprocess.PIPE, check=True).stdout
            assert replay == screen
            native_replays += 1
        checks.append({"name": name, "dimension": dimension,
                       "projective_vectors": count, "affordable_points": len(points),
                       "affordable_span_dimension": dimension - 1,
                       "separating_functional": [1] + [0] * (dimension - 1),
                       "raw_directions_in_selected_basis": raw_coordinates,
                       "contains_zero_section": contains_zero_section,
                       "zero_section_coordinates": zero_coordinates,
                       "contains_all_scalar_multiples_over_extensions": index < 4})
        bases.append((rows, columns, coordinate_map))
        total += count
        affordable += len(points)
        independent += record["independent_rank_vectors"]
    assert total == 5692269 and affordable == 10452 and independent == 10512
    cyclic_change_of_basis = []
    for row, column in zip(bases[1][0], bases[1][1]):
        coordinate = bases[0][2].express(row)
        assert coordinate is not None
        assert reconstruct(bases[0][1], coordinate) == flatten(column)
        cyclic_change_of_basis.append(coordinate)
    assert rank8(cyclic_change_of_basis, 7) == 7

    # A literal occupation-kernel control for the two new, distinct-row maps.
    literal = 0
    low = family.space.outputs
    high = [tuple(2 if a == 1 else a for a in word) for word in low]
    for complement in cases()[-1][1:]:
        rows, columns = raw_cache[tuple(map(tuple, complement))]
        outputs = [[(low[i], 1)] + [(high[j], value) for j, value in enumerate(complement[i]) if value]
                   for i in range(7)]
        for x, y, v, u in product(range(7), repeat=4):
            polynomial = 0
            for (first, a), (second, b) in product(outputs[x], outputs[u]):
                binary = family.space.kernel(first, second, y, v)
                spread = sum(((binary >> degree) & 1) << (3 * degree) for degree in range(13))
                polynomial ^= spread * TABLE[a][b]
            for k, degree in enumerate(family.degrees):
                value = (polynomial >> (3 * degree)) & 7
                assert ((rows[k][7 * x + y] >> (3 * (7 * v + u))) & 7) == value
                assert ((columns[k][7 * x + u] >> (3 * (7 * v + y))) & 7) == value
                literal += 1
    result = {"source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
              "configurations": checks, "projective_vectors_tested": total,
              "affordable_projective_points": affordable,
              "source_independent_rank_vectors": independent,
              "literal_physical_coefficient_entries_checked": literal,
              "deterministic_native_control_replays": native_replays,
              "cyclic_configurations_have_the_same_full_space": True,
              "second_cyclic_basis_in_first": cyclic_change_of_basis,
              "all_five_affordable_spans_are_the_E_coefficient_hyperplane": True,
              "all_recorded_F8_spaces_fail_seed_criterion": True,
              "no_claim_about_ranks_over_field_extensions": True}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "configurations"}), flush=True)


if __name__ == "__main__":
    main()
