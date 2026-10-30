#!/usr/bin/env python3
"""Exact F4 scan for a nonbinary section changing all 24 output rows."""

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import subprocess

from explore_two_even_site_sections import parent_coefficients
from check_modular_crossing import rank_two
from check_symplectic_seed_tensor_fields import integer_gate


def field_scale(low, high, scalar):
    if scalar == 0:
        return 0, 0
    if scalar == 1:
        return low, high
    if scalar == 2:
        return high, low ^ high
    return low ^ high, low


def field_rank(low, high):
    pivots = {}
    for first, second in zip(low, high):
        while first | second:
            pivot = (first | second).bit_length() - 1
            scalar = ((first >> pivot) & 1) + 2 * ((second >> pivot) & 1)
            if pivot not in pivots:
                pivots[pivot] = field_scale(first, second, (0, 1, 3, 2)[scalar])
                break
            a, b = field_scale(*pivots[pivot], scalar)
            first ^= a
            second ^= b
    return len(pivots)


def embed(low, high):
    size = len(low)
    return [a | (b << size) for a, b in zip(low, high)] + [
        b | ((a ^ b) << size) for a, b in zip(low, high)]


def construct():
    parent = parent_coefficients()
    labels = parent["code"]["labels"]
    syndromes = parent["code"]["syndromes"]
    bucket = defaultdict(list)
    for i in range(labels):
        for j in range(labels):
            bucket[syndromes[i] ^ syndromes[j]].append((i, j))
    groups = [bucket[key] for key in sorted(bucket)]
    lifted, original = [], []
    for cut, name in enumerate(("rows", "columns")):
        binary_generators, field_generators = [], []
        for matrix in parent[name]:
            blocks = []
            for group in groups:
                size = len(group)
                low, high = [0] * size, [0] * size
                for i, (x, y) in enumerate(group):
                    for j, (v, u) in enumerate(group):
                        values = []
                        for first in range(2):
                            for second in range(2):
                                if cut == 0:
                                    rowpair, colpair = (2 * x + first, y), (v, 2 * u + second)
                                else:
                                    rowpair, colpair = (2 * x + first, 2 * y + second), (v, u)
                                block, row = parent["positions"][2 * cut][rowpair]
                                other, col = parent["positions"][2 * cut + 1][colpair]
                                assert block == other
                                values.append((matrix[block][row] >> col) & 1)
                        # lambda=omega, coefficients (1+omega,omega).
                        # Their products are omega,1,1,1+omega.
                        low[i] |= (values[1] ^ values[2] ^ values[3]) << j
                        high[i] |= (values[0] ^ values[3]) << j
                blocks.append((low, high))
            field_generators.append(blocks)
            binary_generators.append([embed(low, high) for low, high in blocks])
            binary_generators.append([embed(high, [a ^ b for a, b in zip(low, high)])
                                      for low, high in blocks])
        lifted.append(binary_generators)
        original.append(field_generators)
    for cut in lifted:
        assert all(block == [1 << i for i in range(len(block))] for block in cut[0])
    # Each bit of a matrix entry supplies one exact linear-independence test.
    coordinates = []
    for block in range(len(groups)):
        for row in range(2 * len(groups[block])):
            values = [matrix[block][row] for matrix in lifted[0]]
            for col in range(2 * len(groups[block])):
                coordinates.append(sum(((value >> col) & 1) << i for i, value in enumerate(values)))
    assert rank_two(coordinates) == len(lifted[0]) == 12
    return labels, groups, lifted, original, parent["ambient_entries_checked"]


def native_input(labels, groups, cuts):
    dimension = len(cuts[0])
    # Restriction of scalars doubles matrix ranks.  Pad the doubled cut by
    # zeros to the square dimension (2d)^2 required by the native format.
    result = [f"{2 * labels} {dimension} {2 * integer_gate(labels)}"]
    for matrices in cuts:
        result.append(str(len(groups) + 1))
        result.append(f"1 {2 * labels * labels}")
        result.extend("0" for _ in range(dimension))
        for block, group in enumerate(groups):
            size = 2 * len(group)
            result.append(f"{size} 1")
            for matrix in matrices:
                for row in matrix[block]:
                    result.append(" ".join(f"{(row >> shift) & ((1 << 64) - 1):x}"
                                           for shift in range(0, size, 64)))
    return "\n".join(result) + "\n"


def direct_ranks(generators, mask):
    ranks = []
    for cut in generators:
        total = 0
        for block in range(len(cut[0])):
            size = len(cut[0][block][0])
            low, high = [0] * size, [0] * size
            for i, matrix in enumerate(cut):
                scalar = (mask >> (2 * i)) & 3
                for row, (a, b) in enumerate(zip(*matrix[block])):
                    x, y = field_scale(a, b, scalar)
                    low[row] ^= x
                    high[row] ^= y
            total += field_rank(low, high)
        ranks.append(total)
    return ranks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--screen-output", type=Path)
    args = parser.parse_args()
    labels, groups, lifted, original, contracted = construct()
    encoded_input = native_input(labels, groups, lifted)
    run = subprocess.run([str(args.native), "--field2"], input=encoded_input,
                         text=True, stdout=subprocess.PIPE, check=True)
    if args.screen_output:
        args.screen_output.write_text(run.stdout)
    events = [json.loads(line) for line in run.stdout.splitlines()]
    assert events[0]["field_size"] == 2 and events[0]["dimension"] == 12
    assert events[0]["code_start"] == 1 and events[0]["code_count"] == 4095
    assert events[-1]["event"] == "complete" and events[-1]["projective_tested"] == 4095
    eligible, pivots = [], {}
    for entry in events:
        if entry["event"] != "eligible":
            continue
        mask = entry["index"]
        ranks = direct_ranks(original, mask)
        assert [2 * value for value in ranks] == [entry["row_rank"], entry["column_rank"]]
        assert ranks[0] * ranks[1] < integer_gate(labels) ** 2
        eligible.append({"F4_coefficient_code": mask, "row_rank": ranks[0], "column_rank": ranks[1]})
        vector = mask
        while vector:
            pivot = vector.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = vector
                break
            vector ^= pivots[pivot]
    remainder = 1
    while remainder and remainder.bit_length() - 1 in pivots:
        remainder ^= pivots[remainder.bit_length() - 1]
    assert len(eligible) == events[-1]["eligible"] == 15
    assert len(pivots) == 8 and remainder == 1
    eligible_codes = {entry["F4_coefficient_code"] for entry in eligible}
    assert all(sum((0, 2, 3, 1)[(code >> (2 * i)) & 3] << (2 * i)
                   for i in range(6)) in eligible_codes for code in eligible_codes)
    result = {"labels": labels, "field_size": 4, "section_parameter": "omega, omega^2=omega+1",
              "output_section": "(1+omega) lower + omega upper at every logical label",
              "coefficient_dimension_F4": 6, "restriction_of_scalars_dimension_F2": 12,
              "nonzero_vectors_tested": 4095, "projective_F4_classes_represented": 1365,
              "integer_cost_gate_over_F4": integer_gate(labels),
              "eligible_vectors": eligible, "eligible_projective_F4_classes": 5,
              "eligible_span_dimension_F2": 8, "eligible_span_dimension_F4": 4,
              "crossing_remainder": remainder, "excluded_over_F4": True,
              "literal_parent_entries_checked": contracted,
              "all_eligible_rank_pairs_independently_checked_in_F4": len(eligible),
              "native_lift_rank_multiplier": 2, "native_zero_padding_dimension": 2 * labels * labels,
              "native_input_sha256": hashlib.sha256(encoded_input.encode()).hexdigest(),
              "native_screen_sha256": hashlib.sha256(run.stdout.encode()).hexdigest(),
              "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                  for name in (Path(__file__).name, "explore_two_even_site_sections.py",
                               "screen_weighted_pfaffian_span.cpp", "screen_pfaffian_span_f4.cpp")}}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
