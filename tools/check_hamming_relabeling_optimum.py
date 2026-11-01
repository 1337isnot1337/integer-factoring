#!/usr/bin/env python3
"""Exact F2 relabeling cutoff four via 60 complete slice intersections."""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_hamming_pfaffian import independent_ranks
from check_modular_recombination import (
    crossing,
    decode,
    encode,
    insert,
    ranks,
    represent,
)
from hamming_relabeling import seed_data

ROOT = Path(__file__).resolve().parents[1]
ELL_POSITIONS = (250, 251, 254, 255)
ELL = sum(1 << p for p in ELL_POSITIONS)


def kernel(rows, dimension):
    equations = {}
    for row in rows:
        insert(equations, row, 0)
    result = []
    for j in range(dimension):
        if j in equations:
            continue
        vector = 1 << j
        for pivot, (row, _) in sorted(equations.items()):
            if (row & vector).bit_count() % 2:
                vector ^= 1 << pivot
        assert all((row & vector).bit_count() % 2 == 0 for row in rows)
        result.append(vector)
    assert len(result) + len(equations) == dimension
    return result


def lower_bound(generators):
    forms = [
        sum(((g >> i) & 1) << j for j, g in enumerate(generators)) for i in range(256)
    ]
    histogram, slice_span, systems, checked = Counter(), {}, 0, 0
    for leg in range(4):
        for factor in range(1, 16):
            pivot = (factor & -factor).bit_length() - 1
            equations = []
            for rest in product(range(4), repeat=3):
                positions = []
                for value in range(4):
                    coordinates = list(rest)
                    coordinates.insert(leg, value)
                    position = 0
                    for coordinate in coordinates:
                        position = 4 * position + coordinate
                    positions.append(position)
                for value in range(4):
                    if value == pivot:
                        continue
                    equations.append(
                        forms[positions[value]]
                        ^ (forms[positions[pivot]] if factor >> value & 1 else 0)
                    )
            nullspace = kernel(equations, len(generators))
            histogram[len(nullspace)] += 1
            for coefficient in nullspace:
                vector = 0
                for j, generator in enumerate(generators):
                    if coefficient >> j & 1:
                        vector ^= generator
                assert (vector & ELL).bit_count() % 2 == 0
                insert(slice_span, vector, 0)
                checked += 1
            systems += 1
    target = encode(crossing(4), 4)
    assert (target & ELL).bit_count() % 2 == 1
    assert represent(slice_span, target) is None
    assert histogram == {0: 16, 2: 40, 14: 4}
    assert len(slice_span) == 60
    return {
        "complete_slice_intersection_systems": systems,
        "slice_kernel_basis_vectors_checked": checked,
        "intersection_dimension_histogram": dict(sorted(histogram.items())),
        "slice_span_dimension": len(slice_span),
        "separator_entry_positions": list(ELL_POSITIONS),
        "separator_hex": hex(ELL),
        "separator_on_crossing": 1,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    evidence = ROOT / "evidence/hamming-relabeling-optimum.json"
    saved = json.loads(evidence.read_text())
    _, images, _, _ = seed_data()
    space = {}
    for j, (vector, _) in enumerate(images):
        insert(space, vector, 1 << j)
    generators = [v for v, _ in space.values()]
    assert len(generators) == 102
    total, costs, widths = 0, [], []
    for piece in saved["witness"]:
        vector = int(piece["tensor_hex"], 16)
        mask = int(piece["source_image_mask"], 16)
        assert mask >> len(images) == 0
        replay = 0
        for j, (image, _) in enumerate(images):
            if mask >> j & 1:
                replay ^= image
        assert replay == vector and represent(space, vector) is not None
        row, column = ranks(decode(vector, 4), 4)
        assert (row, column) == independent_ranks(vector, 4)
        assert (row, column) == (piece["row_rank"], piece["column_rank"])
        total ^= vector
        costs.append(row * column)
        widths.append((row - 1).bit_length() + (column - 1).bit_length())
    assert total == encode(crossing(4), 4)
    assert max(costs) == 4 and max(widths) == 2
    result = {
        "status": "exact F2 cutoff four in the 102-dimensional relabeling span",
        "binary_width": 2,
        "span_dimension": 102,
        "exact_minimum_rank_product_cutoff": 4,
        "exact_minimum_maximum_binary_width_sum": 2,
        "witness_piece_count": len(costs),
        "crossing_entries_checked": 256,
        "lower_bound": lower_bound(generators),
        "witness": saved["witness"],
    }
    paths = [Path(__file__), ROOT / "tools/hamming_relabeling.py"]
    result["source_sha256"] = {
        p.name: sha256(p.read_bytes()).hexdigest() for p in paths
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "witness"}))


if __name__ == "__main__":
    main()
