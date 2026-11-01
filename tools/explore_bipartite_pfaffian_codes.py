#!/usr/bin/env python3
"""Exact ternary crossings on disjoint involution-invariant monomial codes.

The operator labels and state labels are different sets of equal size.
Every coefficient is retained. Both unwanted endpoint tensors vanish.
"""

import argparse
from collections import defaultdict
from functools import lru_cache
import hashlib
from itertools import product
import json
from math import isqrt, sqrt
from pathlib import Path
import random
import time

from complement_shell_tensor import multiply
from explore_complement_edge_pencils import (
    all_length_labels, evaluate_local, independent_coefficients,
    minimum_decomposition, syndrome,
)
from check_modular_crossing import rank_two


def mode_complement(word):
    return tuple((0, 2, 1, 3)[mask] for mask in word)


def invariant_partition_space(labels):
    lookup = {word: i for i, word in enumerate(labels)}
    orbits = []
    for i, word in enumerate(labels):
        j = lookup[mode_complement(word)]
        if i <= j:
            orbits.append((i,) if i == j else (i, j))
    target = len(labels) // 2
    if len(labels) % 2:
        raise ValueError("the full label count must be even")
    @lru_cache(None)
    def count(position, remaining):
        if remaining == 0:
            return 1
        if position == len(orbits) or remaining < 0:
            return 0
        return count(position + 1, remaining - len(orbits[position])) + count(position + 1, remaining)

    total = count(1, target - len(orbits[0]))

    def unrank(index):
        if not 0 <= index < total:
            raise ValueError("partition index is out of range")
        # Fix the first orbit to remove exchange A,B duplication.
        chosen = [0]
        remaining = target - len(orbits[0])
        for position in range(1, len(orbits)):
            if remaining == 0:
                break
            included = count(position + 1, remaining - len(orbits[position]))
            if index < included:
                chosen.append(position)
                remaining -= len(orbits[position])
            else:
                index -= included
        assert remaining == 0 and index == 0
        first = tuple(i for k in chosen for i in orbits[k])
        present = set(first)
        return first, tuple(i for i in range(len(labels)) if i not in present)

    return total, unrank


def coefficient_basis(labels, first, second):
    d = len(first)
    assert d == len(second) and not set(first) & set(second)
    a, b = [labels[i] for i in first], [labels[i] for i in second]
    assert set(map(mode_complement, a)) == set(a)
    assert set(map(mode_complement, b)) == set(b)
    length = sum(x.bit_count() for x in a[0])
    base = length + 1
    local = evaluate_local((2, 1 << base, 1, 1, 1, 1))
    groups = defaultdict(list)
    for x, y in product(range(d), repeat=2):
        groups[syndrome(a[x]) ^ syndrome(b[y])].append((x, y))
    keys = sorted(groups)
    positions = {pair: (block, i) for block, key in enumerate(keys)
                 for i, pair in enumerate(groups[key])}
    size = length * base + 1
    rows = [[[0] * len(groups[key]) for key in keys] for _ in range(size)]
    cols = [[[0] * len(groups[key]) for key in keys] for _ in range(size)]
    digest = hashlib.sha256()
    count = 0
    for block, key in enumerate(keys):
        group = groups[key]
        for i, (x, y) in enumerate(group):
            for j, (v, u) in enumerate(group):
                value = 1
                for aa, bb, cc, ee in zip(a[x], b[y], b[u], a[v]):
                    value = multiply(value, local[((aa * 4 + bb) * 4 + cc) * 4 + ee])
                    if not value:
                        break
                assert value & 1 == int(x == v and y == u)
                assert not (value >> length) & 1
                assert not (value >> (length * base)) & 1
                digest.update(value.to_bytes((size + 7) // 8, "little"))
                cblock, ci = positions[x, u]
                other, cj = positions[v, y]
                assert cblock == other
                while value:
                    bit = value & -value
                    degree = bit.bit_length() - 1
                    assert degree % base + degree // base <= length
                    rows[degree][block][i] ^= 1 << j
                    cols[degree][cblock][ci] ^= 1 << cj
                    value ^= bit
                count += 1
    indices = independent_coefficients(rows)
    assert indices[0] == 0
    return {
        "name": "bipartite_pfaffian", "labels": d, "length": length,
        "operator_indices": first, "state_indices": second,
        "ambient_occupations": labels, "basis_dimension": len(indices),
        "basis_monomials_p_q_r": [[i % base, i // base, length - i % base - i // base]
                                    for i in indices],
        "block_sizes": [len(groups[key]) for key in keys],
        "endpoint_coordinates": [1] + [0] * (len(indices) - 1),
        "basis_row_blocks": [[[hex(row) for row in block] for block in rows[i]] for i in indices],
        "basis_column_blocks": [[[hex(row) for row in block] for block in cols[i]] for i in indices],
        "complete_nonzero_syndrome_entries_checked": count,
        "full_dense_entry_count": d ** 4,
        "complete_polynomial_sha256": digest.hexdigest(),
        "both_unwanted_endpoints_zero": True,
    }


def exact_binary_cost(data):
    dimension = data["basis_dimension"]
    if dimension > 16:
        raise ValueError("binary exact search is limited to dimension16")
    matrices = [data["basis_row_blocks"], data["basis_column_blocks"]]
    matrices = [[[[int(row, 16) for row in block] for block in tensor] for tensor in cut]
                for cut in matrices]
    working = [[[0] * size for size in data["block_sizes"]] for _ in range(2)]
    ranks = [[0, 0] for _ in range(1 << dimension)]
    previous = 0
    for index in range(1, 1 << dimension):
        gray = index ^ (index >> 1)
        changed = (gray ^ previous).bit_length() - 1
        for cut in range(2):
            for target, source in zip(working[cut], matrices[cut][changed]):
                for i, row in enumerate(source):
                    target[i] ^= row
            ranks[gray][cut] = sum(rank_two(block) for block in working[cut])
        previous = gray
    upper, witness = minimum_decomposition([sqrt(r * c) for r, c in ranks], target=1)
    lower, _ = minimum_decomposition([isqrt(r * c) for r, c in ranks], target=1)
    return {"minimum_cost_float": upper, "integer_lower_bound": lower,
            "decomposition_masks": witness, "all_binary_combination_ranks": ranks,
            "finite_threshold_float": data["labels"] ** sqrt(2),
            "binary_exclusion_by_integer_comparison": lower ** 70 > data["labels"] ** 99}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", default="2:2,3:3")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--export-directory", type=Path)
    args = parser.parse_args()
    records = []
    if args.export_directory:
        args.export_directory.mkdir(parents=True, exist_ok=True)
    for spec in args.cases.split(","):
        pairs, length = map(int, spec.split(":"))
        labels = tuple(all_length_labels(pairs, length))
        total, unrank = invariant_partition_space(labels)
        selected = range(total)
        if args.limit and total > args.limit:
            selected = sorted(random.Random(args.seed).sample(selected, args.limit))
        best = float("inf")
        for index in selected:
            started = time.monotonic()
            data = coefficient_basis(labels, *unrank(index))
            data["name"] = f"bipartite_P{pairs}_L{length}_partition{index}"
            analysis = exact_binary_cost(data)
            record = {"pairs": pairs, "length": length, "partition_index": index,
                      "partitions_available": total, "labels": data["labels"],
                      "basis_dimension": data["basis_dimension"],
                      "operator_indices": data["operator_indices"],
                      "state_indices": data["state_indices"],
                      "complete_polynomial_sha256": data["complete_polynomial_sha256"],
                      **analysis, "seconds": time.monotonic() - started}
            records.append(record)
            if analysis["minimum_cost_float"] < best:
                best = analysis["minimum_cost_float"]
                if args.export_directory:
                    path = args.export_directory / f"bipartite_P{pairs}_L{length}_best.json"
                    path.write_text(json.dumps(data, sort_keys=True, indent=2) + "\n")
            report = {"status": "exact finite binary coefficient-space searches; no asymptotic improvement asserted",
                      "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      "cases": records}
            temporary = args.output.with_suffix(args.output.suffix + ".tmp")
            temporary.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n")
            temporary.replace(args.output)
            print(json.dumps({k: v for k, v in record.items()
                              if k not in ("all_binary_combination_ranks", "operator_indices", "state_indices")}), flush=True)


if __name__ == "__main__":
    main()
