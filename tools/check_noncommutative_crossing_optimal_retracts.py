#!/usr/bin/env python3
"""Optimal two-vertex-null retracts for a permutation involution."""

from itertools import combinations, product
import argparse
import hashlib
import json
from pathlib import Path

from check_modular_crossing import rank_two


def indices(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def involution(mask, pairs, fixed):
    low = mask & ((1 << pairs) - 1)
    high = (mask >> pairs) & ((1 << pairs) - 1)
    tail = mask >> (2 * pairs)
    return high | (low << pairs) | (tail << (2 * pairs))


def no_fixed_retract(pairs):
    columns, rows = [], []
    for i in range(max(0, pairs - 1)):
        z = (1 << i) | (1 << (pairs - 1))
        a = sum(1 << j for j in range(i + 1, pairs - 1))
        columns.append(a | ((a ^ z) << pairs))
        b = (1 << i) - 1
        rows.append(b | ((b ^ (1 << i)) << pairs))
    return columns, rows


def colored_retract(pairs, fixed, colors=None):
    if not pairs:
        return [], []
    assert fixed > 0
    colors = min(pairs, fixed) if colors is None else colors
    assert 1 <= colors <= min(pairs, fixed)
    classes = [[i for i in range(pairs) if i % colors == c] for c in range(colors)]
    columns, rows = [0] * pairs, [0] * pairs
    for c, group in enumerate(classes):
        for position, i in enumerate(group):
            column = (1 << i) | (1 << (2 * pairs + c))
            row = 1 << i
            for j in group[:position]:
                column ^= (1 << j) | (1 << (pairs + j))
            for j in group[position + 1:]:
                row ^= (1 << j) | (1 << (pairs + j))
            columns[i], rows[i] = column, row
    return columns, rows


def audit(pairs, fixed, colors=None):
    columns, rows = (colored_retract(pairs, fixed, colors) if fixed else
                     no_fixed_retract(pairs))
    dimension = pairs if fixed else max(0, pairs - 1)
    assert len(columns) == len(rows) == dimension
    assert rank_two(columns) == rank_two(rows) == dimension
    shifted = [involution(column, pairs, fixed) for column in columns]
    for i, row in enumerate(rows):
        for j, column in enumerate(columns):
            assert (row & column).bit_count() % 2 == int(i == j)
            assert not ((row & shifted[j]).bit_count() & 1)
    assert all(not ((a & b).bit_count() & 1) for a, b in product(columns, repeat=2))
    assert rank_two(columns + shifted) == 2 * dimension
    return {"swap_orbits": pairs, "fixed_labels": fixed, "labels": dimension,
            "PI_identity": True, "PJI_zero": True, "input_Gram_zero": True,
            "transverse_dimension": 2 * dimension,
            "max_input_terms": max(map(int.bit_count, columns), default=0),
            "max_output_terms": max(map(int.bit_count, rows), default=0),
            "total_input_terms": sum(map(int.bit_count, columns)),
            "total_output_terms": sum(map(int.bit_count, rows))}


def physical_complement(mask, pairs):
    return sum(((mask >> (2 * i)) & 1) << (2 * i + 1)
               | ((mask >> (2 * i + 1)) & 1) << (2 * i)
               for i in range(pairs))


def exterior_case(pairs, length):
    basis = tuple(sum(1 << i for i in subset)
                  for subset in combinations(range(2 * pairs), length))
    index = {mask: i for i, mask in enumerate(basis)}
    representatives, partners, fixed = [], [], []
    for i, mask in enumerate(basis):
        j = index[physical_complement(mask, pairs)]
        if i < j:
            representatives.append(i)
            partners.append(j)
        elif i == j:
            fixed.append(i)
    m, f = len(representatives), len(fixed)
    expected_fixed = len(tuple(combinations(range(pairs), length // 2))) if length % 2 == 0 else 0
    assert f == expected_fixed
    result = audit(m, f)
    columns, rows = colored_retract(m, f) if f else no_fixed_retract(m)
    reorder = representatives + partners + fixed

    def map_vector(vector):
        return sum(1 << reorder[i] for i in indices(vector))

    actual_columns = [map_vector(column) for column in columns]
    actual_rows = [map_vector(row) for row in rows]
    actual_J = [index[physical_complement(mask, pairs)] for mask in basis]
    shifted = [sum(1 << actual_J[i] for i in indices(column))
               for column in actual_columns]
    for i, row in enumerate(actual_rows):
        for j, column in enumerate(actual_columns):
            assert (row & column).bit_count() % 2 == int(i == j)
            assert not ((row & shifted[j]).bit_count() & 1)
    result.update({"physical_pairs": pairs, "exterior_degree": length,
                   "ambient_dimension": len(basis),
                   "simple_representative_pairing_labels": m // 2,
                   "physical_J_matrix_audited": True})
    if pairs <= 4:
        result.update({"basis_mode_masks": list(basis),
                       "J_representative_indices": representatives,
                       "J_partner_indices": partners, "J_fixed_indices": fixed,
                       "input_columns_hex": [hex(a) for a in actual_columns],
                       "output_rows_hex": [hex(a) for a in actual_rows]})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    count = 0
    for pairs in range(41):
        for fixed in range(7):
            audit(pairs, fixed)
            count += 1
    cases = [exterior_case(2, 2), exterior_case(3, 3),
             exterior_case(4, 2), exterior_case(4, 3), exterior_case(4, 4),
             exterior_case(5, 5), exterior_case(6, 6)]
    result = {"status": "Exact PI, PJI, isotropy, and transversality checks passed.",
              "scope": "Optimal dimension for simultaneous p^L and q^L cancellation; no improved rank-sum claim.",
              "abstract_involution_cases": count, "exterior_cases": cases,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
