#!/usr/bin/env python3
"""Independent polynomial and coefficient-rank checks for stable Pfaffian bounds."""

from collections import defaultdict
from itertools import product
import argparse
import hashlib
import json
from math import comb
from pathlib import Path

from check_modular_crossing import rank_two
from complement_shell_tensor import multiply
from explore_complement_edge_pencils import all_length_labels
from pfaffian_grade_types import G_GRADES, H_GRADES, H_STATES, full_grade_source_profile, type_basis
from pfaffian_particle_hole import local_exterior_table, one_particle_matrix, rank_profile
from stable_pfaffian_bounds import (
    HG_INTERTWINER_COLUMNS, centered_row_bounds, centered_score_case, score_case, stable_row_bounds,
)


def verify_HG_intertwiner():
    # Sparse polynomials give an independent arithmetic representation.
    zero, one = frozenset(), frozenset({(0, 0, 0)})
    p, q, r = (frozenset({term}) for term in ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    def add(*values):
        out = set()
        for value in values:
            out.symmetric_difference_update(value)
        return frozenset(out)
    def mul(a, b):
        out = set()
        for x in a:
            for y in b:
                term = tuple(i + j for i, j in zip(x, y))
                if term in out:
                    out.remove(term)
                else:
                    out.add(term)
        return frozenset(out)
    def divide_delta(value):
        remainder, result = set(value), set()
        while any(term[0] for term in remainder):
            term = max(term for term in remainder if term[0])
            quotient = term[0] - 1, term[1], term[2]
            result.add(quotient)
            remainder.symmetric_difference_update({term, (quotient[0], quotient[1], quotient[2] + 1)})
        assert not remainder
        return frozenset(result)
    def kron(a, b):
        return [[mul(a[i // len(b)][j // len(b)], b[i % len(b)][j % len(b)])
                 for j in range(len(a) * len(b))] for i in range(len(a) * len(b))]
    def diagonal(blocks):
        size = sum(map(len, blocks))
        result = [[zero] * size for _ in range(size)]
        offset = 0
        for block in blocks:
            for i, row in enumerate(block):
                for j, value in enumerate(row):
                    result[offset + i][offset + j] = value
            offset += len(block)
        return result
    g, h = ((r, p), (one, one)), ((one, one), (q, add(p, q, r)))
    delta = add(p, r)
    big = kron(h, g)
    pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    middle = [[divide_delta(add(mul(big[a][c], big[b][d]), mul(big[a][d], big[b][c])))
               for c, d in pairs] for a, b in pairs]
    source = kron(middle, big)
    target = diagonal([[[mul(delta, value) for value in row] for row in big]] * 2 + [
        kron(big, [[mul(value, value) for value in row] for row in g]),
        kron(big, [[mul(value, value) for value in row] for row in h])])
    columns = HG_INTERTWINER_COLUMNS
    assert rank_two(list(columns)) == 24
    source_grades = [(a + c, b + d) for a, b in H_GRADES for c, d in G_GRADES]
    target_grades = ([(a + 1, b + 1) for a, b in G_GRADES] * 2
                     + [(a + c, b + d) for a, b in G_GRADES for c, d in ((2, 0), (0, 2))]
                     + [(a + c, b + d) for a, b in G_GRADES for c, d in ((0, 0), (2, 2))])
    for row, col in product(range(24), repeat=2):
        left = add(*(source[row][k] for k in range(24) if columns[col] >> k & 1))
        right = add(*(target[k][col] for k in range(24) if columns[k] >> row & 1))
        assert left == right
        if columns[col] >> row & 1:
            assert source_grades[row] == target_grades[col]
    assert all(value == one for value in big[1])
    return {"intertwiner_rank": 24, "generic_polynomial_equations": 576,
            "middle_exterior_minors_divided_exactly": 36, "grade_preserving": True,
            "intertwiner_columns": list(columns), "constant_ancillary_row_entries": 4}


def direct_augmented_profile(pairs, length, states):
    """Construct the augmented matrices before the stable decomposition."""
    g = one_particle_matrix(6, 3, 2)
    local = local_exterior_table(6, 3, 2)
    h = [[local[((x * 4 + y) * 4 + u) * 4 + v] for v, u in H_STATES] for x, y in H_STATES]
    totals = defaultdict(int)
    for scalars in range(pairs + 1):
        for odd in range(pairs - scalars + 1):
            middle = pairs - scalars - odd
            multiplicity = comb(pairs, scalars) * comb(pairs - scalars, odd) * 2 ** (scalars + odd)
            first = type_basis(odd + 1, middle, length - scalars, states + 1 - scalars)
            second = type_basis(odd + 1, middle, length + 1 - scalars, states - scalars)
            matrices = (g,) * (odd + 1) + (h,) * middle
            rows = {}
            for i, a in enumerate(first):
                for j, b in enumerate(second):
                    value = 1 << (2 * scalars)
                    for matrix, x, y in zip(matrices, a, b):
                        value = multiply(value, matrix[x][y])
                    while value:
                        bit = value & -value
                        degree = bit.bit_length() - 1
                        if degree not in rows:
                            rows[degree] = [0] * len(first)
                        rows[degree][i] ^= 1 << j
                        value ^= bit
            for degree, matrix in rows.items():
                totals[degree] += multiplicity * rank_two(matrix)
    return {degree: rank for degree, rank in totals.items() if rank}


def verify_augmented_ranks():
    cases = []
    for pairs, length, states in ((1, 1, 1), (2, 1, 1), (2, 2, 2), (2, 1, 2),
                                  (3, 2, 2), (3, 3, 3), (3, 2, 3)):
        expected = direct_augmented_profile(pairs, length, states)
        actual = {degree: rank for degree, rank in stable_row_bounds(pairs, length, states)["row_bounds"] if rank}
        assert actual == expected, (pairs, length, states)
        cases.append({"pairs": pairs, "operator_grade": length, "state_grade": states,
                      "exact_augmented_ranks": sorted(expected.items())})
    return cases


def verify_centered_bounds():
    cases = []
    for pairs, length in ((2, 2), (3, 3), (4, 2)):
        stride = 2 * pairs + 3
        actual = rank_profile(all_length_labels(pairs, length), p=3, q=1 | (1 << stride), r=1)
        upper = centered_score_case(pairs, length)
        bounds = {(a, b): (r, c) for a, b, r, c in upper["complete_centered_jet_bounds"]}
        profile = []
        for degree, r, c in actual["profile"]:
            a, b = degree % stride, degree // stride
            ur, uc = bounds.get((a, b), (0, 0))
            assert r <= ur and c <= uc, (pairs, length, a, b, r, c, ur, uc)
            profile.append([a, b, r, c])
        cases.append({"pairs": pairs, "length": length, "exact_centered_jet_ranks": sorted(profile)})
    # An ancillary row of parameter grade one can add one augmented degree.
    # Those extra coefficients do not belong to the original tensor.
    for pairs, length, states in ((2, 1, 2), (3, 1, 3)):
        for column in (0, 2):
            bounds, _ = centered_row_bounds(pairs, length, states, 0, column)
            assert all(a + b <= length + 1 for a, b in bounds)
    return cases


def verify_full_grade_types():
    source = Path(__file__).resolve().parents[1] / "evidence/pfaffian-particle-hole.json"
    cases = []
    for previous in json.loads(source.read_text())["complete_profiles"]:
        if "energy" in previous:
            continue
        pairs, length = previous["complement_pairs"], previous["length"]
        report = full_grade_source_profile(pairs, length)
        assert report["profile"] == previous["profile"]
        cases.append({"pairs": pairs, "length": length,
                      "nonempty_types": len(report["types"]),
                      "largest_type_matrix": report["largest_type_matrix"],
                      "type_entries_evaluated": report["type_entries_evaluated"]})
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {
        "status": "pass",
        "scope": "explicit stable Pfaffian module and complete coefficient-rank upper bounds; no improved factoring exponent",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "module_sha256": hashlib.sha256(Path(__file__).with_name("stable_pfaffian_bounds.py").read_bytes()).hexdigest(),
        "generic_HG_intertwiner": verify_HG_intertwiner(),
        "exact_augmented_profile_checks": verify_augmented_ranks(),
        "centered_jet_actual_rank_controls": verify_centered_bounds(),
        "exact_full_grade_type_comparisons": verify_full_grade_types(),
    }
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
