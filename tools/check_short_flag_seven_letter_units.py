#!/usr/bin/env python3
"""Short-flag rank-one-row-cut units for the seven-letter quadratic inverse."""
import json
from itertools import combinations

from check_quadratic_crossings import binary_rank
from check_quadratic_inverse_paths import (
    compose,
    construct,
    inverse,
    kernel_flag,
    quotient_jet,
    row_cut,
)
from check_block_smith_crossings import smith_at_precision

Q = 7
N = Q * Q
IDENTITY = [1 << i for i in range(N)]


def inverse_row_cut(matrix):
    answer = [0] * N
    for x in range(Q):
        for y in range(Q):
            column = matrix[Q * x + y]
            for u in range(Q):
                for v in range(Q):
                    if column >> (Q * u + v) & 1:
                        answer[Q * v + y] ^= 1 << (Q * x + u)
    return answer


def rank_one(a, b):
    return [a if b >> j & 1 else 0 for j in range(N)]


def nilpotency_index(unit):
    nilpotent = [a ^ b for a, b in zip(IDENTITY, unit)]
    power = [1 << i for i in range(N)]
    for index in range(1, N + 2):
        if all(value == 0 for value in power):
            return index - 1
        power = compose(nilpotent, power)
    return None


def first_jet_multiplicity(unit, one_levels):
    nilpotent = [a ^ b for a, b in zip(IDENTITY, unit)]
    basis, levels = kernel_flag(nilpotent)
    inverse_basis = inverse(basis)
    triangular = compose(compose(inverse_basis, nilpotent), basis)
    constant = row_cut(unit, Q)
    assert binary_rank(constant) == 1
    jets = {}
    for level in range(2, max(levels) + 1):
        selected = [column if levels[j] == level else 0
                    for j, column in enumerate(triangular)]
        physical = compose(compose(basis, selected), inverse_basis)
        jets[level] = quotient_jet(row_cut(compose(physical, unit), Q), constant)
    rows = [0] * (N - 1)
    for level in one_levels:
        rows = [a ^ b for a, b in zip(rows, jets[level])]
    multiplicity = binary_rank(rows)
    physical, _ = construct(Q, basis, levels, triangular, inverse_basis, set(one_levels))
    degree = max(value.bit_length() - 1 for row in physical for value in row)
    return basis, levels, triangular, inverse_basis, physical, degree, multiplicity


def main():
    base = [(1 << (Q * b + a)) ^ ((1 << (Q * (b - 1) + a)) if b else 0)
            for a in range(Q) for b in range(Q)]
    constant = row_cut(base, Q)
    a = constant[0]
    b = 0
    for j, column in enumerate(constant):
        if column:
            b |= 1 << j
    assert all(constant[j] == (a if b >> j & 1 else 0) for j in range(N))
    assert inverse_row_cut(constant) == base
    base_index = nilpotency_index(base)

    witnesses = []
    for which in ("a", "b"):
        for bit in range(N):
            aa = a ^ (1 << bit) if which == "a" else a
            bb = b ^ (1 << bit) if which == "b" else b
            unit = inverse_row_cut(rank_one(aa, bb))
            if binary_rank(unit) != N:
                continue
            index = nilpotency_index(unit)
            if index is None:
                continue
            witnesses.append((index, which, bit))
    witnesses.sort()
    minimum = witnesses[0][0]
    best = witnesses[0]
    aa = a ^ (1 << best[2]) if best[1] == "a" else a
    bb = b ^ (1 << best[2]) if best[1] == "b" else b
    unit = inverse_row_cut(rank_one(aa, bb))
    assert binary_rank(unit) == N
    assert binary_rank(row_cut(unit, Q)) == 1
    assert nilpotency_index(unit) == minimum
    _, levels = kernel_flag([x ^ y for x, y in zip(IDENTITY, unit)])
    maximum_level = max(levels)
    possible = list(range(2, maximum_level + 1))
    candidates = []
    for count in (1, 2):
        if count > len(possible):
            continue
        best_choice = None
        for selection in combinations(possible, count):
            _, _, _, _, physical, degree, multiplicity = first_jet_multiplicity(unit, selection)
            if best_choice is None or multiplicity < best_choice["first_jet_multiplicity"]:
                best_choice = {
                    "one_levels": list(selection),
                    "degree": degree,
                    "first_jet_multiplicity": multiplicity,
                }
        candidates.append(best_choice)
    _, _, _, _, physical, _, _ = first_jet_multiplicity(unit, candidates[0]["one_levels"])
    degree = candidates[0]["degree"]
    zero_orders = smith_at_precision(physical, degree + 1)
    assert zero_orders == (0,) * N
    report = {
        "alphabet": Q,
        "dimension": N,
        "base_nilpotency_index": base_index,
        "single_flips_screened": 2 * N,
        "invertible_nilpotent_flips": len(witnesses),
        "minimum_nilpotency_index": minimum,
        "minimum_witnesses": [
            {"factor": which, "bit": bit} for index, which, bit in witnesses if index == minimum
        ],
        "best_witness": {"factor": best[1], "bit": best[2]},
        "best_kernel_levels": levels,
        "candidates": candidates,
        "degree_one_physical_zero_orders": list(zero_orders),
        "scope": (
            "The D=1 physical matrix has invertible constant term, so its zero cut "
            "has all Smith orders zero. The unit is an infinity-side building block; "
            "the positive zero orders of the crossing profile require auxiliary dimension > 1."
        ),
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
