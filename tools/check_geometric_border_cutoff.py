#!/usr/bin/env python3
"""Exact weighted exclusion through 64 for supplied formal crossing limits."""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from check_geometric_subspace_compression import geometric_complement, split_system, three_system
from check_weighted_cut_rank_budget import (
    column_flattening, crossing, encoded_fraction, exact_above_sqrt2_power,
    exact_below_sqrt2_power, lower_interpolant, matrix_rank, weights,
)
from check_weighted_cut_border_stability import (
    add, at_one, compose, constant, identity, kron, matrix_add,
    reduction, shift,
)


@lru_cache(None)
def formal_complement(d, m):
    assert d >= 2 and m >= 0
    if m == 0:
        return d * d
    if m == 1:
        return d * d // 4
    return geometric_complement(d, m)


@lru_cache(None)
def formal_budget(d, s, minimum_products=0):
    assert 0 <= s <= d and minimum_products >= 0
    return min(formal_complement(d, m) + s * s * m
               for m in range(minimum_products, max(minimum_products, 2 * d - 1) + 1))


def harmonic_bound(d):
    return Fraction(2 * d * d, 2 * d - 1) + 2 * sum(
        (w * formal_budget(d, s) for s, w in enumerate(weights(d), 1)), Fraction())


def scalar_controls():
    closure = knots = one_product = 0
    for d in range(2, 71):
        assert formal_complement(d, 1) == formal_complement(d, 2) == d * d // 4
        assert formal_complement(d, 1) == max(s * (d - s) for s in range(d + 1))
        one_product += 1
        for s in range(d + 1):
            for minimum in (0, 1, 2, d, 2 * d - 1, 2 * d + 2):
                extended = min(formal_complement(d, m) + s * s * m
                               for m in range(minimum, max(minimum, 4 * d + 2) + 1))
                assert formal_budget(d, s, minimum) == extended
                closure += 1
            assert formal_budget(d, s) >= s * d
        for j in range(1, d + 1):
            assert lower_interpolant(d, j * j, 1) == j
            assert lower_interpolant(d, 1, j * j) == j
            knots += 2
    return {"one_product_maxima": one_product, "finite_minimum_closures": closure,
            "exact_interpolation_knots": knots}


def moving_case(d, s, q, As, Ds, pole, prime):
    ell, p, m = d - q, d - s, len(As)
    unit = identity(d, prime)
    # Different unipotent changes force both initial-coordinate subspaces to move.
    G = [[add(unit[i][j], {1: 1} if i == j + 1 else {}, prime)
          for j in range(d)] for i in range(d)]
    Ginv = [[{i - j: (-1) ** (i - j) % prime} if i >= j else {}
             for j in range(d)] for i in range(d)]
    H = [[add(unit[i][j], {2: 1} if i == j + 1 else {}, prime)
          for j in range(d)] for i in range(d)]
    Hinv = [[{2 * (i - j): (-1) ** (i - j) % prime} if i >= j else {}
             for j in range(d)] for i in range(d)]
    assert compose(G, Ginv, prime) == compose(Ginv, G, prime) == unit
    assert compose(H, Hinv, prime) == compose(Hinv, H, prime) == unit
    inclusion = [row[:s] for row in H]
    quotient = Ginv[ell:]
    assert any(entry for row in inclusion[s:] for entry in row)
    assert any(entry for row in quotient for entry in row[:ell])
    assert matrix_rank(reduction(inclusion), prime) == s
    assert matrix_rank(reduction(quotient), prime) == q
    B = [[{} for _ in range(d * d)] for _ in range(d * d)]
    for index, (Ablock, Dblock) in enumerate(zip(As, Ds)):
        P0 = [row + [0] * p for row in Ablock] + [[0] * s + row for row in Dblock]
        assert len(P0) == d and all(len(row) == d for row in P0)
        P = compose(compose(G, constant(P0, prime), prime), Hinv, prime)
        P = shift(P, -pole - index, 1, prime)
        assert all(not value for row in compose(compose(quotient, P, prime), inclusion, prime)
                   for value in row)
        Q = constant([[int(i * d + j == index) for j in range(d)] for i in range(d)], prime)
        B = matrix_add(B, kron(P, Q, prime), prime)
    assert matrix_rank(column_flattening(at_one(B, prime), d), prime) == m
    n = d * d
    perturbation = constant([[(i + 2 * j + 1) % prime for j in range(n)] for i in range(n)], prime)
    T = matrix_add(constant(crossing(d), prime), shift(perturbation, 1, 1, prime), prime)
    A = matrix_add(T, shift(B, 0, -1, prime), prime)
    assert matrix_add(A, B, prime) == T
    assert any(power < 0 for row in A for entry in row for power in entry)
    assert any(power < 0 for row in B for entry in row for power in entry)
    projection, restriction = kron(quotient, unit, prime), kron(inclusion, unit, prime)
    projected_B = compose(compose(projection, B, prime), restriction, prime)
    assert all(not entry for row in projected_B for entry in row)
    projected_A = compose(compose(projection, A, prime), restriction, prime)
    projected_T = compose(compose(projection, T, prime), restriction, prime)
    assert projected_A == projected_T
    assert matrix_rank(reduction(projected_A), prime) == s * q
    # Both unsaturated choices have the correct generic span but lose reduction rank.
    assert matrix_rank(reduction(shift(inclusion, 1, 1, prime)), prime) == 0
    assert matrix_rank(reduction(shift(quotient, 1, 1, prime)), prime) == 0


def moving_controls():
    cases = 0
    systems = []
    for d, s, q, a, b in ((3, 1, 1, 1, 1), (4, 1, 2, 1, 2),
                           (5, 1, 2, 2, 2), (5, 2, 1, 2, 2)):
        As, Ds, _ = split_system(d, s, q, a, b)
        systems.append((d, s, q, As, Ds))
    for d, s, q, a, b in ((3, 1, 1, 0, 0), (4, 1, 2, 0, 0),
                           (4, 2, 1, 0, 0), (5, 2, 2, 0, 0),
                           (4, 1, 1, 1, 0), (4, 1, 1, 0, 1)):
        As, Ds, _ = three_system(d, s, q, a, b)
        systems.append((d, s, q, As, Ds))
    for prime in (2, 3, 5):
        for pole in (1, 3):
            for d, s, q, As, Ds in systems:
                moving_case(d, s, q, As, Ds, pole, prime)
                cases += 1
    return {"moving_input_output_cases": cases, "primes": [2, 3, 5],
            "initial_pole_orders": [1, 3], "unsaturated_negative_controls": 2 * cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert 577 ** 2 - 2 * 408 ** 2 == 1
    assert 1393 ** 2 - 2 * 985 ** 2 == -1
    scalar = scalar_controls()
    moving = moving_controls()
    dimensions = []
    for d in range(2, 65):
        value = harmonic_bound(d)
        assert exact_above_sqrt2_power(value, d)
        dimensions.append({"d": d, "bound": encoded_fraction(value), "above_target": True})
    assert harmonic_bound(64) == Fraction(2347092541, 6538455)
    control = harmonic_bound(65)
    assert control == Fraction(259385858696, 708864255)
    assert exact_below_sqrt2_power(control, 65)
    here = Path(__file__).resolve()
    modules = [here, Path(__import__('check_geometric_subspace_compression').__file__),
               Path(__import__('check_weighted_cut_rank_budget').__file__),
               Path(__import__('check_weighted_cut_border_stability').__file__)]
    result = {"status": "PASS", "scope": "Supplied one-parameter formal crossing limits with Laurent pieces",
              "last_certified_dimension": 64, "scalar_controls": scalar,
              "moving_lattice_controls": moving, "dimensions": dimensions,
              "dimension_65_control": {"bound": encoded_fraction(control), "below_target": True,
                  "meaning": "The scalar certificate stops here; no realizable seed is asserted."},
              "exact_half_complement_imported": False,
              "source_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in modules}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'dimensions'}, indent=2))


if __name__ == '__main__':
    main()
