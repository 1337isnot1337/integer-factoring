#!/usr/bin/env python3
"""Exact adaptive-channel gain beyond a fixed factorization's word rate."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from adaptive_smith_lattices import AdaptiveSmithLattices, InconclusiveSearch
import check_cancellation_module_automata as fixed
import check_exact_module_orbits as exact
import check_minimal_boundary_modules as minimal


def apply(matrix, vector):
    return tuple(sum(value * vector[target] for target, value in row.items())
                 for row in matrix)


def dense(matrix):
    return [[str(row.get(j, Fraction(0))) for j in range(len(matrix))] for row in matrix]


def check(maximum_length=3):
    local = exact.companion_two()
    for (x, y, u, v), matrix in local.items():
        assert all(matrix[i][0].bit_count() % 2 == int(i == 0 and x == v and y == u)
                   for i in range(2))

    point = Fraction(2)
    row = AdaptiveSmithLattices(local, 1, cut=0, side=0, infinity=True)
    column = AdaptiveSmithLattices(local, 1, cut=1, side=0, infinity=False)
    assert (len(row.states), len(column.states)) == (6, 5)
    row_matrix, row_terminal = row.weighted_matrix(point)
    col_matrix, col_terminal = column.weighted_matrix(point)

    row_potential = tuple(map(Fraction, (3, 3, 4, 4, 4, 3)))
    row_rate = Fraction(6)
    assert apply(row_matrix, row_potential) == tuple(row_rate * value for value in row_potential)
    assert all(value <= height / 3 for value, height in zip(row_terminal, row_potential))
    row_prefactor = row_potential[0] / 3
    assert row_prefactor == 1

    col_potential = tuple(map(Fraction, (1000, 911, 730, 1075, 911)))
    col_rate = Fraction(41, 18)
    col_image = apply(col_matrix, col_potential)
    assert all(value <= col_rate * height for value, height in zip(col_image, col_potential))
    col_constant = max(value / height for value, height in zip(col_terminal, col_potential))
    col_prefactor = col_constant * col_potential[0]
    assert col_prefactor == Fraction(1000, 911)
    adaptive_bound = row_rate * col_rate
    assert adaptive_bound == Fraction(41, 3)

    # Rebuild the old, fixed-channel exact orbits independently.  These are
    # full lattices, not finite-precision approximations of their word rate.
    old_orbits, old_matrices, old_records = [], [], []
    for cut, infinity in ((0, True), (1, False)):
        factors = fixed.local_factors(local, cut, 0, (0, 0), infinity)
        representation = minimal.regularize(
            minimal.minimize(factors["channels"], (fixed.ONE, fixed.ZERO),
                             (fixed.ONE, fixed.ZERO)), factors["exponents"])
        orbit = exact.ExactModuleOrbit(representation["channels"])
        matrix = orbit.weighted_matrix(representation["exponents"], point)
        old_orbits.append(orbit)
        old_matrices.append(matrix)
        old_records.append({"states": len(orbit.states),
                            "exponents": representation["exponents"],
                            "transfer": dense(matrix)})
    assert tuple(len(orbit.states) for orbit in old_orbits) == (2, 7)
    old_row_potential = tuple(map(Fraction, (4, 3)))
    assert apply(old_matrices[0], old_row_potential) == tuple(
        row_rate * value for value in old_row_potential)
    classes = (0, 1, 2, 2, 3, 1, 2)
    quotient = exact.lump_matrix(old_matrices[1], classes)
    old_col_potential = tuple(map(Fraction, (1000, 840, 947, 758)))
    old_col_lower = Fraction(16, 7)
    old_col_image = apply(quotient, old_col_potential)
    assert all(value >= old_col_lower * height
               for value, height in zip(old_col_image, old_col_potential))
    fixed_lower = row_rate * old_col_lower
    assert fixed_lower == Fraction(96, 7) > adaptive_bound

    # Literal coefficient ranks use polynomial matrix products of the
    # original crossing; they do not use a lattice or its factorization.
    finite_records = []
    literal_endpoint_entries = 0
    rank_inequalities = 0
    for length in range(1, maximum_length + 1):
        ranks, endpoint_entries = fixed.coefficient_ranks(local, 1, length)
        literal_endpoint_entries += endpoint_entries
        row_words = fixed.weighted_value(row_matrix, row_terminal, length)
        col_words = fixed.weighted_value(col_matrix, col_terminal, length)
        assert row_words <= row_prefactor * row_rate ** length
        assert col_words <= col_prefactor * col_rate ** length
        for degree in range(length + 1):
            assert ranks[0][degree] <= (length - degree + 1) * point ** (-degree) * row_words
            assert ranks[1][degree] <= (degree + 1) * point ** degree * col_words
            assert ranks[0][degree] * ranks[1][degree] <= (
                (length + 1) ** 2 * row_prefactor * col_prefactor * adaptive_bound ** length)
            rank_inequalities += 3
        finite_records.append({
            "length": length, "row_ranks": ranks[0], "column_ranks": ranks[1],
            "row_word_sum": str(row_words), "column_word_sum": str(col_words),
        })

    # Independent scalar and embedded-scalar controls recover the elementary
    # Smith profile; an auxiliary identity creates no artificial advantage.
    scalar_controls = []
    for dimension in (1, 2):
        elementary, degree = fixed.elementary(dimension)
        for cut, infinity in ((0, True), (1, False)):
            control = AdaptiveSmithLattices(
                elementary, degree, cut=cut, side=0, infinity=infinity)
            assert len(control.states) == 1
            assert control.transitions == [((0, 0), (1, 0), (1, 0), (2, 0))]
            scalar_controls.append({"dimension": dimension, "infinity": infinity,
                                    "profile": [order for order, _ in control.transitions[0]]})

    # Full-tail state construction has exact regular membership witnesses.
    # Its checked cases are comparisons, not a general monotonicity theorem.
    capped_controls = []
    for precision in (1, 2, 3):
        capped_row = AdaptiveSmithLattices(local, 1, cut=0, side=0, infinity=True,
                                          precision=precision)
        capped_col = AdaptiveSmithLattices(local, 1, cut=1, side=0, infinity=False,
                                          precision=precision)
        rm, rt = capped_row.weighted_matrix(point)
        cm, ct = capped_col.weighted_matrix(point)
        for length in range(1, maximum_length + 1):
            row_words = fixed.weighted_value(rm, rt, length)
            col_words = fixed.weighted_value(cm, ct, length)
            ranks = finite_records[length - 1]
            for degree in range(length + 1):
                assert ranks["row_ranks"][degree] <= (length - degree + 1) * point ** (-degree) * row_words
                assert ranks["column_ranks"][degree] <= (degree + 1) * point ** degree * col_words
                rank_inequalities += 2
        capped_controls.append({"precision": precision,
                                "states": [len(capped_row.states), len(capped_col.states)]})
    try:
        AdaptiveSmithLattices(local, 1, cut=0, infinity=True, maximum_states=1)
    except InconclusiveSearch:
        pass
    else:
        raise AssertionError("a state-budget failure was accepted as closure")

    return {
        "status": "pass", "field": 2, "point": str(point),
        "adaptive_rate_product_upper": str(adaptive_bound),
        "fixed_exact_word_rate_product_lower": str(fixed_lower),
        "row_rate": str(row_rate), "column_rate_upper": str(col_rate),
        "fixed_column_rate_lower": str(old_col_lower),
        "row_prefactor": str(row_prefactor), "column_prefactor": str(col_prefactor),
        "row_positive_potential": list(map(str, row_potential)),
        "column_positive_potential": list(map(str, col_potential)),
        "fixed_row_positive_potential": list(map(str, old_row_potential)),
        "fixed_column_positive_lower_potential": list(map(str, old_col_potential)),
        "row_transfer": dense(row_matrix), "column_transfer": dense(col_matrix),
        "fixed_column_quotient": dense(quotient),
        "row": row.summary(include_factors=True),
        "column": column.summary(include_factors=True),
        "fixed_orbits": old_records, "finite_records": finite_records,
        "scalar_controls": scalar_controls, "capped_controls": capped_controls,
        "literal_endpoint_entries": literal_endpoint_entries,
        "literal_rank_inequalities": rank_inequalities,
        "state_budget_rejected": True,
        "scope": ("State-dependent local channel refactorization beats the exact "
                  "fixed-factorization grouped-word certificate for one crossing "
                  "at the same point. No improved factoring exponent or global "
                  "optimum over other factorizations is claimed."),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--maximum-length", type=int, default=3)
    args = parser.parse_args()
    assert 1 <= args.maximum_length <= 4
    report = check(args.maximum_length)
    report["checker_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    report["engine_sha256"] = sha256(Path(__file__).with_name("adaptive_smith_lattices.py").read_bytes()).hexdigest()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
