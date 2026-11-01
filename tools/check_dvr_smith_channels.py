#!/usr/bin/env python3
"""Exact complete jet oracle and a strict gain over every constant basis."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import check_cancellation_module_automata as arithmetic
import check_exact_module_orbits as lattices
from dvr_smith_channels import DVRSmithLattices, weighted_combination
from matroid_smith_channels import binary_plane_lattices, capped_lattice


def apply(matrix, vector):
    return tuple(sum(value * vector[target] for target, value in row.items())
                 for row in matrix)


def dense(matrix):
    return [[str(row.get(j, Fraction(0))) for j in range(len(matrix))] for row in matrix]


def literal_complete_actions(graph):
    """Independent full-v-jet enumeration at precision one, before compression."""
    assert graph.precision == 1 and graph.dimension == 2 and graph.alphabet == 2
    indices = {state: i for i, state in enumerate(graph.states)}
    total = 0
    for state_index, factors in enumerate(graph.local_factors):
        bound = graph.jet_bounds[state_index]
        assert bound <= 3, "this independent exhaustive control is deliberately small"
        rank = len(factors)
        actions, cache = set(), {}
        for coefficients in product(range(1 << bound), repeat=rank):
            residue = sum((value & 1) << i for i, value in enumerate(coefficients))
            if not residue:
                continue
            order = min(valuation + (value & -value).bit_length() - 1
                        for value, (valuation, _, _) in zip(coefficients, factors) if value)
            left = weighted_combination(coefficients, factors, order)
            jet = tuple(arithmetic.jet(value, 1) for value in left)
            if jet not in cache:
                rows = tuple(tuple(arithmetic.rat(jet[(i * 4 + physical) * 2 + j]) for j in range(2))
                             for i in range(2) for physical in range(4))
                cache[jet] = indices[capped_lattice(rows, 2, 1)]
            actions.add((residue, cache[jet], order))
            total += 1
        assert actions == set(graph.candidates[state_index])
    return total


def check(maximum_length=3):
    local, point = lattices.companion_two(), Fraction(5)
    states = binary_plane_lattices(2)
    row = DVRSmithLattices(local, 1, seed_states=states, cut=0, infinity=True)
    column = DVRSmithLattices(local, 1, seed_states=states, cut=1, infinity=False)
    assert len(row.states) == len(column.states) == 10
    row_potential = (10000, 5320, 10000, 10000, 5320, 10000, 5320, 10000, 10000, 10000)
    col_potential = (100000, 94634, 63719, 62303, 94634, 62828, 94634, 63719, 54764, 61230)
    row_rate, col_rate = Fraction(266, 25), Fraction(241, 200)
    records, raw_matrices = [], []
    for graph, potential, rate in ((row, row_potential, row_rate),
                                   (column, col_potential, col_rate)):
        matrix, choices = graph.choose(potential, point)
        image = apply(matrix, potential)
        assert all(value <= rate * height for value, height in zip(image, potential))
        certificate = graph.verify(choices)
        terminal = graph.terminal_weights(point)
        prefactor = potential[graph.start] * max(
            value / potential[i] for i, value in enumerate(terminal))
        raw_matrices.append(matrix)
        records.append({"graph": graph.summary(), "policy": choices,
                        "factor_certificate": certificate,
                        "potential": list(potential), "rate": str(rate),
                        "largest_row_ratio": str(max(value / height for value, height in zip(image, potential))),
                        "transfer": dense(matrix), "terminal": list(map(str, terminal)),
                        "prefactor": str(prefactor)})
    assert records[0]["prefactor"] == "1"
    assert records[1]["prefactor"] == "47317/31414"
    upper = row_rate * col_rate
    assert upper == Fraction(32053, 2500) < Fraction(24073, 1875) < Fraction(129, 10)

    # All constant GL_r(F_2) changes, including changes across Smith orders,
    # are covered by this lower potential.  Every stationary policy matrix
    # has image at least the exact rowwise minimum computed here.
    lower_potential = (10000, 9365, 6316, 6168, 9365, 6211, 9365, 6316, 5423, 6062)
    lower_rate = Fraction(603, 500)
    constant_matrix, constant_choices = column.choose(lower_potential, point, constant_only=True)
    lower_image = apply(constant_matrix, lower_potential)
    assert all(value >= lower_rate * height for value, height in zip(lower_image, lower_potential))
    assert col_rate < lower_rate
    constant_certificate = column.verify(constant_choices)

    # Compare the efficient normalized-w enumeration with a direct complete
    # enumeration of all primitive v-jets on the four precision-one states.
    control = DVRSmithLattices(local, 1, precision=1, seed_states=binary_plane_lattices(1),
                              cut=1, infinity=False)
    independent_primitive_jets = literal_complete_actions(control)

    finite_records, literal_entries, rank_checks = [], 0, 0
    for length in range(1, maximum_length + 1):
        ranks, entries = arithmetic.coefficient_ranks(local, 1, length)
        literal_entries += entries
        words = []
        for graph, matrix, record in zip((row, column), raw_matrices, records):
            vector = graph.terminal_weights(point)
            for _ in range(length):
                vector = apply(matrix, vector)
            words.append(vector[graph.start])
            assert words[-1] <= Fraction(record["prefactor"]) * Fraction(record["rate"]) ** length
        for degree in range(length + 1):
            assert ranks[0][degree] <= (length - degree + 1) * point ** (-degree) * words[0]
            assert ranks[1][degree] <= (degree + 1) * point ** degree * words[1]
            rank_checks += 2
        finite_records.append({"length": length, "row_ranks": ranks[0],
                               "column_ranks": ranks[1], "word_sums": list(map(str, words))})
    return {"status": "pass", "point": str(point), "product_upper": str(upper),
            "row": records[0], "column": records[1],
            "constant_basis_column_lower": {
                "rate": str(lower_rate), "positive_potential": list(lower_potential),
                "smallest_row_ratio": str(min(value / height for value, height in zip(lower_image, lower_potential))),
                "minimizing_transfer": dense(constant_matrix),
                "factor_certificate": constant_certificate,
                "scope": "Every constant channel-basis policy on these ten full-tail states."},
            "independent_primitive_jets": independent_primitive_jets,
            "independent_jet_control": control.summary(),
            "literal_endpoint_entries": literal_entries, "literal_rank_inequalities": rank_checks,
            "finite_records": finite_records,
            "scope": ("A complete finite-precision oracle for regular unimodular channel-basis "
                      "changes, and a strict gain over every constant-basis policy in the same "
                      "declared state model. No factoring-exponent improvement is claimed.")}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--maximum-length", type=int, default=3)
    args = parser.parse_args()
    assert 1 <= args.maximum_length <= 4
    report = check(args.maximum_length)
    report["checker_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    report["engine_sha256"] = sha256(Path(__file__).with_name("dvr_smith_channels.py").read_bytes()).hexdigest()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
