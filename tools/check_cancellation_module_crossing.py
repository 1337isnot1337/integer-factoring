#!/usr/bin/env python3
"""An explicit crossing with a strict module/min-plus certificate gap."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import check_cancellation_module_automata as module


def local_crossing():
    local, _ = module.elementary(2)
    for (x, y, u, v), matrix in local.items():
        if v == 0:
            matrix[0][1] ^= 6
        if x == y == 1 and u == v == 0:
            matrix[1][0] ^= 6
    for (x, y, u, v), matrix in local.items():
        for a, b in product(range(2), repeat=2):
            assert matrix[a][b].bit_count() % 2 == int(x == v and y == u and a == b)
    return local


def capped_transfer(factors, cap, s, state_order):
    costs = []
    for choices in factors["channels"]:
        entries = [[cap] * 2 for _ in range(2)]
        for a, b in product(range(2), repeat=2):
            values = [module.rorder(matrix[a][b]) for matrix in choices if matrix[a][b] != module.ZERO]
            entries[a][b] = min(cap, min(values)) if values else cap
        costs.append(entries)
    states = tuple(state_order)
    index = {state: i for i, state in enumerate(states)}
    matrix = [[Fraction(0)] * len(states) for _ in states]
    transitions = []
    for i, state in enumerate(states):
        row = []
        for exponent, entries in zip(factors["exponents"], costs):
            changed = tuple(min(state[a] + entries[a][b] for a in range(2)) for b in range(2))
            increment = min(changed)
            target = tuple(value - increment for value in changed)
            assert target in index
            matrix[i][index[target]] += s ** (exponent - increment)
            row.append((increment, index[target]))
        transitions.append(row)
    initial = index[(0, cap)]
    terminal = tuple(s ** (-state[0]) for state in states)
    return matrix, initial, terminal, transitions


def matrix_value(matrix, initial, terminal, length):
    vector = [Fraction(int(i == initial)) for i in range(len(matrix))]
    for _ in range(length):
        vector = [sum(vector[i] * matrix[i][j] for i in range(len(matrix))) for j in range(len(matrix))]
    return sum(a * b for a, b in zip(vector, terminal))


def check(maximum_word=3):
    local = local_crossing()
    s = Fraction(3, 2)
    row = module.local_factors(local, 0, 1, (0, 0), True)
    column = module.local_factors(local, 1, 1, (0, 0), False)
    assert row["exponents"] == (2, 2, 1, -1)
    assert column["exponents"] == (0, -1, -1, -2)
    row_auto = module.ModuleAutomaton(row["channels"], 2)
    col_auto = module.ModuleAutomaton(column["channels"], 1)
    row_matrix, row_terminal = row_auto.weighted_matrix(row["exponents"], s)
    col_matrix, col_terminal = col_auto.weighted_matrix(column["exponents"], s)
    assert row_matrix == [{0: Fraction(19, 6), 1: Fraction(9, 4)},
                          {0: Fraction(9, 4), 1: Fraction(19, 6)}]
    assert col_matrix == [{0: Fraction(19, 9), 1: Fraction(2, 3)},
                          {1: Fraction(25, 9)}]
    row_beta, col_beta = Fraction(65, 12), Fraction(25, 9)
    assert all(sum(values.values()) == row_beta for values in row_matrix)
    assert all(sum(values.values()) == col_beta for values in col_matrix)
    assert row_terminal == (Fraction(1), Fraction(4, 9))
    assert col_terminal == (Fraction(1), Fraction(1))
    B = row_beta * col_beta
    assert B == Fraction(1625, 108)

    old_row, old_start, old_end, old_transitions = capped_transfer(
        row, 2, s, ((2, 0), (0, 0), (0, 2)))
    expected_old = [
        [Fraction(13, 6), Fraction(1), Fraction(9, 4)],
        [Fraction(9, 4), Fraction(13, 6), Fraction(9, 4)],
        [Fraction(9, 4), Fraction(1), Fraction(13, 6)],
    ]
    assert old_row == expected_old
    old_h, old_lower = (Fraction(4, 5), Fraction(1), Fraction(4, 5)), Fraction(17, 3)
    assert all(sum(value * h for value, h in zip(row_, old_h)) >= old_lower * h
               for row_, h in zip(old_row, old_h))
    assert old_lower > row_beta
    old_col, old_col_start, old_col_end, _ = capped_transfer(
        column, 1, s, ((0, 1), (0, 0)))
    assert old_col == [[Fraction(19, 9), Fraction(2, 3)], [Fraction(0), Fraction(25, 9)]]
    assert all(sum(row_) == col_beta for row_ in old_col)
    old_B_lower = old_lower * col_beta
    assert old_B_lower == Fraction(425, 27) > B

    counts = Counter()
    rank_records = []
    for length in range(1, maximum_word + 1):
        ranks, endpoints = module.coefficient_ranks(local, 2, length)
        counts["literal_endpoint_entries"] += endpoints
        row_scores, col_scores = Counter(), Counter()
        for word in product(range(4), repeat=length):
            row_bound, col_bound = row_auto.word_bound(word), col_auto.word_bound(word)
            assert row_bound >= module.minplus_bound(row_auto, word)
            assert col_bound >= module.minplus_bound(col_auto, word)
            if length <= 2:
                row_actual = module.literal_word_order(row_auto, word)
                col_actual = module.literal_word_order(col_auto, word)
                assert row_actual is None or row_actual >= row_bound
                assert col_actual is None or col_actual >= col_bound
                counts["literal_rational_word_orders"] += 2
            row_scores[sum(row["exponents"][i] for i in word) - row_bound] += 1
            col_scores[sum(column["exponents"][i] for i in word) - col_bound] += 1
            counts["min_plus_word_comparisons"] += 2
        row_value = module.weighted_value(row_matrix, row_terminal, length)
        col_value = module.weighted_value(col_matrix, col_terminal, length)
        assert row_value == sum(count * s ** degree for degree, count in row_scores.items())
        assert col_value == sum(count * s ** degree for degree, count in col_scores.items())
        assert row_value <= matrix_value(old_row, old_start, old_end, length)
        assert col_value == matrix_value(old_col, old_col_start, old_col_end, length)
        assert row_value <= row_beta ** length and col_value <= col_beta ** length
        for k in range(2 * length + 1):
            row_admitted = sum(count for degree, count in row_scores.items() if degree >= k)
            col_admitted = sum(count for degree, count in col_scores.items() if degree >= -k)
            assert ranks[0][k] <= (2 * length - k + 1) * row_admitted
            assert ranks[1][k] <= (k + 1) * col_admitted
            assert ranks[0][k] <= (2 * length - k + 1) * s ** (-k) * row_beta ** length
            assert ranks[1][k] <= (k + 1) * s ** k * col_beta ** length
            assert ranks[0][k] * ranks[1][k] <= (2 * length + 1) ** 2 * B ** length
            counts["literal_coefficient_rank_inequalities"] += 5
        rank_records.append({"length": length, "row": ranks[0], "column": ranks[1]})

    def strings(matrix):
        return [[str(value) for value in row_] for row_ in matrix]

    return {
        "status": "pass", "field": 2, "maximum_word": maximum_word,
        "local_formula": "P(z)I2+z(1+z)([v=0]E01+[x=y=1,u=v=0]E10)",
        "s": str(s), "B_module": str(B), "B_min_plus_lower": str(old_B_lower),
        "row_beta": str(row_beta), "column_beta": str(col_beta),
        "row_automaton": row_auto.summary(), "column_automaton": col_auto.summary(),
        "row_transfer": strings([[row_matrix[i].get(j, Fraction(0)) for j in range(2)] for i in range(2)]),
        "column_transfer": strings([[col_matrix[i].get(j, Fraction(0)) for j in range(2)] for i in range(2)]),
        "old_row_transfer": strings(old_row), "old_row_positive_lower_potential": list(map(str, old_h)),
        "old_row_rate_lower": str(old_lower), "rank_records": rank_records, **dict(counts),
        "scope": "Strict same-factorization, same-cap word-certificate improvement for one explicit crossing; no factoring exponent improvement.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--maximum-word", type=int, default=3)
    args = parser.parse_args()
    report = check(args.maximum_word)
    report["checker_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    report["module_checker_sha256"] = sha256(Path(module.__file__).read_bytes()).hexdigest()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
