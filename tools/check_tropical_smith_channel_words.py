#!/usr/bin/env python3
"""Exact finite checks for deterministic min-plus Smith-channel counting."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
from itertools import product
import json
from pathlib import Path
from random import Random

import check_degree_sensitive_smith_paths as paths
import check_smith_profile_transfer as basic


def update(costs, matrix):
    return tuple(min(costs[a] + matrix[a][b] for a in range(len(costs)) if costs[a] is not None)
                 for b in range(len(costs)))


def normalize(costs):
    minimum = min(costs)
    return minimum, tuple(cost - minimum for cost in costs)


def build_automaton(costs, initial=(0,), terminal=(0,)):
    dimension = len(costs[0])
    cap = max(value for matrix in costs for row in matrix for value in row)
    first = tuple(0 if a in initial else None for a in range(dimension))
    states = []
    indices = {}

    def get_index(state):
        assert min(state) == 0 and max(state) <= cap
        if state not in indices:
            indices[state] = len(states)
            states.append(state)
        return indices[state]

    starts = []
    for symbol, matrix in enumerate(costs):
        increment, state = normalize(update(first, matrix))
        starts.append((symbol, increment, get_index(state)))
    transitions = []
    cursor = 0
    while cursor < len(states):
        next_states = []
        for symbol, matrix in enumerate(costs):
            increment, state = normalize(update(states[cursor], matrix))
            next_states.append((symbol, increment, get_index(state)))
        transitions.append(next_states)
        cursor += 1
    assert len(states) <= (cap + 1)**dimension - cap**dimension
    return {"states": states, "starts": starts, "transitions": transitions,
            "terminal": [min(state[a] for a in terminal) for state in states],
            "initial_support": initial, "terminal_support": terminal, "cap": cap}


def capped_costs(data, cap):
    return [[[cap if entry is None else min(entry, cap) for entry in row] for row in matrix]
            for matrix in data["deltas"]]


def word_cost(costs, symbols, initial, terminal):
    dimension = len(costs[0])
    current = tuple(0 if a in initial else None for a in range(dimension))
    for symbol in symbols:
        current = update(current, costs[symbol])
    return min(current[a] for a in terminal)


def literal_path_cost(costs, symbols, initial, terminal):
    dimension = len(costs[0])
    return min(sum(costs[symbol][a][b] for symbol, a, b in zip(symbols, trajectory, trajectory[1:]))
               for first in initial for last in terminal
               for middle in product(range(dimension), repeat=len(symbols) - 1)
               for trajectory in [(first, *middle, last)])


def run_word(automaton, symbols):
    first_symbol, increment, state = automaton["starts"][symbols[0]]
    assert first_symbol == symbols[0]
    total = increment
    for symbol in symbols[1:]:
        selected, increment, state = automaton["transitions"][state][symbol]
        assert selected == symbol
        total += increment
    return total + automaton["terminal"][state]


def weighted_automaton(automaton, exponents, s):
    size = len(automaton["states"])
    matrix = [[Fraction(0)] * size for _ in range(size)]
    initial = [Fraction(0)] * size
    for symbol, increment, state in automaton["starts"]:
        initial[state] += s ** (exponents[symbol] - increment)
    for state, transitions in enumerate(automaton["transitions"]):
        for symbol, increment, target in transitions:
            matrix[state][target] += s ** (exponents[symbol] - increment)
    terminal = [s ** (-cost) for cost in automaton["terminal"]]
    return matrix, initial, terminal


def weighted_value(matrix, initial, terminal, length):
    vector = initial[:]
    for _ in range(length - 1):
        vector = [sum(vector[a] * matrix[a][b] for a in range(len(matrix))) for b in range(len(matrix))]
    return sum(a * b for a, b in zip(vector, terminal))


def positive_potential(matrix):
    """Generate candidates numerically in rationals; certify every bound exactly."""
    size = len(matrix)
    vector = [Fraction(1)] * size
    best = None
    for _ in range(13):
        image = [sum(entry * vector[b] for b, entry in enumerate(row)) for row in matrix]
        beta = max(value / coordinate for value, coordinate in zip(image, vector))
        assert beta > 0
        assert all(value <= beta * coordinate for value, coordinate in zip(image, vector))
        if best is None or beta < best[0]:
            best = beta, tuple(vector)
        scale = max(image)
        vector = [Fraction(max(1, round(value / scale * 4096)), 4096) for value in image]
    return best


def nondeterministic_value(costs, exponents, s, length, initial=(0,), terminal=(0,)):
    dimension = len(costs[0])
    matrix = [[sum(s ** (exponent - edge[a][b]) for exponent, edge in zip(exponents, costs))
               for b in range(dimension)] for a in range(dimension)]
    vector = [Fraction(int(a in initial)) for a in range(dimension)]
    for _ in range(length):
        vector = [sum(vector[a] * matrix[a][b] for a in range(dimension)) for b in range(dimension)]
    return sum(vector[a] for a in terminal)


def check_synthetic():
    random = Random(551070)
    stats = Counter()
    for dimension, alphabet, cap in ((2, 3, 2), (3, 2, 2), (3, 3, 1)):
        for _ in range(4):
            costs = [[[random.randrange(cap + 1) for _ in range(dimension)] for _ in range(dimension)]
                     for _ in range(alphabet)]
            initial, terminal = (0, dimension - 1), (1,)
            automaton = build_automaton(costs, initial, terminal)
            exponents = [random.randrange(-2, 3) for _ in range(alphabet)]
            for length in range(1, 5):
                for symbols in product(range(alphabet), repeat=length):
                    expected = literal_path_cost(costs, symbols, initial, terminal)
                    assert word_cost(costs, symbols, initial, terminal) == expected
                    assert run_word(automaton, symbols) == expected
                    stats["literal_synthetic_word_checks"] += 1
                for s in (Fraction(1), Fraction(3, 2), Fraction(3)):
                    matrix, start, end = weighted_automaton(automaton, exponents, s)
                    expected = sum(s ** (sum(exponents[symbol] for symbol in symbols)
                                        - literal_path_cost(costs, symbols, initial, terminal))
                                   for symbols in product(range(alphabet), repeat=length))
                    assert weighted_value(matrix, start, end, length) == expected
                    assert expected <= nondeterministic_value(costs, exponents, s, length, initial, terminal)
                    stats["synthetic_weighted_identities"] += 1
    # Nondeterministic paths multiply while the one-symbol language has
    # exactly one word at every length.
    costs = [[[0, 0], [0, 0]]]
    automaton = build_automaton(costs)
    matrix, start, end = weighted_automaton(automaton, [0], Fraction(2))
    assert weighted_value(matrix, start, end, 6) == 1
    assert nondeterministic_value(costs, [0], Fraction(2), 6) == 32
    # Removed common minima and terminal state offsets are both necessary.
    costs = [[[1, 1], [1, 1]]]
    assert run_word(build_automaton(costs), (0, 0, 0)) == 3
    costs = [[[0, 2], [0, 2]]]
    assert run_word(build_automaton(costs, (0,), (1,)), (0,)) == 2
    # Exact diagonal transitions of costs zero and one, with both states
    # initially reachable and off-diagonal edges absent, have no finite
    # normalized-difference bound.
    diagonal_costs = (0, 0)
    for _ in range(20):
        diagonal_costs = (diagonal_costs[0], diagonal_costs[1] + 1)
    assert normalize(diagonal_costs)[1] == (0, 20)
    return {**dict(stats), "strict_path_multiplicity_control": {"word_count": 1, "path_count": 32},
            "common_increment_control": 3, "terminal_offset_control": 2,
            "uncapped_infinite_edges_can_give_unbounded_normalized_differences": [0, 20]}


def check_family(name, local, degree, maximum_word):
    dimension = len(next(iter(local.values())))
    heights = [(0,)] if dimension == 1 else [(0, 0), (0, 1)]
    factors = [paths.local_path_data(local, cut, side, height, infinity)
               for cut, side, height, infinity in product(range(2), range(2), heights, (False, True))]
    rank_data = {length: basic.full_coefficients(local, degree, length)[0]
                 for length in range(1, maximum_word + 1)}
    stats = Counter()
    options = defaultdict(list)
    for data in factors:
        for cap in (0, 1, 2, 4):
            costs = capped_costs(data, cap)
            automaton = build_automaton(costs)
            stats["automata_constructed"] += 1
            stats["maximum_reachable_states"] = max(stats["maximum_reachable_states"], len(automaton["states"]))
            counts_by_length = {}
            for length in range(1, maximum_word + 1):
                counts = Counter()
                for symbols in product(range(4), repeat=length):
                    expected = literal_path_cost(costs, symbols, (0,), (0,))
                    assert run_word(automaton, symbols) == expected
                    assert word_cost(costs, symbols, (0,), (0,)) == expected
                    counts[sum(data["exponents"][symbol] for symbol in symbols) - expected] += 1
                    stats["literal_field_word_checks"] += 1
                counts_by_length[length] = counts
            for s in map(Fraction, ("1", "3/2", "3", "5")):
                matrix, start, end = weighted_automaton(automaton, data["exponents"], s)
                beta, potential = positive_potential(matrix)
                constant = sum(a * b for a, b in zip(start, potential)) * max(
                    value / coordinate for value, coordinate in zip(end, potential))
                options[s, data["cut"], data["infinity"]].append((beta, data, cap, automaton, potential, constant))
                for length in range(1, maximum_word + 1):
                    weighted = weighted_value(matrix, start, end, length)
                    assert weighted == sum(count * s**score for score, count in counts_by_length[length].items())
                    assert weighted <= nondeterministic_value(costs, data["exponents"], s, length)
                    assert weighted <= constant * beta**(length - 1)
                    for k in range(degree * length + 1):
                        threshold = (k if data["infinity"] else -k) - data["boundary_shift"]
                        admitted = sum(count for score, count in counts_by_length[length].items() if score >= threshold)
                        splits = max(0, length * max(data["exponents"]) - threshold + 1)
                        actual = rank_data[length][data["cut"]][k]
                        assert actual <= splits * admitted
                        assert admitted <= s**(-threshold) * weighted
                        assert actual <= splits * s**(-threshold) * constant * beta**(length - 1)
                        stats["coefficient_rank_checks"] += 1
                        stats["rational_word_and_potential_checks"] += 3
    best = None
    for s, cut in product(map(Fraction, ("1", "3/2", "3", "5")), range(2)):
        high = min(options[s, cut, True], key=lambda item: item[0])
        low = min(options[s, 1 - cut, False], key=lambda item: item[0])
        value = high[0] * low[0]
        if best is None or value < best[0]:
            best = value, s, high, low

    def describe(option):
        beta, data, cap, automaton, potential, constant = option
        return {"beta": str(beta), "local": paths.describe(data), "cap": cap,
                "states": automaton["states"], "starts": automaton["starts"],
                "transitions": automaton["transitions"], "terminal_costs": automaton["terminal"],
                "potential": list(map(str, potential)), "boundary_constant": str(constant)}

    return {"name": name, "degree": degree, "auxiliary_dimension": dimension, **dict(stats),
            "best_sampled_certificate": {"B_det": str(best[0]), "s": str(best[1]),
                                         "infinity": describe(best[2]), "zero": describe(best[3])}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--maximum-word", type=int, default=3)
    args = parser.parse_args()
    cases = [("scalar elementary", *basic.elementary()),
             ("auxiliary elementary", *basic.elementary(2)),
             ("ordered pairing", *basic.pairing((1, 6, 9), 0)),
             ("mixed nilpotent loop", *basic.loop(13, 18))]
    results = [check_family(name, local, degree, args.maximum_word) for name, local, degree in cases]
    assert Fraction(results[0]["best_sampled_certificate"]["B_det"]) == Fraction(256, 27)
    report = {
        "status": "pass", "field": 2, "maximum_word": args.maximum_word,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "path_checker_sha256": hashlib.sha256(Path(paths.__file__).read_bytes()).hexdigest(),
        "cases": results, "synthetic_checks_and_controls": check_synthetic(),
        "scope": "Exact deterministic automaton, finite coefficient-rank and rational-potential checks; no improving seed.",
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
