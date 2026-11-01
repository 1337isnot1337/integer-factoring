#!/usr/bin/env python3
"""Exact checks for scalar-twisted factor gathering and channel-walk bounds."""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import isqrt
from pathlib import Path
import random
import time

from check_commuting_halves_matrix_seeds import identity, madd, mmul, ordered


def scale(a, value, p):
    return tuple(value*x % p for x in a)


def crossing_factors(coefficient, cut, p):
    unit, zero = identity(2), (0,)*4
    factors = []
    for a, b in product(range(2), repeat=2):
        # For both cuts, the left input pairs the external coordinates
        # with the two reversed right coordinates of the crossing.
        left = [coefficient if x == b and y == a else zero
                for x, y in product(range(2), repeat=2)]
        right = [unit if x == a and y == b else zero
                 for x, y in product(range(2), repeat=2)]
        factors.append((left, right))
    return factors


def entry(factors, i, indices, cut, p):
    x, y, u, v = indices
    left, right = (2*x+y, 2*u+v) if cut == "row" else (2*x+u, 2*y+v)
    answer = (0,)*4
    for a, b in factors[i]:
        answer = madd(answer, mmul(a[left], b[right], 2, p), p)
    return answer


def graph_and_twists(factors, p):
    states = [(i, a) for i, piece in enumerate(factors) for a in range(len(piece))]
    twists, graph = {}, {}
    commutator_tests = zero_tests = 0
    for s, t in product(states, repeat=2):
        ls, rs = factors[s[0]][s[1]]
        lt, rt = factors[t[0]][t[1]]
        candidates = []
        for value in range(1, p):
            valid = True
            for right, left in product(rs, lt):
                commutator_tests += 1
                if mmul(right, left, 2, p) != scale(mmul(left, right, 2, p), value, p):
                    valid = False
                    break
            if valid:
                candidates.append(value)
        assert candidates
        twists[s, t] = candidates[0]
        left_zero = right_zero = True
        for a, b in product(ls, lt):
            zero_tests += 1
            left_zero &= not any(mmul(a, b, 2, p))
        for a, b in product(rs, rt):
            zero_tests += 1
            right_zero &= not any(mmul(a, b, 2, p))
        graph[s, t] = not (left_zero or right_zero)
    return states, graph, twists, commutator_tests, zero_tests


def walk_types(states, graph, k):
    count = max(i for i, _ in states)+1
    dynamic = defaultdict(int)
    for state in states:
        kind = [0]*count
        kind[state[0]] = 1
        dynamic[state, tuple(kind)] += 1
    for _ in range(1, k):
        following = defaultdict(int)
        for (state, kind), value in dynamic.items():
            for next_state in states:
                if graph[state, next_state]:
                    next_kind = list(kind)
                    next_kind[next_state[0]] += 1
                    following[next_state, tuple(next_kind)] += value
        dynamic = following
    answer = defaultdict(int)
    for (_, kind), value in dynamic.items():
        answer[kind] += value
    return dict(answer)


def rational_certificate(states, graph, weights, h, beta):
    assert all(value > 0 for value in weights+h+[beta])
    for i, state in enumerate(states):
        actual = sum(weights[next_state[0]]*h[j]
                     for j, next_state in enumerate(states) if graph[state, next_state])
        if actual > beta*h[i]:
            return None
    return sum(weights[state[0]]*h[i] for i, state in enumerate(states))/min(h)


def weighted_count(counts, weights):
    answer = Fraction(0)
    for kind, value in counts.items():
        weight = Fraction(value)
        for exponent, base in zip(kind, weights):
            weight *= base**exponent
        answer += weight
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    rng = random.Random(20261030)
    unit, zero = identity(2), (0,)*4
    # A nontrivial scalar twist in F3: Z X = -X Z.
    x, z = (0, 1, 1, 0), (1, 0, 0, 2)
    plus = ([x, zero, zero, zero], [z, zero, zero, zero])
    minus = ([scale(x, 2, 3), zero, zero, zero], [z, zero, zero, zero])
    twisted = {cut: [[plus], [minus], crossing_factors(unit, cut, 3)] for cut in ("row", "col")}
    # Three noncommuting square-zero matrices sum to I in characteristic two.
    nilpotents = [(0, 1, 0, 0), (0, 0, 1, 0), (1, 1, 1, 1)]
    assert madd(madd(nilpotents[0], nilpotents[1], 2), nilpotents[2], 2) == unit
    assert all(not any(mmul(n, n, 2, 2)) for n in nilpotents)
    assert all(any(mmul(nilpotents[i], nilpotents[j], 2, 2))
               for i, j in product(range(3), repeat=2) if i != j)
    nilseed = {cut: [crossing_factors(n, cut, 2) for n in nilpotents] for cut in ("row", "col")}
    crossing_checks = factor_checks = commutator_tests = zero_tests = 0
    dropped_twist_changes_answer = False
    graph_checks = 0
    seed_summaries = []
    for seed, p, name in [(twisted, 3, "scalar_twist"), (nilseed, 2, "square_zero")]:
        graphs = {}
        for cut in ("row", "col"):
            states, graph, phases, a, b = graph_and_twists(seed[cut], p)
            commutator_tests += a
            zero_tests += b
            graphs[cut] = states, graph, phases
        for indices in product(range(2), repeat=4):
            total = zero
            for i in range(3):
                value = entry(seed["row"], i, indices, "row", p)
                assert value == entry(seed["col"], i, indices, "col", p)
                total = madd(total, value, p)
            assert total == scale(unit, int(indices[0] == indices[3] and indices[1] == indices[2]), p)
            crossing_checks += 1
        for k in (2, 3):
            for sample in range(12):
                external = ([(0, 0, 0, 0)]*k if sample == 0 else
                            [tuple(rng.randrange(2) for _ in range(4)) for _ in range(k)])
                for word in product(range(3), repeat=k):
                    actual = ordered([entry(seed["row"], word[j], external[j], "row", p)
                                      for j in range(k)], 2, p)
                    for cut in ("row", "col"):
                        states, graph, phases = graphs[cut]
                        gathered = untwisted = zero
                        alphabets = [[(piece, a) for a in range(len(seed[cut][piece]))] for piece in word]
                        for path in product(*alphabets):
                            if any(not graph[path[j], path[j+1]] for j in range(k-1)):
                                continue
                            lefts, rights = [], []
                            for j, (piece, channel) in enumerate(path):
                                xx, yy, uu, vv = external[j]
                                l, r = (2*xx+yy, 2*uu+vv) if cut == "row" else (2*xx+uu, 2*yy+vv)
                                lefts.append(seed[cut][piece][channel][0][l])
                                rights.append(seed[cut][piece][channel][1][r])
                            value = mmul(ordered(lefts, 2, p), ordered(rights, 2, p), 2, p)
                            phase = 1
                            for j in range(k):
                                for l in range(j+1, k):
                                    phase = phase*phases[path[j], path[l]] % p
                            gathered = madd(gathered, scale(value, phase, p), p)
                            untwisted = madd(untwisted, value, p)
                        assert gathered == actual
                        dropped_twist_changes_answer |= untwisted != actual
                        factor_checks += 1
        states, graph, phases = graphs["row"]
        for k in range(1, 5):
            counts = walk_types(states, graph, k)
            if k <= 3:
                direct = defaultdict(int)
                for path in product(states, repeat=k):
                    if all(graph[path[j], path[j+1]] for j in range(k-1)):
                        kind = tuple(sum(state[0] == i for state in path) for i in range(3))
                        direct[kind] += 1
                assert counts == dict(direct)
            if name == "square_zero":
                assert sum(counts.values()) == 3*2**(k-1)*4**k
                beta = Fraction(8)
            else:
                assert sum(counts.values()) == 6**k
                beta = Fraction(6)
            constant = rational_certificate(states, graph, [Fraction(1)]*3,
                                            [Fraction(1)]*len(states), beta)
            assert constant is not None
            assert sum(counts.values()) <= constant*beta**(k-1)
            graph_checks += 1
        seed_summaries.append({"name": name, "field": p,
                               "unweighted_rank_sum_base": len(states),
                               "certified_walk_base": int(beta),
                               "channels_per_cut": len(states),
                               "allowed_edges_per_cut": sum(graph.values()),
                               "nontrivial_cross_phases": sum(value != 1 for value in phases.values())})
    assert dropped_twist_changes_answer
    states, graph, _ = graphs["row"]
    assert rational_certificate(states, graph, [Fraction(1)]*3,
                                [Fraction(1)]*len(states), Fraction(7)) is None

    # Different row/column channel profiles test the reciprocal weighting
    # and exact Cauchy equality when the graphs are complete.
    row_states = [(0, 0)]+[(1, i) for i in range(4)]
    col_states = [(0, i) for i in range(4)]+[(1, 0)]
    row_graph = {(s, t): True for s, t in product(row_states, repeat=2)}
    col_graph = {(s, t): True for s, t in product(col_states, repeat=2)}
    weights = [Fraction(2), Fraction(1, 2)]
    reciprocal = [1/value for value in weights]
    assert rational_certificate(row_states, row_graph, weights, [Fraction(1)]*5, Fraction(4)) == 4
    assert rational_certificate(col_states, col_graph, reciprocal, [Fraction(1)]*5, Fraction(4)) == 4
    weighted_checks = 0
    for k in range(1, 6):
        rows = walk_types(row_states, row_graph, k)
        cols = walk_types(col_states, col_graph, k)
        rank_sum = 0
        for kind in rows:
            term = isqrt(rows[kind]*cols[kind])
            assert term*term == rows[kind]*cols[kind]
            rank_sum += term
        pr, pc = weighted_count(rows, weights), weighted_count(cols, reciprocal)
        assert pr == pc == rank_sum == 4**k
        assert rank_sum*rank_sum == pr*pc
        weighted_checks += 1

    output = {"status": "passed", "seed_summaries": seed_summaries,
              "single_copy_crossing_entries": crossing_checks,
              "ordered_gathered_factor_checks": factor_checks,
              "all_input_scalar_commutation_checks": commutator_tests,
              "all_input_zero_product_checks": zero_tests,
              "type_walk_count_checks": graph_checks,
              "reciprocal_weight_exact_cauchy_checks": weighted_checks,
              "incorrect_twist_negative_control": dropped_twist_changes_answer,
              "invalid_rational_supervector_rejected": True,
              "scope": "Both actual toy seeds are far above the improvement threshold; the 12-to-8 base reduction only demonstrates certified ordered nilpotence.",
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "seconds": time.monotonic()-started}
    print(json.dumps(output, indent=2, sort_keys=True))
    if args.output:
        args.output.write_text(json.dumps(output, indent=2, sort_keys=True)+"\n")


if __name__ == "__main__":
    main()
