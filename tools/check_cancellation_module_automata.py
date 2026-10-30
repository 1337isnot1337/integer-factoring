#!/usr/bin/env python3
"""Exact cancellation-aware Smith-channel automata over F_2.

The module transition and its verification are new here.  Rational arithmetic,
the elementary/pairing/loop controls, and the ordered coefficient interface are
adapted from PR #30 at 28105165edb8e9a0427468dc833f3c09bc0a7fb7:
check_degree_sensitive_smith_paths.py and check_smith_profile_transfer.py.
No source checkout or third-party package is needed to run this checker.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random


ZERO = (0, 1)
ONE = (1, 1)


def pmul(a, b):
    answer = 0
    while b:
        bit = b & -b
        answer ^= a << (bit.bit_length() - 1)
        b ^= bit
    return answer


def pdivrem(a, b):
    assert b
    quotient = 0
    while a and a.bit_length() >= b.bit_length():
        shift = a.bit_length() - b.bit_length()
        quotient ^= 1 << shift
        a ^= b << shift
    return quotient, a


def pgcd(a, b):
    while b:
        a, b = b, pdivrem(a, b)[1]
    return a


def rat(n, d=1):
    assert d
    if not n:
        return ZERO
    divisor = pgcd(n, d)
    n, nr = pdivrem(n, divisor)
    d, dr = pdivrem(d, divisor)
    assert nr == dr == 0
    return n, d


def radd(a, b):
    return rat(pmul(a[0], b[1]) ^ pmul(b[0], a[1]), pmul(a[1], b[1]))


def rmul(a, b):
    return rat(pmul(a[0], b[0]), pmul(a[1], b[1]))


def rdiv(a, b):
    assert b[0]
    return rat(pmul(a[0], b[1]), pmul(a[1], b[0]))


def rshift(a, shift):
    return rat(a[0] << shift, a[1]) if shift >= 0 else rat(a[0], a[1] << -shift)


def porder(a):
    assert a
    return (a & -a).bit_length() - 1


def rorder(a):
    return None if a == ZERO else porder(a[0]) - porder(a[1])


@lru_cache(maxsize=None)
def jet(a, precision):
    assert precision >= 1
    if a == ZERO:
        return 0
    assert rorder(a) >= 0 and a[1] & 1
    inverse = 1
    for degree in range(1, precision):
        if pmul(a[1], inverse) >> degree & 1:
            inverse ^= 1 << degree
    mask = (1 << precision) - 1
    assert pmul(a[1], inverse) & mask == 1
    return pmul(a[0], inverse) & mask


def binary_basis(rows):
    """Canonical reduced basis with highest-bit pivots."""
    pivots = {}
    for row in rows:
        for column in sorted(pivots, reverse=True):
            if row >> column & 1:
                row ^= pivots[column]
        if row:
            column = row.bit_length() - 1
            for old in list(pivots):
                if pivots[old] >> column & 1:
                    pivots[old] ^= row
            pivots[column] = row
    return tuple(pivots[column] for column in sorted(pivots, reverse=True))


def pack(vector, precision):
    mask = (1 << precision) - 1
    assert all(0 <= entry <= mask for entry in vector)
    return sum(entry << (precision * j) for j, entry in enumerate(vector))


def unpack(row, dimension, precision):
    mask = (1 << precision) - 1
    return tuple((row >> (precision * j)) & mask for j in range(dimension))


def module_basis(generators, dimension, precision):
    """R-span in R^D, R=F_2[tau]/tau^precision."""
    mask = (1 << precision) - 1
    rows = []
    for generator in generators:
        assert len(generator) == dimension
        for shift in range(precision):
            rows.append(pack(tuple((value << shift) & mask for value in generator), precision))
    result = binary_basis(rows)
    # Check stability under multiplication by tau independently.
    for row in result:
        moved = pack(tuple((value << 1) & mask for value in unpack(row, dimension, precision)), precision)
        assert binary_basis((*result, moved)) == result
    return result


def lift_generators(state, dimension, precision):
    """Generators over F_2[[tau]], including the indispensable full tail."""
    answer = list(unpack(row, dimension, precision) for row in state)
    answer.extend(tuple((1 << precision) if j == a else 0 for j in range(dimension))
                  for a in range(dimension))
    return answer


def polynomial_vector_product(vector, matrix, precision):
    dimension = len(vector)
    mask = (1 << precision) - 1
    answer = [0] * dimension
    for a, b in product(range(dimension), repeat=2):
        answer[b] ^= pmul(vector[a], matrix[a][b]) & mask
    return tuple(answer)


class ModuleAutomaton:
    """Finite over-approximation of normalized ordered row-product modules."""

    def __init__(self, channels, precision, initial=None, terminal=None, maximum_states=10000,
                 free_boundaries=False):
        self.channels = tuple(tuple(tuple(tuple(entry for entry in row) for row in matrix)
                                   for matrix in choices) for choices in channels)
        self.precision = precision
        self.dimension = len(self.channels[0][0])
        self.free_boundaries = free_boundaries
        assert precision >= 1 and self.dimension >= 1
        self.initial = tuple(initial or (ONE,) + (ZERO,) * (self.dimension - 1))
        self.terminal = tuple(terminal or (ONE,) + (ZERO,) * (self.dimension - 1))
        self.initial_rows = (tuple(tuple(ONE if j == a else ZERO for j in range(self.dimension))
                                   for a in range(self.dimension)) if free_boundaries else (self.initial,))
        assert any(value != ZERO for value in self.terminal)
        self.deltas = []
        self.jets = []
        for choices in self.channels:
            orders = [rorder(entry) for matrix in choices for row in matrix for entry in row
                      if entry != ZERO]
            assert orders and min(orders) >= 0
            delta = min(orders)
            self.deltas.append(delta)
            available = 2 * precision + delta
            self.jets.append(tuple(tuple(tuple(jet(entry, available) for entry in row)
                                         for row in matrix) for matrix in choices))
        start = module_basis([tuple(jet(entry, precision) for entry in row) for row in self.initial_rows],
                             self.dimension, precision)
        self.states = [start]
        self.indices = {start: 0}
        self.transitions = []
        cursor = 0
        while cursor < len(self.states):
            following = []
            for symbol in range(len(channels)):
                increment, state = self.advance(self.states[cursor], symbol)
                if state not in self.indices:
                    assert len(self.states) < maximum_states, "state budget reached, not an impossibility certificate"
                    self.indices[state] = len(self.states)
                    self.states.append(state)
                following.append((increment, self.indices[state]))
            self.transitions.append(tuple(following))
            cursor += 1
        self.final_costs = tuple(self.terminal_cost(state) for state in self.states)

    def advance(self, state, symbol):
        precision = self.precision
        available = 2 * precision + self.deltas[symbol]
        images = [polynomial_vector_product(generator, matrix, available)
                  for generator in lift_generators(state, self.dimension, precision)
                  for matrix in self.jets[symbol]]
        increment = min(porder(entry) for image in images for entry in image if entry)
        assert increment <= precision + self.deltas[symbol]
        assert available >= precision + increment
        mask = (1 << precision) - 1
        normalized = [tuple((entry >> increment) & mask for entry in image) for image in images]
        state = module_basis(normalized, self.dimension, precision)
        assert state and any(entry & 1 for row in state
                             for entry in unpack(row, self.dimension, precision))
        return increment, state

    def terminal_cost(self, state):
        if self.free_boundaries:
            values = [porder(entry) for generator in lift_generators(state, self.dimension, self.precision)
                      for entry in generator if entry]
            return min(values)
        delta = min(rorder(entry) for entry in self.terminal if entry != ZERO)
        precision = self.precision + delta + 1
        terminal = tuple(jet(entry, precision) for entry in self.terminal)
        values = []
        mask = (1 << precision) - 1
        for generator in lift_generators(state, self.dimension, self.precision):
            value = 0
            for left, right in zip(generator, terminal):
                value ^= pmul(left, right) & mask
            if value:
                values.append(porder(value))
        assert values and min(values) <= self.precision + delta
        return min(values)

    def word_bound(self, word):
        state, total = 0, 0
        for symbol in word:
            increment, state = self.transitions[state][symbol]
            total += increment
        return total + self.final_costs[state]

    def weighted_matrix(self, exponents, s):
        matrix = []
        for transitions in self.transitions:
            row = Counter()
            for exponent, (increment, target) in zip(exponents, transitions):
                row[target] += s ** (exponent - increment)
            matrix.append(dict(row))
        terminal = tuple(s ** (-value) for value in self.final_costs)
        return matrix, terminal

    def summary(self):
        return {"precision": self.precision, "states": [list(state) for state in self.states],
                "transitions": [[list(edge) for edge in row] for row in self.transitions],
                "terminal_costs": list(self.final_costs), "dimension": self.dimension,
                "free_boundaries": self.free_boundaries}


def weighted_value(matrix, terminal, length):
    vector = {0: Fraction(1)}
    for _ in range(length):
        following = Counter()
        for state, coefficient in vector.items():
            for target, weight in matrix[state].items():
                following[target] += coefficient * weight
        vector = dict(following)
    return sum(value * terminal[state] for state, value in vector.items())


def potential(matrix, rounds=32):
    vector = tuple(Fraction(1) for _ in matrix)
    best = None
    for _ in range(rounds):
        image = tuple(sum(weight * vector[target] for target, weight in row.items()) for row in matrix)
        beta = max(value / entry for value, entry in zip(image, vector))
        assert beta > 0
        assert all(value <= beta * entry for value, entry in zip(image, vector))
        if best is None or beta < best[0]:
            best = beta, vector
        scale = max(image)
        vector = tuple(Fraction(max(1, round(value / scale * 65536)), 65536) for value in image)
    return best


def minplus_bound(automaton, word):
    precision, dimension = automaton.precision, automaton.dimension
    costs = ([0] * dimension if automaton.free_boundaries else
             [min(precision, rorder(entry)) if entry != ZERO else precision
              for entry in automaton.initial])
    for symbol in word:
        entry_costs = [[precision] * dimension for _ in range(dimension)]
        for a, b in product(range(dimension), repeat=2):
            available = [rorder(matrix[a][b]) for matrix in automaton.channels[symbol]
                         if matrix[a][b] != ZERO]
            if available:
                entry_costs[a][b] = min(precision, min(available))
        costs = [min(costs[a] + entry_costs[a][b] for a in range(dimension))
                 for b in range(dimension)]
    return (min(costs) if automaton.free_boundaries else
            min(cost + rorder(entry) for cost, entry in zip(costs, automaton.terminal) if entry != ZERO))


def rational_vector_product(vector, matrix):
    dimension = len(vector)
    answer = [ZERO] * dimension
    for a, b in product(range(dimension), repeat=2):
        answer[b] = radd(answer[b], rmul(vector[a], matrix[a][b]))
    return tuple(answer)


def literal_word_order(automaton, word):
    """Independent exact rational products, enumerating physical half-labels."""
    nonzero_orders = []
    for initial in automaton.initial_rows:
        for choices in product(*(range(len(automaton.channels[symbol])) for symbol in word)):
            vector = initial
            for symbol, physical in zip(word, choices):
                vector = rational_vector_product(vector, automaton.channels[symbol][physical])
            if automaton.free_boundaries:
                nonzero_orders.extend(rorder(entry) for entry in vector if entry != ZERO)
            else:
                result = ZERO
                for left, right in zip(vector, automaton.terminal):
                    result = radd(result, rmul(left, right))
                if result != ZERO:
                    nonzero_orders.append(rorder(result))
    return min(nonzero_orders) if nonzero_orders else None


def schur_factors(matrix):
    """Exact rational local rank-one factors, with complete reconstruction."""
    work = [[rat(entry) for entry in row] for row in matrix]
    factors = []
    for _ in range(len(matrix[0])):
        choices = [(rorder(value), a, b) for a, row in enumerate(work)
                   for b, value in enumerate(row) if value != ZERO]
        assert choices
        order, a, b = min(choices)
        pivot = work[a][b]
        left = [rshift(row[b], -order) for row in work]
        right = [rdiv(value, pivot) for value in work[a]]
        assert all(value == ZERO or rorder(value) >= 0 for value in left + right)
        factors.append((order, left, right))
        for i, j in product(range(len(work)), range(len(work[0]))):
            work[i][j] = radd(work[i][j], rshift(rmul(left[i], right[j]), order))
    assert all(value == ZERO for row in work for value in row)
    for i, j in product(range(len(matrix)), range(len(matrix[0]))):
        value = ZERO
        for order, left, right in factors:
            value = radd(value, rshift(rmul(left[i], right[j]), order))
        assert value == rat(matrix[i][j])
    return factors


def local_factors(local, cut, side, heights, infinity):
    dimension = len(heights)
    entries = [[set() for _ in range(4)] for _ in range(4 * dimension * dimension)]
    direction = -1 if infinity else 1
    for (x, y, u, v), value in local.items():
        left, right = (2 * x + y, 2 * u + v) if cut == 0 else (2 * x + u, 2 * y + v)
        for a, b in product(range(dimension), repeat=2):
            physical, other = (left, right) if side == 0 else (right, left)
            entries[(physical * dimension + a) * dimension + b][other] = {
                direction * k + heights[a] - heights[b]
                for k in range(value[a][b].bit_length()) if value[a][b] >> k & 1}
    common = min(k for row in entries for entry in row for k in entry)
    matrix = [[sum(1 << (k - common) for k in entry) for entry in row] for row in entries]
    channels, exponents = [], []
    for order, left, right in schur_factors(matrix):
        scalar_order = min(rorder(value) for value in right if value != ZERO)
        exponents.append(-common - order - scalar_order)
        channels.append(tuple(tuple(tuple(left[(physical * dimension + a) * dimension + b]
                                          for b in range(dimension)) for a in range(dimension))
                              for physical in range(4)))
    assert heights[0] == 0  # The control boundaries stay e_0 on both sides.
    return {"channels": tuple(channels), "exponents": tuple(exponents), "cut": cut,
            "side": side, "heights": tuple(heights), "infinity": infinity}


def elementary(dimension=1):
    local = {}
    for x, y, u, v in product(range(2), repeat=4):
        identity, crossing = int(x == u and y == v), int(x == v and y == u)
        value = identity | ((identity ^ crossing) << 1)
        local[x, y, u, v] = [[value if a == b else 0 for b in range(dimension)]
                             for a in range(dimension)]
    return local, 1


def pairing():
    local = {}
    for x, y, u, v in product(range(2), repeat=4):
        R, C, E = int(x == y and u == v), int(x == u and y == v), int(x == v and y == u)
        value = [[0] * 2 for _ in range(2)]
        for selected, power, matrix in zip((R, R ^ C ^ E, C), (0, 0, 1), (1, 6, 9)):
            if selected:
                for a, b in product(range(2), repeat=2):
                    value[a][b] ^= (int(a == b) ^ (3 if matrix >> (2 * a + b) & 1 else 0)) << power
        local[x, y, u, v] = value
    return local, 2


def mixed_loop():
    local, _ = elementary(2)
    first, second = 13, 18
    for (x, y, u, v), value in local.items():
        for a, b in product(range(2), repeat=2):
            row, column, swapped = 2 * (2 * x + y) + a, 2 * (2 * u + v) + b, 2 * (2 * v + u) + b
            if first >> row & 1:
                if second >> column & 1:
                    value[a][b] ^= 6
                if ((second >> column) ^ (second >> swapped)) & 1:
                    value[a][b] ^= 12
    return local, 3


def coefficient_ranks(local, degree, length):
    labels, dimension = 1 << length, len(next(iter(local.values())))
    matrices = [[[0] * (labels * labels) for _ in range(degree * length + 1)] for _ in range(2)]
    endpoints = 0
    for x, y, u, v in product(range(labels), repeat=4):
        state = tuple(int(j == 0) for j in range(dimension))
        for position in range(length):
            physical = tuple(value >> position & 1 for value in (x, y, u, v))
            following = [0] * dimension
            for a, b in product(range(dimension), repeat=2):
                following[b] ^= pmul(state[a], local[physical][a][b])
            state = tuple(following)
        value = state[0]
        assert value.bit_count() % 2 == int(x == v and y == u)
        endpoints += 1
        for k in range(degree * length + 1):
            if value >> k & 1:
                matrices[0][k][x * labels + y] |= 1 << (u * labels + v)
                matrices[1][k][x * labels + u] |= 1 << (y * labels + v)
    return [[len(binary_basis(matrix)) for matrix in cut] for cut in matrices], endpoints


def verify_synthetic():
    counts = Counter()
    # Dense entry support hides U^2=0; min-plus sees a zero-cost path forever.
    nilpotent = ((ONE, ONE), (ONE, ONE))
    strict = []
    for precision in (1, 2, 3):
        automaton = ModuleAutomaton(((nilpotent,),), precision)
        for length in range(1, 9):
            word = (0,) * length
            assert automaton.word_bound(word) == (length - 1) * precision
            assert minplus_bound(automaton, word) == 0
            actual = literal_word_order(automaton, word)
            assert actual == (0 if length == 1 else None)
            counts["nilpotent_cancellation_controls"] += 1
        matrix, end = automaton.weighted_matrix((0,), Fraction(2))
        assert weighted_value(matrix, end, 8) == Fraction(1, 2 ** (7 * precision))
        strict.append({"precision": precision, "module_word_weight_at_length_8": str(weighted_value(matrix, end, 8)),
                       "min_plus_word_weight_at_length_8": "1", "automaton": automaton.summary()})

    random = Random(680234)
    for case in range(4):
        channels = tuple(tuple(tuple(tuple(rat(random.randrange(8), random.choice((1, 3, 5)))
                                           for _ in range(2)) for _ in range(2))
                               for _ in range(2)) for _ in range(2))
        previous_bounds = {}
        for precision in (1, 2, 3):
            automaton = ModuleAutomaton(channels, precision)
            counts["random_automata"] += 1
            counts["maximum_random_states"] = max(counts["maximum_random_states"], len(automaton.states))
            for length in range(1, 4):
                for word in product(range(2), repeat=length):
                    lower = automaton.word_bound(word)
                    assert lower >= minplus_bound(automaton, word)
                    actual = literal_word_order(automaton, word)
                    assert actual is None or actual >= lower
                    assert min(lower, precision) == (precision if actual is None else min(actual, precision))
                    assert lower >= previous_bounds.get(word, 0)
                    previous_bounds[word] = lower
                    counts["literal_rational_word_checks"] += 1
                    counts["precision_monotonicity_and_exactness_checks"] += 2

        # Free auxiliary boundaries give a submultiplicative word pressure.
        previous_counts = {}
        exponents, s = (1, -1), Fraction(3, 2)
        for precision in (1, 2, 3):
            automaton = ModuleAutomaton(channels, precision, free_boundaries=True)
            transfer, terminal = automaton.weighted_matrix(exponents, s)
            exact_counts, approximate_counts = {}, {}
            for length in range(1, 4):
                exact = Fraction(0)
                approximate = Fraction(0)
                for word in product(range(2), repeat=length):
                    actual = literal_word_order(automaton, word)
                    lower = automaton.word_bound(word)
                    exponent = sum(exponents[symbol] for symbol in word)
                    assert min(lower, precision) == (precision if actual is None else min(actual, precision))
                    assert lower >= minplus_bound(automaton, word)
                    exact += Fraction(0) if actual is None else s ** (exponent - actual)
                    approximate += s ** (exponent - lower)
                    counts["literal_free_boundary_word_checks"] += 1
                assert approximate == weighted_value(transfer, terminal, length)
                assert exact <= approximate <= previous_counts.get(length, approximate)
                previous_counts[length] = approximate
                exact_counts[length], approximate_counts[length] = exact, approximate
            for a, b in ((1, 1), (1, 2), (2, 1)):
                assert exact_counts[a + b] <= exact_counts[a] * exact_counts[b]
                assert approximate_counts[a + b] <= approximate_counts[a] * approximate_counts[b]
                counts["free_boundary_submultiplicativity_checks"] += 2

    # The full tail is mandatory.  Both initial rows agree modulo tau^2,
    # but multiplication and division reveal a hidden high coefficient.
    vector = (ONE, rat(1 << 2))
    hidden = ((ZERO, ZERO), (ONE, ZERO))
    assert jet(vector[1], 2) == 0
    assert rational_vector_product((ONE, ZERO), hidden) == (ZERO, ZERO)
    assert rational_vector_product(vector, hidden) == (rat(1 << 2), ZERO)
    tail_control = ModuleAutomaton(((hidden,),), 2, initial=vector)
    assert tail_control.word_bound((0,)) == literal_word_order(tail_control, (0,)) == 2
    return {**dict(counts), "strict_cancellation_examples": strict,
            "discarded_high_coefficients_reappear_after_division": True}


def check_family(name, local, degree, maximum_word, precisions):
    dimension = len(next(iter(local.values())))
    heights = [(0,)] if dimension == 1 else [(0, 0), (0, 1), (0, -1)]
    ranks = {}
    counts = Counter()
    for length in range(1, maximum_word + 1):
        ranks[length], endpoints = coefficient_ranks(local, degree, length)
        counts["ordered_endpoint_entries"] += endpoints
    candidates = {}
    for cut, side, height, infinity in product(range(2), range(2), heights, (False, True)):
        factors = local_factors(local, cut, side, height, infinity)
        for precision in precisions:
            automaton = ModuleAutomaton(factors["channels"], precision)
            counts["automata_constructed"] += 1
            counts["maximum_reachable_states"] = max(counts["maximum_reachable_states"], len(automaton.states))
            distributions = {}
            for length in range(1, maximum_word + 1):
                scores = Counter()
                for word in product(range(4), repeat=length):
                    lower = automaton.word_bound(word)
                    assert lower >= minplus_bound(automaton, word)
                    if length <= 2:
                        actual = literal_word_order(automaton, word)
                        assert actual is None or actual >= lower
                        counts["literal_channel_product_checks"] += 1
                    scores[sum(factors["exponents"][symbol] for symbol in word) - lower] += 1
                    counts["termwise_min_plus_comparisons"] += 1
                distributions[length] = scores
            for s in map(Fraction, ("1", "3/2", "2", "3", "5")):
                matrix, terminal = automaton.weighted_matrix(factors["exponents"], s)
                beta, h = potential(matrix)
                constant = h[0] * max(value / bound for value, bound in zip(terminal, h))
                key = (s, cut, infinity)
                option = (beta, precision, factors, automaton, h, constant)
                if key not in candidates or beta < candidates[key][0]:
                    candidates[key] = option
                for length, scores in distributions.items():
                    weighted = weighted_value(matrix, terminal, length)
                    assert weighted == sum(count * s ** score for score, count in scores.items())
                    assert weighted <= constant * beta ** length
                    for k in range(degree * length + 1):
                        threshold = k if infinity else -k
                        admitted = sum(count for score, count in scores.items() if score >= threshold)
                        splits = max(0, length * max(factors["exponents"]) - threshold + 1)
                        actual = ranks[length][cut][k]
                        assert actual <= splits * admitted, (name, cut, side, height, infinity, precision, length, k)
                        assert admitted <= s ** (-threshold) * weighted
                        counts["coefficient_rank_checks"] += 1
                        counts["rational_weighted_checks"] += 2
    best = None
    for s, cut in product(map(Fraction, ("1", "3/2", "2", "3", "5")), range(2)):
        high, low = candidates[s, cut, True], candidates[s, 1 - cut, False]
        value = high[0] * low[0]
        if best is None or value < best[0]:
            best = value, s, high, low

    def describe(option):
        beta, precision, factors, automaton, h, constant = option
        return {"beta": str(beta), "precision": precision,
                "factors": {key: factors[key] for key in ("exponents", "cut", "side", "heights", "infinity")},
                "potential": list(map(str, h)), "boundary_constant": str(constant),
                "automaton": automaton.summary()}

    return {"name": name, "degree": degree, "auxiliary_dimension": dimension, **dict(counts),
            "best_sampled_certificate": {"B_module": str(best[0]), "s": str(best[1]),
                                         "infinity": describe(best[2]), "zero": describe(best[3])}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--maximum-word", type=int, default=2)
    parser.add_argument("--precisions", default="1,2,3")
    args = parser.parse_args()
    precisions = tuple(map(int, args.precisions.split(",")))
    assert args.maximum_word >= 1 and all(value >= 1 for value in precisions)
    cases = [("scalar elementary", *elementary()), ("auxiliary elementary", *elementary(2)),
             ("ordered pairing", *pairing()), ("mixed nilpotent loop", *mixed_loop())]
    results = [check_family(name, local, degree, args.maximum_word, precisions)
               for name, local, degree in cases]
    assert Fraction(results[0]["best_sampled_certificate"]["B_module"]) == Fraction(256, 27)
    report = {"status": "pass", "field": 2, "maximum_word": args.maximum_word,
              "precisions": precisions, "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "cases": results, "synthetic_controls": verify_synthetic(),
              "scope": "Exact finite module automata, literal rational-product and coefficient-rank checks; no improved factoring exponent."}
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
