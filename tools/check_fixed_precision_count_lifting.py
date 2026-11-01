#!/usr/bin/env python3
"""Exact checks for fixed-precision counting and local tensor lifts.

All arithmetic is exact and uses the Python standard library.  The local
Smith decompositions are recomputed by power-series elimination.  Word
tensor coefficients are then compared with literal nilpotent expansions.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
from itertools import combinations, product
import json
from math import comb, gcd
from pathlib import Path

from check_unitary_orbit_recovery import (
    compare_trace,
    factor_reference,
    perfect_power,
    recover_nonzero,
)


LABELS = tuple(product(range(2), repeat=4))


def check_digits(limit: int) -> dict[str, int]:
    stats: Counter[str] = Counter()
    for prime, precision in ((2, 5), (3, 4), (5, 3), (7, 3)):
        for count in range(limit + 1):
            digits = [comb(count, prime**j) % prime if prime**j <= count else 0
                      for j in range(precision)]
            assert digits == [(count // prime**j) % prime for j in range(precision)]
            assert sum(digit * prime**j for j, digit in enumerate(digits)) == count % prime**precision
            stats["count_digit_instances"] += 1
            stats["individual_digits"] += precision
    # Enumerate all finite acceptance sets, then count sorted tuples by
    # literal tuple enumeration rather than binomial coefficients.
    for mask in range(1 << 6):
        accepted = tuple(i for i in range(6) if (mask >> i) & 1)
        for length in range(1, 5):
            found = 0
            for values in product(accepted, repeat=length):
                states = [compare_trace(left, right, 3)
                          for left, right in zip(values, values[1:])]
                assert all(trace[0] == 0 and len(trace) == 4 for _, trace in states)
                if all(state == -1 for state, _ in states):
                    found += 1
            assert found == (comb(len(accepted), length) if length <= len(accepted) else 0)
            stats["sorted_tuple_instances"] += 1
    # Exhaustive traces: the scan has exactly one possible next state.
    for left, right in product(range(16), repeat=2):
        result, trace = compare_trace(left, right, 4)
        assert result == (left > right) - (left < right)
        for index in range(4):
            a, b = (left >> (3 - index)) & 1, (right >> (3 - index)) & 1
            allowed = [state for state in (-1, 0, 1)
                       if state == (trace[index] if trace[index] else (a > b) - (a < b))]
            assert allowed == [trace[index + 1]]
        stats["unique_comparison_traces"] += 1
    return dict(stats)


def check_roots(limit: int) -> dict[str, int]:
    roots = {}
    stats: Counter[str] = Counter()
    for value in range(3, limit + 1, 2):
        roots[value] = tuple(z for z in range(value) if (z * z - 1) % value == 0)
        omega = len(set(factor_reference(value)))
        assert len(roots[value]) == 2**omega
        reduced = tuple((z,) for z in roots[value] if z not in (1, value - 1))
        assert len(reduced) % 4 == (0 if omega == 1 else 2)
        lifted = (len(reduced) % 2) + 2 * (comb(len(reduced), 2) % 2)
        assert lifted == len(reduced) % 4
        stats["odd_root_inputs"] += 1
        stats["residues_examined"] += value

    for value in range(2, limit + 1):
        output = []
        current = value
        while current % 2 == 0:
            output.append(2)
            current //= 2
        pending = [(current, 1)] if current > 1 else []
        while pending:
            current, multiplicity = pending.pop()
            power = perfect_power(current)
            if power is not None:
                base, exponent = power
                pending.append((base, multiplicity * exponent))
                continue
            witnesses = tuple((z,) for z in roots[current] if z not in (1, current - 1))
            if len(witnesses) % 4 == 0:
                output.extend([current] * multiplicity)
            else:
                (z,), queries = recover_nonzero(witnesses, current.bit_length(), 1, 4)
                divisor = gcd(current, z - 1)
                assert 1 < divisor < current
                pending.extend(((divisor, multiplicity), (current // divisor, multiplicity)))
                stats["prefix_queries"] += queries
        assert sorted(output) == factor_reference(value)
        stats["complete_factorizations"] += 1
    return dict(stats)


def field_factor(matrix: list[list[int]], prime: int):
    """Exact sum-of-outer-products elimination over the prime field."""
    work = [[entry % prime for entry in row] for row in matrix]
    factors = []
    while True:
        pivot = next(((i, j) for i in range(len(work)) for j in range(len(work[0]))
                      if work[i][j]), None)
        if pivot is None:
            break
        row, column = pivot
        inverse = pow(work[row][column], -1, prime)
        left = [work[i][column] * inverse % prime for i in range(len(work))]
        right = work[row][:]
        factors.append((left, right))
        work = [[(entry - left[i] * right[j]) % prime for j, entry in enumerate(values)]
                for i, values in enumerate(work)]
    assert all(sum(left[i] * right[j] for left, right in factors) % prime == matrix[i][j] % prime
               for i in range(len(matrix)) for j in range(len(matrix[0])))
    return factors


def indices(label, cut):
    x, y, u, v = label
    return (2 * x + y, 2 * u + v) if cut == "row" else (2 * x + u, 2 * y + v)


def scalar_seed_replay(maximum_word: int) -> dict[str, object]:
    stats: Counter[str] = Counter()
    profiles = []
    for prime, precision in ((2, 2), (2, 3), (3, 2), (5, 2)):
        modulus = prime**precision
        crossing = [int(x == v and y == u) for x, y, u, v in LABELS]
        identity = [int(x == u and y == v) for x, y, u, v in LABELS]
        pieces = [identity, [(e - i) % prime for e, i in zip(crossing, identity)]]
        endpoint = [(e - a - b) // prime for e, a, b in zip(crossing, *pieces)]
        decompositions = {}
        ranks = {}
        for cut in ("row", "column"):
            by_piece = []
            for piece in pieces:
                matrix = [[0] * 4 for _ in range(4)]
                for label, entry in zip(LABELS, piece):
                    i, j = indices(label, cut)
                    matrix[i][j] = entry
                factors = field_factor(matrix, prime)
                ranks.setdefault(cut, []).append(len(factors))
                good = []
                residual = []
                for label, entry in zip(LABELS, piece):
                    i, j = indices(label, cut)
                    lifted_product = sum(left[i] * right[j] for left, right in factors)
                    assert (entry - lifted_product) % prime == 0
                    good.append(lifted_product % modulus)
                    residual.append(((entry - lifted_product) // prime) % modulus)
                by_piece.append((good, residual))
            decompositions[cut] = by_piece
        profiles.append({"prime": prime, "precision": precision, "row": ranks["row"], "column": ranks["column"]})

        for word_length in range(1, maximum_word + 1):
            for labels in product(range(16), repeat=word_length):
                direct = {(0, 0, 0): 1}
                for label in labels:
                    updated = defaultdict(int)
                    for counts, coefficient in direct.items():
                        for symbol, value in enumerate((pieces[0][label], pieces[1][label], prime * endpoint[label])):
                            new = list(counts)
                            new[symbol] += 1
                            if new[2] < precision:
                                updated[tuple(new)] = (updated[tuple(new)] + coefficient * value) % modulus
                    direct = dict(updated)
                assert sum(direct.values()) % modulus == product_value(crossing[label] for label in labels) % modulus
                for cut in ("row", "column"):
                    expanded = {((0, 0, 0), 0): 1}
                    for label in labels:
                        updated = defaultdict(int)
                        options = [(2, 1, prime * endpoint[label])]
                        for symbol, (good, residual) in enumerate(decompositions[cut]):
                            options.extend(((symbol, 0, good[label]), (symbol, 1, prime * residual[label])))
                        for (counts, defects), coefficient in expanded.items():
                            for symbol, extra, value in options:
                                if defects + extra >= precision:
                                    continue
                                new = list(counts)
                                new[symbol] += 1
                                key = (tuple(new), defects + extra)
                                updated[key] = (updated[key] + coefficient * value) % modulus
                        expanded = dict(updated)
                    regrouped = defaultdict(int)
                    for (counts, _), coefficient in expanded.items():
                        regrouped[counts] = (regrouped[counts] + coefficient) % modulus
                    assert all(regrouped[counts] == coefficient for counts, coefficient in direct.items())
                    assert all(counts in direct or coefficient == 0 for counts, coefficient in regrouped.items())
                    stats["literal_seed_word_cut_instances"] += 1
                    stats["grouped_seed_coefficients"] += len(direct)
    return {**dict(stats), "local_profiles": profiles}


def product_value(values):
    result = 1
    for value in values:
        result *= value
    return result


def series_mul(left, right, modulus, precision):
    result = [0] * precision
    for i, a in enumerate(left[:precision]):
        if a:
            for j, b in enumerate(right[:precision - i]):
                result[i + j] = (result[i + j] + a * b) % modulus
    return result


def series_inverse(series, prime, precision):
    inverse = [0] * precision
    inverse[0] = pow(series[0], -1, prime)
    for degree in range(1, precision):
        inverse[degree] = -inverse[0] * sum(series[j] * inverse[degree - j]
                                         for j in range(1, min(degree + 1, len(series)))) % prime
    assert series_mul(series, inverse, prime, precision) == [1] + [0] * (precision - 1)
    return inverse


def series_factors(matrix, prime, precision):
    work = [[series[:] for series in row] for row in matrix]
    result = []
    while True:
        candidates = [(next((k for k, value in enumerate(series) if value), precision), i, j)
                      for i, row in enumerate(work) for j, series in enumerate(row)]
        order, row, column = min(candidates)
        if order == precision:
            break
        unit = work[row][column][order:]
        inverse = series_inverse(unit, prime, precision - order)
        left = [series_mul(work[i][column][order:], inverse, prime, precision - order) + [0] * order
                for i in range(len(work))]
        right = [series[order:] + [0] * order for series in work[row]]
        result.append((order, left, right))
        for i in range(len(work)):
            for j in range(len(work[0])):
                contribution = [0] * order + series_mul(left[i], right[j], prime, precision - order)
                work[i][j] = [(a - b) % prime for a, b in zip(work[i][j], contribution)]
    for i, row in enumerate(matrix):
        for j, series in enumerate(row):
            reconstructed = [0] * precision
            for order, left, right in result:
                contribution = [0] * order + series_mul(left[i], right[j], prime, precision - order)
                reconstructed = [(a + b) % prime for a, b in zip(reconstructed, contribution)]
            assert reconstructed == series
    return result


def smith_local_lift(cut, infinity, ring_precision, series_precision):
    prime = 2
    modulus = prime**ring_precision
    actual = [[None] * 4 for _ in range(4)]
    for label in LABELS:
        x, y, u, v = label
        crossing = int(x == v and y == u)
        identity = int(x == u and y == v)
        slope = crossing ^ identity
        # Correct the lifted constant coefficient, preserving reduction.
        constant = (crossing - slope) % modulus
        assert (constant + slope) % modulus == crossing
        assert constant % prime == identity
        coefficients = [constant, slope]
        if infinity:
            coefficients.reverse()
        coefficients += [0] * (series_precision - 2)
        i, j = indices(label, cut)
        actual[i][j] = coefficients
    field = [[[coefficient % prime for coefficient in series] for series in row] for row in actual]
    factors = series_factors(field, prime, series_precision)
    orders = sorted(order for order, _, _ in factors)
    expected = [0, 1, 1, 2] if infinity or cut == "column" else [0, 0, 0, 0]
    assert orders == expected, (cut, infinity, orders)
    good = [[None] * 4 for _ in range(4)]
    error = [[None] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            series = [0] * series_precision
            for order, left, right in factors:
                contribution = [0] * order + series_mul(left[i], right[j], modulus, series_precision - order)
                series = [(a + b) % modulus for a, b in zip(series, contribution)]
            good[i][j] = series
            assert all((a - b) % prime == 0 for a, b in zip(actual[i][j], series))
            error[i][j] = [((a - b) // prime) % modulus for a, b in zip(actual[i][j], series)]
            assert [(a + prime * b) % modulus for a, b in zip(series, error[i][j])] == actual[i][j]
    return actual, good, error, orders


def smith_replay(maximum_word: int) -> dict[str, object]:
    stats: Counter[str] = Counter()
    spectra = []
    # The certificate uses row at infinity and column at zero.
    for ring_precision in (2, 3):
        modulus = 2**ring_precision
        for cut, infinity in (("row", True), ("column", False)):
            for word_length in range(1, maximum_word + 1):
                precision = max(4, word_length + 1)
                actual, good, error, orders = smith_local_lift(cut, infinity, ring_precision, precision)
                spectra.append({"ring_precision": ring_precision, "cut": cut,
                                "at_infinity": infinity, "word_length": word_length, "orders": orders})
                for labels in product(range(16), repeat=word_length):
                    target = [1] + [0] * word_length
                    expanded = [[1] + [0] * word_length] + [[0] * (word_length + 1)
                                                                         for _ in range(ring_precision - 1)]
                    for label_id in labels:
                        i, j = indices(LABELS[label_id], cut)
                        target = series_mul(target, actual[i][j], modulus, word_length + 1)
                        next_expanded = []
                        for defects in range(ring_precision):
                            value = series_mul(expanded[defects], good[i][j], modulus, word_length + 1)
                            if defects:
                                added = series_mul(expanded[defects - 1], error[i][j], modulus, word_length + 1)
                                value = [(a + 2 * b) % modulus for a, b in zip(value, added)]
                            next_expanded.append(value)
                        expanded = next_expanded
                    reconstructed = [sum(series[degree] for series in expanded) % modulus
                                     for degree in range(word_length + 1)]
                    assert reconstructed == target
                    stats["literal_smith_word_cut_instances"] += 1
                    stats["smith_coefficient_comparisons"] += word_length + 1
    return {**dict(stats), "recomputed_profiles": spectra}


def negative_controls():
    # A zero field matrix can have an arbitrarily large nonzero diagonal
    # lift.  Locality, rather than an arbitrary global rank assertion, is
    # essential to the theorem.
    lifted_zero = [[2 * int(i == j) for j in range(4)] for i in range(4)]
    assert all(value % 2 == 0 for row in lifted_zero for value in row)
    assert len({tuple(sum(lifted_zero[i][j] * vector[j] for j in range(4)) % 4 for i in range(4))
                for vector in product(range(4), repeat=4)}) == 16
    # Ordered distinct pairs count twice the desired binomial witness set.
    assert 3 * 2 % 2 == 0 and comb(3, 2) % 2 == 1
    # An uncorrected binary pencil does not have the exact Z/4 endpoint.
    label = (0, 1, 0, 1)
    x, y, u, v = label
    crossing, identity = int(x == v and y == u), int(x == u and y == v)
    assert crossing == 0 and (identity + (crossing ^ identity)) % 4 == 2
    # At precision three, two-error terms cannot be dropped.
    assert (1 + 2) ** 2 % 8 == 1 and (1 + 2 + 2) % 8 == 5
    return {
        "arbitrary_zero_matrix_lift": lifted_zero,
        "ordered_pairs_fail": {"count": 3, "ordered_mod2": 0, "unordered_mod2": 1},
        "endpoint_correction_required": {"label": label, "uncorrected_mod4": 2, "desired": 0},
        "two_defects_required_mod8": {"exact": 1, "one_defect_only": 5},
        "even_N_root_formula_fails": {"N": 8, "roots": [1, 3, 5, 7]},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--count-limit", type=int, default=4096)
    parser.add_argument("--root-limit", type=int, default=2048)
    parser.add_argument("--maximum-word", type=int, default=3)
    args = parser.parse_args()
    report = {
        "status": "pass",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "parameters": {"count_limit": args.count_limit, "root_limit": args.root_limit,
                       "maximum_word": args.maximum_word},
        "digit_and_comparison_checks": check_digits(args.count_limit),
        "root_recovery": check_roots(args.root_limit),
        "finite_seed_lifts": scalar_seed_replay(args.maximum_word),
        "smith_profile_lifts": smith_replay(args.maximum_word),
        "negative_controls": negative_controls(),
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
