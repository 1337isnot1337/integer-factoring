#!/usr/bin/env python3
"""Independent exact-rational checks of degree-sensitive Smith paths.

The models and full word-rank routine are shared with the existing Smith
checker.  All local factorization, jet, minor, path, and potential checks
here are independent and use rational functions over F_2.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
from itertools import permutations, product
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import check_smith_profile_transfer as basic


ZERO = (0, 1)
ONE = (1, 1)


def poly_mul(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        a <<= 1
        b >>= 1
    return out


def poly_divrem(a, b):
    assert b
    quotient = 0
    while a and a.bit_length() >= b.bit_length():
        shift = a.bit_length() - b.bit_length()
        quotient ^= 1 << shift
        a ^= b << shift
    return quotient, a


def poly_gcd(a, b):
    while b:
        a, b = b, poly_divrem(a, b)[1]
    return a


def rational(numerator, denominator=1):
    assert denominator
    if not numerator:
        return ZERO
    divisor = poly_gcd(numerator, denominator)
    n, nr = poly_divrem(numerator, divisor)
    d, dr = poly_divrem(denominator, divisor)
    assert nr == dr == 0
    return n, d


def add(a, b):
    return rational(poly_mul(a[0], b[1]) ^ poly_mul(b[0], a[1]), poly_mul(a[1], b[1]))


def mul(a, b):
    return rational(poly_mul(a[0], b[0]), poly_mul(a[1], b[1]))


def divide(a, b):
    assert b[0]
    return rational(poly_mul(a[0], b[1]), poly_mul(a[1], b[0]))


def shifted(a, degree):
    return rational(a[0] << degree, a[1]) if degree >= 0 else rational(a[0], a[1] << -degree)


def poly_order(a):
    assert a
    return (a & -a).bit_length() - 1


def order(a):
    return None if a == ZERO else poly_order(a[0]) - poly_order(a[1])


def jet(a, precision):
    assert precision > 0
    if a == ZERO:
        return 0
    assert order(a) >= 0 and a[1] & 1
    inverse = 1
    for degree in range(1, precision):
        if poly_mul(a[1], inverse) >> degree & 1:
            inverse ^= 1 << degree
    mask = (1 << precision) - 1
    assert poly_mul(a[1], inverse) & mask == 1
    return poly_mul(a[0], inverse) & mask


def determinant(matrix):
    if not matrix:
        return 1
    answer = 0
    for permutation in permutations(range(len(matrix))):
        term = 1
        for row, column in enumerate(permutation):
            term = poly_mul(term, matrix[row][column])
        answer ^= term  # All determinant signs coincide in characteristic two.
    return answer


def schur_factors(matrix, degree):
    """Raw complete Schur elimination, with exact bordered-minor checks."""
    original = matrix
    work = [[rational(entry) for entry in row] for row in matrix]
    rows, columns = [], []
    factors = []
    minor_checks = 0
    for step in range(len(matrix[0])):
        denominator = determinant([[original[i][j] for j in columns] for i in rows])
        assert denominator
        for i, j in product(range(len(work)), range(len(work[0]))):
            if i in rows or j in columns:
                assert work[i][j] == ZERO
                continue
            bordered = [[original[a][b] for b in columns + [j]] for a in rows + [i]]
            expected = rational(determinant(bordered), denominator)
            assert work[i][j] == expected
            if expected != ZERO:
                assert order(expected) <= (step + 1) * degree
            minor_checks += 1
        candidates = [(order(entry), i, j) for i, row in enumerate(work)
                      for j, entry in enumerate(row) if entry != ZERO]
        assert candidates
        valuation, pi, pj = min(candidates)
        pivot = work[pi][pj]
        left = [shifted(work[i][pj], -valuation) for i in range(len(work))]
        right = [divide(entry, pivot) for entry in work[pi]]
        assert all(entry == ZERO or order(entry) >= 0 for entry in left + right)
        factors.append((valuation, left, right))
        for i, j in product(range(len(work)), range(len(work[0]))):
            work[i][j] = add(work[i][j], shifted(mul(left[i], right[j]), valuation))
        rows.append(pi)
        columns.append(pj)
    assert all(entry == ZERO for row in work for entry in row)
    for i, j in product(range(len(matrix)), range(len(matrix[0]))):
        reconstructed = ZERO
        for valuation, left, right in factors:
            reconstructed = add(reconstructed, shifted(mul(left[i], right[j]), valuation))
        assert reconstructed == rational(matrix[i][j])
    return factors, minor_checks


def normalized_matrix(local, cut, side, heights, infinity):
    auxiliary = len(heights)
    entries = [[set() for _ in range(4)] for _ in range(4 * auxiliary * auxiliary)]
    direction = -1 if infinity else 1
    for (x, y, u, v), value in local.items():
        left, right = (2 * x + y, 2 * u + v) if cut == 0 else (2 * x + u, 2 * y + v)
        for a, b in product(range(auxiliary), repeat=2):
            physical, other = (left, right) if side == 0 else (right, left)
            entries[(physical * auxiliary + a) * auxiliary + b][other] = {
                direction * k + heights[a] - heights[b]
                for k in range(value[a][b].bit_length()) if value[a][b] >> k & 1
            }
    support = [k for row in entries for entry in row for k in entry]
    minimum, maximum = min(support), max(support)
    matrix = [[sum(1 << (k - minimum) for k in entry) for entry in row] for row in entries]
    return matrix, minimum, maximum - minimum


def local_path_data(local, cut, side, heights, infinity):
    matrix, common, degree = normalized_matrix(local, cut, side, heights, infinity)
    factors, minor_checks = schur_factors(matrix, degree)
    precision = 4 * degree + 1
    auxiliary = len(heights)
    exponents = []
    deltas = []
    finite_deltas = []
    jet_checks = exact_zeros = 0
    for pivot_order, left, right in factors:
        known = precision - pivot_order
        assert known > 0
        restored = common + pivot_order
        exponents.append(-restored)
        left_jets = [jet(entry, known) for entry in left]
        right_jets = [jet(entry, known) for entry in right]
        for entry, truncated in zip(left + right, left_jets + right_jets):
            bound = poly_order(truncated) if truncated else known
            assert entry == ZERO or order(entry) >= bound
            # Complete degree precision proves that every zero jet is exact.
            assert bool(truncated) == (entry != ZERO)
            if truncated:
                assert bound == order(entry)
            else:
                exact_zeros += 1
            jet_checks += 1
        scalar_order = min(order(entry) for entry in right if entry != ZERO)
        edge_costs = [[None] * auxiliary for _ in range(auxiliary)]
        finite_costs = [[0] * auxiliary for _ in range(auxiliary)]
        for a, b in product(range(auxiliary), repeat=2):
            indices = [(physical * auxiliary + a) * auxiliary + b for physical in range(4)]
            possible = [order(left[index]) for index in indices if left[index] != ZERO]
            edge_costs[a][b] = min(possible) + scalar_order if possible else None
            finite_costs[a][b] = min(poly_order(left_jets[index]) if left_jets[index] else known
                                      for index in indices) + scalar_order
            if edge_costs[a][b] is not None:
                assert edge_costs[a][b] == finite_costs[a][b]
        deltas.append(edge_costs)
        finite_deltas.append(finite_costs)
    assert len(factors) == 4
    return {
        "cut": cut, "side": side, "heights": heights, "infinity": infinity,
        "exponents": exponents, "deltas": deltas, "finite_deltas": finite_deltas,
        "boundary_shift": abs(heights[0]), "precision": precision, "degree": degree,
        "minor_checks": minor_checks, "jet_checks": jet_checks, "exact_zeros": exact_zeros,
    }


def transition_matrix(data, s, finite=False):
    auxiliary = len(data["heights"])
    matrix = [[Fraction(0) for _ in range(auxiliary)] for _ in range(auxiliary)]
    for exponent, delta in zip(data["exponents"], data["finite_deltas"] if finite else data["deltas"]):
        for a, b in product(range(auxiliary), repeat=2):
            if delta[a][b] is not None:
                matrix[a][b] += s ** (exponent - delta[a][b])
    return matrix


def best_potential(matrix):
    dimension = len(matrix)
    potentials = [(Fraction(1),)] if dimension == 1 else [
        (Fraction(1), ratio) for ratio in map(Fraction, ("1/9", "1/3", "1/2", "1", "3/2", "2", "3", "9"))
    ]
    answer = None
    for vector in potentials:
        beta = max(sum(matrix[a][b] * vector[b] for b in range(dimension)) / vector[a]
                   for a in range(dimension))
        beta = max(beta, Fraction(1, 10**6))
        assert all(sum(matrix[a][b] * vector[b] for b in range(dimension)) <= beta * vector[a]
                   for a in range(dimension))
        if answer is None or beta < answer[0]:
            answer = beta, vector
    return answer


def weighted_path_value(matrix, length):
    vector = [Fraction(1)] + [Fraction(0)] * (len(matrix) - 1)
    for _ in range(length):
        vector = [sum(vector[a] * matrix[a][b] for a in range(len(matrix))) for b in range(len(matrix))]
    return vector[0]


def path_scores(data, length):
    dimension = len(data["heights"])
    counts = {(0, 0): 1}
    for _ in range(length):
        following = defaultdict(int)
        for (state, score), count in counts.items():
            for exponent, delta in zip(data["exponents"], data["deltas"]):
                for target in range(dimension):
                    if delta[state][target] is not None:
                        following[target, score + exponent - delta[state][target]] += count
        counts = dict(following)
    return {score: count for (state, score), count in counts.items() if state == 0}


def describe(data):
    return {key: data[key] for key in ("cut", "side", "heights", "infinity", "exponents", "deltas", "precision")}


def check_family(name, local, degree, maximum_word):
    dimension = len(next(iter(local.values())))
    heights_list = [(0,)] if dimension == 1 else [(0, 0), (0, 1), (0, -1), (1, -1)]
    factors = [local_path_data(local, cut, side, heights, infinity)
               for cut, side, heights, infinity in product(range(2), range(2), heights_list, (False, True))]
    counts = Counter()
    for data in factors:
        for field in ("minor_checks", "jet_checks", "exact_zeros"):
            counts[field] += data[field]
    ranks_by_length = {}
    for length in range(1, maximum_word + 1):
        ranks_by_length[length], endpoints = basic.full_coefficients(local, degree, length)
        counts["ordered_endpoint_entries"] += endpoints
    best = None
    for s in map(Fraction, ("1", "3/2", "3", "5")):
        potential_data = []
        for data in factors:
            matrix = transition_matrix(data, s)
            finite = transition_matrix(data, s, finite=True)
            assert all(matrix[a][b] <= finite[a][b] for a, b in product(range(dimension), repeat=2))
            beta, vector = best_potential(matrix)
            potential_data.append((data, beta, vector, matrix))
            for length in range(1, maximum_word + 1):
                score_counts = path_scores(data, length)
                weighted = weighted_path_value(matrix, length)
                assert weighted == sum(count * s**score for score, count in score_counts.items())
                assert weighted <= beta**length  # Both boundary supports are the first coordinate.
                for k in range(degree * length + 1):
                    threshold = (k if data["infinity"] else -k) - data["boundary_shift"]
                    admitted = sum(count for score, count in score_counts.items() if score >= threshold)
                    split = max(0, length * max(data["exponents"]) - threshold + 1)
                    actual = ranks_by_length[length][data["cut"]][k]
                    assert actual <= split * admitted, (name, data, length, k, actual, split, admitted)
                    assert admitted <= s ** (-threshold) * weighted
                    assert actual <= split * s ** (-threshold) * beta**length
                    counts["coefficient_rank_checks"] += 1
                    counts["rational_path_and_potential_checks"] += 3
        for first_cut in range(2):
            high = [entry for entry in potential_data if entry[0]["infinity"] and entry[0]["cut"] == first_cut]
            low = [entry for entry in potential_data if not entry[0]["infinity"] and entry[0]["cut"] == 1 - first_cut]
            first = min(high, key=lambda entry: entry[1])
            second = min(low, key=lambda entry: entry[1])
            B = first[1] * second[1]
            if best is None or B < best[0]:
                best = (B, s, first, second)
    B, s, first, second = best
    return {
        "name": name, "polynomial_degree": degree, "auxiliary_dimension": dimension,
        **dict(counts), "best_sampled_certificate": {
            "B_path": str(B), "s": str(s),
            "infinity": {"data": describe(first[0]), "beta": str(first[1]),
                         "potential": list(map(str, first[2])), "matrix": [[str(x) for x in row] for row in first[3]]},
            "zero": {"data": describe(second[0]), "beta": str(second[1]),
                     "potential": list(map(str, second[2])), "matrix": [[str(x) for x in row] for row in second[3]]},
        },
    }


def negative_controls():
    # A truncated zero alone gives no exact-zero certificate.
    later_term = rational(1 << 7, 1 ^ (1 << 3))
    assert jet(later_term, 7) == 0 and order(later_term) == 7
    # Division by the pivot reduces the known precision.
    original = rational(1 << 9)
    quotient = shifted(original, -4)
    assert jet(original, 9) == 0 and jet(quotient, 5) == 0 and order(quotient) == 5
    # Original polynomial degree, not merely observed low coefficients,
    # is needed in the bordered-minor zero certificate.
    assert determinant([[1, 0], [0, 1 << 7]]) == 1 << 7
    # A reducible Perron root need not admit a positive potential at equality.
    matrix = [[Fraction(1), Fraction(1)], [Fraction(0), Fraction(1)]]
    assert weighted_path_value(matrix, 4) == 1
    assert all(matrix[0][0] * 1 + matrix[0][1] * ratio > 1
               for ratio in map(Fraction, ("1/9", "1/2", "1", "2")))
    # Dropping a Laurent boundary pole changes coefficient indices.
    assert mul(shifted(ONE, -2), shifted(ONE, 2)) == ONE
    assert shifted(ONE, 2) != ONE
    return {
        "zero_jet_not_zero_series": {"numerator": 1 << 7, "denominator": 9, "jet_precision": 7},
        "pivot_precision_loss": {"input_precision": 9, "pivot_order": 4, "valid_quotient_precision": 5},
        "uncertified_degree_invalidates_exact_zero": True,
        "positive_potential_may_need_strict_perron_slack": True,
        "boundary_poles_require_restoration": True,
    }


def rational_matrix_product(left, right):
    result = [[ZERO] * 2 for _ in range(2)]
    for a, b, c in product(range(2), repeat=3):
        result[a][c] = add(result[a][c], mul(left[a][b], right[b][c]))
    return result


def native_crosscheck(probe):
    """Compare the separate C++ jets with literal rational gauge products."""
    cases = [("elementary", *basic.elementary(2)),
             ("pairing", *basic.pairing((1, 6, 9), 0)),
             ("loop", *basic.loop(13, 18))]
    choices = [(1, 0, 0, 0, (0, 0, 0, 0), 9),
               (1, 0, 1, 0, (0, 0, 0, 0), 9),
               (1, 1, 0, 0, (0, 0, 0, 0), 9),
               (1, 0, 2, 3, (0, 0, 0, 0), 9),
               (1, 1, 1, 1, (1, -1, 0, 0), 13),
               (2, 1, 0, 0, (0, 1, -1, 0), 6)]
    totals = Counter()
    with TemporaryDirectory(prefix="smith-path-probe-") as directory:
        for name, local, degree in cases:
            source = Path(directory) / f"{name}.txt"
            source.write_text("\n".join(str(entry) for labels in sorted(local)
                                        for row in local[labels] for entry in row) + "\n")
            for ramification, a, b, p, physical, basis in choices:
                command = [str(probe), str(source), str(degree), str(ramification), str(a), str(b), str(p),
                           *(str(value) for value in physical), "0", str(basis)]
                native = json.loads(subprocess.check_output(command, text=True))
                S = [[rational((basis >> (2 * i + j)) & 1) for j in range(2)] for i in range(2)]
                Sinv = [[S[1][1], S[0][1]], [S[1][0], S[0][0]]]
                H = [[shifted(ONE, a), rational(p)], [ZERO, shifted(ONE, b)]]
                Hinv = [[shifted(ONE, -a), shifted(rational(p), -a - b)], [ZERO, shifted(ONE, -b)]]
                gauge = rational_matrix_product(H, S)
                inverse = rational_matrix_product(Sinv, Hinv)
                assert rational_matrix_product(gauge, inverse) == [[ONE, ZERO], [ZERO, ONE]]
                for cut, side, infinity in product(range(2), range(2), (False, True)):
                    entries = [[set() for _ in range(4)] for _ in range(16)]
                    for (x, y, u, v), coefficients in local.items():
                        raw = [[ZERO] * 2 for _ in range(2)]
                        for i, j in product(range(2), repeat=2):
                            for k in range(coefficients[i][j].bit_length()):
                                if coefficients[i][j] >> k & 1:
                                    raw[i][j] = add(raw[i][j], shifted(ONE, (-1 if infinity else 1) * k * ramification))
                        changed = rational_matrix_product(rational_matrix_product(gauge, raw), inverse)
                        exponent = sum((2 * bit - 1) * weight for bit, weight in zip((x, y, u, v), physical))
                        left, right = (2 * x + y, 2 * u + v) if cut == 0 else (2 * x + u, 2 * y + v)
                        for i, j in product(range(2), repeat=2):
                            value = shifted(changed[i][j], exponent)
                            assert value[1] & (value[1] - 1) == 0
                            shift = value[1].bit_length() - 1
                            entry = {k - shift for k in range(value[0].bit_length()) if value[0] >> k & 1}
                            row, column = (left, right) if side == 0 else (right, left)
                            entries[(row * 2 + i) * 2 + j][column] = entry
                    support = [value for row in entries for entry in row for value in entry]
                    common, highest = min(support), max(support)
                    matrix = [[sum(1 << (value - common) for value in entry) for entry in row] for row in entries]
                    factors, minor_checks = schur_factors(matrix, highest - common)
                    base = [-(common + valuation) for valuation, _, _ in factors]
                    kappa = []
                    for _, left, _ in factors:
                        values = []
                        for i, j in product(range(2), repeat=2):
                            possible = [order(left[(physical_index * 2 + i) * 2 + j]) for physical_index in range(4)
                                        if left[(physical_index * 2 + i) * 2 + j] != ZERO]
                            values.append(min(possible) if possible else None)
                        kappa.append(values)
                    found = native[4 * cut + 2 * side + int(infinity)]
                    assert found["base_exponents"] == base, (name, command, found, base)
                    assert found["kappa"] == kappa, (name, command, found, kappa)
                    totals["native_factors_compared"] += 1
                    totals["native_edge_orders_compared"] += 16
                    totals["additional_bordered_minor_checks"] += minor_checks
                totals["native_model_gauge_instances"] += 1
    return dict(totals)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--maximum-word", type=int, default=3)
    parser.add_argument("--native-probe", type=Path,
                        help="optional prebuilt smith_walk_probe executable from the separate native screen")
    args = parser.parse_args()
    cases = [("scalar elementary", *basic.elementary()),
             ("auxiliary elementary", *basic.elementary(2)),
             ("ordered pairing", *basic.pairing((1, 6, 9), 0)),
             ("mixed nilpotent loop", *basic.loop(13, 18))]
    results = [check_family(name, local, degree, args.maximum_word) for name, local, degree in cases]
    assert Fraction(results[0]["best_sampled_certificate"]["B_path"]) == Fraction(256, 27)
    assert sum(item["exact_zeros"] for item in results) > 0
    report = {
        "status": "pass", "field": 2, "maximum_word": args.maximum_word,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "base_checker_sha256": hashlib.sha256(Path(basic.__file__).read_bytes()).hexdigest(),
        "cases": results, "negative_controls": negative_controls(),
        "scope": "Exact rational-function, jet, bordered-minor, path, and finite-rank checks; no improving seed.",
    }
    if args.native_probe is not None:
        report["independent_native_crosscheck"] = native_crosscheck(args.native_probe.resolve())
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
