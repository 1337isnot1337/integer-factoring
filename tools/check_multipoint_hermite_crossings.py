#!/usr/bin/env python3
"""Exact Hermite identities, rank controls, a witness, and content accounting.

This checks a reusable multipoint construction and a rho < 19/20 example.
It does not improve the known factoring exponent one half.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path


def rank(matrix, prime):
    if not matrix:
        return 0
    work = [[value % prime for value in row] for row in matrix]
    count = 0
    for column in range(len(work[0])):
        pivot = next((i for i in range(count, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[count], work[pivot] = work[pivot], work[count]
        inverse = pow(work[count][column], -1, prime)
        work[count] = [value * inverse % prime for value in work[count]]
        for i in range(count + 1, len(work)):
            factor = work[i][column]
            if factor:
                work[i] = [(a - factor * b) % prime
                           for a, b in zip(work[i], work[count])]
        count += 1
        if count == len(work):
            break
    return count


def solve(matrix, vector, prime):
    size = len(vector)
    work = [[value % prime for value in row] + [vector[i] % prime]
            for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next(i for i in range(column, size) if work[i][column])
        work[column], work[pivot] = work[pivot], work[column]
        inverse = pow(work[column][column], -1, prime)
        work[column] = [value * inverse % prime for value in work[column]]
        for i in range(size):
            if i == column:
                continue
            factor = work[i][column]
            if factor:
                work[i] = [(a - factor * b) % prime
                           for a, b in zip(work[i], work[column])]
    return [row[-1] for row in work]


def multiply(first, second, prime):
    result = [0] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            result[i + j] = (result[i + j] + a * b) % prime
    return result


def local_polynomial(value, point, degree, prime):
    padded = list(value) + [0] * (degree + 1 - len(value))
    if point is None:
        return padded[::-1]
    result = [0] * (degree + 1)
    for j, coefficient in enumerate(padded):
        for k in range(j + 1):
            result[k] = (result[k] + coefficient * comb(j, k)
                         * pow(point, j - k, prime)) % prime
    return result


def flatten(table, alphabet, cut, coefficient=0):
    pairs = list(product(range(alphabet), repeat=2))
    result = []
    for a, b in pairs:
        row = []
        for c, d in pairs:
            key = (a, b, c, d) if cut == 0 else (a, c, b, d)
            value = table.get(key, ())
            row.append(value[coefficient] if coefficient < len(value) else 0)
        result.append(row)
    return result


def local_profiles(seed, alphabet, degree, point, prime):
    table = {key: local_polynomial(value, point, degree, prime)
             for key, value in seed.items()}
    profiles = []
    size = alphabet * alphabet
    for cut in (0, 1):
        coefficients = [flatten(table, alphabet, cut, k) for k in range(degree + 1)]
        previous_rank = 0
        previous_increment = 0
        orders = []
        for length in range(1, degree * size + 2):
            toeplitz = [[0] * (length * size) for _ in range(length * size)]
            for i in range(length):
                for j in range(i + 1):
                    if i - j > degree:
                        continue
                    block = coefficients[i - j]
                    for r in range(size):
                        toeplitz[i * size + r][j * size:(j + 1) * size] = block[r]
            current_rank = rank(toeplitz, prime)
            increment = current_rank - previous_rank
            assert previous_increment <= increment <= size
            orders.extend([length - 1] * (increment - previous_increment))
            previous_rank, previous_increment = current_rank, increment
            if increment == size:
                break
        assert len(orders) == size
        profiles.append(orders)
    return profiles


def interpolation_weights(degree, points, multiplicities, prime):
    assert sum(multiplicities) == degree + 1
    rows = []
    for point, count in zip(points, multiplicities):
        for k in range(count):
            rows.append([int(j == degree - k) if point is None else
                         (comb(j, k) * pow(point, j - k, prime) % prime if j >= k else 0)
                         for j in range(degree + 1)])
    # Evaluation at z=1 is the all-one row in the monomial basis.
    weights = solve(list(map(list, zip(*rows))), [1] * (degree + 1), prime)
    assert [sum(weights[i] * rows[i][j] for i in range(degree + 1)) % prime
            for j in range(degree + 1)] == [1] * (degree + 1)
    return weights


def tensor_power(seed, alphabet, copies, prime):
    words = list(product(range(alphabet), repeat=copies))
    result = {}
    for labels in product(range(len(words)), repeat=4):
        value = [1]
        for position in range(copies):
            key = tuple(words[label][position] for label in labels)
            factor = seed.get(key, ())
            if not factor:
                value = []
                break
            value = multiply(value, factor, prime)
        if value and any(value):
            result[labels] = value
    return result


def low_order_word_count(orders, copies, cutoff):
    counts = Counter({0: 1})
    for _ in range(copies):
        updated = Counter()
        for total, count in counts.items():
            for order in orders:
                if total + order <= cutoff:
                    updated[total + order] += count
        counts = updated
    return sum(counts.values())


def check_power(seed, alphabet, degree, copies, points, multiplicities, prime):
    table = tensor_power(seed, alphabet, copies, prime)
    full_degree = degree * copies
    weights = interpolation_weights(full_degree, points, multiplicities, prime)
    jets = []
    rank_checks = 0
    for point, count in zip(points, multiplicities):
        local = {key: local_polynomial(value, point, full_degree, prime)
                 for key, value in table.items()}
        profiles = local_profiles(seed, alphabet, degree, point, prime)
        for k in range(count):
            jets.append({key: value[k] for key, value in local.items() if value[k]})
            for cut in (0, 1):
                actual = rank(flatten(local, alphabet ** copies, cut, k), prime)
                bound = (k + 1) * low_order_word_count(profiles[cut], copies, k)
                assert actual <= bound, (prime, point, copies, cut, k, actual, bound)
                rank_checks += 1
    dimension = alphabet ** copies
    for key in product(range(dimension), repeat=4):
        x, y, u, v = key
        expected = int(x == v and y == u)
        assert sum(table.get(key, ())) % prime == expected
        assert sum(weight * jet.get(key, 0) for weight, jet in zip(weights, jets)) % prime == expected
    return {'prime': prime, 'alphabet': alphabet, 'copies': copies,
            'points': points, 'multiplicities': multiplicities,
            'weights': weights, 'complete_entries': dimension ** 4,
            'rank_checks': rank_checks}


def affine_seed(prime):
    # K(z)=(1-z) I + z Swap, with a literal swap endpoint at one.
    seed = {}
    for x, y, u, v in product(range(2), repeat=4):
        identity, crossing = int(x == u and y == v), int(x == v and y == u)
        value = [identity, (crossing - identity) % prime]
        if any(value):
            seed[x, y, u, v] = value
    return seed


def mixing_seed():
    prime = alphabet = 3
    seed = {}
    for x, y, u, v in product(range(alphabet), repeat=4):
        if (x + y - u - v) % prime:
            continue
        a, b = (x - u) % prime, (u - v) % prime
        value = [(1 - a * a) % prime, (-2 * a * b) % prime, (-b * b) % prime]
        if any(value):
            seed[x, y, u, v] = value
    return seed


def content_ledger_control(seed, expected):
    # 2*z*(z-2) equals one at z=1 over F_3.  Degree-five homogenization
    # contributes one further common zero at infinity beyond its degree four.
    scaled = {key: multiply(value, [0, 2, 2], 3)
              for key, value in seed.items()}
    slopes = {0: Fraction(15, 14), 2: Fraction(13, 28),
              None: Fraction(13, 28)}
    scales = {0: (Fraction(1), Fraction(9, 5)),
              2: (Fraction(5, 2), Fraction(1)),
              None: (Fraction(5, 2), Fraction(1))}
    shifted = {}
    for point, profiles in expected.items():
        actual = local_profiles(scaled, 3, 5, point, 3)
        assert actual == [[order + 1 for order in profile]
                          for profile in profiles]
        shifted[str(point)] = actual
        s, r = scales[point]
        before = (sum(s ** -order for order in profiles[0]) ** 28
                  * sum(r ** -order for order in profiles[1]) ** 28
                  * (s * r) ** int(28 * slopes[point]))
        after = (sum(s ** -order for order in actual[0]) ** 28
                 * sum(r ** -order for order in actual[1]) ** 28
                 * (s * r) ** int(28 * (slopes[point] + 1)))
        assert before == after
    primitive_totals = [sum(Fraction(sum(profiles[cut]), 9)
                            for profiles in expected.values())
                        for cut in (0, 1)]
    shifted_totals = [sum(Fraction(sum(profiles[cut]), 9)
                          for profiles in shifted.values())
                      for cut in (0, 1)]
    assert primitive_totals == [2, 2]
    assert shifted_totals == [5, 5]
    assert sum(slopes.values()) == 2
    assert sum(value + 1 for value in slopes.values()) == 5
    strict_floor_controls = {}
    for alphabet in (2, 3, 4, 5):
        n = alphabet * alphabet
        assert 4 ** n * n ** n > (n + 1) ** (n + 1)
        strict_floor_controls[str(alphabet)] = True
    return scaled, {
        'primitive_degree': 2, 'homogeneous_degree': 5,
        'common_homogeneous_content_degree': 3,
        'selected_content_orders': [1, 1, 1],
        'primitive_mean_order_totals': [str(value) for value in primitive_totals],
        'shifted_mean_order_totals': [str(value) for value in shifted_totals],
        'shifted_profiles': shifted,
        'exact_pressure_redistribution_checks': 3,
        'strict_three_quarter_floor_integer_controls': strict_floor_controls,
        'scope': 'Finite controls for the analytic homogeneous-content proof; not a finite verification of all alphabets.'}


def check():
    cases = []
    for copies in range(1, 5):
        counts = [(copies + 2) // 2, (copies + 1) // 2]
        cases.append(check_power(affine_seed(2), 2, 1, copies, [0, None], counts, 2))
    for prime in (3, 5):
        cases.append(check_power(affine_seed(prime), 2, 1, 3, [0, 2, None], [2, 1, 1], prime))
    seed = mixing_seed()
    expected = {
        0: [[0] * 9, [0, 1, 1, 1, 1, 2, 4, 4, 4]],
        2: [[0] * 3 + [1] * 3 + [2] * 3, [0] * 9],
        None: [[0] * 3 + [1] * 3 + [2] * 3, [0] * 9],
    }
    for point, orders in expected.items():
        assert local_profiles(seed, 3, 2, point, 3) == orders
    scaled_seed, content_control = content_ledger_control(seed, expected)
    cases.append(check_power(seed, 3, 2, 1, [0, 2, None], [1, 1, 1], 3))
    cases.append(check_power(seed, 3, 2, 2, [0, 2, None], [2, 2, 1], 3))
    cases.append(check_power(scaled_seed, 3, 5, 1, [0, 2, None], [2, 2, 2], 3))
    s0, s1 = Fraction(9, 5), Fraction(5, 2)
    b0 = 9 * (1 + 4 / s0 + 1 / s0 ** 2 + 3 / s0 ** 4)
    b1 = 27 * (1 + 1 / s1 + 1 / s1 ** 2)
    left0, left1, right = b0 ** 140 * s0 ** 150, b1 ** 140 * s1 ** 65, 3 ** 532
    assert left0 < right and left1 < right
    return {'status': 'pass', 'cases': cases,
            'complete_entries_checked': sum(case['complete_entries'] for case in cases),
            'jet_rank_checks': sum(case['rank_checks'] for case in cases),
            'ternary_mixing_profiles': {str(point): value for point, value in expected.items()},
            'homogeneous_content_control': content_control,
            'rational_certificate': {'rho_upper': '19/20',
                                     'jet_slopes': ['15/14', '13/28', '13/28'],
                                     'sum_slopes': '2', 'scales': ['9/5', '5/2', '5/2'],
                                     'common_power': 140, 'right_power_of_three': 532,
                                     'all_integer_comparisons_pass': True,
                                     'left_numerator_bits': [left0.numerator.bit_length(), left1.numerator.bit_length()],
                                     'left_denominator_bits': [left0.denominator.bit_length(), left1.denominator.bit_length()]},
            'scope': 'Exact multipoint interpolation, finite rank controls, and homogeneous-content ledger controls. The displayed rho is weaker than half-rate; the unweighted same-place criterion has a universal rho > 3/4 floor.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = check()
    report['source_sha256'] = {Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}
    rendered = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
