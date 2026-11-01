#!/usr/bin/env python3
"""Exact Pauli/Walsh rank interface and 64-label atom-dictionary exclusion.

Only integers and rational arithmetic enter proof decisions. The dictionary
allows both Fourier orientations and any transformations preserving cut ranks.
It does not include arbitrary linear recombinations as single atoms.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
from itertools import product
import json
from pathlib import Path
import random


def matrix_rank(rows, p):
    matrix = [list(row) for row in rows]
    count = 0
    for col in range(len(matrix[0])):
        pivot = next((i for i in range(count, len(matrix)) if matrix[i][col] % p), None)
        if pivot is None:
            continue
        matrix[count], matrix[pivot] = matrix[pivot], matrix[count]
        scale = pow(matrix[count][col] % p, -1, p)
        matrix[count] = [x * scale % p for x in matrix[count]]
        for i in range(count + 1, len(matrix)):
            scale = matrix[i][col] % p
            if scale:
                matrix[i] = [(x - scale * y) % p for x, y in zip(matrix[i], matrix[count])]
        count += 1
        if count == len(matrix):
            break
    return count


def walsh(values):
    result = list(values)
    half = 1
    while half < len(result):
        for start in range(0, len(result), 2 * half):
            for i in range(start, start + half):
                a, b = result[i], result[i + half]
                result[i], result[i + half] = a + b, a - b
        half *= 2
    return result


def pauli_operator(coefficients, m, p):
    d = 1 << m
    matrix = [[0] * (d * d) for _ in range(d * d)]
    for x in range(d):
        for z in range(d):
            coefficient = coefficients[x * d + z] % p
            if not coefficient:
                continue
            for i in range(d):
                for j in range(d):
                    sign = -1 if (z & (i ^ j)).bit_count() & 1 else 1
                    row, col = (i ^ x) * d + (j ^ x), i * d + j
                    matrix[row][col] = (matrix[row][col] + coefficient * sign) % p
    return matrix


def other_cut(matrix, d):
    return [[matrix[x * d + y][u * d + v] for y in range(d) for v in range(d)]
            for x in range(d) for u in range(d)]


def swap(d):
    return [[int(x == v and y == u) for u in range(d) for v in range(d)]
            for x in range(d) for y in range(d)]


def interface_checks():
    rng = random.Random(60241)
    ranks = crossing_entries = eigenvectors = 0
    for p in (3, 5, 7):
        for m in (1, 2, 3):
            d, n = 1 << m, 1 << (2 * m)
            cases = (list(product(range(p), repeat=4)) if m == 1 and p == 3
                     else [tuple(rng.randrange(p) for _ in range(n)) for _ in range(12)])
            for coefficients in cases:
                matrix = pauli_operator(coefficients, m, p)
                spectrum = walsh(coefficients)
                assert matrix_rank(matrix, p) == sum(a % p != 0 for a in spectrum)
                assert matrix_rank(other_cut(matrix, d), p) == sum(a % p != 0 for a in coefficients)
                ranks += 2
                for u, v in ((0, 0), (d - 1, d - 1), (d - 1, 0), (0, d - 1)):
                    vector = [(-1 if (v & i).bit_count() & 1 else 1) if j == i ^ u else 0
                              for i in range(d) for j in range(d)]
                    value = spectrum[v * d + u] % p
                    assert [sum(a * b for a, b in zip(row, vector)) % p for row in matrix] == [value * a % p for a in vector]
                    eigenvectors += 1
            crossing = [(-1 if (x & z).bit_count() & 1 else 1) * pow(d, -1, p) % p
                        for x in range(d) for z in range(d)]
            assert pauli_operator(crossing, m, p) == swap(d)
            crossing_entries += n * n
    return dict(fields=[3, 5, 7], label_bits=[1, 2, 3],
                cut_rank_tests=ranks, bell_eigenvector_tests=eigenvectors,
                crossing_identity_entries=crossing_entries)


def dictionary_checks():
    m, d, p = 6, 64, 3
    kinds = defaultdict(list)
    for x in range(d):
        for z in range(d):
            both = (x & z).bit_count()
            first, second = (x & ~z).bit_count(), (z & ~x).bit_count()
            kinds[m - both - first - second, first, second, both].append(x * d + z)
    assert len(kinds) == 84 and sum(map(len, kinds.values())) == d * d
    table = Counter()
    rows = []
    eligible = []
    # This upper rational exponent is used only to enlarge the admissible
    # individual-atom list. Rejecting the enlarged list is therefore sound.
    assert 577 * 577 - 2 * 408 * 408 == 1
    for kind, indices in sorted(kinds.items()):
        function = [0] * (d * d)
        for i in indices:
            function[i] = 1
        spectrum = walsh(function)
        row = sum(x % p != 0 for x in spectrum)
        col = len(indices)
        assert row * col >= d * d
        transformed_twice = walsh([x % p for x in spectrum])
        assert all((x - d * d * y) % p == 0 for x, y in zip(transformed_twice, function))
        table[row, col] += 1
        cost = row * col
        permitted = cost ** 204 < d ** 577
        if permitted:
            eligible.append((row, col))
        rows.append(dict(type=kind, row_rank=row, column_rank=col,
                         rank_product=cost, passes_enlarged_individual_filter=permitted))
    expected = {(4096, 1): 4, (2688, 6): 12, (2688, 15): 12,
                (2048, 20): 6, (2892, 30): 12, (2688, 60): 24,
                (1440, 90): 4, (3972, 120): 4, (720, 180): 6}
    assert dict(table) == expected
    assert len(eligible) == 46
    ratio = min(Fraction(max(r, c), min(r, c)) for r, c in eligible)
    assert ratio == Fraction(482, 5)
    assert all(r * c < 129600 for r, c in eligible)
    assert 129600 ** 204 > d ** 577
    # From sum min(r_i,c_i)>=d, total S is at least d*sqrt(ratio).
    # Prove that lower bound exceeds d^(577/408)>d^sqrt(2).
    assert (d * d * ratio.numerator) ** 204 > d ** 577 * ratio.denominator ** 204
    return dict(field=p, label_bits=m, labels=d, all_types=len(kinds),
                enlarged_filter_types=len(eligible),
                minimum_large_to_small_ratio={'numerator':ratio.numerator,'denominator':ratio.denominator},
                exact_pell_identity=1, all_cut_preserving_basis_changes_allowed=True,
                fourier_duals_included=True, arbitrary_extension_scalars_allowed=True,
                rank_sum_criterion_excluded=True,
                table=[dict(row_rank=r, column_rank=c, types=count, rank_product=r*c)
                       for (r,c), count in sorted(table.items())],
                types=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = dict(status='PASS',
                  scope='64-label Pauli type-atom dictionary exclusion; no asymptotic factoring claim.',
                  interface=interface_checks(), dictionary=dictionary_checks(),
                  checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({**result, 'dictionary':{key:value for key,value in result['dictionary'].items() if key!='types'}}, indent=2))


if __name__ == '__main__':
    main()
