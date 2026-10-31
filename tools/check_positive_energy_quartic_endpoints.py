#!/usr/bin/env python3
"""Exact projected endpoints for positive-energy quartic perturbations.

Checks all four-mode local entries and every two-copy fixed-mode type.
Also checks termwise denominator cancellation for f=1/[w(1+w)].
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_schur_downset_crossing import sparse_kernel
from weighted_smith_oracle import polynomial_exact_divide, polynomial_multiply


ROOT = Path(__file__).resolve().parent.parent


def local_polynomials(n=4, ambient=4):
    kernel = sparse_kernel(tuple(range(n)), ambient)
    rows = list(range(2 * n, 4 * n))
    columns = list(range(2 * n))
    matrix = [[kernel[i][j] for j in columns] for i in rows]

    @lru_cache(None)
    def minor(row_mask, column_mask):
        if row_mask.bit_count() != column_mask.bit_count():
            return 0
        if not row_mask:
            return 1
        bit = row_mask & -row_mask
        i = bit.bit_length() - 1
        remaining_rows = row_mask ^ bit
        result = 0
        columns_left = column_mask
        while columns_left:
            bit = columns_left & -columns_left
            j = bit.bit_length() - 1
            columns_left ^= bit
            if matrix[i][j]:
                result ^= polynomial_multiply(matrix[i][j], minor(remaining_rows, column_mask ^ bit))
        return result

    q = 1 << n
    gaussian = []
    for x, y, u, v in product(range(q), repeat=4):
        gaussian.append(minor(v | (y << n), u | (x << n)))
    return gaussian, minor


def quartic_polynomials(minor, quartic, n=4):
    a, b, c, d = quartic
    q = 1 << n
    values = []
    for x, y, u, v in product(range(q), repeat=4):
        if (x >> a & 1) and (y >> b & 1) and (u >> c & 1) and (v >> d & 1):
            row_mask = (v ^ (1 << d)) | ((y ^ (1 << b)) << n)
            column_mask = (u ^ (1 << c)) | ((x ^ (1 << a)) << n)
            values.append(minor(row_mask, column_mask))
        else:
            values.append(0)
    return values


def index(x, y, u, v, q=16):
    return ((x * q + y) * q + u) * q + v


def fixed_mode_words(n=4):
    result = defaultdict(list)
    for first, second in product(range(1 << n), repeat=2):
        counts = tuple(((first >> a) & 1) + ((second >> a) & 1) for a in range(n))
        result[counts].append((first, second))
    return result


def endpoint_control(gaussian, quartic_values, quartic, word_types, n=4):
    a, b, c, d = quartic
    delta = b + d - a - c
    assert delta > 0
    energy = [sum(a for a in range(n) if state >> a & 1) for state in range(1 << n)]
    endpoint = [(first ^ second).bit_count() & 1 for first, second in zip(gaussian, quartic_values)]
    positive_entries = 0
    for x, y, u, v in product(range(1 << n), repeat=4):
        value = endpoint[index(x, y, u, v)]
        difference = energy[y] + energy[v] - energy[x] - energy[u]
        if value:
            assert difference >= 0
            positive_entries += difference > 0
        if difference == 0:
            assert value == int(x == v and y == u)
    checks = 0
    for words in word_types.values():
        for x, y, u, v in product(words, repeat=4):
            value = endpoint[index(x[0], y[0], u[0], v[0])]
            value &= endpoint[index(x[1], y[1], u[1], v[1])]
            assert value == int(x == v and y == u)
            checks += 1
    return {'quartic': quartic, 'index_energy_increase': delta,
            'local_entries_checked': 16 ** 4,
            'strictly_positive_endpoint_entries': positive_entries,
            'two_copy_type_count': len(word_types), 'typed_entries_checked': checks}


def rational_control(gaussian, quartic_values, quartic, word_types):
    a, b, c, d = quartic
    delta = b + d - a - c
    assert delta > 3
    # One or two quartic factors require at least two or three internal
    # Gaussian edges. Check the stronger divisibility before canceling
    # the corresponding powers of c=1+w.
    c_powers = (1, 3, 5, 15)
    counts = [0, 0, 0]
    typed_checks = 0
    for words in word_types.values():
        for x, y, u, v in product(words, repeat=4):
            first = index(x[0], y[0], u[0], v[0])
            second = index(x[1], y[1], u[1], v[1])
            terms = ((1, polynomial_multiply(quartic_values[first], gaussian[second])),
                     (1, polynomial_multiply(gaussian[first], quartic_values[second])),
                     (2, polynomial_multiply(quartic_values[first], quartic_values[second])))
            endpoint = (gaussian[first].bit_count() & 1) * (gaussian[second].bit_count() & 1)
            for number, value in terms:
                if not value:
                    continue
                polynomial_exact_divide(value, c_powers[number + 1])
                quotient = polynomial_exact_divide(value, c_powers[number])
                assert quotient.bit_count() % 2 == 0
                endpoint ^= quotient.bit_count() & 1
                counts[number] += 1
            assert endpoint == int(x == v and y == u)
            typed_checks += 1
    return {'quartic': quartic, 'index_energy_increase': delta,
            'coefficient': '1/[w(1+w)]',
            'nonzero_one_quartic_terms_checked': counts[1],
            'nonzero_two_quartic_terms_checked': counts[2],
            'typed_entries_checked': typed_checks,
            'termwise_denominator_cancellation': True,
            'all_perturbed_endpoint_terms_vanish': True}


def check():
    gaussian, minor = local_polynomials()
    words = fixed_mode_words()
    quartics = [(1, 2, c, 3) for c in range(4)] + [(0, 0, 0, 3), (0, 1, 0, 3)]
    records = []
    for quartic in quartics:
        values = quartic_polynomials(minor, quartic)
        records.append(endpoint_control(gaussian, values, quartic, words))
    rational_quartic = (0, 1, 0, 3)
    rational = rational_control(gaussian, quartic_polynomials(minor, rational_quartic),
                                rational_quartic, words)
    # Zero energy increase is not enough: this diagonal endpoint is canceled.
    invalid = quartic_polynomials(minor, (0, 0, 0, 0))
    witness = index(1, 1, 1, 1)
    assert (gaussian[witness] ^ invalid[witness]).bit_count() % 2 == 0
    return {
        'status': 'exact projected-endpoint and rational-cancellation controls; no claimed exponent',
        'field': 2, 'modes': [0, 1, 2, 3], 'ambient': 4,
        'bare_w_minus_two_controls': records, 'rational_control': rational,
        'zero_energy_negative_control': {'quartic': [0, 0, 0, 0], 'labels': [1, 1, 1, 1],
                                         'actual_endpoint': 0, 'required_endpoint': 1},
        'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest()
                          for name in ('tools/check_positive_energy_quartic_endpoints.py',
                                       'tools/check_schur_downset_crossing.py',
                                       'tools/weighted_smith_oracle.py')},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    rendered = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end='')


if __name__ == '__main__':
    main()
