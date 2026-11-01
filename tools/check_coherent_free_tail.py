#!/usr/bin/env python3
"""Check the exact regular free-tail reduction and determinant identity.

Finite polynomial controls support an arbitrary dyadic-tail theorem in
research/coherent-free-tail-reduction.md. No finite scan substitutes for
that proof, and no factoring exponent improvement is claimed.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from check_block_smith_crossings import smith_profile
from check_cancellation_module_automata import pmul
from coherent_hot_compounds_fast import profiles
from search_coherent_hot_retractions import _atomic_write


ROOT = Path(__file__).resolve().parents[1]


def transpose(a):
    return list(map(list, zip(*a)))


def product(a, b):
    out = [[0] * len(b[0]) for _ in a]
    for i, row in enumerate(a):
        for j, value in enumerate(row):
            if value:
                for h, other in enumerate(b[j]):
                    if other:
                        out[i][h] ^= pmul(value, other)
    return out


def hot_matrix(modes):
    return [[(3 << (2 * ((i % modes) | (j % modes)) + 1))
             if i // modes == j // modes else
             ((1 << (2 * (i % modes) + 1)) if i % modes == j % modes else 0)
             for j in range(2 * modes)] for i in range(2 * modes)]


def tail_change(k, m):
    """Rows ordered as core, aggregate tail, then the tail complement."""
    assert k & (k - 1) == m & (m - 1) == 0 and m >= 2 * k
    n = m // k - 1
    rows = []
    for side in (0, 1):
        for c in range(k):
            row = [0] * (2 * m)
            row[side * m + c] = 1
            rows.append(row)
    for side in (0, 1):
        for c in range(k):
            row = [0] * (2 * m)
            for h in range(1, n + 1):
                row[side * m + c + k * h] = 1
            rows.append(row)
    for side in (0, 1):
        for c in range(k):
            for i in range(1, n):
                row = [0] * (2 * m)
                row[side * m + c + k] = 1 << (2 * k * i)
                row[side * m + c + k * (i + 1)] = 1
                rows.append(row)
    assert len(rows) == 2 * m
    return rows


def active_matrix(k, m):
    b = hot_matrix(k)
    alpha = sum(1 << (2 * k * h) for h in range(1, m // k))
    a = [[0] * (4 * k) for _ in range(4 * k)]
    for i in range(2 * k):
        for j in range(2 * k):
            a[i][j] = b[i][j]
            a[2 * k + i][2 * k + j] = pmul(alpha, b[i][j])
            if i // k == j // k:
                a[i][2 * k + j] = a[2 * k + i][j] = pmul(alpha, b[i][j])
    return a


def append_full_tail(vectors, k, m):
    return [sum(1 << (word | (tail << k))
                for word in range(1 << k) if vector >> word & 1)
            for tail in range(1 << (m - k)) for vector in vectors]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    counts, reductions, moments = Counter(), [], []
    for k in (1, 2, 4, 8):
        for m in (2 * k, 4 * k, 8 * k):
            if m > 32:
                continue
            change = tail_change(k, m)
            transformed = product(product(change, hot_matrix(m)), transpose(change))
            active = active_matrix(k, m)
            for i in range(4 * k):
                for j in range(4 * k):
                    assert transformed[i][j] == active[i][j]
                    counts['active_polynomial_entries'] += 1
                for j in range(4 * k, 2 * m):
                    assert transformed[i][j] == transformed[j][i] == 0
                    counts['vanishing_cross_block_entries'] += 2
            n = m // k - 1
            g = sum(1 << (2 * k * h) for h in range(n))
            assert g & 1
            for i in range(n):
                row = change[2 * k] if i == 0 else change[4 * k + i - 1]
                value = 0
                for h in range(n):
                    value ^= pmul(row[k * (h + 1)], 1 << (2 * k * h))
                assert value == (g if i == 0 else 0)
                counts['tail_vector_identities'] += 1
            degree = max(value.bit_length() - 1 for row in active for value in row)
            current, precision = smith_profile(active, degree)
            expected = list(range(1, 4 * k, 2)) + list(range(2 * m + 1, 2 * m + 4 * k, 2))
            assert list(current) == expected
            counts['active_one_particle_orders'] += len(current)
            spectator = [row[4 * k:] for row in transformed[4 * k:]]
            expected_spectator = list(range(4 * k + 1, 2 * m, 2)) + list(range(2 * m + 4 * k + 1, 4 * m, 2))
            if spectator:
                degree = max(value.bit_length() - 1 for row in spectator for value in row)
                current_spectator, spectator_precision = smith_profile(spectator, degree)
                assert list(current_spectator) == expected_spectator
                counts['spectator_one_particle_orders'] += len(current_spectator)
            else:
                spectator_precision = 0
            assert sorted(expected + expected_spectator) == list(range(1, 4 * m, 2))
            reductions.append({'core_modes': k, 'total_modes': m,
                               'active_orders': expected, 'spectator_orders': expected_spectator,
                               'active_precision': precision, 'spectator_precision': spectator_precision})
    for k, words, m in ((1, [0], 4), (1, [0, 1], 4),
                         (2, [0, 1], 4), (2, [0, 1, 2], 4),
                         (2, [0, 1, 2, 3], 4)):
        p = [1 << word for word in words]
        core = profiles(p, p, k)
        extended = append_full_tail(p, k, m)
        full = profiles(extended, extended, m)
        q, tail = len(p), m - k
        occupation_sum = sum(word.bit_count() for word in words)
        row = {'core_modes': k, 'total_modes': m, 'core_words': words,
               'core_q': q, 'full_q': len(extended), 'cuts': []}
        for cut in ('zero', 'infinity'):
            observed = sum(full[cut + '_orders_t'])
            predicted = 4 ** tail * (sum(core[cut + '_orders_t'])
                                    + 4 * tail * q * occupation_sum + 2 * m * tail * q * q)
            assert observed == predicted
            counts['complete_compressed_profiles'] += 2
            counts['complete_compressed_orders'] += q * q + len(extended) ** 2
            row['cuts'].append({'cut': cut, 'core_orders': core[cut + '_orders_t'],
                                'full_orders': full[cut + '_orders_t'],
                                'order_sum': observed})
        moments.append(row)
    source = ROOT / 'evidence/coherent-hot-q32-retractions.json'
    records = json.loads(source.read_text())['records']
    coherent, coordinate = records[55], records[56]
    assert coherent['heights'] == coordinate['heights'] and coherent['charges'] == coordinate['charges']
    difference = sum(sum(coherent[cut + '_orders_t']) - sum(coordinate[cut + '_orders_t'])
                     for cut in ('zero', 'infinity'))
    assert difference == 540
    names = (Path(__file__).name, 'check_block_smith_crossings.py',
             'check_cancellation_module_automata.py', 'coherent_hot_compounds.py',
             'coherent_hot_compounds_fast.py', 'coherent_hot_compounds_fast.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp')
    report = {'scope': __doc__, 'completed': True, 'counts': dict(counts),
              'regular_reductions': reductions, 'determinant_controls': moments,
              'q32_normalized_two_cut_moment_difference': [135, 256],
              'q32_catalogue_sha256': sha256(source.read_bytes()).hexdigest(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'infinite_tail_claim_proved_in_note': True,
              'finite_checks_alone_prove_uniform_theorem': False,
              'factoring_exponent_improved': False}
    _atomic_write(args.output, report)
    print(json.dumps({'completed': True, 'counts': dict(counts)}), flush=True)


if __name__ == '__main__':
    main()
