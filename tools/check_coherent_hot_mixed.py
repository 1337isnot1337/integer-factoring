#!/usr/bin/env python3
"""Literal mixed-occupation Pfaffian replay and complete parity profiles.

A five-mode physical retraction mixes the equal-height words 15 and 24,
which have four and two particles. Its two cuts are reconstructed directly
from the interval-16 alternating kernel and compared with the native
compound reduction. The raw polynomial infinity restoration is checked.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_block_smith_crossings import reverse_polynomial, smith_profile
from check_cancellation_module_automata import pmul
from check_coherent_q32_pressure import displayed_maps
from coherent_hot_compounds import height, support
from coherent_hot_mixed import profiles
from eight_letter_reflection_retraction import halve_exponents, shift_polynomial
from reflection_dual_rail_gaussian import gauged_kernel


def mixed_seventeen_maps(transposed=False):
    base = [1 << word for word in range(15)] + [1 << 16]
    p = base + [(1 << 15) ^ (1 << 24)]
    l = base + [1 << 24]
    return (l, p) if transposed else (p, l)


def selected_pfaffian():
    original_zero = gauged_kernel(16, 'zero')
    original_infinity = gauged_kernel(16, 'infinity')
    indices = [8 * leg + 7 - hot for leg in range(4) for hot in range(5)]
    zero = [[halve_exponents(original_zero[i][j]) for j in indices] for i in indices]
    infinity = [[halve_exponents(original_infinity[i][j]) for j in indices] for i in indices]

    def boundary_swap(index):
        leg, mode = divmod(index, 5)
        return 5 * (0, 2, 1, 3)[leg] + mode

    assert all(infinity[i][j] == zero[boundary_swap(i)][boundary_swap(j)]
               for i, j in product(range(20), repeat=2))

    @lru_cache(None)
    def pfaffian(mask):
        if not mask:
            return 1
        if mask.bit_count() % 2:
            return 0
        first = mask & -mask
        i = first.bit_length() - 1
        rest = mask ^ first
        result, remaining = 0, rest
        while remaining:
            chosen = remaining & -remaining
            remaining ^= chosen
            j = chosen.bit_length() - 1
            if zero[i][j]:
                result ^= pmul(zero[i][j], pfaffian(rest ^ chosen))
        return result

    return pfaffian


def literal_profiles(p, l, pfaffian):
    q = len(p)
    heights = [height(support(value)[0]) for value in p]
    ps, ls = tuple(map(support, p)), tuple(map(support, l))
    zero = [[0] * (q * q) for _ in range(q * q)]
    infinity = [[0] * (q * q) for _ in range(q * q)]
    raw_entries, infinity_entries = [], []
    for x, y, u, v in product(range(q), repeat=4):
        at_zero, at_infinity = 0, 0
        for a, b, c, d in product(ls[x], ls[y], ps[u], ps[v]):
            at_zero ^= pfaffian(a | (b << 5) | (c << 10) | (d << 15))
            at_infinity ^= pfaffian(a | (c << 5) | (b << 10) | (d << 15))
        zero[q * x + y][q * u + v] = at_zero
        infinity[q * x + u][q * y + v] = at_infinity
        raw = shift_polynomial(at_zero, 64 - heights[x] - heights[y])
        assert raw.bit_count() % 2 == int(x == v and y == u)
        raw_entries.append(raw)
        infinity_entries.append((u, v, at_infinity))
    degree = max(value.bit_length() - 1 for value in raw_entries)
    for raw, (u, v, normalized) in zip(raw_entries, infinity_entries):
        restored = shift_polynomial(reverse_polynomial(raw, degree),
                                    heights[u] + heights[v] + 64 - degree)
        assert restored == normalized
    cut_profiles = []
    parity = [(heights[a] + heights[b]) % 2 for a, b in product(range(q), repeat=2)]
    for matrix in (zero, infinity):
        assert all(matrix[i][j] == 0 for i, j in product(range(q * q), repeat=2)
                   if parity[i] != parity[j])
        orders = []
        for sector in (0, 1):
            indices = [i for i, value in enumerate(parity) if value == sector]
            block = [[matrix[i][j] for j in indices] for i in indices]
            block_degree = max(value.bit_length() - 1 for row in block for value in row)
            profile, _ = smith_profile(block, block_degree)
            orders.extend(profile)
        assert len(orders) == q * q
        cut_profiles.append(sorted(orders))
    return {'zero_orders_t': cut_profiles[0], 'infinity_orders_t': cut_profiles[1],
            'literal_raw_degree': degree, 'literal_endpoint_entries': q ** 4,
            'literal_infinity_restoration_entries': q ** 4}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    counts, records = Counter(), []
    pfaffian = selected_pfaffian()
    for transposed in (False, True):
        p, l = mixed_seventeen_maps(transposed)
        computed = profiles(p, l, 5)
        literal = literal_profiles(p, l, pfaffian)
        for key in ('zero_orders_t', 'infinity_orders_t'):
            assert computed[key] == literal[key]
            counts['mixed_complete_literal_cut_comparisons'] += 1
            counts['mixed_complete_smith_orders'] += len(computed[key])
        counts['literal_endpoint_entries'] += literal['literal_endpoint_entries']
        counts['literal_infinity_restoration_entries'] += literal['literal_infinity_restoration_entries']
        records.append({'transposed': transposed, 'profile': computed, 'literal': literal})
    root = Path(__file__).resolve().parents[1]
    old = json.loads((root / 'evidence/coherent-q32-pressure.json').read_text())
    for coherent in (False, True):
        result = profiles(*displayed_maps(coherent), 8)
        expected = old['physical_maps'][int(coherent)]
        for key in ('zero_orders_t', 'infinity_orders_t'):
            assert result[key] == expected[key]
            counts['occupation_catalogue_cut_comparisons'] += 1
        records.append({'q32_coherent': coherent, 'profile': result})
    names = (Path(__file__).name, 'coherent_hot_mixed.py', 'coherent_hot_mixed.cpp',
             'coherent_hot_compounds.py', 'coherent_hot_compounds_fast.cpp',
             'check_coherent_q32_pressure.py', 'check_block_smith_crossings.py',
             'check_cancellation_module_automata.py', 'eight_letter_reflection_retraction.py',
             'reflection_dual_rail_gaussian.py', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp')
    report = {'status': 'pass', 'scope': __doc__, 'counts': dict(counts), 'records': records,
              'selected_pfaffians_evaluated': pfaffian.cache_info().currsize,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}, 'factoring_exponent_improved': False}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': 'pass', 'counts': dict(counts),
                      'selected_pfaffians_evaluated': report['selected_pfaffians_evaluated']}, indent=2))


if __name__ == '__main__':
    main()
