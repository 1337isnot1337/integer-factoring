#!/usr/bin/env python3
"""Prune the complete height-15, 55-letter constant-retraction family.

The literal coordinate core of heights below 15 supplies a field-independent
pressure floor 27/38. No deletion closure or orthogonality is imposed on the
larger candidate. Exact rational intervals certify every natural scale.
"""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import time

import tensor_hot_compounds_f16 as field
from check_tensor_coupled_quadratic_core_floor import homogeneous_evaluate


def coefficients(values):
    counts = Counter(values)
    return [counts[i] for i in range(max(counts) + 1)]


def recompute(words):
    import tensor_literal_hot_pfaffians_f16 as literal
    maps = [{word: 1} for word in words]
    coherent = field.profiles(maps, maps, 7)
    assert coherent['zero_orders_t'] == coherent['infinity_orders_t']
    orders, blocks = [], []
    for parity in (0, 1):
        pairs = [(i, j) for i, u in enumerate(words) for j, v in enumerate(words)
                 if (u.bit_count() + v.bit_count()) % 2 == parity]
        rows = [field.tensor_vector(maps[i], maps[j], 7) for i, j in pairs]
        columns = [field.tensor_vector(maps[j], maps[i], 7) for i, j in pairs]
        result = literal.compound_block(7, rows, columns, order_cap=128, cache_budget=500_000)
        orders.extend(result['orders_t'])
        blocks.append({k: v for k, v in result.items() if k != 'orders_t'})
    assert max(orders) < 128
    assert sorted(orders) == coherent['zero_orders_t']
    return {'zero_orders_t': coherent['zero_orders_t'],
            'infinity_orders_t': coherent['infinity_orders_t']}, blocks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recompute-profiles', action='store_true')
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/tensor-literal-height15-pressure-gate.json'))
    args = parser.parse_args()
    started = time.monotonic()
    words = sorted((w for w in range(128) if field.height(w) < 15),
                   key=lambda w: (field.height(w), w))
    assert len(words) == 24
    ambient = [w for w in range(1024) if field.height(w) <= 20 and w != 128]
    assert len(ambient) == 55
    if args.recompute_profiles:
        profile, blocks = recompute(words)
        recomputed = 3 * 24 * 24
    else:
        prior = json.loads(args.output.read_text())
        assert prior['passed'] and prior['protected_words'] == words
        profile, blocks, recomputed = prior['protected_profiles'], prior['literal_blocks'], 0
    assert len(profile['zero_orders_t']) == len(profile['infinity_orders_t']) == 576
    assert profile['zero_orders_t'] == profile['infinity_orders_t']
    Z, A = coefficients(profile['zero_orders_t']), coefficients(field.height(w) for w in ambient)
    assert Z[0] == Z[1] == 1 and A[:3] == [1, 1, 0] and sum(A) == 55

    # For 0<s<=4/5, A(s^38)<=1+2s^38. The binomial upper bound
    # (1+x)^27 <= 1+27x/(1-27x) then proves the required comparison.
    split = Fraction(4, 5)
    assert 53 * split**76 < 1 and 54 * split**38 < 1
    margin = 19 * (1 - 54 * split**38) - 54 * split**11
    assert margin > 0

    @lru_cache(None)
    def left(s):
        n, d = s.numerator, s.denominator
        return (homogeneous_evaluate(Z, n**27, d**27)**19,
                d**(27 * (len(Z) - 1) * 19))

    @lru_cache(None)
    def right(s):
        n, d = s.numerator, s.denominator
        return (homogeneous_evaluate(A, n**38, d**38)**27,
                d**(38 * (len(A) - 1) * 27))

    intervals, pending = [], [(split, Fraction(1), 0)]
    while pending:
        low, high, depth = pending.pop()
        ln, ld = left(low)
        rn, rd = right(high)
        if ln * rd > rn * ld:
            intervals.append((low, high, depth))
        else:
            assert depth < 23, ('unresolved interval', low, high)
            middle = (low + high) / 2
            pending.extend([(middle, high, depth+1), (low, middle, depth+1)])
    intervals.sort()
    assert intervals[0][0] == split and intervals[-1][1] == 1
    assert all(a[1] == b[0] for a, b in zip(intervals, intervals[1:]))
    names = [Path(__file__).name, 'check_tensor_coupled_quadratic_core_floor.py',
             'tensor_literal_hot_pfaffians_f16.py', 'tensor_literal_hot_pfaffians_f16.cpp',
             'tensor_hot_compounds_f16_fast.cpp', 'tensor_hot_compounds_f16.py',
             'tensor_hot_compounds_f16.cpp', 'tensor_orthogonal_scaling.py',
             'coherent_hot_compounds_f4.cpp', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp']
    report = {'passed': True, 'scope': __doc__, 'rho_lower_bound': [27, 38],
              'protected_words': words, 'protected_profiles': profile, 'literal_blocks': blocks,
              'physical_heights': sorted(field.height(w) for w in ambient),
              'profile_polynomial': Z, 'height_polynomial': A,
              'polynomial_inequality': 'Z(s^27)^19 > A(s^38)^27 for 0<s<=1',
              'small_interval': [[0, 1], [4, 5]],
              'small_interval_margin': [margin.numerator, margin.denominator],
              'rational_intervals': [[[a.numerator, a.denominator], [b.numerator, b.denominator]]
                                     for a, b, _ in intervals],
              'maximum_dyadic_depth': max(d for _, _, d in intervals),
              'smith_orders_recomputed_this_run': recomputed,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'seconds': time.monotonic() - started}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'rho_lower_bound': [27, 38],
                      'rational_intervals': len(intervals),
                      'maximum_dyadic_depth': report['maximum_dyadic_depth'],
                      'smith_orders_recomputed': recomputed, 'seconds': report['seconds']}), flush=True)


if __name__ == '__main__':
    main()
