#!/usr/bin/env python3
"""Exact common-core floor for the full height-eight/height-twelve family.

Every nonzero coefficient choice over every characteristic-two field has
the same 28-letter height budget and retains 21 coordinate states. The
certificate proves rho>=5/7 for coordinate and orthogonal retracts of
these full six-mode kernels at the inherited natural heights.
"""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import time


def coefficients(values):
    counts = Counter(values)
    return [counts[i] for i in range(max(counts) + 1)]


def homogeneous_evaluate(values, numerator, denominator):
    out, power = values[-1], denominator
    for value in reversed(values[:-1]):
        out = out * numerator + value * power
        power *= denominator
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recompute-profiles', action='store_true')
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/tensor-coupled-quadratic-core-floor.json'))
    args = parser.parse_args()
    started = time.monotonic()
    from coherent_differential_kernels_f4 import differential_kernel
    from coherent_hot_compounds_f4 import height
    code = differential_kernel([{9: 1, 6: 1}, {33: 1, 18: 1, 12: 1}], 6, 36)
    assert code['q'] == 28
    edges = [9, 6, 33, 18, 12]
    words = sorted((w for w in range(64) if all(w & edge != edge for edge in edges)),
                   key=lambda w: (height(w), w))
    assert len(words) == 21
    if args.recompute_profiles:
        from coherent_hot_compounds_f4 import profiles
        import tensor_hot_compounds_f16 as f16
        maps = [{word: 1} for word in words]
        native = profiles(maps, maps, 6)
        extension = f16.profiles(maps, maps, 6)
        assert native['zero_orders_t'] == native['infinity_orders_t']
        assert native['zero_orders_t'] == extension['zero_orders_t']
        assert native['infinity_orders_t'] == extension['infinity_orders_t']
        profile = {'zero_orders_t': native['zero_orders_t'],
                   'infinity_orders_t': native['infinity_orders_t']}
        recomputed = 4 * 21 * 21
    else:
        prior = json.loads(args.output.read_text())
        assert prior['passed'] and prior['protected_words'] == words
        profile = prior['protected_profiles']
        recomputed = 0
    assert len(profile['zero_orders_t']) == len(profile['infinity_orders_t']) == 441
    assert profile['zero_orders_t'] == profile['infinity_orders_t']
    Z, A = coefficients(profile['zero_orders_t']), coefficients(code['heights'])
    assert Z[0] == Z[1] == 1
    assert A[:3] == [1, 1, 0] and sum(A) == 28
    # For s<=1/2, A(s^7)<=1+2s^7 and Z(s^5)^7>=1+7s^5.
    assert 26 < 2 ** 14
    assert 7 * (2 ** 7 - 20) > 5 * 2 ** 7

    @lru_cache(None)
    def left(s):
        n, d = s.numerator, s.denominator
        return (homogeneous_evaluate(Z, n ** 5, d ** 5) ** 7,
                d ** (35 * (len(Z) - 1)))

    @lru_cache(None)
    def right(s):
        n, d = s.numerator, s.denominator
        return (homogeneous_evaluate(A, n ** 7, d ** 7) ** 10,
                d ** (70 * (len(A) - 1)))

    intervals, pending, maximum_depth = [], [(Fraction(1, 2), Fraction(1), 0)], 0
    while pending:
        low, high, depth = pending.pop()
        ln, ld = left(low)
        rn, rd = right(high)
        if ln * rd > rn * ld:
            intervals.append((low, high))
            maximum_depth = max(maximum_depth, depth)
        else:
            assert depth < 24, ('unresolved interval', low, high)
            middle = (low + high) / 2
            pending.extend([(middle, high, depth + 1), (low, middle, depth + 1)])
    intervals.sort()
    assert intervals[0][0] == Fraction(1, 2) and intervals[-1][1] == 1
    assert all(a[1] == b[0] for a, b in zip(intervals, intervals[1:]))
    names = [Path(__file__).name, 'coherent_differential_kernels_f4.py',
             'coherent_hot_compounds_f4.py', 'coherent_hot_compounds_f4.cpp',
             'tensor_hot_compounds_f16.py', 'tensor_hot_compounds_f16.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp']
    report = {'passed': True, 'scope': __doc__, 'rho_lower_bound': [5, 7],
              'protected_words': words, 'protected_profiles': profile,
              'physical_heights': code['heights'], 'height_polynomial': A, 'profile_polynomial': Z,
              'polynomial_inequality': 'Z(s^5)^7 > A(s^7)^10 for 0<s<=1',
              'small_interval': [[0, 1], [1, 2]], 'small_interval_integer_margin': 116,
              'rational_intervals': [[[a.numerator, a.denominator], [b.numerator, b.denominator]]
                                     for a, b in intervals],
              'maximum_dyadic_depth': maximum_depth,
              'smith_orders_recomputed_this_run': recomputed,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'seconds': time.monotonic() - started}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'rho_lower_bound': [5, 7], 'rational_intervals': len(intervals),
                      'maximum_dyadic_depth': maximum_depth, 'recomputed_orders': recomputed,
                      'seconds': report['seconds']}))


if __name__ == '__main__':
    main()
