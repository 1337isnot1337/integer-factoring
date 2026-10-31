"""Optional finite controls for the native F4 search tool (stdlib + C++17).

The reference weighted profiles use exact Fraction weights and independent
polynomial Bareiss elimination. This finite replay validates implementation
examples; it does not make floating searches into mathematical certificates.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import random

from native_weighted_smith_f4 import Oracle, PairedOracle
from weighted_smith_oracle_f4 import WeightedSmithOracle, polynomial_multiply


def unpack(value):
    return [sum(((value >> (2 * k + plane)) & 1) << k
                for k in range((value.bit_length() + 1) // 2)) for plane in range(2)]


def pack(pair):
    return sum((((pair[0] >> k) & 1) + 2 * ((pair[1] >> k) & 1)) << (2 * k)
               for k in range(max(pair).bit_length()))


def exact_profile(values, q, heights, scale, cut):
    answer = []
    for side, current in enumerate((cut, 1 - cut)):
        n = q * q
        matrix = [[0] * n for _ in range(n)]
        for index, pair in enumerate(values):
            x, y, u, v = [(index // q ** leg) % q for leg in range(4)]
            row, col = (q * x + u, q * y + v) if current else (q * x + y, q * u + v)
            matrix[row][col] = pack(pair)
        rowlegs, collegs = ((0, 2), (1, 3)) if current else ((0, 1), (2, 3))
        rr = [sum(heights[side][q * leg + label]
                  for leg, label in zip(rowlegs, divmod(i, q))) for i in range(n)]
        cc = [sum(heights[side][q * leg + label]
                  for leg, label in zip(collegs, divmod(i, q))) for i in range(n)]
        profile, _ = WeightedSmithOracle(matrix, infinity=not side).profile(rr, cc, scale)
        answer.extend(profile)
    return answer


def check():
    rng = random.Random(0xF4C0)
    count, maximum = 0, 0.
    scale = Fraction(3, 7)
    for _ in range(12):
        values = []
        for index in range(16):
            x, y, u, v = [(index >> leg) & 1 for leg in range(4)]
            polynomial = polynomial_multiply(rng.randrange(1 << 10), 5) ^ int(x == v and y == u)
            values.append(unpack(polynomial))
        heights = [[Fraction(0) if i % 2 == 0 else Fraction(rng.randrange(-20, 21), 17)
                    for i in range(8)] for _ in range(2)]
        params = [float(h) for side in heights for i, h in enumerate(side) if i % 2] + [float(scale)]
        oracle = Oracle(values, 2)
        for cut in (0, 1):
            _, _, observed = oracle.evaluate(params, [.5, .5], cut, True)
            expected = exact_profile(values, 2, heights, scale, cut)
            error = max(abs(a - float(b)) for a, b in zip(observed, expected))
            assert error < 1e-10, (cut, error)
            maximum = max(maximum, error)
            count += len(expected)
            prefix, _, rows, cols = oracle.prefixes(params, cut)
            assert max(abs(a - b) for a, b in zip(prefix, observed)) < 1e-12
            assert all(rows[i].bit_count() == cols[i].bit_count() == i % 4 + 1 for i in range(8))

    # A generic point avoids ties; compare every returned gradient coordinate.
    params = [rng.uniform(-.25, .25) for _ in range(8)] + [.31]
    _, gradient = oracle.evaluate(params, [.37, .63])
    errors = []
    for i in range(9):
        plus, minus = params[:], params[:]
        plus[i] += 1e-6
        minus[i] -= 1e-6
        difference = (oracle.evaluate(plus, [.37, .63])[0] - oracle.evaluate(minus, [.37, .63])[0]) / 2e-6
        errors.append(abs(difference - gradient[i]))
    assert max(errors) < 1e-6

    # A scalar denominator shifts every infinity order, with the opposite
    # pressure/scale-gradient change. It does not establish admissibility.
    base_cost, base_gradient, base_orders = oracle.evaluate(params, [.37, .63], want_spectra=True)
    restored = Oracle(values, 2, denominator_degree=2)
    cost, gradient, orders = restored.evaluate(params, [.37, .63], want_spectra=True)
    assert abs(cost - base_cost + 2 * params[-1]) < 1e-12
    assert abs(gradient[-1] - base_gradient[-1] + 2) < 1e-12
    assert all(abs(orders[i] - base_orders[i] - (2 * params[-1] if i < 4 else 0)) < 1e-12 for i in range(8))
    try:
        oracle.evaluate(params, [.5, .5])
        raise AssertionError('Inactive oracle was silently reused.')
    except RuntimeError:
        pass

    # The full-support matrix is invertible at w=1, but these gauges hide
    # pivots beyond the finite window. The engine must fail closed.
    dense = []
    for index in range(16):
        x, y, u, v = [(index >> leg) & 1 for leg in range(4)]
        dense.append([3 ^ int(x == v and y == u), 0])
    unresolved = Oracle(dense, 2)
    bad = [0.] * 9
    bad[0], bad[-1] = 600., 1.
    try:
        unresolved.evaluate(bad, [.5, .5])
        raise AssertionError('An unresolved finite-window profile was accepted.')
    except ArithmeticError:
        pass

    # Exact F4 conjugation/product control: (1+omega*c)(1+omega^2*c)
    # equals 1+w+w^2. The shared-core construction has q=8.
    one_copy = []
    for index in range(4 ** 4):
        x, y, u, v = [(index // 4 ** leg) % 4 for leg in range(4)]
        one_copy.append([1, 3] if x == v and y == u else [0, 0])
    paired = PairedOracle(one_copy, 1, denominator_degree=0)
    paired_values = []
    for index in range(8 ** 4):
        x, y, u, v = [(index // 8 ** leg) % 8 for leg in range(4)]
        expected = [7, 0] if x == v and y == u else [0, 0]
        assert paired.entry(x, y, u, v) == expected
        paired_values.append(expected)
    paired_params = [0.] * paired.parameters
    paired_params[-1] = .5
    _, _, paired_orders = paired.evaluate(paired_params, [1 / 8] * 8, want_spectra=True)
    generic = Oracle(paired_values, 8)
    _, _, generic_orders = generic.evaluate(paired_params, [1 / 8] * 8, want_spectra=True)
    assert paired_orders == generic_orders

    crossing = []
    for index in range(16):
        x, y, u, v = [(index >> leg) & 1 for leg in range(4)]
        crossing.append([int(x == v and y == u), 0])
    oracle = Oracle(crossing, 2)
    rho, probabilities, _, _ = oracle.types([0.] * 8 + [1.])
    assert abs(rho - 1) < 1e-12 and probabilities == [.5, .5]
    return {
        'scope': 'finite implementation controls; floating search still requires exact replay',
        'exact_fraction_profile_orders': count,
        'maximum_profile_error': maximum,
        'gradient_coordinates': 9,
        'maximum_gradient_error': max(errors),
        'paired_entries_checked': 8 ** 4,
        'paired_generic_orders_compared': len(paired_orders),
        'prefix_masks_checked': count,
        'denominator_restoration': True,
        'inactive_seed_rejected': True,
        'unresolved_precision_rejected': True,
        'constant_crossing_type_pressure': rho,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('evidence/native-weighted-smith-f4.json'))
    args = parser.parse_args()
    result = check()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
