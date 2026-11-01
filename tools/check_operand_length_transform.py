#!/usr/bin/env python3
"""Exact controls for operand-length-conditioned multiplication networks.

Full-support first transform stages are checked against ordinary integer
convolution, including asymmetric supports longer than the former half
length. Large parameter checks use integers, not actual enormous tables.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import random

import check_adaptive_transform_constructor as adaptive
from check_adaptive_transform_constructor import Parameters, TransformModel, ceil_div


def support_count(length, residue, radix):
    return max(0, (length - 1 - residue) // radix + 1)


@dataclass(frozen=True)
class LengthParameters(Parameters):
    @property
    def width_denominator(self):
        return self.a * self.a - self.a - 1

    @property
    def digit_bits(self):
        return ceil_div(self.n + 1, self.width_denominator)

    def verify(self, totals):
        p, n, a, d, w, t, m = (self.groups, self.n, self.a, self.digits,
                                self.digit_bits, self.log_digits, self.modulus_bits)
        length = a * a
        assert length == 2 * d and a >= 1024 and d >= 1 << 19
        assert 2 * d * d <= n
        assert (w - 1) * self.width_denominator < n + 1 <= w * self.width_denominator
        assert 1000 * d * w < 501 * n
        assert m % d == 0 and 2 * w + t <= m <= 2 * w + t + d - 1
        if p == 2:
            assert 4 * d * d <= n
        else:
            assert 50 * d * d < 19 * n
        wire, ordered, settings = self.volumes()
        if p == 2:
            assert 10 * ordered < 221 * n and 10 * settings < 121 * n
        else:
            assert ordered < 25 * n and settings < 14 * n
        assert 2 * wire <= ordered
        width, local, cluster_degree, cluster_settings = self.table_bounds()
        assert p * p * width <= 49 * n
        assert p * local <= 56 * n
        assert p * (local + 1) <= 57 * n
        assert p * ((p - 1) * width + local) <= 105 * n
        g, f, c = self.clusters, self.residue_width, self.carry_width
        residual_degree = 14 * a * f + 2 * f + 4 * c
        residual_settings = 2 * a * (w + 3 * f) + 2 * (f + c)
        assert g * cluster_degree <= ordered + g * residual_degree
        assert g * cluster_settings <= settings + g * residual_settings
        assert 1000 * p * g * residual_degree < 2627 * ordered
        assert 1000 * p * g * residual_settings < 2627 * settings
        supports = []
        for ell in sorted({1, 2, n // 2, n - 1, n}):
            dp, dq = ceil_div(ell, w), ceil_div(n - ell + 1, w)
            assert dp + dq <= length - a and min(dp, dq) <= d
            assert ceil_div(dp, a) + ceil_div(dq, a) <= a
            # The counts are piecewise constant in r, with breakpoints at
            # the two support remainders; these points cover every value.
            residues = {0, a - 1, dp % a, dq % a}
            residues.update(r - 1 for r in tuple(residues) if r)
            assert all(support_count(dp, r, a) + support_count(dq, r, a) <= a
                       for r in residues)
            rp, rq = ell - w * (dp - 1), n - ell + 1 - w * (dq - 1)
            assert 1 <= rp <= w and 1 <= rq <= w
            assert (ell - 1) + (n - ell + 1) == n
            totals['large_support_cases'] += 1
            supports.append({'ell': str(ell), 'left_digits': dp, 'right_digits': dq,
                             'left_top_width': str(rp), 'right_top_width': str(rq)})
        totals['large_parameter_cases'] += 1
        return {'P': p, 'n': str(n), 'A': a, 'D': d, 'digit_bits': str(w),
                'modulus_bits': str(m), 'half_transform_capacity': str(d * w),
                'ordered_bits_per_input_bit': str(Fraction(ordered, n)),
                'settings_bits_per_input_bit': str(Fraction(settings, n)),
                'clique_volume_per_input_bit': str(Fraction(p * p * width, n)),
                'scaled_entry_bits_per_input_bit': str(Fraction(p * (local + 1), n)),
                'supports': supports}


class LengthTransform(TransformModel):
    def __init__(self, h, n, ell):
        a = 2 * h
        denominator = a * a - a - 1
        if not 1 <= ell <= n or denominator < 1:
            raise ValueError('invalid target length, operand length, or transform shape')
        super().__init__(h, ceil_div(n + 1, denominator))
        self.n, self.ell = n, ell
        self.capacity = self.length * self.digit_bits
        self.left_digits = ceil_div(ell, self.digit_bits)
        self.right_bits = n - ell + 1
        self.right_digits = ceil_div(self.right_bits, self.digit_bits)
        assert self.left_digits + self.right_digits <= self.length - self.a
        assert min(self.left_digits, self.right_digits) <= self.digits
        assert self.exact_recovery

    def admits(self, left, right):
        return (1 << (self.ell - 1)) <= left < (1 << self.ell) and 0 <= right < (1 << self.right_bits)

    def digits_of(self, value):
        if not 0 <= value < 1 << self.capacity:
            raise ValueError('integer outside the full transform capacity')
        return [(value >> (i * self.digit_bits)) & (self.base - 1)
                for i in range(self.length)]

    def fa_layer(self, digits):
        # This intentionally uses A terms, replacing the old h-term input.
        a, modulus, length = self.a, self.modulus, self.length
        return [[self.powers[(r * k1) % length]
                 * sum(digits[r + a * i] * self.powers[(a * i * k1) % length]
                       for i in range(a)) % modulus
                 for k1 in range(a)] for r in range(a)]

    def ordinary_convolution(self, left, right):
        assert self.admits(left, right)
        lhs, rhs = self.digits_of(left), self.digits_of(right)
        result = [0] * self.length
        for i in range(self.left_digits):
            for j in range(self.right_digits):
                result[i + j] += lhs[i] * rhs[j]
        return result

    def verify_pair(self, left, right, totals):
        assert self.admits(left, right)
        actual = self.coefficients(left, right)
        expected = self.ordinary_convolution(left, right)
        assert actual == expected
        assert max(expected) <= min(self.left_digits, self.right_digits) * (self.base - 1) ** 2
        value, carries, valid = self.carry(actual)
        assert valid and value == left * right
        assert carries[0] == carries[-1] == 0
        assert all(0 <= carry <= self.carry_bound for carry in carries)
        assert all(support_count(self.left_digits, r, self.a)
                   + support_count(self.right_digits, r, self.a) <= self.a
                   for r in range(self.a))
        totals['exact_product_assignments'] += 1
        totals['exact_convolution_coefficients'] += self.length
        if max(self.left_digits, self.right_digits) > self.digits:
            totals['asymmetric_beyond_half_length_assignments'] += 1
        return value


def prefix_and_recovery_controls(totals):
    records = []
    for n in range(2, 10):
        counts = Counter()
        pairs = 0
        for ell in range(1, n + 1):
            model = LengthTransform(2, n, ell)
            assert not model.admits(0, 1)
            assert not model.admits(1 << ell, 1)
            assert not model.admits(1 << (ell - 1), 1 << model.right_bits)
            for left in range(1 << (ell - 1), 1 << ell):
                for right in range(1 << model.right_bits):
                    value = model.verify_pair(left, right, totals)
                    pairs += 1
                    if value.bit_length() == n:
                        for k in range(n + 1):
                            counts[value, ell, k, left >> k] += 1
        for value in range(1 << (n - 1), 1 << n):
            divisors = [d for d in range(1, value + 1) if value % d == 0]
            for k in range(n + 1):
                for prefix in range(1 << (n - k)):
                    actual = 0
                    for ell in range(1, n + 1):
                        count = counts[value, ell, k, prefix]
                        assert count == sum(d.bit_length() == ell and d >> k == prefix for d in divisors)
                        actual += count
                        totals['length_conditioned_prefix_queries'] += 1
                    assert actual == sum(d >> k == prefix for d in divisors)
                    totals['summed_prefix_queries'] += 1
            prefix = 0
            for k in reversed(range(n)):
                candidate = 2 * prefix
                count = sum(counts[value, ell, k, candidate] for ell in range(1, n + 1))
                count -= int(1 >> k == candidate)
                assert count >= 0
                prefix = candidate if count else candidate + 1
            assert prefix == min(d for d in divisors if d >= 2)
            assert all(prefix % d for d in range(2, prefix))
            totals['smallest_prime_factor_recoveries'] += 1
            totals['prime_target_recoveries'] += int(len(divisors) == 2)
        records.append({'target_bit_length': n, 'h': 2, 'conditioned_classes': n,
                        'enumerated_input_pairs': pairs,
                        'target_range': [1 << (n - 1), (1 << n) - 1]})
    # Without the exact leading bit, p=q=3 is counted in both ell=2 and
    # ell=3 when N=9, although its true first-factor bit length is two.
    duplicated = [ell for ell in range(1, 5) if 3 < 1 << ell and 3 < 1 << (4 - ell + 1)]
    assert duplicated == [2, 3]
    return records, {'target': 9, 'divisor': 3,
                     'classes_without_leading_bit_check': duplicated,
                     'unique_correct_class': 2, 'duplicate_negative_control_passed': True}


def sampled_support_controls(totals):
    rng = random.Random(82120261031)
    records = []
    for h, n in ((1, 7), (2, 9), (2, 11), (2, 37), (4, 100), (4, 257), (8, 500)):
        for ell in sorted({1, 2, n // 2, n - 1, n}):
            model = LengthTransform(h, n, ell)
            minimum, left_max, right_max = 1 << (ell - 1), (1 << ell) - 1, (1 << model.right_bits) - 1
            pairs = {(minimum, 0), (left_max, right_max), (left_max, 1)}
            pairs.update((minimum | rng.getrandbits(ell - 1), rng.getrandbits(model.right_bits))
                         for _ in range(4))
            for left, right in pairs:
                model.verify_pair(left, right, totals)
            records.append({'h': h, 'n': n, 'ell': ell, 'digit_bits': model.digit_bits,
                            'left_digits': model.left_digits, 'right_digits': model.right_digits,
                            'left_top_width': ell - model.digit_bits * (model.left_digits - 1),
                            'right_top_width': model.right_bits - model.digit_bits * (model.right_digits - 1),
                            'sampled_input_pairs': len(pairs)})
    return records


def large_parameter_controls(totals):
    records = []
    assert 50 * 2 ** 62 < 19 * 3 ** 40
    assert 1 + Fraction(221, 10) * (1 + Fraction(1, 2) + Fraction(2627, 4000)) < 49
    assert Fraction(242, 10) * (1 + Fraction(2627, 2000)) < 56
    assert 1 + 25 * (1 + Fraction(1, 3) + Fraction(2627, 9000)) < 49
    assert 28 * (1 + Fraction(2627, 3000)) < 56
    for groups in range(2, 129):
        threshold = groups ** 40
        initial = LengthParameters(threshold, groups)
        denominator, d, t = initial.width_denominator, initial.digits, initial.log_digits
        near = denominator * ceil_div(threshold + 1, denominator) - 1
        candidates = {threshold, threshold + 1, near, near + 1,
                      near + denominator - 1, near + denominator, 2 * threshold, 11 * threshold}
        block = ceil_div(2 * initial.digit_bits + t, d)
        for modulus_block in (block, block + 1):
            for difference in (-1, 1):
                width = (modulus_block * d + difference - t) // 2
                assert 2 * width + t == modulus_block * d + difference
                for n in ((width - 1) * denominator, width * denominator - 1):
                    if n >= threshold:
                        candidates.add(n)
        for n in sorted(candidates):
            record = LengthParameters(n, groups).verify(totals)
            if groups in (2, 3, 4, 8, 16, 64, 128):
                records.append(record)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    totals = Counter()
    prefix_records, negative = prefix_and_recovery_controls(totals)
    result = {'schema_version': 1, 'status': 'pass',
              'claim': 'Exact operand-length-conditioned constructor with volume <=49n and entry bits <=57n/P; factoring exponent remains one half.',
              'source_sha256': {path.name: sha256(path.read_bytes()).hexdigest() for path in
                                (Path(__file__), Path(adaptive.__file__))},
              'prefix_and_recovery_models': prefix_records,
              'highest_bit_negative_control': negative,
              'asymmetric_support_models': sampled_support_controls(totals),
              'large_parameter_controls': large_parameter_controls(totals),
              'totals': dict(sorted(totals.items()))}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['totals'], sort_keys=True))


if __name__ == '__main__':
    main()
