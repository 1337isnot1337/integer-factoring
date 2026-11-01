#!/usr/bin/env python3
"""Exact controls for the adaptive multiplication-network constructor.

The transform is implemented from its separate FA, MID, and IB formulas.
Its outputs are compared with independent ordinary integer convolution.
Large-cutoff checks use symbolic integer parameter sizes, not enormous
written tables. No new factoring exponent is asserted.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
import json
from pathlib import Path
import random


def ceil_div(a: int, b: int) -> int:
    return (a + b - 1) // b


def next_power_two(value: int) -> int:
    return 1 << (value - 1).bit_length()


@dataclass(frozen=True)
class Parameters:
    n: int
    groups: int

    def __post_init__(self) -> None:
        if self.groups < 2 or self.n < self.groups ** 40:
            raise ValueError("requires P>=2 and P^40<=n")

    @property
    def h(self) -> int:
        return next_power_two(self.groups ** 9)

    @property
    def a(self) -> int:
        return 2 * self.h

    @property
    def digits(self) -> int:
        return 2 * self.h ** 2

    @property
    def digit_bits(self) -> int:
        return ceil_div(self.n, self.digits)

    @property
    def log_digits(self) -> int:
        return self.digits.bit_length() - 1

    @property
    def modulus_bits(self) -> int:
        return self.digits * ceil_div(2 * self.digit_bits + self.log_digits, self.digits)

    @property
    def residue_width(self) -> int:
        return self.modulus_bits + 1

    @property
    def carry_width(self) -> int:
        return self.digit_bits + self.log_digits

    @property
    def capacity(self) -> int:
        return self.digits * self.digit_bits

    @property
    def clusters(self) -> int:
        return 3 * self.groups ** 8

    def volumes(self) -> tuple[int, int, int]:
        a, w, f, c = self.a, self.digit_bits, self.residue_width, self.carry_width
        wire_bits = 4 * a * a * f + (a * a - 1) * c
        ordered_bound = a * a * (8 * f + 2 * c)
        setting_bound = a * a * (w + 4 * f + c)
        return wire_bits, ordered_bound, setting_bound

    def table_bounds(self) -> tuple[int, int, int, int]:
        a, w, f, c, g = self.a, self.digit_bits, self.residue_width, self.carry_width, self.clusters
        ca, ct = a // g + 2, a * a // g + 2
        degree = 7 * ca * a * f + ct * (f + 2 * c)
        settings = ca * a * (w + 3 * f) + ct * (f + c)
        ordered = self.volumes()[1]
        width = max(1, (self.groups * ordered + g * degree) // self.groups ** 3)
        local_bits = settings * (2 * g // self.groups)
        return width, local_bits, degree, settings

    def verify(self, totals: Counter) -> dict:
        p, n, h, d, w, t, m = (self.groups, self.n, self.h, self.digits,
                               self.digit_bits, self.log_digits, self.modulus_bits)
        assert p ** 9 <= h < 2 * p ** 9
        assert d == 1 << t
        assert 2 * p ** 18 <= d < 8 * p ** 18
        assert 2 * d * d <= n
        assert h >= 512 and d >= 1 << 19
        assert n <= d * w < n + d
        assert 2 * w + t <= m <= 2 * w + t + d - 1
        assert m % d == 0
        wire_bits, ordered, settings = self.volumes()
        assert ordered < 45 * n
        assert settings < 25 * n
        assert 2 * wire_bits <= ordered
        width, local_bits, degree, cluster_bits = self.table_bounds()
        assert p * p * width <= 99 * n
        assert p * local_bits <= 116 * n
        assert p * (local_bits + 1) <= 117 * n
        assert p * ((p - 1) * width + local_bits) <= 215 * n
        g, a, f, c = self.clusters, self.a, self.residue_width, self.carry_width
        r = 14 * a * f + 2 * f + 4 * c
        rw = 2 * a * (w + 3 * f) + 2 * (f + c)
        assert g * degree <= ordered + g * r
        assert g * cluster_bits <= settings + g * rw
        assert 1000 * p * g * r < 2627 * ordered
        assert 1000 * p * g * rw < 2627 * settings
        totals["parameter_cases"] += 1
        return {
            "P": p, "n": str(n), "h": h, "D": d,
            "digit_bits": str(w), "modulus_bits": str(m),
            "capacity_padding": str(d * w - n),
            "residue_width": str(f), "carry_width": str(c),
            "wire_bits": str(wire_bits),
            "wire_bits_per_input_bit": str(Fraction(wire_bits, n)),
            "ordered_volume_bound": str(ordered),
            "setting_volume_bound": str(settings),
            "clique_width": str(width),
            "clique_volume_per_input_bit": str(Fraction(p * p * width, n)),
            "entry_bits": str(local_bits + 1),
            "scaled_entry_bits_per_input_bit": str(Fraction(p * (local_bits + 1), n)),
        }


class TransformModel:
    """A finite model; optional smaller m is used only for negative controls."""

    def __init__(self, h: int, digit_bits: int, modulus_bits: int | None = None):
        if h < 1 or h & (h - 1) or digit_bits < 1:
            raise ValueError("power-of-two h and positive digit width required")
        self.h, self.a, self.digits = h, 2 * h, 2 * h * h
        self.length, self.digit_bits = self.a * self.a, digit_bits
        self.log_digits = self.digits.bit_length() - 1
        self.modulus_bits = (self.digits * ceil_div(2 * digit_bits + self.log_digits, self.digits)
                             if modulus_bits is None else modulus_bits)
        if self.modulus_bits < 1 or self.modulus_bits % self.digits:
            raise ValueError("m must be a positive multiple of D")
        self.base = 1 << digit_bits
        self.modulus = (1 << self.modulus_bits) + 1
        self.root = 1 << (self.modulus_bits // self.digits)
        self.inverse_root = pow(self.root, self.length - 1, self.modulus)
        self.inverse_length = pow(2, (2 * self.modulus_bits - 1) * (self.log_digits + 1), self.modulus)
        self.carry_bound = self.digits * (self.base - 1)
        self.carry_width = digit_bits + self.log_digits
        self.capacity = self.digits * digit_bits
        self.exact_recovery = self.digits * (self.base - 1) ** 2 < 1 << self.modulus_bits
        self.powers = [pow(self.root, j, self.modulus) for j in range(self.length)]
        self.inverse_powers = [pow(self.inverse_root, j, self.modulus) for j in range(self.length)]
        assert self.powers[self.digits] == self.modulus - 1
        assert pow(self.root, self.length, self.modulus) == 1
        assert self.inverse_length * self.length % self.modulus == 1
        assert self.carry_bound + 1 <= 1 << self.carry_width

    def digits_of(self, value: int) -> list[int]:
        if not 0 <= value < 1 << self.capacity:
            raise ValueError("integer outside input capacity")
        return [(value >> (i * self.digit_bits)) & (self.base - 1)
                for i in range(self.digits)]

    def fa_layer(self, digits: list[int]) -> list[list[int]]:
        a, modulus, length = self.a, self.modulus, self.length
        result = []
        for r in range(a):
            row = []
            for k1 in range(a):
                value = sum(digits[r + a * i] * self.powers[(a * i * k1) % length]
                            for i in range(self.h))
                row.append(value * self.powers[(r * k1) % length] % modulus)
            result.append(row)
        return result

    def coefficients(self, left: int, right: int) -> list[int]:
        a, modulus, length = self.a, self.modulus, self.length
        fa_left = self.fa_layer(self.digits_of(left))
        fa_right = self.fa_layer(self.digits_of(right))
        middle = []
        for k1 in range(a):
            products = []
            for k2 in range(a):
                hat_left = sum(fa_left[r][k1] * self.powers[(a * r * k2) % length]
                               for r in range(a)) % modulus
                hat_right = sum(fa_right[r][k1] * self.powers[(a * r * k2) % length]
                                for r in range(a)) % modulus
                products.append(hat_left * hat_right % modulus)
            middle.append([
                self.inverse_powers[(k1 * j1) % length]
                * sum(products[k2] * self.inverse_powers[(a * k2 * j1) % length]
                      for k2 in range(a)) % modulus
                for j1 in range(a)
            ])
        output = [0] * length
        for j1 in range(a):
            for j2 in range(a):
                output[j1 + a * j2] = (
                    self.inverse_length
                    * sum(middle[k1][j1] * self.inverse_powers[(a * k1 * j2) % length]
                          for k1 in range(a)) % modulus)
        return output

    def ordinary_convolution(self, left: int, right: int) -> list[int]:
        lhs, rhs = self.digits_of(left), self.digits_of(right)
        result = [0] * self.length
        for i, x in enumerate(lhs):
            for j, y in enumerate(rhs):
                result[i + j] += x * y
        return result

    def carry(self, coefficients: list[int]) -> tuple[int, list[int], bool]:
        incoming, value, place, valid = 0, 0, 1, True
        carries = [0]
        for coefficient in coefficients:
            outgoing, digit = divmod(coefficient + incoming, self.base)
            valid = valid and outgoing <= self.carry_bound
            value += digit * place
            place *= self.base
            incoming = outgoing
            carries.append(incoming)
        return value, carries, valid and incoming == 0


def verify_pair(model: TransformModel, left: int, right: int, totals: Counter) -> None:
    actual = model.coefficients(left, right)
    expected = model.ordinary_convolution(left, right)
    assert actual == [value % model.modulus for value in expected]
    assert model.exact_recovery
    assert actual == expected
    assert max(expected) <= model.digits * (model.base - 1) ** 2 < 1 << model.modulus_bits
    recovered, carries, valid = model.carry(actual)
    assert valid and recovered == left * right
    assert all(0 <= value <= model.carry_bound for value in carries)
    assert carries[0] == carries[-1] == 0
    totals["exact_product_assignments"] += 1
    totals["exact_convolution_coefficients"] += model.length


def small_prefix_checks(totals: Counter) -> list[dict]:
    records = []
    for digit_bits in (1, 2, 3):
        model = TransformModel(1, digit_bits)
        limit = 1 << model.capacity
        counts = Counter()
        for left, right in product(range(limit), repeat=2):
            verify_pair(model, left, right, totals)
            value = left * right
            if 2 <= value < limit:
                for k in range(value.bit_length() + 1):
                    counts[value, k, left >> k] += 1
        for value in range(2, limit):
            divisors = [d for d in range(1, value + 1) if value % d == 0]
            n = value.bit_length()
            for k in range(n + 1):
                for prefix in range(1 << (n - k)):
                    assert counts[value, k, prefix] == sum(d >> k == prefix for d in divisors)
                    totals["divisor_prefix_queries"] += 1

            # The existing prefix-search interface excludes the two trivial
            # divisors and determines a nontrivial factor one bit at a time.
            if len(divisors) > 2:
                prefix = 0
                for k in reversed(range(n)):
                    candidate = 2 * prefix
                    count = (counts[value, k, candidate]
                             - int(1 >> k == candidate) - int(value >> k == candidate))
                    prefix = candidate if count else candidate + 1
                assert 1 < prefix < value and value % prefix == 0
                totals["prefix_factor_recoveries"] += 1
        records.append({"h": 1, "digit_bits": digit_bits,
                        "D": model.digits, "modulus_bits": model.modulus_bits,
                        "enumerated_factor_pairs": limit * limit,
                        "target_range": [2, limit - 1]})
    return records


def sampled_transform_checks(totals: Counter) -> list[dict]:
    rng = random.Random(20261031)
    records = []
    # Widths include both sides of modulus-ceiling transitions and examples
    # with the source-style D-(log2(D)+1)/2 width from the simpler redesign.
    for h, bits in ((1, 4), (2, 1), (2, 2), (2, 6), (2, 7),
                    (4, 13), (4, 14), (4, 29), (8, 124)):
        model = TransformModel(h, bits)
        limit = 1 << model.capacity
        pairs = {(0, 0), (0, limit - 1), (1, limit - 1),
                 (limit - 1, limit - 1), (model.base - 1, model.base - 1)}
        pairs.update((rng.randrange(limit), rng.randrange(limit)) for _ in range(8))
        for left, right in pairs:
            verify_pair(model, left, right, totals)
        for exponent in range(1, model.length):
            assert sum(model.powers[(exponent * j) % model.length]
                       for j in range(model.length)) % model.modulus == 0
            totals["principal_root_sums"] += 1
        records.append({"h": h, "D": model.digits, "digit_bits": bits,
                        "modulus_bits": model.modulus_bits,
                        "carry_width": model.carry_width,
                        "sampled_factor_pairs": len(pairs)})
    return records


def alias_negative_control(totals: Counter) -> dict:
    model = TransformModel(2, 6, modulus_bits=8)
    assert not model.exact_recovery
    left = right = 63
    coefficients = model.coefficients(left, right)
    expected = model.ordinary_convolution(left, right)
    assert expected[0] == 3969 and coefficients[0] == 114
    assert coefficients == [x % model.modulus for x in expected]
    target, carries, valid = model.carry(coefficients)
    assert valid and target == 114 and left * right != target
    assert target < 1 << model.capacity
    totals["false_assignment_negative_controls"] += 1
    # On arbitrary coefficient settings the bounded-message convention must
    # reject overflow even though its clipped outgoing message is legal.
    good_model = TransformModel(2, 6)
    outgoing = (good_model.modulus - 1 + good_model.carry_bound) // good_model.base
    assert outgoing > good_model.carry_bound
    clipped_message = outgoing if outgoing <= good_model.carry_bound else 0
    assert clipped_message < 1 << good_model.carry_width
    totals["carry_overflow_negative_controls"] += 1
    return {"D": model.digits, "digit_bits": model.digit_bits,
            "invalid_smaller_modulus_bits": model.modulus_bits,
            "modulus": model.modulus, "left_factor": left, "right_factor": right,
            "true_product": left * right, "falsely_accepted_target": target,
            "canonical_coefficients": coefficients, "carries": carries,
            "valid_modulus_overflow_is_rejected": True}


def parameter_checks(totals: Counter) -> list[dict]:
    records = []
    for groups in range(2, 129):
        threshold = groups ** 40
        initial = Parameters(threshold, groups)
        d = initial.digits
        near = d * ceil_div(threshold, d)
        candidates = {threshold, threshold + 1, near, near + 1,
                      near + d - 1, near + d, 2 * threshold, 11 * threshold}
        # Put 2w+t immediately before, at, and after a multiple of D when
        # parity permits it. D is even and t is odd in this construction,
        # so exact divisibility is impossible; +/-1 are the fragile cases.
        modulus_block = ceil_div(2 * initial.digit_bits + initial.log_digits, d)
        for block in (modulus_block, modulus_block + 1):
            for difference in (-1, 1):
                numerator = block * d + difference - initial.log_digits
                assert numerator % 2 == 0
                width = numerator // 2
                if width >= 1:
                    for n in (d * (width - 1) + 1, d * width):
                        if n >= threshold:
                            candidates.add(n)
        for n in sorted(candidates):
            record = Parameters(n, groups).verify(totals)
            if groups in (2, 3, 4, 8, 16, 64, 128):
                records.append(record)
    # The generic integer ceiling itself must also handle equality, even
    # though the chosen square transform shape makes t odd.
    for digits in (2, 4, 8, 16, 32):
        for target in range(1, 5 * digits + 1):
            rounded = digits * ceil_div(target, digits)
            assert target <= rounded <= target + digits - 1
            assert rounded % digits == 0
            if target % digits == 0:
                assert rounded == target
            totals["modulus_ceiling_controls"] += 1
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    totals = Counter()
    report = {
        "result": "adaptive exact-transform written top constructor",
        "scope": "smaller explicit network and table bounds; factoring exponent remains one half",
        "prefix_models": small_prefix_checks(totals),
        "transform_models": sampled_transform_checks(totals),
        "negative_controls": alias_negative_control(totals),
        "large_cutoff_parameter_cases": parameter_checks(totals),
        "summary": dict(sorted(totals.items())),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
