#!/usr/bin/env python3
"""Certified divisor maxima below arbitrary binary bounds.

This successor adds a computed final radius and rigorous partial-benefit
pruning to the adaptive oracle. The practical implementation uses ordinary
integers and polynomial space; the note proves a logarithmic-space variant.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import sys
from pathlib import Path

from adaptive_divisor_preprocessing import AdaptiveMaximum, BenefitSearch
from nicolas_threshold_centers import reduced_log


@dataclass
class CertifiedMaximum(AdaptiveMaximum):
    exclusive_bound: int


class CertifiedSearch(BenefitSearch):
    def __init__(self, exclusive_bound: int):
        if exclusive_bound < 3:
            raise ValueError("internal search requires an exclusive bound at least three")
        super().__init__((exclusive_bound - 1).bit_length())
        self.limit = exclusive_bound
        self.log_cache = {}
        self.increment_cache = {}

    def logarithm(self, value: int):
        # Instance-owned caches are released with this query.
        if value not in self.log_cache:
            exponent = value.bit_length() - 1
            shift = max(0, exponent - self.precision - 2)
            top = value >> shift
            estimate = exponent * self.log_two + reduced_log(
                top, 1 << (exponent - shift), self.precision)
            error = (exponent + 1) * (8 * self.precision + 5)
            self.log_cache[value] = estimate, error
        return self.log_cache[value]

    def increment_log(self, exponent: int):
        if exponent not in self.increment_cache:
            self.increment_cache[exponent] = (
                reduced_log(exponent + 1, exponent, self.precision),
                8 * self.precision + 5)
        return self.increment_cache[exponent]

    def required_radius_upper(self, best_divisors: int):
        ln, en = self.logarithm(self.superior_value)
        ld, ed = self.logarithm(best_divisors)
        lt, et = self.logarithm(self.superior_divisors)
        lx, ex = self.logarithm(self.limit)
        numerator = (self.numerator * (lx + ex - ln)
                     + self.denominator * (lt + et - ld))
        return numerator, self.denominator * self.scale

    def benefit_prefixes(self, max_height, max_level):
        """Certified lower costs of exact and still-active prime exponents."""
        reference = []
        remaining = self.superior_value
        for i in range(max_height):
            p = self.prime(i)
            exponent = 0
            while remaining % p == 0:
                remaining //= p
                exponent += 1
            reference.append(exponent)
        assert remaining == 1
        exact, active = [], []
        for level in range(max_level + 1):
            exact_prefix = [0]
            active_prefix = [0]
            le, ee = self.logarithm(level + 1)
            for i, optimal in enumerate(reference):
                lo, eo = self.logarithm(optimal + 1)
                lp, ep = self.logarithm(self.prime(i))
                difference = level - optimal
                lower = (self.denominator * (lo - le - ee)
                         + self.numerator * difference *
                         (lp if difference >= 0 else lp + ep))
                cost = max(0, lower)
                exact_prefix.append(exact_prefix[-1] + cost)
                # If the unconstrained optimum is still allowed, zero is
                # the minimum possible future cost for this active prime.
                active_prefix.append(active_prefix[-1] +
                                     (cost if optimal < level else 0))
            exact.append(exact_prefix)
            active.append(active_prefix)
        return exact, active

    def enumerate_window(self, windows, best_value, best_divisors,
                         best_exponents, radius_numerator, radius_denominator):
        if any(lo > hi for lo, hi in windows):
            return best_value, best_divisors, best_exponents, 0, 0
        max_height = windows[0][1]
        primorial = [1]
        for i in range(max_height):
            primorial.append(primorial[-1] * self.prime(i))
        exact_cost, active_cost = self.benefit_prefixes(max_height, len(windows))
        right = radius_numerator * self.denominator * self.scale
        heights = []
        leaves = nodes = 0

        def visit(level, previous, value, divisors, next_prime, frozen_cost):
            nonlocal best_value, best_divisors, best_exponents, leaves, nodes
            nodes += 1
            if level > len(windows):
                leaves += 1
                if divisors > best_divisors or (divisors == best_divisors and value < best_value):
                    best_value, best_divisors = value, divisors
                    best_exponents = tuple(sum(h > i for h in heights)
                                           for i in range(heights[0] if heights else 0))
                return
            lower, upper = windows[level - 1]
            upper = min(upper, previous)
            if next_prime is not None:
                power = (level + 1) // 2
                cap = 0
                while cap < upper and self.prime(cap) ** power < next_prime:
                    cap += 1
                upper = min(upper, cap)
            for height in range(lower, upper + 1):
                product = value * primorial[height]
                if product >= self.limit:
                    break
                if level < len(windows) and height < windows[level][0]:
                    continue
                if level == 1:
                    # Every prime outside the chosen support has exponent zero.
                    frozen_new = exact_cost[0][-1] - exact_cost[0][height]
                else:
                    # Primes which just left the active set have final exponent
                    # level-1; all earlier frozen exponents remain unchanged.
                    frozen_new = (frozen_cost + exact_cost[level - 1][previous]
                                  - exact_cost[level - 1][height])
                lower_benefit = frozen_new + active_cost[level][height]
                if lower_benefit * radius_denominator > right:
                    continue
                divisor_power = level ** height
                assert divisors % divisor_power == 0
                updated = (divisors // divisor_power) * ((level + 1) ** height)
                heights.append(height)
                visit(level + 1, height, product, updated,
                      self.prime(height) if level == 1 else next_prime, frozen_new)
                heights.pop()

        visit(1, self.n - 1, 1, 1, None, 0)
        return best_value, best_divisors, best_exponents, leaves, nodes

    def run(self):
        initial = self.select_epsilon()
        value, divisors, exponents = self.superior_value, self.superior_divisors, initial
        # The coarse initial radius is unchanged. A fine dyadic grid permits
        # a certified final jump without a factor-two radius overshoot.
        radius_denominator = (1 << ((self.n + 1).bit_length())) * (1 << 16)
        radius_numerator = 1 << 16
        stages = []
        total_leaves = total_nodes = 0
        while True:
            windows = self.height_windows(radius_numerator, radius_denominator)
            value, divisors, exponents, leaves, nodes = self.enumerate_window(
                windows, value, divisors, exponents,
                radius_numerator, radius_denominator)
            total_leaves += leaves
            total_nodes += nodes
            upper, denominator = self.required_radius_upper(divisors)
            certified = upper * radius_denominator <= radius_numerator * denominator
            box_size = 1
            for lo, hi in windows:
                box_size *= max(0, hi - lo + 1)
            stages.append({"radius_numerator": radius_numerator,
                           "radius_denominator": radius_denominator,
                           "height_windows": windows, "box_size": box_size,
                           "visited_feasible_profiles": leaves, "search_nodes": nodes,
                           "best_divisors": divisors, "certified": certified})
            if certified:
                break
            # Round the required upper bound upward, with one extra grid unit.
            # If this step is smaller than doubling, the next stage must finish.
            final_jump = ((upper * radius_denominator + denominator - 1)
                          // denominator + 1)
            radius_numerator = min(2 * radius_numerator, final_jump)
        return CertifiedMaximum(
            self.n, value, divisors, (divisors - 1).bit_length(), exponents,
            self.numerator, self.denominator, self.superior_value,
            self.superior_divisors, stages, total_leaves, total_nodes,
            self.fallbacks, self.limit)


def maximum_below(exclusive_bound: int) -> CertifiedMaximum:
    if exclusive_bound < 2:
        raise ValueError("the exclusive bound must be at least two")
    if exclusive_bound == 2:
        return CertifiedMaximum(1, 1, 1, 0, (), 1, 1, 1, 1, [], 1, 1, 0, 2)
    return CertifiedSearch(exclusive_bound).run()


def is_highly_composite(value: int) -> bool:
    """Strict divisor-record convention; 1 is highly composite."""
    if value < 1:
        return False
    return maximum_below(value + 1).least_witness == value


def main():
    # This CLI intentionally accepts and emits arbitrary-size exact integers.
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--bits', type=int, help='use the exclusive bound 2^n')
    group.add_argument('--below', type=int, help='an arbitrary exclusive integer bound')
    group.add_argument('--is-record', type=int, help='test strict highly-composite membership')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.is_record is not None:
        result = {'integer': args.is_record, 'highly_composite': is_highly_composite(args.is_record)}
    else:
        if args.bits is not None and args.bits < 1:
            parser.error('--bits must be positive')
        bound = 1 << args.bits if args.bits is not None else args.below
        result = asdict(maximum_below(bound))
    text = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
