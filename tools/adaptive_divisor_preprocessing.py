#!/usr/bin/env python3
"""Certified adaptive-benefit search for exact maximum divisor counts.

The algorithm contains no numerical Nicolas constants. Correctness is
certified by a computed benefit bound at every input. The polynomial-time
asymptotic proof uses Nicolas's benefit and gap theorems separately.
This practical Python version uses ordinary integers and polynomial space;
the logarithmic-workspace implementation is described in the research note.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from functools import lru_cache
import json
from pathlib import Path

from nicolas_threshold_centers import reduced_log


@dataclass
class AdaptiveMaximum:
    input_bits: int
    least_witness: int
    maximum_divisors: int
    divisor_parameter_bits: int
    exponents: tuple[int, ...]
    epsilon_numerator: int
    epsilon_denominator: int
    superior_integer: int
    superior_divisors: int
    stages: list[dict]
    visited_profiles: int
    total_search_nodes: int
    exact_power_fallbacks: int


class BenefitSearch:
    def __init__(self, n: int):
        if n < 2:
            raise ValueError("internal search requires n>=2")
        self.n = n
        self.limit = 1 << n
        self.bound = 4 * n * n
        self.prime_log_bound = (self.bound - 1).bit_length()
        self.max_levels = min(n - 1, 2 * self.prime_log_bound - 2)
        q_min = max((n + 2) ** 2,
                    self.prime_log_bound * (self.prime_log_bound + 1))
        self.denominator = 1 << (q_min - 1).bit_length()
        self.precision = 32 * ((n + 2).bit_length() + 10)
        self.scale = 1 << self.precision
        self.log_two = reduced_log(2, 1, self.precision)
        self.primes = [2]
        self.fallbacks = 0

    def prime(self, index: int) -> int:
        while len(self.primes) <= index:
            candidate = self.primes[-1] + 1
            while True:
                valid = True
                for p in self.primes:
                    if p * p > candidate:
                        break
                    if candidate % p == 0:
                        valid = False
                        break
                if valid:
                    self.primes.append(candidate)
                    break
                candidate += 1
        return self.primes[index]

    @lru_cache(maxsize=None)
    def logarithm(self, value: int) -> tuple[int, int]:
        """Lower logarithm and certified error in scaled integer units.

        Only the leading precision+3 bits are used. This is also the
        representation used by the logarithmic-space proof with a virtual
        product-bit oracle.
        """
        if value < 1:
            raise ValueError("logarithm requires a positive integer")
        exponent = value.bit_length() - 1
        shift = max(0, exponent - self.precision - 2)
        top = value >> shift
        estimate = exponent * self.log_two + reduced_log(
            top, 1 << (exponent - shift), self.precision)
        error = (exponent + 1) * (8 * self.precision + 5)
        return estimate, error

    @lru_cache(maxsize=None)
    def increment_log(self, exponent: int) -> tuple[int, int]:
        """log((e+1)/e), with e>=1."""
        return (reduced_log(exponent + 1, exponent, self.precision),
                8 * self.precision + 5)

    def profitable_increment(self, p: int, exponent: int, numerator: int) -> bool:
        """Whether ((e+1)/e)^Q >= p^numerator, exactly."""
        lp, ep = self.logarithm(p)
        lr, er = self.increment_log(exponent)
        lo = self.denominator * lr - numerator * (lp + ep)
        hi = self.denominator * (lr + er) - numerator * lp
        if lo >= 0:
            return True
        if hi < 0:
            return False
        self.fallbacks += 1
        return ((exponent + 1) ** self.denominator
                >= exponent ** self.denominator * p ** numerator)

    def superior(self, numerator: int, stop_at_limit: bool = True):
        """Largest superior maximizer for epsilon=numerator/Q."""
        value = divisors = 1
        exponents = []
        index = 0
        while True:
            p = self.prime(index)
            if not self.profitable_increment(p, 1, numerator):
                break
            # epsilon>=1/(L+1), hence exponent < 2(L+1).
            lo, hi = 1, 2 * (self.prime_log_bound + 2) + 1
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if self.profitable_increment(p, mid, numerator):
                    lo = mid
                else:
                    hi = mid
            exponent = lo
            value *= p ** exponent
            divisors *= exponent + 1
            exponents.append(exponent)
            if stop_at_limit and value >= self.limit:
                return value, divisors, tuple(exponents)
            index += 1
        return value, divisors, tuple(exponents)

    def select_epsilon(self):
        q = self.denominator
        lower = (q + self.prime_log_bound) // (self.prime_log_bound + 1)
        upper = q
        assert self.superior(lower)[0] >= self.limit
        assert self.superior(upper)[0] < self.limit
        while upper - lower > 1:
            mid = (lower + upper) // 2
            if self.superior(mid)[0] >= self.limit:
                lower = mid
            else:
                upper = mid
        self.numerator = upper
        self.superior_value, self.superior_divisors, initial = self.superior(upper)
        assert self.superior_value < self.limit
        assert self.superior(upper - 1)[0] >= self.limit
        return initial

    def height_windows(self, radius_numerator: int, radius_denominator: int):
        windows = []
        right = radius_numerator * self.denominator * self.scale
        for level in range(1, self.max_levels + 1):
            lr, er = self.increment_log(level)
            mandatory = possible = 0
            index = 0
            while True:
                p = self.prime(index)
                if p > self.bound:
                    break
                lp, ep = self.logarithm(p)
                lo = self.numerator * lp - self.denominator * (lr + er)
                hi = self.numerator * (lp + ep) - self.denominator * lr
                if lo * radius_denominator > right:
                    # The true marginal cost is increasing with p, so all
                    # subsequent primes are forbidden as well.
                    break
                possible += 1
                if hi * radius_denominator < -right:
                    mandatory += 1
                index += 1
            windows.append((mandatory, min(self.n - 1, possible)))
        return windows

    def required_radius_upper(self, best_divisors: int):
        ln, en = self.logarithm(self.superior_value)
        ld, ed = self.logarithm(best_divisors)
        lt, et = self.logarithm(self.superior_divisors)
        lx = self.n * self.log_two
        ex = self.n * (8 * self.precision + 5)
        # Positive terms are rounded upward, negative terms downward.
        numerator = (self.numerator * (lx + ex - ln)
                     + self.denominator * (lt + et - ld))
        return numerator, self.denominator * self.scale

    def enumerate_window(self, windows, best_value, best_divisors, best_exponents):
        if any(lo > hi for lo, hi in windows):
            return best_value, best_divisors, best_exponents, 0, 0
        max_height = windows[0][1]
        primorial = [1]
        for i in range(max_height):
            primorial.append(primorial[-1] * self.prime(i))
        heights = []
        leaves = nodes = 0

        def visit(level, previous, value, divisors, next_prime):
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
            # Necessary exchange bound for every strict divisor record.
            if next_prime is not None:
                power = (level + 1) // 2
                cap = 0
                while cap < upper and self.prime(cap) ** power < next_prime:
                    cap += 1
                upper = min(upper, cap)
            if lower > upper:
                return
            for height in range(lower, upper + 1):
                product = value * primorial[height]
                if product >= self.limit:
                    break
                if level < len(windows) and height < windows[level][0]:
                    continue
                updated = (divisors // (level ** height)) * ((level + 1) ** height)
                assert divisors % (level ** height) == 0
                heights.append(height)
                visit(level + 1, height, product, updated,
                      self.prime(height) if level == 1 else next_prime)
                heights.pop()

        visit(1, self.n - 1, 1, 1, None)
        return best_value, best_divisors, best_exponents, leaves, nodes

    def run(self):
        initial = self.select_epsilon()
        value, divisors, exponents = self.superior_value, self.superior_divisors, initial
        radius_denominator = 1 << ((self.n + 1).bit_length())
        radius_numerator = 1
        stages = []
        total_leaves = total_nodes = 0
        while True:
            windows = self.height_windows(radius_numerator, radius_denominator)
            value, divisors, exponents, leaves, nodes = self.enumerate_window(
                windows, value, divisors, exponents)
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
            radius_numerator *= 2
        return AdaptiveMaximum(
            self.n, value, divisors, (divisors - 1).bit_length(), exponents,
            self.numerator, self.denominator, self.superior_value,
            self.superior_divisors, stages, total_leaves, total_nodes, self.fallbacks)


def adaptive_maximum(n: int) -> AdaptiveMaximum:
    if n < 1:
        raise ValueError("input bit bound must be positive")
    if n == 1:
        return AdaptiveMaximum(1, 1, 1, 0, (), 1, 1, 1, 1, [], 1, 1, 0)
    return BenefitSearch(n).run()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_bits", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = asdict(adaptive_maximum(args.input_bits))
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
