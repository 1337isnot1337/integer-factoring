#!/usr/bin/env python3
"""Streamed polynomial-space factors for the original F_2 quotient family.

The factor routines have no memoization or exponentially large tables.
The independent source-reference imports occur only inside verify().
"""

from itertools import combinations
import argparse
import json


def binomial_bits(exponent: int) -> int:
    result = 1
    while exponent:
        bit = exponent & -exponent
        result ^= result << bit
        exponent ^= bit
    return result


def multiply_bits(left: int, right: int, degree: int) -> int:
    mask = (1 << (degree + 1)) - 1
    result = 0
    while left:
        bit = left & -left
        result ^= right << (bit.bit_length() - 1)
        left ^= bit
    return result & mask


def streamed_permanent(matrix: list[list[int]], degree: int) -> int:
    """Ryser's formula over F_2[S], retaining at most degree+1 bits."""
    size = len(matrix)
    if not size:
        return 1
    assert all(len(row) == size for row in matrix)
    value = 0
    for subset in range(1 << size):
        term = 1
        for row in matrix:
            row_sum = 0
            for column in range(size):
                if (subset >> column) & 1:
                    row_sum ^= row[column]
            term = multiply_bits(term, row_sum, degree)
            if not term:
                break
        value ^= term
    return value


def unaffected_subsets(operator, output, input_):
    count = len(output) + len(input_) - len(operator)
    if count < 0 or count % 2:
        return
    count //= 2
    common = tuple(sorted(set(output) & set(input_)))
    if count > len(common):
        return
    yield from combinations(common, count)


def streamed_top_profile(label, output, input_) -> int:
    """Bit i is the original aFam angle-i factor, over F_2."""
    h = sum(label)
    length = len(label)
    polynomial = 0
    for unchanged in unaffected_subsets(label, output, input_):
        unchanged = set(unchanged)
        columns = [(True, k) for k in output if k not in unchanged]
        columns += [(False, k) for k in input_ if k not in unchanged]
        matrix = [[1 << (((n & k) - 1) if creation else (n & ~k))
                   for creation, k in columns] for n in label]
        polynomial ^= streamed_permanent(matrix, h - length)
    return multiply_bits(polynomial, binomial_bits(h + length), 2 * h)


def streamed_bottom_factor(label, output, input_) -> int:
    value = 0
    for unchanged in unaffected_subsets(label, output, input_):
        unchanged = set(unchanged)
        columns = [(True, k) for k in output if k not in unchanged]
        columns += [(False, k) for k in input_ if k not in unchanged]
        matrix = [[int(((n - 1) & ~(k - 1) == 0) if creation
                       else ((n - 1) & k == 0))
                   for creation, k in columns] for n in label]
        value ^= streamed_permanent(matrix, 0)
    return value


def verify():
    from check_modular_crossing import modular_angle, modular_operator, regular_partitions
    from math import comb
    states = tuple(mu for h in range(17) for mu in regular_partitions(2, h))
    top_entries = bottom_entries = folded_entries = 0
    max_particles = 0
    for h in range(1, 17):
        for label in regular_partitions(2, h):
            max_particles = max(max_particles, len(label))
            for output in states:
                for input_ in states:
                    profile = streamed_top_profile(label, output, input_)
                    assert profile.bit_length() <= 2 * h + 1
                    for index in range(2 * h + 1):
                        assert ((profile >> index) & 1) == modular_angle(
                            2, label, index, output, input_)
                        top_entries += 1
                    assert streamed_bottom_factor(label, output, input_) == modular_operator(
                        2, label, output, input_)
                    bottom_entries += 1
    for n in range(1, 16, 2):
        for k in range(1, 32, 2):
            creator = annihilator = 0
            for j in range(1, n + 1, 2):
                if comb(n, j) % 2:
                    if k >= j and comb(k - 1, j - 1) % 2:
                        creator ^= binomial_bits(j)
                    if comb(k + j - 1, j - 1) % 2:
                        annihilator ^= binomial_bits(j)
            assert creator == ((1 << ((n & k) - 1)) ^ (1 << (n & k)))
            assert annihilator == ((1 << (n & ~k)) ^ (1 << ((n & ~k) + 1)))
            folded_entries += 2
    # Omitting the common (1+S)^h angle factor changes an actual entry.
    assert streamed_top_profile((1,), (1,), ()) == 0b101
    assert binomial_bits(1) != 0b101
    return {
        "status": "pass",
        "scope": "finite exact checks of streamed local factors; space bound is the accompanying proof",
        "top_angle_entries": top_entries,
        "bottom_factor_entries": bottom_entries,
        "folded_one_particle_entries": folded_entries,
        "largest_tested_operator_particle_count": max_particles,
        "source_revision": "16a5ce41304fcdaea4e0be7c132e1c361dc26da5",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    args = parser.parse_args()
    encoded = json.dumps(verify(), indent=2, sort_keys=True) + "\n"
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(encoded)
    print(encoded, end="")
