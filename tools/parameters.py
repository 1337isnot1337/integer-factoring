"""Exact parameter schedule for the proposed improved recursion.

This is a parameter/cost-model utility, not an implementation of factoring.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log, log2
import json


@dataclass(frozen=True)
class Level:
    d: int
    k: int
    slack_denominator: int
    logarithm_bound: int
    grid_side: int


def level_parameters(d: int, k: int) -> Level:
    if d < 0 or k < 48:
        raise ValueError("require d >= 0 and k >= 48")
    denominator = (d + 48) ** 2
    logarithm = k.bit_length()
    grid = 12 * denominator ** 2 * k ** 2 * logarithm
    return Level(d, k, denominator, logarithm, grid)


def schedule(level_count: int) -> list[Level]:
    if level_count < 0:
        raise ValueError("level_count must be nonnegative")
    levels = []
    k = 48
    for d in range(level_count):
        level = level_parameters(d, k)
        levels.append(level)
        k = level.grid_side ** 2
    return levels


def admissible_depth(input_bits: int) -> int:
    """Find the maximal depth with K_depth^40 <= input_bits.

    Returns zero below the first threshold, matching the finite-input
    fallback convention in the mathematical note.
    """
    if input_bits < 0:
        raise ValueError("input_bits must be nonnegative")
    d, k = 0, 48
    while True:
        next_k = level_parameters(d, k).grid_side ** 2
        if next_k ** 40 > input_bits:
            return d
        d, k = d + 1, next_k


def summary(level_count: int = 5) -> dict:
    rho = (8 + log2(7)) / 12
    return {
        "kind": "recurrence_parameters_not_a_factoring_benchmark",
        "rho": rho,
        "alpha": -log(rho) / log(4),
        "rational_rho": "281/312",
        "rational_alpha": log(312 / 281) / log(4),
        "levels": [
            {
                "d": level.d,
                "K_bit_length": level.k.bit_length(),
                "D": level.slack_denominator,
                "B": level.logarithm_bound,
                "G_bit_length": level.grid_side.bit_length(),
                "log2_input_length_threshold": 40 * log2(level.k),
            }
            for level in schedule(level_count)
        ],
    }


if __name__ == "__main__":
    print(json.dumps(summary(), indent=2))
