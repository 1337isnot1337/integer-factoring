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


def _natural(name: str, value: int, minimum: int = 0) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer at least {minimum}")
    return value


@dataclass(frozen=True)
class WidthStep:
    level: Level
    input_width: int
    ordered_total_bound: int
    output_width: int

    @property
    def input_volume(self) -> int:
        return self.level.grid_side ** 4 * self.input_width

    @property
    def output_volume(self) -> int:
        return self.level.k ** 2 * self.output_width


def width_step(level: Level, width: int) -> WidthStep:
    """One branch-independent padding step, using only integer arithmetic.

    The input clique has G^2 vertices; the output has K vertices. This
    computes equation (8), without constructing tables or a grouping.
    """
    _natural("width", width)
    k = _natural("K", level.k, 1)
    d = _natural("D", level.slack_denominator, 1)
    g = _natural("G", level.grid_side, 1)
    total = 281 * g ** 4 * width // 312 + g ** 3 * width + 3 * g ** 4
    return WidthStep(level, width, total, (d * total + g ** 4 * width) // (d * k ** 2))


def width_schedule(levels: list[Level], top_width: int) -> list[WidthStep]:
    """Propagate common widths from the top level down to the base clique.

    ``levels`` is in the bottom-up order returned by ``schedule``; the
    returned steps are in evaluation order, from top to bottom.
    """
    _natural("top_width", top_width)
    steps = []
    width = top_width
    for level in reversed(levels):
        step = width_step(level, width)
        steps.append(step)
        width = step.output_width
    return steps


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
