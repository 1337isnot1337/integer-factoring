#!/usr/bin/env python3
"""Compare bounded-precision interval jets against exact integer Newton sums."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import random

from interval_unit_jets import (QuadraticJets, interval_elementaries,
                                unit_progression_jet as exact_unit_progression)
from modular_interval_jets import interval_elementaries_mod, unit_progression_jet


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = random.Random(44251)
    for _ in range(1000):
        lower = rng.randrange(0, 100)
        upper = lower + rng.randrange(0, 80)
        degree = rng.randrange(0, 30)
        precision = rng.randrange(1, 40)
        expected = tuple(value % (1 << precision) for value in interval_elementaries(lower, upper, degree))
        assert interval_elementaries_mod(lower, upper, degree, precision) == expected
    for _ in range(1000):
        precision = rng.randrange(1, 18)
        degree = rng.randrange(0, 7)
        ring = QuadraticJets(precision, degree)
        unit = [tuple(rng.randrange(1 << precision) for _ in range(2)) for _ in range(degree + 1)]
        if not (unit[0][0] % 2 or unit[0][1] % 2):
            unit[0] = (unit[0][0] + 1, unit[0][1])
        unit = tuple(unit)
        step = tuple(tuple(2 * rng.randrange(1 << (precision - 1)) for _ in range(2))
                     for _ in range(degree + 1))
        if rng.randrange(2):
            step = ((0, 0),) + step[1:]
        lower = rng.randrange(0, 10)
        upper = lower + rng.randrange(0, 50)
        actual = unit_progression_jet(ring, (lower, upper), unit, step)
        assert actual == exact_unit_progression(ring, (lower, upper), unit, step)
        literal = ring.constant(1)
        for j in range(lower, upper):
            literal = ring.multiply(literal, ring.add(unit, ring.scale(step, j)))
        assert actual == literal
    huge_lower, huge_upper, huge_degree, huge_precision = (1 << 300) + 17, (1 << 301) + 91, 80, 83
    huge = interval_elementaries_mod(huge_lower, huge_upper, huge_degree, huge_precision)
    exact = interval_elementaries(huge_lower, huge_upper, huge_degree)
    assert huge == tuple(value % (1 << huge_precision) for value in exact)
    result = {"claim": "Precision-bounded integer power sums and exact Newton division",
              "elementary_comparisons": 1000, "jet_comparisons": 1000,
              "literal_product_comparisons": 1000,
              "huge_control": {"endpoint_bits": huge_upper.bit_length(), "degree": huge_degree,
                               "precision": huge_precision,
                               "working_precision": huge_precision + huge_degree - huge_degree.bit_count()},
              "scope": "Arithmetic speed and bit-size improvement; no factoring-exponent claim."}
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
