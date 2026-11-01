#!/usr/bin/env python3
"""Exact finite-difference identities and dense repeated-pole contractions."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, factorial
from pathlib import Path
import random

from check_pole_cauchy_conditioning import expand, inverse, mod_pair
from check_pole_smith_stabilization import Quad, matrix_at
from pole_cauchy_repeated_poles import inverse_repeated_pole_contraction

ROOT = Path(__file__).resolve().parents[1]


def controls():
    rng = random.Random(2026110110)
    identity_count = dense_count = approximation_count = 0
    examples = []
    for order in range(1, 8):
        for _ in range(40):
            x, c = rng.randrange(1, 20, 2), rng.randrange(-20, 0)
            step = 1 << rng.randrange(6, 11)
            value = sum(((-1) ** (order - 1 - j) * comb(order - 1, j)
                         * Fraction(1, x - c - j * step) for j in range(order)), Fraction())
            value /= step ** (order - 1) * factorial(order - 1)
            product = Fraction(1)
            for j in range(order):
                product /= x - c - j * step
            assert value == product
            identity_count += 1
    for index in range(36):
        size = rng.randrange(1, 5)
        left = tuple((a, 1) for a in sorted(rng.sample(list(range(3, 20, 2)), size)))
        right = ((rng.randrange(1, 6, 2), size),)
        aa, bb = expand(left), expand(right)
        n = rng.randrange(31, 200, 2)
        k = max(abs(n - a * b) for a in aa for b in bb).bit_length() + 1
        c = (-1, 0, max(aa) + 2, aa[0] + (1 << 20))[index % 4]
        d = (-3, max(bb) + 2, 0)[index % 3]
        r, s = 1 + index % 3, 1 + (index // 3) % 3
        dense = inverse(matrix_at(n, aa, bb, k))
        expected = sum((dense[j][i] / ((a - c) ** r * (b - d) ** s)
                        for i, a in enumerate(aa) for j, b in enumerate(bb)), Quad())
        result = inverse_repeated_pole_contraction(n, left, right, k, c, d, r, s, 5)
        scale = result["scale_exponent"]
        assert tuple(result["residue"]) == mod_pair(expected * (1 << (-scale)),
                                                    result["working_precision"])
        assert result["working_precision"] + scale == 5
        assert result["simple_pole_calls"] == r * s
        if r != 1 or s != 1:
            step = 1 << result["finite_difference_step_exponent"]
            row_weights, column_weights = [], []
            for nodes, pole, order, weights in ((aa, c, r, row_weights), (bb, d, s, column_weights)):
                for node in nodes:
                    weight = Fraction(1)
                    for j in range(order):
                        weight /= node - pole - j * step
                    weights.append(weight)
            approximant = sum((dense[j][i] * row_weights[i] * column_weights[j]
                               for i in range(size) for j in range(size)), Quad())
            assert (approximant - expected).valuation() >= 5
            approximation_count += 1
        dense_count += 1
        if len(examples) < 12:
            examples.append({"n": n, "left": left, "right": right, "k": k,
                             "row_pole": c, "column_pole": d, "row_order": r, "column_order": s,
                             "step_exponent": result["finite_difference_step_exponent"],
                             "division_valuation": result["finite_difference_denominator_valuation"],
                             "working_precision": result["working_precision"],
                             "scale_exponent": scale})
    rejected = 0
    for r, s, p in ((0, 1, 5), (1, -1, 5), (1.5, 1, 5), (1, 1, 0)):
        try:
            inverse_repeated_pole_contraction(101, ((3, 2),), ((3, 2),), 10, -1, -1, r, s, p)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("An invalid order or precision was accepted")
    return {"exact_univariate_finite_difference_identities": identity_count,
            "dense_repeated_pole_contractions": dense_count,
            "independent_product_approximation_bounds": approximation_count,
            "invalid_contract_rejections": rejected, "examples": examples}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"claim": "Higher external integer poles with fully charged finite-difference precision",
              "controls": controls(), "factoring_exponent_improved": False,
              "complexity_scope": "Polynomial in numerical pole orders, encoded domains and poles, and precision",
              "source_sha256": {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, "pole_cauchy_repeated_poles.py",
                                             "pole_cauchy_rational_poles.py", "pole_cauchy_domain_precision.py")}}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["controls"], indent=2))


if __name__ == "__main__":
    main()
