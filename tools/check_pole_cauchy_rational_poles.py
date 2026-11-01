#!/usr/bin/env python3
"""Exact controls for nonunit rational weights on implicit Cauchy inverses."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import random
from time import perf_counter

from check_pole_cauchy_conditioning import expand, inverse, mod_pair
from check_pole_smith_stabilization import Quad, matrix_at
from interval_unit_jets import QuadraticJets
from pole_cauchy_rational_poles import affine_quadratic_product, affine_valuation
from pole_cauchy_rational_poles import inverse_two_pole_contraction

ROOT = Path(__file__).resolve().parents[1]


def controls():
    rng = random.Random(202611019)
    product_checks = clipped_checks = 0
    for _ in range(1000):
        first, step, length = rng.randrange(-500, 501), rng.randrange(-30, 31), rng.randrange(17)
        shift = None if rng.randrange(4) == 0 else rng.randrange(11)
        if shift is None and any(first + step * i == 0 for i in range(length)):
            continue
        precision = rng.randrange(1, 14)
        ring = QuadraticJets(precision, 0)
        expected = Quad(1)
        for i in range(length):
            expected *= Quad(first + step * i, 0 if shift is None else 1 << shift)
            if shift is not None and (first + step * i) % (1 << shift) == 0:
                clipped_checks += 1
        valuation, unit = affine_quadratic_product(ring, first, step, length, shift)
        assert valuation == affine_valuation(first, step, length, shift) == expected.valuation()
        assert unit[0] == mod_pair(expected / (1 << valuation), precision)
        product_checks += 1
    matrix_checks = contractions = external_high_valuation_controls = zero_product_controls = 0
    examples = []
    for index in range(80):
        size = rng.randrange(1, 6)
        if index % 3:
            left = ((rng.randrange(3, 10, 2), size),)
            right = ((rng.randrange(3, 10, 2), size),)
        else:
            left = tuple((a, 1) for a in sorted(rng.sample(list(range(3, 30, 2)), size)))
            right = tuple((b, 1) for b in sorted(rng.sample(list(range(3, 30, 2)), size)))
        aa, bb = expand(left), expand(right)
        n = rng.randrange(31, 300, 2)
        k = max(abs(n - a * b) for a in aa for b in bb).bit_length() + 1
        dense = inverse(matrix_at(n, aa, bb, k))
        t, z = Quad(Fraction(0), Fraction(1 << k)), Quad(Fraction(n), Fraction(1 << k))
        c_high = n * pow(bb[0], -1, 1 << (k + 3)) % (1 << (k + 3)) + (1 << (k + 3))
        d_high = n * pow(aa[0], -1, 1 << (k + 4)) % (1 << (k + 4)) + (1 << (k + 4))
        pairs = ((-1, max(bb) + 2), (0, -3), (1, n), (c_high, d_high))
        for c, d in pairs:
            if c in aa or d in bb:
                continue
            expected = sum((dense[j][i] / ((a - c) * (b - d))
                            for i, a in enumerate(aa) for j, b in enumerate(bb)), Quad())
            quotient = Quad(1)
            for a in aa:
                quotient *= (z - a * d) / (a - c)
            for b in bb:
                quotient *= (z - b * c) / (b - d)
            for _ in aa:
                quotient /= z
            assert expected == z / (t * (c * d - z)) * (1 - quotient)
            actual = inverse_two_pole_contraction(n, left, right, k, c, d, 7)
            scale = actual["scale_exponent"]
            assert scale <= 0
            assert tuple(actual["residue"]) == mod_pair(expected * (1 << (-scale)), actual["working_precision"])
            assert actual["working_precision"] + scale == 7
            assert actual["product_valuation"] == quotient.valuation()
            external_high_valuation_controls += (c, d) == (c_high, d_high)
            zero_product_controls += c == 0 or d == 0
            contractions += 1
            if len(examples) < 8:
                examples.append({"n": n, "left": left, "right": right, "k": k,
                                 "row_pole": c, "column_pole": d,
                                 "product_valuation": actual["product_valuation"],
                                 "working_precision": actual["working_precision"],
                                 "scale_exponent": scale})
        matrix_checks += 1
    rejections = 0
    for c, d, precision in ((3, -1, 5), (-1, 3, 5), (-1, -1, 0)):
        try:
            inverse_two_pole_contraction(101, ((3, 2),), ((3, 2),), 10, c, d, precision)
        except ValueError:
            rejections += 1
        else:
            raise AssertionError("An invalid pole or precision was accepted")
    return {"literal_affine_product_checks": product_checks, "clipped_factor_controls": clipped_checks,
            "dense_matrices": matrix_checks, "exact_rational_identities": contractions,
            "complete_inverse_contractions": contractions,
            "external_differences_beyond_matrix_displacement": external_high_valuation_controls,
            "zero_pole_controls": zero_product_controls, "invalid_contract_rejections": rejections,
            "examples": examples}


def large_control():
    n, size = (1 << 401) + 45, 1 << 200
    left, right = ((1, size),), ((3, size),)
    k = max(abs(n - 3), abs(n - (2 * size - 1) * (2 * size + 1))).bit_length() + 12
    c, d = -1, 2 * size + 3
    start = perf_counter()
    result = inverse_two_pole_contraction(n, left, right, k, c, d, 12)
    elapsed = perf_counter() - start
    higher = inverse_two_pole_contraction(n, left, right, k, c, d, 16)
    assert result["scale_exponent"] == higher["scale_exponent"]
    assert result["residue"] == [value % (1 << result["working_precision"]) for value in higher["residue"]]
    assert result["row_denominator_bound"] == result["column_denominator_bound"] == 201
    return {"n": n, "input_bits": n.bit_length(), "domain_size": size, "left": left,
            "right": right, "k": k, "row_pole": c, "column_pole": d,
            "result": result, "higher_precision_result": higher,
            "lower_precision_seconds": elapsed, "enumerated_matrix_nodes": 0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--small-only", action="store_true")
    args = parser.parse_args()
    checked = controls()
    print(json.dumps({"controls": checked}, indent=2), flush=True)
    large = None if args.small_only else large_control()
    report = {"claim": "Complete implicit inverse contractions with two external integer poles",
              "controls": checked, "large_control": large,
              "source_sha256": {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, "pole_cauchy_rational_poles.py",
                                             "pole_cauchy_domain_precision.py", "modular_interval_jets.py")},
              "factoring_exponent_improved": False,
              "scope": "Separable rational poles may be even on odd nodes; the joint product-pole trace remains open."}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    if large is not None:
        print(json.dumps({"input_bits": large["input_bits"], "domain_size": large["domain_size"],
                          "product_valuation": large["result"]["product_valuation"],
                          "working_precision": large["result"]["working_precision"],
                          "absolute_precision": large["result"]["absolute_precision"],
                          "lower_precision_seconds": large["lower_precision_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
