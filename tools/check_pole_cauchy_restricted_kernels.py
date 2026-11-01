#!/usr/bin/env python3
"""Dense and implicit controls for residue-restricted regular joint kernels."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import random
from time import perf_counter

from check_pole_cauchy_analytic_kernels import evaluate, catalogue
from check_pole_cauchy_conditioning import expand, inverse, mod_pair
from check_pole_smith_stabilization import Quad, matrix_at
from pole_cauchy_restricted_kernels import inverse_restricted_analytic_kernel, polynomial_times_selector
from pole_cauchy_residue_blocks import selector_rounding_coefficients

ROOT = Path(__file__).resolve().parents[1]


def controls():
    rng = random.Random(2026110116)
    polynomial_checks = dense_checks = 0
    for _ in range(60):
        precision = rng.randrange(1, 5)
        degree = rng.randrange(precision)
        coefficients = [rng.randrange(-4, 5) for _ in range(degree + 1)]
        coefficients[0] = 1
        bits, residue = rng.randrange(5), rng.randrange(8)
        residue %= 1 << bits
        pole = Quad(Fraction(residue), Fraction(-(1 << bits)))
        actual = polynomial_times_selector(coefficients, residue, bits, precision)
        for a in (1, 3, 9, 13):
            inverse_value = 1 / (a - pole)
            power, value = Quad(Fraction(1)), Quad()
            for coefficient in actual:
                power *= inverse_value
                value += Quad(*(Fraction(x) for x in coefficient)) * power
            chi = ((1 << bits) * inverse_value)
            chi = chi * chi * chi
            selector, power = Quad(), Quad(Fraction(1))
            rounding = selector_rounding_coefficients(precision)
            for exponent in range(2 * precision):
                selector += rounding.get(exponent, 0) * power
                power *= chi
            polynomial = sum(value * a ** d for d, value in enumerate(coefficients))
            assert value == polynomial * selector
            polynomial_checks += 1
    circuits = catalogue(rng)
    examples = []
    for index in range(30):
        size = rng.randrange(1, 4)
        left = ((rng.randrange(1, 8, 2), size),)
        right = tuple((b, 1) for b in sorted(rng.sample(list(range(1, 16, 2)), size)))
        aa, bb = expand(left), expand(right)
        n = rng.randrange(31, 180, 2)
        k = max(abs(n - a * b) for a in aa for b in bb).bit_length() + 1
        precision, denominator_bits = 1 + index % 3, index % 2
        row_bits, column_bits = index % 4, (index // 4) % 4
        row_residue, column_residue = rng.randrange(8), rng.randrange(8)
        circuit = circuits[index % len(circuits)]
        actual = inverse_restricted_analytic_kernel(n, left, right, k, circuit,
                    row_residue, row_bits, column_residue, column_bits,
                    precision, denominator_bits=denominator_bits)
        dense = inverse(matrix_at(n, aa, bb, k))
        modulus = 1 << (precision + denominator_bits)
        expected = sum((dense[j][i] * evaluate(circuit, a, b, modulus)
                        for i, a in enumerate(aa) for j, b in enumerate(bb)
                        if (a - row_residue) % (1 << row_bits) == 0
                        and (b - column_residue) % (1 << column_bits) == 0), Quad())
        expected *= 1 << actual["inverse_normalization_exponent"]
        assert tuple(actual["residue"]) == mod_pair(expected, precision + denominator_bits)
        assert actual["working_precision"] + actual["scale_exponent"] == precision
        dense_checks += 1
        if index < 10:
            examples.append({"n": n, "size": size, "absolute_precision": precision,
                             "external_denominator_bits": denominator_bits,
                             "compiler_total_degree": actual["compiler"]["total_degree"],
                             "pole_polynomial_contractions": actual["pole_polynomial_contractions"],
                             "residue": actual["residue"]})
    # Force nonzero restricted values, including an external dyadic division.
    # Random residue classes above can be empty or vanish at low precision.
    nonzero = []
    circuit = [('var', 0), ('var', 1), ('const', 2), ('mul', 1, 2), ('add', 0, 3)]
    nodes = [1, 3, 5]
    dense = inverse(matrix_at(15, nodes, nodes, 5))
    for a, b in ((1, 1), (3, 5), (5, 3)):
        for division in (0, 1):
            actual = inverse_restricted_analytic_kernel(
                15, ((1, 3),), ((1, 3),), 5, circuit, a, 3, b, 3,
                5, denominator_bits=division)
            expected = dense[nodes.index(b)][nodes.index(a)] * (a + 2 * b)
            expected *= 1 << actual['inverse_normalization_exponent']
            residue = mod_pair(expected, 5 + division)
            assert tuple(actual['residue']) == residue != (0, 0)
            assert actual['scale_exponent'] == -division
            nonzero.append({'row': a, 'column': b, 'denominator_bits': division,
                            'residue': list(residue), 'absolute_precision': 5})
    return {"exact_polynomial_selector_products": polynomial_checks,
            "dense_restricted_joint_kernel_contractions": dense_checks, "examples": examples,
            "nonzero_singleton_kernel_controls": len(nonzero),
            "nonzero_singleton_examples": nonzero}


def large_control():
    n, size = (1 << 401) + 45, (1 << 200) + 13
    left, right = ((1, size),), ((3, size),)
    k = max(abs(n - 3), abs(n - (2 * size - 1) * (2 * size + 1))).bit_length() + 12
    circuit = [('var', 0), ('var', 1), ('mul', 0, 1), ('const', 2),
               ('mul', 2, 3), ('const', 1), ('add', 4, 5),
               ('pow', 6, -((1 << 512) + (1 << 257) + 123))]
    started = perf_counter()
    lower = inverse_restricted_analytic_kernel(n, left, right, k, circuit, 5, 7, 11, 9, 2)
    elapsed = perf_counter() - started
    higher = inverse_restricted_analytic_kernel(n, left, right, k, circuit, 5, 7, 11, 9, 4)
    assert lower["residue"] == [x % 4 for x in higher["residue"]]
    return {"n": n, "input_bits": n.bit_length(), "domain_size": size, "k": k,
            "circuit": circuit, "lower_result": lower, "higher_result": higher,
            "lower_precision_seconds": elapsed, "enumerated_matrix_nodes": 0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"claim": "Regular joint kernels on residue blocks of a normalized implicit Cauchy inverse",
              "controls": controls(), "factoring_exponent_improved": False}
    if not args.small_only:
        report["large_control"] = large_control()
    report["source_sha256"] = {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
        for name in (Path(__file__).name, "pole_cauchy_restricted_kernels.py",
                     "pole_cauchy_pole_polynomials.py", "pole_cauchy_analytic_kernels.py")}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["controls"], indent=2))


if __name__ == "__main__":
    main()
