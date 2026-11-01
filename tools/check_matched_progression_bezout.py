#!/usr/bin/env python3
"""Exact independent controls for the deeply matched modular Bezout compiler."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb
from pathlib import Path
from time import perf_counter

from balanced_bezout_modular import _linear, _scalar_coefficients, _unit_product
from check_balanced_bezout_modular import literal_bezout, multiply, product, residue
from interval_unit_jets import QuadraticJets
from matched_progression_bezout import (
    matched_newton_coefficient, matched_progression_bezout, reflected_root_incidence,
    _factorial_order, _ratio_jet,
)

ROOT = Path(__file__).resolve().parents[1]


def direct_difference(m, r, h, a, eta):
    c, lam = Fraction(m * eta, a), Fraction(1 << h, a)
    nodes = [Fraction(i) for i in range(m)] + [(j + c) / (1 - lam * j) for j in range(r)]
    values = [Fraction(0)] * m + [Fraction(1)] * r
    for width in range(1, len(nodes)):
        values = [(values[j + 1] - values[j]) / (nodes[j + width] - nodes[j])
                  for j in range(len(values) - 1)]
    return values[0]


def controls():
    newton = series = pairs = coefficients = incidence = 0
    examples = []
    for t in (1, 2, 3, 4):
        m = 1 << t
        for h in (2 * t + 1, 2 * t + 2):
            a = (1 << h) * m + 1
            for p in (1, 4, 8):
                for r in sorted(set([1, m // 2, m - 1, m])):
                    result = matched_newton_coefficient(m, t, h, a, 3, r, p)
                    scale = h * (m - r) + 2 * _factorial_order(m - 1)
                    expected = residue((1 << scale) * direct_difference(m, r, h, a, 3), p)
                    assert result["coefficient"] == expected
                    newton += 1
                    series += int(result["fallback_nodes"] == 0)
    for t in (1, 2, 3):
        m, h = 1 << t, 2 * t + 1
        step, a = 1 << h, (1 << h) * m + 1
        aa, bb = [a + step * i for i in range(m)], [a - step * i for i in range(m)]
        for eta in (-1, 1, 3):
            for p in (1, 4, 9):
                result = matched_progression_bezout(m, h, a, eta, p)
                u, v = literal_bezout(result["n"], aa, bb, result["bezout_exponent"])
                expanded = []
                for key, poly in (("u_reversed_coefficients", u), ("v_reversed_coefficients", v)):
                    moved = [sum(poly[j] * comb(j, i) * (a - 1) ** (j - i)
                                 for j in range(i, m)) for i in range(m)]
                    expected = [residue(value, p) for value in reversed(moved)]
                    actual = result[key] + [0] * (m - len(result[key]))
                    assert expected == actual
                    expanded.append(list(reversed(actual)))
                    coefficients += m
                pp = product([[-step * i, 1] for i in range(m)])
                qq = product([[result["n"] - a * b, -b] for b in bb])
                left, right = multiply(expanded[0], pp), multiply(expanded[1], qq)
                combined = [residue(x + y, p) for x, y in zip(left, right)]
                assert combined == [(1 << result["bezout_exponent"]) % (1 << p)] + [0] * (2 * m - 1)
                if len(examples) < 3:
                    examples.append(result)
                pairs += 1
        for eta in [0, 1, 3] + [-(step // m) * (i * i + extra) for i in range(m) for extra in (0, 1)]:
            n = a * a + step * m * eta
            if n <= 0:
                continue
            actual = reflected_root_incidence(m, h, a, eta)
            matches = [(i, j) for i, x in enumerate(aa) for j, y in enumerate(bb) if x * y == n]
            assert len(matches) == actual["gcd_degree"]
            assert not matches or matches == [(actual["matched_index"], actual["matched_index"])]
            incidence += 1
    return {"exact_newton_coefficients": newton, "perturbation_path_controls": series,
            "complete_rational_bezout_pairs": pairs, "complete_coefficient_comparisons": coefficients,
            "complete_polynomial_identities": pairs, "nonconstant_incidence_controls": incidence,
            "examples": examples}


def implicit_controls():
    results = []
    for t, p in ((32, 70), (64, 74), (200, 12)):
        m, h = 1 << t, 2 * t + 1
        step, a = 1 << h, (1 << h) * m + 1
        started = perf_counter()
        lower = matched_progression_bezout(m, h, a, 3, p)
        seconds = perf_counter() - started
        higher = matched_progression_bezout(m, h, a, 3, p + 3)
        modulus = 1 << p
        for key in ("u_reversed_coefficients", "v_reversed_coefficients"):
            length = max(len(lower[key]), len(higher[key]))
            first = lower[key] + [0] * (length - len(lower[key]))
            second = [x % modulus for x in higher[key]] + [0] * (length - len(higher[key]))
            assert first == second
        assert lower["bezout_exponent"] == higher["bezout_exponent"]
        degree = min(m, (p - 1) // h)
        ring = QuadraticJets(p, degree)
        # Complete reversed input products: every omitted coefficient has valuation >=p.
        pp = _scalar_coefficients(_unit_product(
            ring, (0, m), _linear(ring, 1), _linear(ring, 0, -step)))
        qq = _scalar_coefficients(_unit_product(
            ring, (0, m), _linear(ring, -a, step * m * 3), _linear(ring, step, step * a)))
        first = multiply(lower["u_reversed_coefficients"], pp)
        second = multiply(lower["v_reversed_coefficients"], qq)
        size = max(len(first), len(second))
        first += [0] * (size - len(first))
        second += [0] * (size - len(second))
        assert all((x + y) % modulus == 0 for x, y in zip(first, second))
        baseline = _ratio_jet(m, m, 0, t, a, 3, p, 0)[0]
        changed = baseline != lower["u_reversed_coefficients"][0]
        if t < 200:
            assert changed
        assert lower["enumerated_nodes"] == higher["enumerated_nodes"] == 0
        results.append({"lower_seconds": seconds, "lower": lower, "higher": higher,
                        "translated_grid_leading_coefficient": baseline,
                        "reciprocal_deformation_changes_output": changed,
                        "complete_modular_identity_verified": True})
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {"claim": "Primitive modular Bezout pairs for deeply matched reflected progressions",
              "controls": controls(), "factoring_exponent_improved": False}
    if not args.small_only:
        result["implicit_controls"] = implicit_controls()
    result["source_sha256"] = {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
        for name in (Path(__file__).name, "matched_progression_bezout.py", "balanced_bezout_modular.py",
                     "interval_unit_jets.py", "pole_cauchy_conditioning.py")}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result["controls"].items() if isinstance(value, int)}, indent=2))
    if not args.small_only:
        print(json.dumps([{"dimension_bits": x["lower"]["parameters"]["count"].bit_length(),
                           "N_bits": x["lower"]["n"].bit_length(),
                           "precision": x["lower"]["parameters"]["precision"],
                           "working_precision": x["lower"]["maximum_working_precision"],
                           "u_terms": len(x["lower"]["u_terms"]),
                           "v_terms": len(x["lower"]["v_terms"] ),
                           "deformation_changed": x["reciprocal_deformation_changes_output"],
                           "lower_seconds": x["lower_seconds"]} for x in result["implicit_controls"]], indent=2))


if __name__ == "__main__":
    main()
