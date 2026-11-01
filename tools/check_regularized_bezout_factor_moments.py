#!/usr/bin/env python3
"""Finite exact provider controls for the regularized Bezout moment interface.

The dense provider enumerates its domains. It is a finite control, not the
missing compressed global-jet algorithm.
"""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import prod
from pathlib import Path

from check_balanced_bezout_modular import residue
from check_pole_cauchy_conditioning import mod_pair
from check_pole_smith_stabilization import Quad, W
from regularized_bezout_factor_moments import (
    factor_moments_from_bezout_jet, regularized_bezout_parameters,
)

ROOT = Path(__file__).resolve().parents[1]


def multiply(first, second):
    result = [Quad() for _ in range(len(first) + len(second) - 1)]
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            result[i + j] += a * b
    return result


def product(factors):
    result = [Quad.of(1)]
    for factor in factors:
        result = multiply(result, [Quad.of(x) for x in factor])
    return result


def quotient(numerator, denominator):
    work = list(numerator)
    answer = [Quad() for _ in range(len(numerator) - len(denominator) + 1)]
    for i in range(len(answer) - 1, -1, -1):
        answer[i] = work[i + len(denominator) - 1] / denominator[-1]
        for j, coefficient in enumerate(denominator):
            work[i + j] -= answer[i] * coefficient
    assert not any(work)
    return answer


def dense_regularized_jet(n, first_a, rows, first_b, columns, precision, *, extra_padding=0):
    """Enumerative exact provider, solely for tests with at most32 nodes per axis."""
    if max(rows, columns) > 32:
        raise ValueError("The exact dense provider is restricted to finite controls")
    contract = regularized_bezout_parameters(n, first_a, rows, first_b, columns, precision)
    contract["padding_exponent"] += extra_padding
    contract["clearing_exponent"] += extra_padding
    aa = [first_a + 2 * i for i in range(rows)]
    bb = [first_b + 2 * j for j in range(columns)]
    z = n + (1 << contract["padding_exponent"]) * W
    scale = (1 << contract["clearing_exponent"]) * W
    pp = product([[1 - a, 1] for a in aa])
    qq = product([[z - b, -b] for b in bb])
    vv = [Quad() for _ in aa]
    for a in aa:
        denominator = prod((z - a * b for b in bb), start=Quad.of(1))
        denominator *= prod(a - other for other in aa if other != a)
        basis = product([[1 - other, 1] for other in aa if other != a])
        vv = [old + scale * term / denominator for old, term in zip(vv, basis)]
    numerator = [-value for value in multiply(vv, qq)]
    numerator[0] += scale
    uu = quotient(numerator, pp)
    v_pairs = [mod_pair(value, precision) for value in reversed(vv)]
    u_pairs = [mod_pair(value, precision) for value in reversed(uu)]
    assert all(imaginary == 0 for real, imaginary in v_pairs + u_pairs)
    return {"contract": contract, "rows": aa, "columns": bb, "p": pp, "q": qq,
            "u": uu, "v": vv, "u_reverse_mod": [x[0] for x in u_pairs],
            "v_reverse_mod": [x[0] for x in v_pairs],
            "v_jet": [real for real, imaginary in v_pairs[:min(rows, precision)]],
            "enumerated_domain_nodes": rows + columns}


def controls():
    pairs = coefficients = moments = padding = 0
    examples = []
    for index in range(60):
        rows, columns = 1 << (index % 4), 1 << ((index // 4) % 4)
        n = (9, 15, 21, 25, 35, 45, 49, 63, 77, 105)[index % 10]
        first_a, first_b = 1 + 2 * (index % 3), 1 + 2 * ((index // 3) % 3)
        precision = (1, 3, 6)[index % 3]
        data = dense_regularized_jet(n, first_a, rows, first_b, columns, precision)
        contract, aa, bb = data["contract"], data["rows"], data["columns"]
        matched = [(i, j, a, b) for i, a in enumerate(aa) for j, b in enumerate(bb) if a * b == n]
        weights = []
        for i, j, a, b in matched:
            alpha = prod(a - other for other in aa if other != a) // (1 << contract["row_derivative_order"])
            beta = prod(b - other for other in bb if other != b) // (1 << contract["column_derivative_order"])
            assert alpha % 2 and beta % 2
            weights.append((a, Fraction(1, a ** (columns - 1) * alpha * beta)))
        predicted = [Quad() for _ in aa]
        for a, weight in weights:
            basis = product([[1 - other, 1] for other in aa if other != a])
            predicted = [old + weight * value for old, value in zip(predicted, basis)]
        assert [mod_pair(x, precision) for x in predicted] == [mod_pair(x, precision) for x in data["v"]]
        for key in ("u", "v"):
            for depth, value in enumerate(reversed(data[key])):
                actual = value.valuation()
                assert actual is None or actual >= depth
                coefficients += 1
        count = max(6, 2 * len(matched) + 2)
        result = factor_moments_from_bezout_jet(first_a, rows, precision, data["v_jet"], count)
        expected = [residue(sum((weight * a ** j for a, weight in weights), Fraction()), precision)
                    for j in range(count)]
        assert result["ordinary_moments"] == expected
        assert result["enumerated_domain_nodes"] == 0
        moments += count
        identity, other = multiply(data["u"], data["p"]), multiply(data["v"], data["q"])
        assert all(not (x + y) for x, y in zip(identity[1:], other[1:]))
        assert identity[0] + other[0] == (1 << contract["clearing_exponent"]) * W
        if index % 10 == 0:
            higher = dense_regularized_jet(n, first_a, rows, first_b, columns, precision, extra_padding=2)
            assert data["v_jet"] == higher["v_jet"]
            padding += 1
        if len(examples) < 6 or (matched and rows != columns and len(examples) < 10):
            examples.append({"contract": contract, "factor_pairs": [[a, b] for _, _, a, b in matched],
                             "v_jet": data["v_jet"], "ordinary_moments": expected})
        pairs += 1
    cancellation = dense_regularized_jet(15, 1, 4, 1, 4, 1)
    assert cancellation["v_jet"] == [0]
    assert all(value == 0 for value in cancellation["u_reverse_mod"] + cancellation["v_reverse_mod"])
    return {"complete_regularized_pairs": pairs, "complete_coefficient_bounds": coefficients,
            "factor_support_moments": moments, "independent_padding_comparisons": padding,
            "nonempty_zero_at_low_precision_controls": 1,
            "examples": examples}


def reciprocal_carry_controls():
    identities = carries = 0
    for t in range(6):
        count = 1 << t
        factorial = prod(range(1, count + 1))
        factorial_unit = factorial >> (count - 1)
        for a in (1, 3, 5, 11):
            for first_b in (1, 5):
                for n in (3, 15, 35, 77):
                    row_product = prod(n - a * (first_b + 2 * j) for j in range(count))
                    scale = 1 << (2 * count - 1)
                    assert row_product % scale == 0
                    normalized = row_product // scale
                    q = (Fraction(n, a) - first_b) / 2
                    binomial = prod((q - j for j in range(count)), start=Fraction(1)) / factorial
                    assert normalized == a ** count * factorial_unit * binomial
                    q_mod = residue(q, t + 1)
                    assert normalized % 2 == (q_mod >> t) & 1
                    closest = q_mod % count
                    difference = n - a * (first_b + 2 * closest)
                    if difference:
                        order = (abs(difference) & -abs(difference)).bit_length() - 1
                        actual = (abs(normalized) & -abs(normalized)).bit_length() - 1
                        assert actual == order - t - 1
                    else:
                        assert normalized == 0
                    identities += 1
                    carries += 1
    return {"normalized_row_product_identities": identities, "reciprocal_carry_values": carries}


def implicit_postprocessor_control():
    from balanced_bezout_modular import _inverse, _linear, _multiply, _scalar_coefficients, _unit_product
    from interval_unit_jets import QuadraticJets
    count, precision, length = 1 << 300, 17, 12
    ring = QuadraticJets(precision, length - 1)
    p_reverse = _scalar_coefficients(_unit_product(
        ring, (0, count), _linear(ring, 1), _linear(ring, 0, -2)))
    # Prescribed support is used only to test the supplied-jet postprocessor.
    # It is not a provider for an unknown input N.
    support = ((15, 3), (41, 5))
    shifted = [sum(weight * (a - 1) ** j for a, weight in support) % (1 << precision)
               for j in range(length)]
    supplied = _multiply(p_reverse, shifted, 1 << precision, length)
    result = factor_moments_from_bezout_jet(1, count, precision, supplied, length)
    assert result["ordinary_moments"] == [sum(weight * a ** j for a, weight in support) % (1 << precision)
                                          for j in range(length)]
    return {"dimension_bits": count.bit_length(), "precision": precision, "moment_count": length,
            "supplied_jet": supplied, "prescribed_support": support, "result": result,
            "scope": "Known-support synthetic input tests only the postprocessor; no global provider is claimed."}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"claim": "Fixed-scale global Bezout coefficients contain factor-supported moments",
              "controls": controls(), "row_product_controls": reciprocal_carry_controls(),
              "implicit_postprocessor_control": implicit_postprocessor_control(),
              "fast_global_provider_implemented": False, "factoring_exponent_improved": False}
    report["source_sha256"] = {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
        for name in (Path(__file__).name, "regularized_bezout_factor_moments.py",
                     "balanced_bezout_modular.py", "interval_unit_jets.py", "check_pole_smith_stabilization.py")}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report["controls"].items() if isinstance(value, int)}, indent=2))
    print(json.dumps(report["row_product_controls"], indent=2))
    print(json.dumps({key: report["implicit_postprocessor_control"][key]
                      for key in ("dimension_bits", "precision", "moment_count", "scope")}, indent=2))


if __name__ == "__main__":
    main()
