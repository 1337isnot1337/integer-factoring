#!/usr/bin/env python3
"""Independent rational and coefficientwise checks of the modular compiler."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, factorial, prod
from pathlib import Path
import random
from time import perf_counter

from balanced_bezout_modular import balanced_bezout_modular, normalized_parity_difference
from interval_unit_jets import QuadraticJets, interval_elementaries, unit_progression_jet

ROOT = Path(__file__).resolve().parents[1]


def multiply(a, b):
    result = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] += x * y
    return result


def product(factors):
    result = [Fraction(1)]
    for factor in factors:
        result = multiply(result, factor)
    return result


def quotient(numerator, denominator):
    work = list(numerator)
    result = [Fraction(0)] * (len(numerator) - len(denominator) + 1)
    for i in range(len(result) - 1, -1, -1):
        value = work[i + len(denominator) - 1] / denominator[-1]
        result[i] = value
        for j, coefficient in enumerate(denominator):
            work[i + j] -= value * coefficient
    assert all(value == 0 for value in work)
    return result


def residue(value, precision):
    value = Fraction(value)
    assert value.denominator % 2
    return value.numerator * pow(value.denominator, -1, 1 << precision) % (1 << precision)


def translated(polynomial):
    return [sum(polynomial[j] * comb(j, i) for j in range(i, len(polynomial)))
            for i in range(len(polynomial))]


def literal_bezout(n, aa, bb, exponent):
    """Lagrange at A, independent of parity and Newton jets."""
    m = len(aa)
    p = product([[-a, 1] for a in aa])
    q = product([[n, -b] for b in bb])
    v = [Fraction(0)] * m
    for a in aa:
        basis = product([[-other, 1] for other in aa if other != a])
        denominator = prod(n - a * b for b in bb) * prod(a - other for other in aa if other != a)
        v = [x + y / denominator for x, y in zip(v, basis)]
    numerator = [-value for value in multiply(v, q)]
    numerator[0] += 1
    u = quotient(numerator, p)
    scale = 1 << exponent
    return translated([scale * value for value in u]), translated([scale * value for value in v])


def literal_difference(parameters, order, precision):
    n, first_a, step_a, first_b, step_b, height = parameters
    nodes = [Fraction(first_a + step_a * (i // 2) - 1, 1 << height) if i % 2 == 0
             else (Fraction(n, first_b + step_b * (i // 2)) - 1) / (1 << height)
             for i in range(order + 1)]
    values = [Fraction(i % 2) for i in range(order + 1)]
    for width in range(1, order + 1):
        values = [(values[i + 1] - values[i]) / (nodes[i + width] - nodes[i])
                  for i in range(len(values) - 1)]
    normalized = (-1) ** (order + 1) * factorial(order) * values[0] / (1 << (order - 1))
    return residue(normalized, precision)


def controls():
    rng = random.Random(2026110141)
    cases = differences = coefficient_checks = identities = rejected = 0
    examples = []
    for index in range(36):
        height = 1 + index % 3
        stride = 1 << (height + 1)
        first_a, first_b = 1 + stride * rng.randrange(6), 1 + stride * rng.randrange(6)
        step_a, step_b = stride * rng.randrange(1, 8, 2), stride * rng.randrange(1, 8, 2)
        n = 1 + (1 << height) + stride * rng.randrange(20)
        count, precision = 1 + index % 8, 1 + index % 7
        result = balanced_bezout_modular(n, first_a, step_a, first_b, step_b,
                                         count, height, precision)
        aa = [first_a + step_a * i for i in range(count)]
        bb = [first_b + step_b * i for i in range(count)]
        u, v = literal_bezout(n, aa, bb, result["bezout_exponent"])
        expected_u = [residue(value, precision) for value in reversed(u)]
        expected_v = [residue(value, precision) for value in reversed(v)]
        actual_u = result["u_reversed_coefficients"] + [0] * (count - len(result["u_reversed_coefficients"]))
        actual_v = result["v_reversed_coefficients"] + [0] * (count - len(result["v_reversed_coefficients"]))
        assert actual_u == expected_u and actual_v == expected_v
        coefficient_checks += 2 * count
        p = product([[1 - a, 1] for a in aa])
        q = product([[n - b, -b] for b in bb])
        first = multiply(list(reversed(actual_u)), p)
        second = multiply(list(reversed(actual_v)), q)
        combined = [residue(a + b, precision) for a, b in zip(first, second)]
        assert combined == [(1 << result["bezout_exponent"]) % (1 << precision)] + [0] * (2 * count - 1)
        identities += 1
        parameters = (n, first_a, step_a, first_b, step_b, height)
        for order in range(1, min(2 * count, 7)):
            actual = normalized_parity_difference(*parameters, order, precision)
            assert actual["unit"] == literal_difference(parameters, order, precision)
            differences += 1
        assert len(result["retained_newton_coefficients"]) <= 2 * precision
        assert len(result["u_reversed_coefficients"]) <= 4 * precision
        assert len(result["v_reversed_coefficients"]) <= 4 * precision
        if len(examples) < 4:
            examples.append(result)
        cases += 1
    invalid = [(3, 1, 4, 1, 4, 2, 1, 0), (3, 1, 2, 1, 4, 2, 1, 4),
               (3, 3, 4, 1, 4, 2, 1, 4), (5, 1, 4, 1, 4, 2, 1, 4),
               (3, 1, 4, 1, 4, 0, 1, 4)]
    for parameters in invalid:
        try:
            balanced_bezout_modular(*parameters)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("An invalid domain contract was accepted")
    return {"dense_rational_cases": cases, "complete_coefficient_comparisons": coefficient_checks,
            "literal_parity_divided_differences": differences,
            "complete_polynomial_identities": identities,
            "invalid_contract_rejections": rejected, "examples": examples}


def _complete_reversals(parameters, precision):
    """Independent short reverse products using scaled elementary sums."""
    n, aa, da, bb, db, m, h = (parameters[key] for key in
                               ("n", "first_a", "step_a", "first_b", "step_b", "count", "height"))
    modulus = 1 << precision
    degree = min(m, (precision - 1) // h)
    ring = QuadraticJets(precision, degree)
    elementary = interval_elementaries(0, m, degree)
    p = []
    for j in range(degree + 1):
        value = sum(elementary[i] * comb(m - i, j - i)
                    * (aa - 1) ** (j - i) * da ** i for i in range(j + 1))
        p.append((-1) ** j * value % modulus)
    unit, step = list(ring.constant(-bb)), list(ring.constant(-db))
    if degree:
        unit[1], step[1] = ((n - bb) % modulus, 0), (-db % modulus, 0)
    q_jet = unit_progression_jet(ring, (0, m), tuple(unit), tuple(step))
    assert all(b == 0 for _, b in q_jet)
    return p, [a for a, _ in q_jet]


def implicit_controls():
    records = []
    for count, height, precision in [((1 << 200) + 13, 1, 5), (1 << 200, 1, 4),
                                      ((1 << 128) + 9, 2, 5)]:
        stride = 1 << (height + 1)
        n = (1 << 401) + 1 + (1 << height)
        parameters = (n, 1, stride, 1 + 2 * stride, 3 * stride, count, height)
        started = perf_counter()
        lower = balanced_bezout_modular(*parameters, precision)
        seconds = perf_counter() - started
        higher = balanced_bezout_modular(*parameters, precision + 2)
        modulus = 1 << precision
        for key in ("u_reversed_coefficients", "v_reversed_coefficients"):
            low, high = lower[key], [x % modulus for x in higher[key]]
            assert low + [0] * (len(high) - len(low)) == high
        assert lower["bezout_exponent"] == higher["bezout_exponent"]
        p, q = _complete_reversals(lower["parameters"], precision)
        first = multiply(lower["u_reversed_coefficients"], p)
        second = multiply(lower["v_reversed_coefficients"], q)
        size = max(len(first), len(second))
        first += [0] * (size - len(first))
        second += [0] * (size - len(second))
        assert all((a + b) % modulus == 0 for a, b in zip(first, second))
        assert lower["bezout_exponent"] >= precision
        assert lower["u_reversed_coefficients"][0] % 2 == 1
        assert lower["v_reversed_coefficients"][0] % 2 == 1
        records.append({"lower_seconds": seconds, "lower": lower, "higher": higher,
                        "complete_modular_identity_verified": True,
                        "input_nodes_enumerated": 0})
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {"claim": "Complete modular optimal Bezout polynomials on separated arithmetic progressions",
              "controls": controls(), "factoring_exponent_improved": False}
    if not args.small_only:
        result["implicit_controls"] = implicit_controls()
    result["source_sha256"] = {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
                               for name in (Path(__file__).name, "balanced_bezout_modular.py",
                                            "interval_unit_jets.py", "pole_cauchy_conditioning.py")}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result["controls"].items() if isinstance(value, int)}, indent=2))
    if not args.small_only:
        print(json.dumps([{"domain_size_bits": item["lower"]["parameters"]["count"].bit_length(),
                           "exponent_bits": item["lower"]["bezout_exponent"].bit_length(),
                           "maximum_working_precision": item["lower"]["maximum_working_precision"],
                           "u_terms": len(item["lower"]["u_terms"]),
                           "v_terms": len(item["lower"]["v_terms"]),
                           "lower_seconds": item["lower_seconds"]}
                          for item in result["implicit_controls"]], indent=2))


if __name__ == "__main__":
    main()
