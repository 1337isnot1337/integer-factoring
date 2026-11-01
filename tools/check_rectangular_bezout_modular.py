#!/usr/bin/env python3
"""Independent exact controls of complete rectangular modular Bezout pairs."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
import random
from time import perf_counter

from check_balanced_bezout_modular import literal_bezout, multiply, product, residue
from interval_unit_jets import QuadraticJets
from rectangular_bezout_modular import (
    rectangular_bezout_modular, low_surplus_binomials, unit_progression_general,
    _linear, _scalar_coefficients,
)

ROOT = Path(__file__).resolve().parents[1]


def literal_surplus(nodes, surplus):
    a = len(nodes) - 1
    values = [product_scalar(node - j for j in range(a + surplus)) / factorial(a + surplus)
              for node in nodes]
    for width in range(1, len(nodes)):
        values = [(values[j + 1] - values[j]) / (nodes[j + width] - nodes[j])
                  for j in range(len(values) - 1)]
    return factorial(a) * values[0]


def product_scalar(values):
    result = Fraction(1)
    for value in values:
        result *= value
    return result


def controls():
    rng = random.Random(2026110153)
    pairs = coefficients = surpluses = unit_products = rejected = 0
    examples = []
    cases = [(4, 2, 1), (5, 2, 3), (1, 1, 4), (1, 7, 5), (7, 1, 5)]
    cases += [(rng.randrange(1, 10), rng.randrange(1, 10), rng.randrange(1, 8)) for _ in range(59)]
    for index, (m, s, p) in enumerate(cases):
        a, b = 1 + 4 * rng.randrange(6), 1 + 4 * rng.randrange(6)
        da, db = 4 * rng.randrange(1, 8, 2), 4 * rng.randrange(1, 8, 2)
        n = 3 + 4 * rng.randrange(30)
        result = rectangular_bezout_modular(n, a, da, m, b, db, s, p)
        aa, bb = [a + da * i for i in range(m)], [b + db * j for j in range(s)]
        u, v = literal_bezout(n, aa, bb, result["bezout_exponent"])
        for key, expected, size in (("u_reversed_coefficients", u, s), ("v_reversed_coefficients", v, m)):
            actual = result[key] + [0] * (size - len(result[key]))
            assert actual == [residue(value, p) for value in reversed(expected)]
            coefficients += size
        pp = product([[1 - x, 1] for x in aa])
        qq = product([[n - y, -y] for y in bb])
        first, second = multiply(u, pp), multiply(v, qq)
        assert [residue(x + y, p) for x, y in zip(first, second)] == [
            (1 << result["bezout_exponent"]) % (1 << p)] + [0] * (m + s - 1)
        if index < 24:
            for reciprocal, first_node, step, count in ((False, a, da, m), (True, b, db, s)):
                actual = low_surplus_binomials(n, first_node, step, count, p, reciprocal=reciprocal)
                nodes = [(Fraction(n, first_node + step * j) - 3) / 4 if reciprocal
                         else Fraction(first_node + step * j - 1, 4) for j in range(count)]
                for surplus, value in enumerate(actual["coefficients"]):
                    assert value == residue(literal_surplus(nodes, surplus), p)
                    surpluses += 1
        if len(examples) < 5:
            examples.append(result)
        pairs += 1
    for index in range(96):
        p, degree, count = 1 + index % 8, index % 5, index % 11
        ring = QuadraticJets(p, degree)
        unit = tuple((rng.randrange(ring.modulus), 0) for _ in range(degree + 1))
        unit = ((unit[0][0] | 1, 0),) + unit[1:]
        step = tuple((rng.randrange(ring.modulus), 0) for _ in range(degree + 1))
        step = (((step[0][0] & ~1) if index % 3 else 0, 0),) + step[1:]
        expected = ring.constant(1)
        for i in range(count):
            expected = ring.multiply(expected, ring.add(unit, ring.scale(step, i)))
        assert unit_progression_general(ring, count, unit, step) == expected
        unit_products += 1
    invalid = [(3, 1, 4, 0, 1, 4, 2, 4), (3, 1, 2, 2, 1, 4, 2, 4),
               (5, 1, 4, 2, 1, 4, 2, 4), (3, 3, 4, 2, 1, 4, 2, 4),
               (3, 1, 4, 2, 1, 4, 2, 0)]
    for values in invalid:
        try:
            rectangular_bezout_modular(*values)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("An invalid rectangle contract was accepted")
    assert examples[0]["v_reversed_coefficients"][0] == 0
    assert [x["index"] for x in examples[1]["newton_records"]] == [4, 2]
    return {"complete_rational_rectangles": pairs, "complete_coefficient_comparisons": coefficients,
            "low_surplus_rational_checks": surpluses, "general_increment_product_checks": unit_products,
            "invalid_contract_rejections": rejected, "examples": examples}


def implicit_controls():
    records = []
    for rows, columns, p in [((1 << 200) + 13, (1 << 173) + 5, 5),
                              (1 << 130, 1 << 200, 4), (1, (1 << 256) + 7, 5)]:
        n, a, da, b, db = (1 << 401) + 3, 1, 12, 9, 20
        started = perf_counter()
        lower = rectangular_bezout_modular(n, a, da, rows, b, db, columns, p)
        seconds = perf_counter() - started
        higher = rectangular_bezout_modular(n, a, da, rows, b, db, columns, p + 2)
        modulus = 1 << p
        for key in ("u_reversed_coefficients", "v_reversed_coefficients"):
            length = max(len(lower[key]), len(higher[key]))
            first = lower[key] + [0] * (length - len(lower[key]))
            second = [x % modulus for x in higher[key]] + [0] * (length - len(higher[key]))
            assert first == second
        assert lower["bezout_exponent"] == higher["bezout_exponent"]
        ring = QuadraticJets(p, p - 1)
        pp = _scalar_coefficients(unit_progression_general(
            ring, rows, _linear(ring, 1, 1 - a), _linear(ring, 0, -da)))
        qq = _scalar_coefficients(unit_progression_general(
            ring, columns, _linear(ring, -b, n - b), _linear(ring, -db, -db)))
        first = multiply(lower["u_reversed_coefficients"], pp)
        second = multiply(lower["v_reversed_coefficients"], qq)
        length = max(len(first), len(second))
        first += [0] * (length - len(first))
        second += [0] * (length - len(second))
        assert all((x + y) % modulus == 0 for x, y in zip(first, second))
        assert lower["bezout_exponent"] >= p
        assert lower["enumerated_progression_nodes"] == higher["enumerated_progression_nodes"] == 0
        records.append({"lower_seconds": seconds, "lower": lower, "higher": higher,
                        "complete_modular_identity_verified": True})
    assert records[1]["lower"]["scale_and_window"]["defect"] >= 100
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {"claim": "Complete primitive modular Bezout pairs on arbitrary separated arithmetic rectangles",
              "controls": controls(), "factoring_exponent_improved": False}
    if not args.small_only:
        result["implicit_controls"] = implicit_controls()
    result["source_sha256"] = {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
        for name in (Path(__file__).name, "rectangular_bezout_modular.py", "rectangular_bezout_scale.py",
                     "interval_unit_jets.py", "pole_cauchy_conditioning.py")}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result["controls"].items() if isinstance(value, int)}, indent=2))
    if not args.small_only:
        print(json.dumps([{"row_bits": x["lower"]["parameters"]["rows"].bit_length(),
                           "column_bits": x["lower"]["parameters"]["columns"].bit_length(),
                           "defect": x["lower"]["scale_and_window"]["defect"],
                           "maximum_working_precision": x["lower"]["maximum_working_precision"],
                           "u_terms": len(x["lower"]["u_terms"]), "v_terms": len(x["lower"]["v_terms"]),
                           "lower_seconds": x["lower_seconds"]} for x in result["implicit_controls"]], indent=2))


if __name__ == "__main__":
    main()
