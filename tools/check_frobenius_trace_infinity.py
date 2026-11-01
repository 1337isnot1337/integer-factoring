#!/usr/bin/env python3
"""Exact controls for the degree-free rational-star infinity correction."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import random

from check_frobenius_rational_trace import BinaryField, field_of_degree, rational_value, rational_trace
from frobenius_trace_infinity_correction import degree_free_trace, zero_safe_binomial_trace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = random.Random(0x1AF1A17)
    counts, queries, models = Counter(), {}, []
    for w in range(1, 5):
        field = field_of_degree(2 * w)
        base = 1 << w
        subfield = [value for value in range(field.size) if field.power(value, base) == value]
        outside = [value for value in range(field.size) if field.power(value, base) != value]
        models.append({"w": w, "field_degree": field.degree,
                       "irreducible_polynomial": field.polynomial})
        for trial in range(100):
            denominator_count = trial % 7
            numerator_count = trial % 13
            root_pool = [rng.choice(outside) for _ in range(1 + trial % 3)]
            denominators = []
            for j in range(denominator_count):
                exponent = rng.randrange(w)
                slope = rng.randrange(1, field.size)
                root = root_pool[j % len(root_pool)]
                constant = field.mul(slope, field.power(root, 1 << exponent))
                denominators.append((constant, slope, exponent))
            numerators = [(rng.randrange(field.size), rng.randrange(field.size),
                           rng.randrange(w)) for _ in range(numerator_count)]
            if numerators and trial % 5 == 0:
                a, _, e = numerators[0]
                numerators[0] = (a, 0, e)
            if numerators and denominators and trial % 11 == 0:
                numerators[0] = denominators[0]
                counts["literal_pole_cancellations"] += 1
            counts["cases_with_zero_leading_numerator_factor"] += any(
                not b for _, b, _ in numerators)
            scalar = rng.randrange(field.size)
            actual = degree_free_trace(field, w, scalar, numerators, denominators, queries)
            full = nonzero = 0
            for value in subfield:
                result = rational_value(field, scalar, numerators, denominators, value)
                full ^= result
                if value:
                    nonzero ^= result
            assert actual == (full, nonzero)
            counts["complete_rational_traces"] += 1
            counts["enumerated_field_points"] += len(subfield)
            formal_degree = (sum(1 << e for _, _, e in numerators)
                             - sum(1 << e for _, _, e in denominators))
            counts["positive_formal_degree_cases"] += formal_degree > 0
            counts["denominator_free_cases"] += not denominators
            if formal_degree <= 0:
                assert actual == rational_trace(field, w, scalar, numerators, denominators)
                counts["original_proper_trace_comparisons"] += 1

        for repetitions in (1, 2, 5):
            numerators = [(0, 1, exponent) for _ in range(repetitions)
                          for exponent in range(w)]
            assert degree_free_trace(field, w, 1, numerators, [], queries) == (1, 1)
            counts["positive_multiple_monomial_controls"] += 1
        for scalar in range(field.size):
            assert degree_free_trace(field, w, scalar, [], [], queries) == (0, scalar)
            counts["constant_polynomial_controls"] += 1
        try:
            degree_free_trace(field, w, 1, [(0, 1, 0)], [(1, 1, 0)])
        except ValueError:
            counts["pole_on_subfield_rejections"] += 1
        else:
            raise AssertionError("A genuine pole on the summation field was accepted.")

    # Formal degree beyond several whole field periods, with both a
    # nontrivial denominator and a nonzero polynomial-part contribution.
    field = field_of_degree(12)
    w, base = 6, 64
    subfield = [value for value in range(field.size) if field.power(value, base) == value]
    outside = [value for value in range(field.size) if field.power(value, base) != value]
    for trial in range(12):
        root = rng.choice(outside)
        denominators = [(root, 1, 0)]
        numerators = [(rng.randrange(field.size), rng.randrange(1, field.size),
                       rng.randrange(w)) for _ in range(24)]
        actual = degree_free_trace(field, w, 1, numerators, denominators, queries)
        full = nonzero = 0
        for value in subfield:
            result = rational_value(field, 1, numerators, denominators, value)
            full ^= result
            if value:
                nonzero ^= result
        assert actual == (full, nonzero)
        counts["larger_polynomial_part_traces"] += 1
        counts["enumerated_field_points"] += len(subfield)

    zero_safe_models = []
    for degree, w in ((1, 1), (2, 1), (2, 2), (3, 1), (3, 3),
                      (4, 2), (4, 4), (6, 3), (8, 4)):
        # The earlier irreducibility helper is used only at degree>=2;
        # the linear polynomial X+1 gives the explicit prime field here.
        field = BinaryField(3) if degree == 1 else field_of_degree(degree)
        base = 1 << w
        subfield = [value for value in range(field.size) if field.power(value, base) == value]
        zero_safe_models.append({"field_degree": degree, "subfield_degree": w})
        for trial in range(50):
            numerators = [(rng.randrange(field.size), rng.randrange(field.size), rng.randrange(w))
                          for _ in range(trial % 7)]
            denominators = [(rng.randrange(field.size), rng.randrange(field.size), rng.randrange(w))
                            for _ in range(trial % 6)]
            if denominators and trial % 5 == 0:
                denominators[0] = (0, 0, 0)
            scalar = rng.randrange(field.size)
            actual = zero_safe_binomial_trace(
                w, degree, scalar, numerators, denominators, multiply=field.mul)
            full = nonzero = 0
            for value in subfield:
                term = scalar
                for a, b, exponent in numerators:
                    term = field.mul(term, a ^ field.mul(b, field.power(value, 1 << exponent)))
                for a, b, exponent in denominators:
                    factor = a ^ field.mul(b, field.power(value, 1 << exponent))
                    if not factor:
                        counts["zero_denominator_evaluations"] += 1
                    term = field.mul(term, field.inverse(factor) if factor else 0)
                full ^= term
                if value:
                    nonzero ^= term
            assert actual == (full, nonzero)
            counts["zero_safe_complete_traces"] += 1
            counts["zero_safe_reference_points"] += len(subfield)

    # The same circuit is a polynomial extension on an algebra with
    # nilpotents. Compare that stated power function, not a fictitious
    # inverse of a zero divisor.
    field = field_of_degree(2)
    def dual_multiply(first, second):
        a, b, c, d = first & 3, first >> 2, second & 3, second >> 2
        return field.mul(a, c) | ((field.mul(a, d) ^ field.mul(b, c)) << 2)
    def dual_power(value, exponent):
        result = 1
        for _ in range(exponent):
            result = dual_multiply(result, value)
        return result
    for trial in range(40):
        numerators = [(rng.randrange(16), rng.randrange(16), rng.randrange(2))
                      for _ in range(trial % 4)]
        denominators = [(rng.randrange(16), rng.randrange(16), rng.randrange(2))
                        for _ in range(trial % 3)]
        scalar = rng.randrange(16)
        actual = zero_safe_binomial_trace(
            2, 2, scalar, numerators, denominators, multiply=dual_multiply)
        full = nonzero = 0
        for value in range(4):
            term = scalar
            for a, b, exponent in numerators:
                term = dual_multiply(term, a ^ dual_multiply(b, field.power(value, 1 << exponent)))
            for a, b, exponent in denominators:
                factor = a ^ dual_multiply(b, field.power(value, 1 << exponent))
                term = dual_multiply(term, dual_power(factor, 5))
            full ^= term
            if value:
                nonzero ^= term
        assert actual == (full, nonzero)
        counts["polynomial_zero_divisor_extension_traces"] += 1

    schedules = []
    for coefficient in range(4):
        calls = Counter()
        def counted_multiply(first, second):
            calls["multiplications"] += 1
            return field.mul(first, second)
        def counted_add(first, second):
            calls["additions"] += 1
            return first ^ second
        zero_safe_binomial_trace(
            2, 2, coefficient, [(coefficient, 1, 0), (3, coefficient, 1)],
            [(coefficient, coefficient, 0), (1, coefficient, 1)],
            multiply=counted_multiply, add=counted_add)
        schedules.append(dict(calls))
    assert all(row == schedules[0] for row in schedules)
    counts["zero_safe_fixed_gate_schedules"] = len(schedules)

    names = (Path(__file__).name, "frobenius_trace_infinity_correction.py",
             "check_frobenius_rational_trace.py", "check_carry_fourier_elimination.py",
             "fermat_torus_binomial_contraction.py")
    report = {
        "status": "passed",
        "scope": "ordinary degree-free rational traces and the separate all-coefficient zero-safe binomial circuit",
        "counts": dict(counts), "coefficient_query_counts": queries, "models": models,
        "zero_safe_models": zero_safe_models, "zero_safe_fixed_gate_counts": schedules[0],
        "source_sha256": {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in names},
    }
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "source_sha256"}))


if __name__ == "__main__":
    main()
