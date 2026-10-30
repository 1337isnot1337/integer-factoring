#!/usr/bin/env python3
"""Exact negative controls for finite-mode versions of Appendix B's operators.

This checker uses sparse integer polynomials and literal creation/derivative
operations. It does not test the asymptotic factoring bound or formalize Lean.
The all-h mismatch and grading proofs are in half-recurrence-barrier.md.
"""

from __future__ import annotations

import json
from math import comb, factorial


# A monomial is the sorted tuple of mode numbers of its variables.
Polynomial = dict[tuple[int, ...], int]
LinearField = dict[int, int]


def binomial(n: int, k: int) -> int:
    """The integer polynomial binomial coefficient, including negative n."""
    if k < 0:
        return 0
    if n >= 0:
        return comb(n, k) if k <= n else 0
    return (-1) ** k * comb(k - n - 1, k)


def add_scaled(target: Polynomial, source: Polynomial, scale: int) -> None:
    for monomial, value in source.items():
        target[monomial] = target.get(monomial, 0) + scale * value
        if target[monomial] == 0:
            del target[monomial]


def nabla(mode: int, poly: Polynomial, cutoff: int) -> Polynomial:
    """nabla_-k = t_k times; nabla_k = k d/dt_k; omitted modes are zero."""
    if mode == 0 or abs(mode) > cutoff:
        return {}
    result: Polynomial = {}
    if mode < 0:
        for monomial, coefficient in poly.items():
            result[tuple(sorted(monomial + (-mode,)))] = coefficient
    else:
        for monomial, coefficient in poly.items():
            multiplicity = monomial.count(mode)
            if multiplicity:
                reduced = list(monomial)
                reduced.remove(mode)
                key = tuple(reduced)
                result[key] = result.get(key, 0) + mode * multiplicity * coefficient
    return result


def op_linear(field: LinearField, r: int, poly: Polynomial, cutoff: int) -> Polynomial:
    """Literal Appendix B operator for a linear polynomial sum_j field[j] t_j."""
    coefficient = sum(value * binomial(r - 1, j - 1)
                      for j, value in field.items())
    return {key: coefficient * value for key, value in nabla(-r, poly, cutoff).items()
            if coefficient * value}


def q(h: int, s: int, a: LinearField, b: LinearField,
      poly: Polynomial, cutoff: int, reverse: bool = False) -> Polynomial:
    result: Polynomial = {}
    for i in range(2 * h + 1):
        angle = {j: coefficient * binomial(j + h, i)
                 for j, coefficient in a.items()}
        r = s - h + i
        if reverse:
            term = op_linear(b, -r, op_linear(angle, r, poly, cutoff), cutoff)
        else:
            term = op_linear(angle, r, op_linear(b, -r, poly, cutoff), cutoff)
        add_scaled(result, term, (-1) ** i)
    return result


def difference(left: Polynomial, right: Polynomial) -> Polynomial:
    result = dict(left)
    add_scaled(result, right, -1)
    return result


def verify() -> dict:
    general_cases = 0
    for h in range(1, 33):
        expected = (h + 1) * comb(2 * h, h - 1)
        for vector in ({(): 1}, {(1, 1, 1): 1}, {(h,): 1}):
            defect = difference(q(h, h, {1: 1}, {h: 1}, vector, h),
                                q(h, h, {1: 1}, {h: 1}, vector, h, True))
            assert defect == {key: expected * value for key, value in vector.items()}
            # Exactly the omitted next mode repairs this particular identity.
            repaired = difference(q(h, h, {1: 1}, {h: 1}, vector, h + 1),
                                  q(h, h, {1: 1}, {h: 1}, vector, h + 1, True))
            assert repaired == {}
            general_cases += 1

    actual_family = []
    for h in range(1, 13):
        # The actual a^mu, b^mu for the one-part partition mu=(h).
        a = {j: (-1) ** (h - 1) * comb(h, j) for j in range(1, h + 1)}
        b = {h: factorial(h - 1)}
        vector = {(h,): 1}
        truncated_sum: Polynomial = {}
        full_sum: Polynomial = {}
        for s in range(-2 * h, 2 * h + 1):
            add_scaled(truncated_sum, q(h, s, a, b, vector, h), 1)
            add_scaled(full_sum, q(h, s, a, b, vector, 3 * h), 1)
            assert q(h, s, a, b, vector, 3 * h) == q(h, s, a, b, vector, 3 * h, True)
        expected_sum = {(h,): factorial(h)}
        assert full_sum == expected_sum
        assert truncated_sum != expected_sum
        actual_family.append({
            "h": h,
            "cutoff": h,
            "expected_crossing_coefficient": str(factorial(h)),
            "truncated_crossing_coefficient": str(truncated_sum.get((h,), 0)),
        })

    assert actual_family[0]["truncated_crossing_coefficient"] == "0"
    return {
        "scope": "finite exact negative controls; the all-h obstruction has a written proof",
        "general_order_swap_cases": general_cases,
        "repaired_order_swap_cases": general_cases,
        "actual_one_part_family_cases": actual_family,
        "full_operator_checks": "cutoff 3h contains every mode used in these linear-field tests",
        "asymptotic_or_factoring_runtime_test": False,
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
