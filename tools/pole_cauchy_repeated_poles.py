#!/usr/bin/env python3
"""Higher external poles by precision-certified finite differences.

The running time is polynomial in the numerical pole orders, not merely
their binary lengths. Every division by a power of two is charged.
"""

from math import comb, factorial

from interval_unit_jets import QuadraticJets
from pole_cauchy_domain_precision import inverse_precision_certificate
from pole_cauchy_rational_poles import inverse_two_pole_contraction, pole_denominator_bound


def inverse_repeated_pole_contraction(n, left, right, k, row_pole, column_pole,
                                      row_order, column_order, absolute_precision):
    """sum_ab (M^-1)[b,a]/((a-c)^r (b-d)^s) modulo absolute 2^p."""
    if not isinstance(absolute_precision, int) or absolute_precision < 1:
        raise ValueError("Require positive absolute precision")
    if any(not isinstance(order, int) or order < 1 for order in (row_order, column_order)):
        raise ValueError("Require positive integer pole orders")
    r, s, p = row_order, column_order, absolute_precision
    if r == s == 1:
        result = inverse_two_pole_contraction(n, left, right, k, row_pole, column_pole, p)
        return {**result, "row_order": r, "column_order": s,
                "finite_difference_step_exponent": 0,
                "finite_difference_denominator_valuation": 0,
                "simple_pole_calls": 1, "approximation_loss_bound": 0,
                "scope": "External repeated integer poles; numerical pole orders are charged."}
    certificate = inverse_precision_certificate(n, left, right, k)
    loss = certificate["inverse_loss_bound"]
    u = pole_denominator_bound(left, row_pole)
    v = pole_denominator_bound(right, column_pole)
    approximation_loss = loss + r * u + s * v + max(u, v)
    step_exponent = p + approximation_loss
    assert step_exponent > max(u, v)
    step = 1 << step_exponent
    factorial_product = factorial(r - 1) * factorial(s - 1)
    factorial_valuation = ((r - 1) - (r - 1).bit_count()
                           + (s - 1) - (s - 1).bit_count())
    denominator_valuation = (r + s - 2) * step_exponent + factorial_valuation
    input_precision = p + denominator_valuation
    evaluations = []
    for j in range(r):
        for ell in range(s):
            value = inverse_two_pole_contraction(
                n, left, right, k, row_pole + j * step, column_pole + ell * step,
                input_precision)
            assert value["row_denominator_bound"] == u
            assert value["column_denominator_bound"] == v
            coefficient = (-1) ** (r + s - 2 - j - ell) * comb(r - 1, j) * comb(s - 1, ell)
            evaluations.append((coefficient, value))
    scale_floor = min(value["scale_exponent"] for _, value in evaluations)
    working_precision = input_precision - scale_floor
    ring = QuadraticJets(working_precision, 0)
    numerator = ring.constant(0)
    for coefficient, value in evaluations:
        shifted_coefficient = coefficient * (1 << (value["scale_exponent"] - scale_floor))
        numerator = ring.add(numerator,
                             ring.scale(ring.constant(*value["residue"]), shifted_coefficient))
    odd_denominator = factorial_product >> factorial_valuation
    result = ring.multiply(numerator, ring.inverse(ring.constant(odd_denominator)))
    return {"residue": list(result[0]),
            "scale_exponent": scale_floor - denominator_valuation,
            "working_precision": working_precision, "absolute_precision": p,
            "row_pole": row_pole, "column_pole": column_pole,
            "row_order": r, "column_order": s,
            "row_denominator_bound": u, "column_denominator_bound": v,
            "finite_difference_step_exponent": step_exponent,
            "finite_difference_denominator_valuation": denominator_valuation,
            "factorial_denominator_valuation": factorial_valuation,
            "simple_pole_absolute_precision": input_precision,
            "maximum_simple_pole_working_precision": max(value["working_precision"]
                                                          for _, value in evaluations),
            "simple_pole_calls": len(evaluations),
            "approximation_loss_bound": approximation_loss,
            "precision_certificate": certificate, "enumerated_matrix_nodes": 0,
            "scope": "External repeated integer poles; numerical pole orders are charged."}
