#!/usr/bin/env python3
"""Arbitrary-precision residue-block sums via polynomially sized selectors."""

from math import comb

from interval_unit_jets import QuadraticJets
from pole_cauchy_quadratic_selectors import inverse_quadratic_pole_contraction
from pole_cauchy_rational_poles import divisibility_count


def selector_rounding_coefficients(precision):
    """Hermite polynomial P: P=0 mod x^p and P=1 mod (1-x)^p."""
    if not isinstance(precision, int) or precision < 1:
        raise ValueError("Require positive selector precision")
    p = precision
    return {p + j: (-1) ** j * comb(2 * p - 1, p + j) * comb(p + j - 1, j)
            for j in range(p)}


def inverse_residue_block(n, left, right, k, row_residue, row_bits,
                          column_residue, column_bits, precision):
    """Sum of 2^L M^-1 over a pair of residue classes, modulo 2^precision.

The implementation uses p^2 repeated-pole contractions with orders at
most 6p-3. Thus p-bit selector lifting never expands a 2^p-th power.
"""
    if any(not isinstance(bits, int) or bits < 0 for bits in (row_bits, column_bits)):
        raise ValueError("Require nonnegative residue precisions")
    if any(not isinstance(residue, int) for residue in (row_residue, column_residue)):
        raise ValueError("Require integer residues")
    coefficients = selector_rounding_coefficients(precision)
    row_residue %= 1 << row_bits
    column_residue %= 1 << column_bits
    c, d = (row_residue, -(1 << row_bits)), (column_residue, -(1 << column_bits))
    ring = QuadraticJets(precision, 0)
    result = ring.constant(0)
    budgets = []
    certificate = None
    for j, row_coefficient in coefficients.items():
        for ell, column_coefficient in coefficients.items():
            contraction = inverse_quadratic_pole_contraction(
                n, left, right, k, c, d, 3 * j, 3 * ell, precision)
            certificate = contraction["precision_certificate"]
            normalization = certificate["inverse_loss_bound"] + 3 * j * row_bits + 3 * ell * column_bits
            scale = contraction["scale_exponent"] + normalization
            if scale >= precision:
                contribution = ring.constant(0)
            elif scale >= 0:
                contribution = ring.scale(ring.constant(*contraction["residue"]), 1 << scale)
            else:
                assert all(value % (1 << (-scale)) == 0 for value in contraction["residue"])
                assert contraction["working_precision"] + scale >= precision
                contribution = ring.constant(*(value >> (-scale) for value in contraction["residue"]))
            result = ring.add(result, ring.scale(contribution, row_coefficient * column_coefficient))
            budgets.append({"row_order": 3 * j, "column_order": 3 * ell,
                            "working_precision": contraction["working_precision"],
                            "jet_shift_exponent": contraction["jet_shift_exponent"],
                            "scale_exponent": contraction["scale_exponent"],
                            "integral_normalization_exponent": normalization})
    return {"residue": list(result[0]), "precision": precision,
            "inverse_normalization_exponent": certificate["inverse_loss_bound"],
            "row_residue": row_residue, "row_bits": row_bits,
            "column_residue": column_residue, "column_bits": column_bits,
            "selected_rows": sum(divisibility_count(first - row_residue, 2, length, row_bits)
                                  for first, length in left),
            "selected_columns": sum(divisibility_count(first - column_residue, 2, length, column_bits)
                                     for first, length in right),
            "selector_degree": 2 * precision - 1,
            "maximum_pole_order": 6 * precision - 3,
            "repeated_pole_contractions": len(budgets), "contraction_budgets": budgets,
            "precision_certificate": certificate, "enumerated_matrix_nodes": 0,
            "scope": "Full-precision residue-block sums of the normalized inverse; no implicit rank evaluator."}
