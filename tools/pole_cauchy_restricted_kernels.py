#!/usr/bin/env python3
"""Regular joint kernels restricted to residue blocks of normalized inverses."""

from math import comb

from interval_unit_jets import QuadraticJets
from pole_cauchy_analytic_kernels import compile_regular_odd_kernel
from pole_cauchy_domain_precision import inverse_precision_certificate
from pole_cauchy_pole_polynomials import inverse_pole_polynomial_contraction
from pole_cauchy_quadratic_selectors import _pair_product, _pole
from pole_cauchy_residue_blocks import selector_rounding_coefficients


def polynomial_times_selector(coefficients, residue, bits, precision):
    """Proper pole coefficients for P(a)*P_q(chi_(residue,bits)(a))."""
    coefficients = [_pole(value) for value in coefficients]
    while coefficients and coefficients[-1] == (0, 0):
        coefficients.pop()
    if not coefficients:
        return []
    if len(coefficients) > 3 * precision:
        raise ValueError("The polynomial degree must be less than three times selector precision")
    pole = residue, -(1 << bits)
    powers = [(1, 0)]
    for _ in range(len(coefficients) - 1):
        powers.append(_pair_product(powers[-1], pole))
    result = [(0, 0)] * (6 * precision - 3)
    for selector_degree, selector_coefficient in selector_rounding_coefficients(precision).items():
        common = selector_coefficient << (3 * bits * selector_degree)
        for degree, coefficient in enumerate(coefficients):
            if coefficient == (0, 0):
                continue
            for ell in range(degree + 1):
                index = 3 * selector_degree - ell - 1
                value = _pair_product(coefficient, powers[degree - ell])
                multiplier = common * comb(degree, ell)
                result[index] = tuple(a + multiplier * b for a, b in zip(result[index], value))
    return result


def inverse_restricted_analytic_kernel(n, left, right, k, circuit,
                                       row_residue, row_bits, column_residue, column_bits,
                                       absolute_precision, *, denominator_bits=0):
    """Sum selected 2^L M^-1 entries weighted by K(a,b)/2^denominator_bits."""
    if not isinstance(absolute_precision, int) or absolute_precision < 1:
        raise ValueError("Require positive absolute precision")
    if not isinstance(denominator_bits, int) or denominator_bits < 0:
        raise ValueError("Require a nonnegative external dyadic denominator")
    if any(not isinstance(bits, int) or bits < 0 for bits in (row_bits, column_bits)):
        raise ValueError("Require nonnegative residue precisions")
    if any(not isinstance(residue, int) for residue in (row_residue, column_residue)):
        raise ValueError("Require integer residues")
    row_residue %= 1 << row_bits
    column_residue %= 1 << column_bits
    precision = absolute_precision + denominator_bits
    certificate = inverse_precision_certificate(n, left, right, k)
    loss = certificate["inverse_loss_bound"]
    compiled = compile_regular_odd_kernel(circuit, precision)
    groups = {}
    for d, e, value in compiled["terms"]:
        groups.setdefault(d, {})[e] = value
    ring = QuadraticJets(precision, 0)
    answer, budgets = ring.constant(0), []
    for d, polynomial in sorted(groups.items()):
        row = polynomial_times_selector([0] * d + [1], row_residue, row_bits, precision)
        row = [tuple(value << loss for value in coefficient) for coefficient in row]
        column_coefficients = [polynomial.get(e, 0) for e in range(1 + max(polynomial))]
        column = polynomial_times_selector(column_coefficients, column_residue, column_bits, precision)
        result = inverse_pole_polynomial_contraction(
            n, left, right, k, (row_residue, -(1 << row_bits)),
            (column_residue, -(1 << column_bits)), row, column, precision)
        scale = result["scale_exponent"]
        if scale >= precision:
            contribution = ring.constant(0)
        elif scale >= 0:
            contribution = ring.scale(ring.constant(*result["residue"]), 1 << scale)
        else:
            assert all(value % (1 << (-scale)) == 0 for value in result["residue"])
            contribution = ring.constant(*(value >> (-scale) for value in result["residue"]))
        answer = ring.add(answer, contribution)
        budgets.append({"row_monomial_degree": d,
                        "column_polynomial_degree": max(polynomial),
                        "row_pole_degree": result["row_degree"],
                        "column_pole_degree": result["column_degree"],
                        "working_precision": result["working_precision"],
                        "scale_exponent": scale,
                        "assembled_product_jets": result["assembled_product_jets"]})
    return {"residue": list(answer[0]), "working_precision": precision,
            "scale_exponent": -denominator_bits, "absolute_precision": absolute_precision,
            "inverse_normalization_exponent": loss,
            "row_residue": row_residue, "row_bits": row_bits,
            "column_residue": column_residue, "column_bits": column_bits,
            "compiler": compiled, "compilation_precision": precision,
            "selector_precision": precision, "pole_polynomial_contractions": len(budgets),
            "contraction_budgets": budgets, "precision_certificate": certificate,
            "enumerated_matrix_nodes": 0,
            "scope": "Regular joint kernels restricted to normalized inverse residue blocks; no singular joint kernel."}
