#!/usr/bin/env python3
"""Full inverse contractions against arbitrary explicit separated rational weights."""

from interval_unit_jets import QuadraticJets
from pole_cauchy_domain_precision import inverse_precision_certificate
from pole_cauchy_pole_polynomials import inverse_pole_polynomial_contraction
from rational_weight_valuation_partition import polynomial_valuation_partition, compile_partition_weight


def inverse_rational_weight_contraction(n, left, right, k, row_numerator, row_denominator,
                                        column_numerator, column_denominator, absolute_precision,
                                        *, normalize_inverse=False):
    """sum_ab (M^-1)[b,a] A(a)/P(a) B(b)/Q(b), with all divisions charged."""
    if not isinstance(absolute_precision, int) or absolute_precision < 1:
        raise ValueError("Require positive absolute precision")
    certificate = inverse_precision_certificate(n, left, right, k)
    normalization = certificate["inverse_loss_bound"] if normalize_inverse else 0
    effective_loss = certificate["inverse_loss_bound"] - normalization
    row_partition = polynomial_valuation_partition(left, row_denominator)
    column_partition = polynomial_valuation_partition(right, column_denominator)
    u = row_partition["maximum_denominator_valuation"]
    v = column_partition["maximum_denominator_valuation"]
    precision = absolute_precision + effective_loss + u + v
    row = compile_partition_weight(row_partition, row_numerator, precision)
    column = compile_partition_weight(column_partition, column_numerator, precision)
    contributions, budgets = [], []
    for row_weight in row["leaf_weights"]:
        row_coefficients = [tuple(value << normalization for value in coefficient)
                            for coefficient in row_weight["coefficients"]]
        for column_weight in column["leaf_weights"]:
            division = row_weight["denominator_bits"] + column_weight["denominator_bits"]
            value = inverse_pole_polynomial_contraction(
                n, left, right, k, row_weight["pole"], column_weight["pole"],
                row_coefficients, column_weight["coefficients"], absolute_precision + division)
            scale = value["scale_exponent"] - division
            if value["working_precision"] and scale < absolute_precision:
                contributions.append((scale, value["residue"]))
            budgets.append({"row_residue": row_weight["residue"], "row_bits": row_weight["bits"],
                            "column_residue": column_weight["residue"], "column_bits": column_weight["bits"],
                            "external_denominator_bits": division,
                            "working_precision": value["working_precision"],
                            "scale_after_external_division": scale,
                            "assembled_product_jets": value["assembled_product_jets"]})
    scale_floor = min([0] + [scale for scale, _ in contributions])
    ring = QuadraticJets(absolute_precision - scale_floor, 0)
    result = ring.constant(0)
    for scale, residue in contributions:
        result = ring.add(result, ring.scale(ring.constant(*residue), 1 << (scale - scale_floor)))
    return {"residue": list(result[0]), "scale_exponent": scale_floor,
            "working_precision": ring.precision, "absolute_precision": absolute_precision,
            "local_approximation_precision": precision,
            "inverse_normalization_exponent": normalization,
            "effective_inverse_loss_bound": effective_loss,
            "row_denominator_valuation_bound": u, "column_denominator_valuation_bound": v,
            "row_partition": row_partition, "column_partition": column_partition,
            "row_nonzero_leaf_weights": len(row["leaf_weights"]),
            "column_nonzero_leaf_weights": len(column["leaf_weights"]),
            "row_maximum_pole_order": row["maximum_pole_order"],
            "column_maximum_pole_order": column["maximum_pole_order"],
            "pole_polynomial_contractions": len(budgets), "contraction_budgets": budgets,
            "precision_certificate": certificate, "enumerated_matrix_nodes": 0,
            "scope": "Explicit univariate rational weights with no domain poles; no implicit exponential-degree or joint denominator input."}
