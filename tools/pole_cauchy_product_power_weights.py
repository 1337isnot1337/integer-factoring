#!/usr/bin/env python3
"""Scaled inverse contractions for products of binary-encoded polynomial powers.

Every base polynomial is explicit and nonzero on its odd domain. Signed
exponents are charged by encoding length. The result has a binary scale;
it does not promise absolute precision after an exponentially large division.
"""

from math import comb

from pole_cauchy_domain_precision import inverse_precision_certificate
from pole_cauchy_pole_polynomials import inverse_pole_polynomial_contraction
from pole_cauchy_quadratic_selectors import _pair_product
from pole_cauchy_residue_blocks import selector_rounding_coefficients
from rational_weight_valuation_partition import polynomial_affine, polynomial_valuation_partition


def _contains(outer, inner):
    return (outer["bits"] <= inner["bits"]
            and (inner["residue"] - outer["residue"]) % (1 << outer["bits"]) == 0)


def product_power_partition(domain, factors):
    """Common regular leaves for product P_i(x)^e_i, including signed e_i.

The common refinement has at most the sum of the input leaf counts. Factor
zeros are rejected even for a positive exponent; the contract avoids them.
"""
    records = []
    for polynomial, exponent in factors:
        if not isinstance(exponent, int):
            raise ValueError("Power exponents must be encoded integers")
        if exponent == 0:
            continue
        partition = polynomial_valuation_partition(domain, polynomial)
        records.append({"polynomial": partition["denominator"], "exponent": exponent,
                        "partition": partition})
    # This also validates an empty product's domain.
    baseline = polynomial_valuation_partition(domain, [1])
    sources = [record["partition"] for record in records] or [baseline]
    candidates = {}
    for source in sources:
        for leaf in source["leaves"]:
            candidates[leaf["residue"], leaf["bits"]] = leaf
    # Dyadic residue classes are either nested or disjoint. Every point's
    # finest source leaf is a maximal member, so these leaves cover the domain.
    finest = [leaf for leaf in candidates.values()
              if not any(other["bits"] > leaf["bits"] and _contains(leaf, other)
                         for other in candidates.values())]
    finest.sort(key=lambda leaf: (leaf["bits"], leaf["residue"]))
    leaves = []
    for leaf in finest:
        residue, bits = leaf["residue"], leaf["bits"]
        valuation, units, factor_valuations = 0, [], []
        for record in records:
            parents = [parent for parent in record["partition"]["leaves"]
                       if _contains(parent, leaf)]
            if len(parents) != 1:
                raise ArithmeticError("The common refinement has no unique containing leaf")
            order = parents[0]["valuation"]
            shifted = polynomial_affine(record["polynomial"], residue, 1 << bits)
            if any(value % (1 << order) for value in shifted):
                raise ArithmeticError("The source valuation does not divide a refined polynomial")
            unit = [value >> order for value in shifted]
            if unit[0] % 2 == 0 or any(value % (1 << i) for i, value in enumerate(unit[1:], 1)):
                raise ArithmeticError("Refinement lost the regular unit coefficient bound")
            valuation += record["exponent"] * order
            units.append(unit)
            factor_valuations.append(order)
        leaves.append({"residue": residue, "bits": bits, "valuation": valuation,
                       "factor_valuations": factor_valuations, "unit_polynomials": units,
                       "domain_nodes": leaf["domain_nodes"]})
    if sum(leaf["domain_nodes"] for leaf in leaves) != baseline["domain_size"]:
        raise ArithmeticError("The common refinement does not cover the domain exactly")
    leaf_bound = sum(len(source["leaves"]) for source in sources)
    if len(leaves) > leaf_bound:
        raise ArithmeticError("The common refinement exceeds its leaf budget")
    minimum = min(leaf["valuation"] for leaf in leaves)
    return {"factors": [{"polynomial": record["polynomial"], "exponent": record["exponent"]}
                         for record in records],
            "factor_partitions": [record["partition"] for record in records],
            "leaves": leaves, "minimum_valuation": minimum,
            "maximum_valuation": max(leaf["valuation"] for leaf in leaves),
            "domain_size": baseline["domain_size"], "source_leaf_bound": leaf_bound,
            "power_encoding_bits": sum(abs(record["exponent"]).bit_length() for record in records),
            "enumerated_domain_nodes": 0}


def _multiply(a, b, modulus):
    size = len(a)
    result = [0] * size
    for i, value in enumerate(a):
        if value:
            for j in range(size - i):
                result[i + j] = (result[i + j] + value * b[j]) % modulus
    return result


def _inverse(a, modulus):
    constant = pow(a[0], -1, modulus)
    result = [constant]
    for i in range(1, len(a)):
        result.append(-constant * sum(a[j] * result[i - j] for j in range(1, i + 1)) % modulus)
    return result


def _power(a, exponent, modulus):
    if exponent < 0:
        a, exponent = _inverse(a, modulus), -exponent
    result = [1] + [0] * (len(a) - 1)
    while exponent:
        if exponent & 1:
            result = _multiply(result, a, modulus)
        exponent >>= 1
        if exponent:
            a = _multiply(a, a, modulus)
    return result


def _pole_coefficients(series, residue, bits, precision):
    """Hermite selector times the local polynomial in (a-residue)/2^bits."""
    negative_w_powers = [(1, 0)]
    for _ in range(len(series) - 1):
        negative_w_powers.append(_pair_product(negative_w_powers[-1], (0, -1)))
    coefficients = [(0, 0)] * (6 * precision - 3)
    for j, selector in selector_rounding_coefficients(precision).items():
        for i, value in enumerate(series):
            if not value:
                continue
            for ell in range(i + 1):
                order = 3 * j - ell
                scalar = (selector * value * comb(i, ell)) << (bits * order)
                coefficients[order - 1] = tuple(a + scalar * b for a, b in zip(
                    coefficients[order - 1], negative_w_powers[i - ell]))
    while coefficients and coefficients[-1] == (0, 0):
        coefficients.pop()
    return coefficients


def compile_normalized_product_power(partition, precision):
    """Approximate 2^-mu product P_i(a)^e_i modulo 2^precision everywhere.

mu is the exact minimum valuation on the domain, supplied by the common
partition. Valuation gaps >= precision are skipped before any binary shift.
"""
    if not isinstance(precision, int) or precision < 1:
        raise ValueError("Require positive normalized output precision")
    minimum = partition["minimum_valuation"]
    compiled, skipped = [], 0
    for leaf in partition["leaves"]:
        gap = leaf["valuation"] - minimum
        if gap >= precision:
            skipped += 1
            continue
        local_precision = precision - gap
        modulus = 1 << local_precision
        series = [1] + [0] * (local_precision - 1)
        for record, unit in zip(partition["factors"], leaf["unit_polynomials"]):
            base = [(unit[i] if i < len(unit) else 0) % modulus for i in range(local_precision)]
            series = _multiply(series, _power(base, record["exponent"], modulus), modulus)
        series = [value << gap for value in series]
        if any(value % (1 << i) for i, value in enumerate(series)):
            raise ArithmeticError("A local normalized power lost coefficient divisibility")
        coefficients = _pole_coefficients(series, leaf["residue"], leaf["bits"], precision)
        if coefficients:
            compiled.append({"pole": [leaf["residue"], -(1 << leaf["bits"])],
                             "residue": leaf["residue"], "bits": leaf["bits"],
                             "valuation_gap": gap, "local_unit_precision": local_precision,
                             "local_quotient_coefficients": series,
                             "coefficients": coefficients, "domain_nodes": leaf["domain_nodes"]})
    return {"precision": precision, "normalization_exponent": minimum,
            "leaf_weights": compiled, "skipped_valuation_leaves": skipped,
            "maximum_pole_order": max((len(leaf["coefficients"]) for leaf in compiled), default=0),
            "enumerated_domain_nodes": 0}


def inverse_product_power_contraction(n, left, right, k, row_factors, column_factors, precision):
    """Return R and binary sigma with original contraction = 2^sigma(R+O(2^p)).

R is the complete contraction of the integral normalized weights against
2^L M^-1. sigma=mu_row+mu_column-L is a lower scale, not an exact valuation.
The numerical p is charged; sigma and the power exponents use binary encoding.
"""
    if not isinstance(precision, int) or precision < 1:
        raise ValueError("Require positive normalized output precision")
    certificate = inverse_precision_certificate(n, left, right, k)
    loss = certificate["inverse_loss_bound"]
    row_partition = product_power_partition(left, row_factors)
    column_partition = product_power_partition(right, column_factors)
    row = compile_normalized_product_power(row_partition, precision)
    column = compile_normalized_product_power(column_partition, precision)
    modulus, residue, budgets = 1 << precision, [0, 0], []
    for row_leaf in row["leaf_weights"]:
        row_coefficients = [tuple(value << loss for value in pair)
                            for pair in row_leaf["coefficients"]]
        for column_leaf in column["leaf_weights"]:
            value = inverse_pole_polynomial_contraction(
                n, left, right, k, row_leaf["pole"], column_leaf["pole"],
                row_coefficients, column_leaf["coefficients"], precision)
            scale = value["scale_exponent"]
            if scale >= precision:
                contribution = [0, 0]
            elif scale >= 0:
                contribution = [(x << scale) % modulus for x in value["residue"]]
            else:
                if any(x % (1 << -scale) for x in value["residue"]):
                    raise ArithmeticError("An integral normalized leaf contraction has a nonintegral residue")
                contribution = [(x >> -scale) % modulus for x in value["residue"]]
            residue = [(a + b) % modulus for a, b in zip(residue, contribution)]
            budgets.append({"row_residue": row_leaf["residue"], "row_bits": row_leaf["bits"],
                            "column_residue": column_leaf["residue"], "column_bits": column_leaf["bits"],
                            "working_precision": value["working_precision"],
                            "scale_before_integral_extraction": scale,
                            "assembled_product_jets": value["assembled_product_jets"]})
    return {"residue": residue, "precision": precision,
            "scale_exponent": row_partition["minimum_valuation"] + column_partition["minimum_valuation"] - loss,
            "row_normalization_exponent": row_partition["minimum_valuation"],
            "column_normalization_exponent": column_partition["minimum_valuation"],
            "inverse_normalization_exponent": loss,
            "row_partition": row_partition, "column_partition": column_partition,
            "row_nonzero_leaf_weights": len(row["leaf_weights"]),
            "column_nonzero_leaf_weights": len(column["leaf_weights"]),
            "row_skipped_valuation_leaves": row["skipped_valuation_leaves"],
            "column_skipped_valuation_leaves": column["skipped_valuation_leaves"],
            "maximum_pole_order": max(row["maximum_pole_order"], column["maximum_pole_order"]),
            "pole_polynomial_contractions": len(budgets), "contraction_budgets": budgets,
            "precision_certificate": certificate, "enumerated_matrix_nodes": 0,
            "scope": "Scaled product-power weights with explicit nonvanishing polynomial bases; no singular joint pole or factor-count evaluator."}
