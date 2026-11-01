#!/usr/bin/env python3
"""Polynomial-size valuation partitions for explicit univariate denominators."""

from math import comb

from pole_cauchy_quadratic_selectors import _pair_product
from pole_cauchy_rational_poles import divisibility_count
from pole_cauchy_residue_blocks import selector_rounding_coefficients


def _trim(coefficients):
    coefficients = list(coefficients)
    if any(not isinstance(value, int) for value in coefficients):
        raise ValueError("Polynomial coefficients must be explicit integers")
    while coefficients and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def polynomial_affine(coefficients, offset, scale):
    """Exact coefficients of P(offset+scale*X)."""
    coefficients = _trim(coefficients)
    result = [0] * len(coefficients)
    for degree, value in enumerate(coefficients):
        for j in range(degree + 1):
            result[j] += value * comb(degree, j) * offset ** (degree - j) * scale ** j
    return _trim(result)


def _content(coefficients):
    nonzero = [value for value in coefficients if value]
    if not nonzero:
        raise ValueError("The denominator polynomial is zero")
    return min((value & -value).bit_length() - 1 for value in nonzero)


def _regular_unit(coefficients):
    return coefficients[0] % 2 and all(value % (1 << i) == 0 for i, value in enumerate(coefficients[1:], 1))


def _mod_two_degree(coefficients):
    return max(i for i, value in enumerate(coefficients) if value % 2)


def _count(domain, residue, bits):
    return sum(divisibility_count(first - residue, 2, length, bits) for first, length in domain)


def _first_match(domain, residue, bits):
    modulus = 1 << (bits - 1)
    for first, length in domain:
        offset = ((residue - first) // 2) % modulus
        if offset < length:
            return first + 2 * offset
    raise AssertionError("A nonempty residue class had no representative")


def polynomial_value(coefficients, argument):
    result = 0
    for value in reversed(coefficients):
        result = result * argument + value
    return result


def polynomial_valuation_partition(domain, denominator):
    """Disjoint nonempty leaves where P(a)=2^s U((a-alpha)/2^h).

Every leaf U has odd constant and coefficient_i divisible by 2^i.
Actual denominator zeros are rejected with integer witnesses.
"""
    denominator = _trim(denominator)
    if not denominator:
        raise ValueError("The denominator polynomial is zero")
    previous, size = 0, 0
    for first, length in domain:
        if (not isinstance(first, int) or not isinstance(length, int) or first <= previous
                or first % 2 == 0 or length < 1):
            raise ValueError("Use ordered disjoint positive odd progressions")
        previous = first + 2 * (length - 1)
        size += length
    if not size:
        raise ValueError("Use a nonempty odd domain")
    absolute_bound = polynomial_value([abs(value) for value in denominator], previous)
    value_bits = absolute_bound.bit_length()
    initial = polynomial_affine(denominator, 1, 2)
    initial_content = _content(initial)
    initial = [value >> initial_content for value in initial]
    stack = [(1, 1, initial_content, initial)]
    leaves, width_budget = [], {}
    processed = degree_checks = maximum_depth = 0
    while stack:
        residue, bits, valuation, polynomial = stack.pop()
        count = _count(domain, residue, bits)
        if not count:
            continue
        if valuation >= value_bits:
            witness = _first_match(domain, residue, bits)
            assert polynomial_value(denominator, witness) == 0
            raise ValueError(f"The denominator vanishes at domain node {witness}")
        processed += 1
        degree = _mod_two_degree(polynomial)
        width_budget[bits] = width_budget.get(bits, 0) + max(1, degree)
        maximum_depth = max(maximum_depth, bits)
        if _regular_unit(polynomial):
            leaves.append({"residue": residue, "bits": bits, "valuation": valuation,
                           "unit_polynomial": polynomial, "domain_nodes": count})
            continue
        for parity in (0, 1):
            child_residue, child_bits = residue + (parity << bits), bits + 1
            child_count = _count(domain, child_residue, child_bits)
            if not child_count:
                continue
            child = polynomial_affine(polynomial, parity, 2)
            removed = _content(child)
            if removed == 0:
                assert _regular_unit(child)
                leaves.append({"residue": child_residue, "bits": child_bits,
                               "valuation": valuation, "unit_polynomial": child,
                               "domain_nodes": child_count})
            else:
                multiplicity = min(j for j, value in enumerate(child) if (value >> j) % 2)
                normalized = [value >> removed for value in child]
                assert 1 <= removed <= multiplicity
                assert _mod_two_degree(normalized) <= multiplicity
                degree_checks += 1
                stack.append((child_residue, child_bits, valuation + removed, normalized))
    degree_bound = max(1, len(denominator) - 1)
    assert max(width_budget.values()) <= degree_bound
    assert processed <= degree_bound * (value_bits + 1)
    assert sum(leaf["domain_nodes"] for leaf in leaves) == size
    return {"denominator": denominator, "domain_size": size,
            "value_absolute_bound": absolute_bound, "value_bit_bound": value_bits,
            "maximum_denominator_valuation": max(leaf["valuation"] for leaf in leaves),
            "processed_primitive_nodes": processed, "root_multiplicity_degree_checks": degree_checks,
            "maximum_active_degree_budget": max(width_budget.values()),
            "active_degree_bound": degree_bound, "maximum_processed_residue_bits": maximum_depth,
            "leaves": leaves, "enumerated_domain_nodes": 0}


def compile_partition_weight(partition, numerator, precision):
    """Approximate numerator/P by external quadratic pole polynomials.

The pointwise absolute error precision is precision-u, with u the
maximum denominator valuation. Each leaf keeps its explicit2^s divisor.
"""
    if not isinstance(precision, int) or precision < 1:
        raise ValueError("Require positive approximation precision")
    numerator = _trim(numerator)
    modulus = 1 << precision
    rounding = selector_rounding_coefficients(precision)
    negative_w_powers = [(1, 0)]
    for _ in range(precision - 1):
        negative_w_powers.append(_pair_product(negative_w_powers[-1], (0, -1)))
    compiled = []
    for leaf in partition["leaves"]:
        residue, bits, valuation = leaf["residue"], leaf["bits"], leaf["valuation"]
        unit = leaf["unit_polynomial"]
        local_numerator = polynomial_affine(numerator, residue, 1 << bits)
        inverse_constant = pow(unit[0], -1, modulus)
        series = []
        for i in range(precision):
            value = local_numerator[i] if i < len(local_numerator) else 0
            value -= sum(unit[j] * series[i - j] for j in range(1, min(i + 1, len(unit))))
            value = value * inverse_constant % modulus
            assert value % (1 << i) == 0
            series.append(value)
        coefficients = [(0, 0)] * (6 * precision - 3)
        for j, selector_coefficient in rounding.items():
            for i, series_coefficient in enumerate(series):
                if not series_coefficient:
                    continue
                for ell in range(i + 1):
                    exponent = 3 * j - ell
                    assert exponent >= 2 * precision + 1
                    scalar = (selector_coefficient * series_coefficient * comb(i, ell)) << (bits * exponent)
                    value = negative_w_powers[i - ell]
                    coefficients[exponent - 1] = tuple(a + scalar * b
                        for a, b in zip(coefficients[exponent - 1], value))
        while coefficients and coefficients[-1] == (0, 0):
            coefficients.pop()
        if coefficients:
            compiled.append({"pole": [residue, -(1 << bits)], "coefficients": coefficients,
                             "denominator_bits": valuation, "residue": residue, "bits": bits,
                             "local_quotient_coefficients": series,
                             "domain_nodes": leaf["domain_nodes"]})
    return {"numerator": numerator, "precision": precision,
            "pointwise_absolute_precision": precision - partition["maximum_denominator_valuation"],
            "leaf_weights": compiled, "omitted_zero_jets": len(partition["leaves"]) - len(compiled),
            "maximum_pole_order": max((len(leaf["coefficients"]) for leaf in compiled), default=0),
            "enumerated_domain_nodes": 0}
