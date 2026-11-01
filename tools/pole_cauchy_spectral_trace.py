#!/usr/bin/env python3
"""Degree-linear spectral projection and an explicitly charged trace interface."""

from pole_cauchy_residue_blocks import selector_rounding_coefficients


def spectral_trace_contract(size, precision=None):
    """Only the listed traces are inputs; this does not evaluate them."""
    if not isinstance(size, int) or size < 1:
        raise ValueError("Require a positive matrix size")
    precision = size.bit_length() if precision is None else precision
    if not isinstance(precision, int) or precision < size.bit_length():
        raise ValueError("The trace precision must separate every count through the matrix size")
    coefficients = selector_rounding_coefficients(precision)
    return {"matrix_size": size, "precision": precision,
            "trace_coefficients": coefficients, "maximum_matrix_power": 2 * precision - 1,
            "number_of_trace_queries": precision,
            "scope": "Supplied traces of S=M M^T; no implicit trace evaluator is included."}


def padded_trace_contract(size, difference_valuation_bound, precision=None):
    """Pad the displacement so a supplied first trace is the integer count."""
    contract = spectral_trace_contract(size, precision)
    if not isinstance(difference_valuation_bound, int) or difference_valuation_bound < 0:
        raise ValueError("Require a nonnegative valuation bound for all nonzero differences")
    precision = contract["precision"]
    return {"matrix_size": size, "precision": precision,
            "nonzero_difference_valuation_bound": difference_valuation_bound,
            "one_trace_displacement": difference_valuation_bound + (precision + 1) // 2,
            "entry_sum_displacement": difference_valuation_bound + precision,
            "count_sign": 1,
            "scope": "A supplied Tr(M M^T), or the more padded entry sum; neither aggregation is evaluated."}


def count_from_supplied_spectral_traces(size, traces, precision=None):
    contract = spectral_trace_contract(size, precision)
    modulus = 1 << contract["precision"]
    answer = [0, 0]
    for degree, coefficient in contract["trace_coefficients"].items():
        value = traces[degree]
        if isinstance(value, int):
            value = value, 0
        if len(value) != 2 or any(not isinstance(x, int) for x in value):
            raise ValueError("A supplied trace must have two integer quadratic coordinates")
        answer = [(a + coefficient * b) % modulus for a, b in zip(answer, value)]
    if answer[1] or answer[0] > size:
        raise ValueError("The supplied traces do not give an admissible idempotent rank")
    return {**contract, "count": answer[0], "projector_trace_residue": answer,
            "scope": "The count conclusion assumes the supplied traces satisfy the stated matrix contract."}


def leading_inverse_partner(n, size, node):
    """One node of the leading inverse partial involution, without its expansion."""
    if not isinstance(size, int) or size < 1 or size & (size - 1):
        raise ValueError("Use a power-of-two dimension")
    if not isinstance(n, int) or n < 3 or n % 2 == 0:
        raise ValueError("Use a positive odd factoring input")
    if not isinstance(node, int) or node % 2 == 0 or not 1 <= node < 2 * size:
        raise ValueError("The node must lie in the complete odd residue domain")
    modulus = 2 * size
    partner = n * pow(node, -1, modulus) % modulus
    difference = n - node * partner
    bad = difference != 0 and (difference // modulus) % 2 != 0
    return {"node": node, "partner": partner, "retained": bad,
            "residue_modulus": modulus, "retained_difference_valuation": modulus.bit_length() - 1,
            "scope": "One leading inverse edge; enumerating all retained edges costs the numerical size."}
