#!/usr/bin/env python3
"""Convert a supplied global Bezout jet to exact factor-support moments.

This module does not compute the supplied global jet. Its normalization,
truncated interval product, and moment conversion have polynomial bit cost.
"""

from math import comb

from balanced_bezout_modular import _inverse, _linear, _multiply, _scalar_coefficients, _unit_product
from interval_unit_jets import QuadraticJets


def regularized_bezout_parameters(n, first_a, rows, first_b, columns, precision):
    if (any(not isinstance(x, int) for x in (n, first_a, rows, first_b, columns, precision))
            or min(n, first_a, rows, first_b, columns, precision) < 1
            or n % 2 != 1 or first_a % 2 != 1 or first_b % 2 != 1
            or rows & (rows - 1) or columns & (columns - 1)):
        raise ValueError("Use positive odd N and first nodes, two power-of-two counts, and p>=1")
    row_bits, column_bits = rows.bit_length() - 1, columns.bit_length() - 1
    last_a, last_b = first_a + 2 * (rows - 1), first_b + 2 * (columns - 1)
    height = max(1, row_bits + 1, column_bits + 1, abs(n - first_a * first_b).bit_length(),
                 abs(n - last_a * last_b).bit_length())
    padding = height + precision
    row_order, column_order = 2 * rows - 2 - row_bits, 2 * columns - 2 - column_bits
    return {"n": n, "first_a": first_a, "first_b": first_b, "rows": rows, "columns": columns,
            "precision": precision, "height_bound": height, "padding_exponent": padding,
            "row_derivative_order": row_order, "column_derivative_order": column_order,
            "clearing_exponent": padding + row_order + column_order,
            "unit_multiplier": "omega", "reverse_v_degree": min(rows - 1, precision - 1),
            "reverse_u_degree": min(columns - 1, precision - 1),
            "provider_obligation": "Supply the globally correct cleared V jet; this function computes only its contract."}


def factor_moments_from_bezout_jet(first_a, count, precision, reversed_v_coefficients, moment_count):
    """Convert omega*2^(K+E_A+E_B)*V0's supplied reverse jet to a^j moments.

    K>=H+p is part of the provider contract. The theorem then proves that
    the entire supplied jet is real modulo2^p and supported on true factors.
    Missing coefficients within the required prefix are rejected, not guessed.
    """
    if (any(not isinstance(x, int) for x in (first_a, count, precision, moment_count))
            or min(first_a, count, precision) < 1 or moment_count < 0
            or first_a % 2 != 1 or count & (count - 1)):
        raise ValueError("Use an odd positive first node, a power-of-two count, p>=1 and nonnegative output length")
    coefficients = list(reversed_v_coefficients)
    needed = min(count, precision, moment_count)
    if len(coefficients) < needed or any(not isinstance(value, int) for value in coefficients):
        raise ValueError("Supply every required real coefficient of the proven global jet")
    if not moment_count:
        return {"ordinary_moments": [], "shifted_moments": [], "enumerated_domain_nodes": 0}
    modulus = 1 << precision
    size = min(precision, moment_count)
    ring = QuadraticJets(precision, size - 1)
    p_reverse = _scalar_coefficients(_unit_product(
        ring, (0, count), _linear(ring, 1, 1 - first_a), _linear(ring, 0, -2)))
    numerator = [value % modulus for value in coefficients[:needed]]
    numerator += [0] * (size - len(numerator))
    shifted = _multiply(numerator, _inverse(p_reverse, modulus, size), modulus, size)
    shifted += [0] * (moment_count - size)
    for j, value in enumerate(shifted):
        if value % (1 << min(j, precision)):
            raise ArithmeticError("The supplied jet violates the necessary shifted-moment divisibility")
    ordinary = [sum(comb(j, k) * shifted[k] for k in range(min(j + 1, size))) % modulus
                for j in range(moment_count)]
    return {"ordinary_moments": ordinary, "shifted_moments": shifted,
            "precision": precision, "supplied_coefficients_used": needed,
            "interval_product_degree": size - 1, "enumerated_domain_nodes": 0,
            "scope": "Polynomial conversion of a supplied global Bezout jet; no global jet provider is implemented."}
