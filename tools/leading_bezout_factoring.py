#!/usr/bin/env python3
"""Factoring reduction to one supplied global regularized Bezout coefficient.

No coefficient provider is implemented here. The oracle supplied to the
wrapper must meet the exact contract returned by leading_coefficient_contract.
"""

from __future__ import annotations

from regularized_bezout_factor_moments import regularized_bezout_parameters


def leading_coefficient_contract(n: int, first_a: int) -> dict:
    if (not isinstance(n, int) or n < 3 or n % 2 != 1
            or not isinstance(first_a, int) or first_a < 3
            or first_a > n or first_a % 2 != 1):
        raise ValueError("Use odd N>=3 and an odd row threshold between 3 and N")
    bits = n.bit_length()
    precision = 2 * bits + 1
    size_bits = bits + precision
    count = 1 << size_bits
    contract = regularized_bezout_parameters(n, first_a, count, 3, count, precision)
    return {**contract, "input_bits": bits, "domain_size_bits": size_bits,
            "requested_coefficient_degree": count - 1,
            "requested_coefficients": 1,
            "coefficient_variable": "X (equivalently X-1 at leading degree)",
            "support": "a divides N, a>=first_a, N/a>=3",
            "decision": "The supplied real leading coefficient is nonzero modulo 2^precision",
            "provider_obligation": "Supply [X^(m-1)] omega*2^(K+2E)*V_t, where V_t is the reduced inverse of Q_t modulo P. No provider is implemented."}


def coefficient_has_factor(contract: dict, coefficient: int) -> bool:
    if not isinstance(coefficient, int):
        raise ValueError("Supply a real integer coefficient residue")
    return coefficient % (1 << contract["precision"]) != 0


def smallest_prime_from_leading_coefficient(n: int, coefficient_oracle) -> dict:
    """Conditional wrapper; all global coefficient queries go to the caller."""
    if not isinstance(n, int) or n < 2:
        raise ValueError("Require an integer N>=2")
    if n % 2 == 0:
        return {"smallest_prime_factor": 2, "coefficient_queries": 0,
                "scope": "Even input handled without a coefficient oracle"}
    queries = 0

    def query(threshold):
        nonlocal queries
        contract = leading_coefficient_contract(n, threshold)
        coefficient = coefficient_oracle(contract)
        queries += 1
        return coefficient_has_factor(contract, coefficient)

    if not query(3):
        return {"smallest_prime_factor": n, "coefficient_queries": queries,
                "scope": "Conditional on the supplied global coefficient oracle"}
    last_possible = n // 3
    if last_possible < 3:
        raise ValueError("The supplied coefficient contradicts the factor-support contract")
    low, high = 0, (last_possible - 3) // 2 + 1
    while low + 1 < high:
        middle = (low + high) // 2
        if query(3 + 2 * middle):
            low = middle
        else:
            high = middle
    largest_divisor = 3 + 2 * low
    if n % largest_divisor:
        raise ValueError("The supplied coefficients do not locate an integer divisor")
    return {"smallest_prime_factor": n // largest_divisor,
            "largest_proper_divisor": largest_divisor,
            "coefficient_queries": queries,
            "scope": "Conditional on the supplied global coefficient oracle"}
