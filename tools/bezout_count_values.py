#!/usr/bin/env python3
"""Fixed-field count-value reduction and a restricted exact parity provider.

The general value oracle is supplied by the caller. The parity provider is
implemented for two ordinary odd intervals with bounded symmetric difference.
Its full boundary enumeration cost is checked before work begins.
"""

from __future__ import annotations

from math import isqrt

from regularized_bezout_factor_moments import regularized_bezout_parameters


def count_value_contract(n: int, first_a: int, precision: int) -> dict:
    if n < 3 or n % 2 != 1 or first_a < 3 or first_a > n or first_a % 2 != 1 or precision < 1:
        raise ValueError("Use odd N>=3, an odd threshold between 3 and N, and positive precision")
    size_exponent = n.bit_length() + precision
    count = 1 << size_exponent
    contract = regularized_bezout_parameters(n, first_a, count, 3, count, precision)
    return {**contract, "size_exponent": size_exponent, "evaluation_point": 0,
            "requested_values": 1,
            "provider_obligation": "Supply omega*2^(K+2E)*V_t(0) modulo 2^precision. This general value provider is not implemented.",
            "count_relation": "The supplied value is an explicitly known odd unit times the number of proper factor rows at least first_a."}


def value_from_reversed_jet(count: int, precision: int, coefficients) -> int:
    needed = min(count, precision)
    if count < 1 or precision < 1 or len(coefficients) < needed:
        raise ValueError("Supply every required top coefficient of the globally correct cleared inverse")
    modulus = 1 << precision
    answer = 0
    for depth, coefficient in enumerate(coefficients[:needed]):
        if not isinstance(coefficient, int) or coefficient % (1 << depth):
            raise ValueError("A supplied coefficient violates its necessary divisibility bound")
        answer += (-1 if (count - 1 - depth) % 2 else 1) * coefficient
    return answer % modulus


def smallest_prime_from_count_values(n: int, value_oracle, *, trial_bits: int) -> dict:
    if n < 2 or trial_bits < 1:
        raise ValueError("Require N>=2 and at least one roughness bit")
    if n % 2 == 0:
        return {"smallest_prime_factor": 2, "value_queries": 0, "trial_divisions": 0}
    bound, root = 1 << trial_bits, isqrt(n)
    divisions = 0
    for candidate in range(3, min(bound, root + 1), 2):
        divisions += 1
        if n % candidate == 0:
            return {"smallest_prime_factor": candidate, "value_queries": 0,
                    "trial_divisions": divisions, "scope": "Direct initial trial division"}
    if bound > root:
        return {"smallest_prime_factor": n, "value_queries": 0,
                "trial_divisions": divisions, "scope": "Trial division certified primality"}
    precision = (n.bit_length() + trial_bits - 1) // trial_bits + 1
    queries = 0

    def query(threshold):
        nonlocal queries
        contract = count_value_contract(n, threshold, precision)
        value = value_oracle(contract)
        if not isinstance(value, int):
            raise ValueError("Supply a real integer value residue")
        queries += 1
        return value % (1 << precision) != 0

    if not query(3):
        return {"smallest_prime_factor": n, "value_queries": queries,
                "trial_divisions": divisions, "precision": precision,
                "scope": "Conditional on the supplied normalized value provider"}
    low, high = 0, (n // 3 - 3) // 2 + 1
    while low + 1 < high:
        middle = (low + high) // 2
        if query(3 + 2 * middle):
            low = middle
        else:
            high = middle
    largest = 3 + 2 * low
    if n % largest:
        raise ValueError("The supplied values do not locate an integer divisor")
    return {"smallest_prime_factor": n // largest, "largest_proper_divisor": largest,
            "value_queries": queries, "trial_divisions": divisions, "precision": precision,
            "scope": "Conditional on the supplied normalized value provider"}


def smallest_prime_from_count_value_windows(n: int, value_oracle, *, trial_bits: int = 1) -> dict:
    """Prime-only windows require logarithmic precision over the fixed field."""
    if n < 2 or trial_bits < 1:
        raise ValueError("Require N>=2 and a positive roughness parameter")
    if n % 2 == 0:
        return {"smallest_prime_factor": 2, "value_queries": 0, "trial_divisions": 0}
    bound, root = 1 << trial_bits, isqrt(n)
    divisions = 0
    for candidate in range(3, min(bound, root + 1), 2):
        divisions += 1
        if n % candidate == 0:
            return {"smallest_prime_factor": candidate, "value_queries": 0,
                    "trial_divisions": divisions, "scope": "Direct initial trial division"}
    if bound > root:
        return {"smallest_prime_factor": n, "value_queries": 0,
                "trial_divisions": divisions, "scope": "Trial division certified primality"}
    prime_count_bound = (n.bit_length() + trial_bits - 1) // trial_bits
    precision = prime_count_bound.bit_length()
    modulus, cache, windows = 1 << precision, {}, 0

    def value(threshold):
        if threshold not in cache:
            supplied = value_oracle(count_value_contract(n, threshold, precision))
            if not isinstance(supplied, int):
                raise ValueError("Supply a real integer normalized value residue")
            cache[threshold] = supplied % modulus
        return cache[threshold]

    lower = bound
    while lower <= root:
        upper = min(lower * lower, root + 1)
        left = lower if lower % 2 else lower + 1
        right = upper if upper % 2 else upper + 1
        windows += 1
        if left < right and (value(left) - value(right)) % modulus:
            while right - left > 2:
                middle = left + 2 * ((right - left) // 4)
                if (value(left) - value(middle)) % modulus:
                    right = middle
                else:
                    left = middle
            if n % left:
                raise ValueError("The supplied values do not locate an integer divisor")
            return {"smallest_prime_factor": left, "value_queries": len(cache),
                    "trial_divisions": divisions, "precision": precision,
                    "prime_only_windows": windows,
                    "scope": "Conditional on the supplied normalized value provider"}
        lower = upper
    return {"smallest_prime_factor": n, "value_queries": len(cache),
            "trial_divisions": divisions, "precision": precision,
            "prime_only_windows": windows,
            "scope": "Conditional on the supplied normalized value provider"}


def _difference(first, last, common_first, common_last):
    if common_first > common_last:
        return [(first, last)]
    return [(lo, hi) for lo, hi in ((first, common_first - 2), (common_last + 2, last))
            if lo <= hi]


def _trim(segments, lower, upper):
    answer = []
    for lo, hi in segments:
        start, stop = max(lo, lower), min(hi, upper)
        start += (start - lo) % 2
        stop -= (stop - lo) % 2
        if start <= stop:
            answer.append((start, stop))
    return answer


def boundary_parity_provider(n: int, first_a: int, rows: int, first_b: int,
                             columns: int, *, boundary_limit: int) -> dict:
    """An exact precision-one provider; it visits only the symmetric difference."""
    contract = regularized_bezout_parameters(n, first_a, rows, first_b, columns, 1)
    last_a, last_b = first_a + 2 * (rows - 1), first_b + 2 * (columns - 1)
    common_first, common_last = max(first_a, first_b), min(last_a, last_b)
    a_only = _difference(first_a, last_a, common_first, common_last)
    b_only = _difference(first_b, last_b, common_first, common_last)
    raw_boundary = sum((hi - lo) // 2 + 1 for lo, hi in a_only + b_only)
    a_only = _trim(a_only, (n + last_b - 1) // last_b, n // first_b)
    b_only = (_trim(b_only, (n + common_last - 1) // common_last, n // common_first)
              if common_first <= common_last else [])
    boundary = sum((hi - lo) // 2 + 1 for lo, hi in a_only + b_only)
    if boundary_limit < 0 or boundary > boundary_limit:
        raise ValueError("The explicit symmetric-difference boundary exceeds the supplied work limit")

    def belongs(value, first, last):
        return first <= value <= last and (value - first) % 2 == 0

    root = isqrt(n)
    diagonal = root * root == n and belongs(root, common_first, common_last)
    parity = int(diagonal)
    witnesses = [[root, root]] if diagonal else []
    for lo, hi in a_only:
        for a in range(lo, hi + 1, 2):
            if n % a == 0 and belongs(n // a, first_b, last_b):
                parity ^= 1
                witnesses.append([a, n // a])
    for lo, hi in b_only:
        for b in range(lo, hi + 1, 2):
            if n % b == 0 and belongs(n // b, common_first, common_last):
                parity ^= 1
                witnesses.append([n // b, b])
    return {"contract": contract, "leading_coefficient_mod2": parity,
            "value_at_zero_mod2": parity,
            "u_shifted_monomial": {"degree": columns - 1, "coefficient_mod2": parity},
            "v_shifted_monomial": {"degree": rows - 1, "coefficient_mod2": parity},
            "raw_boundary_nodes": raw_boundary,
            "visited_boundary_nodes": boundary, "interior_nodes_enumerated": 0,
            "uncancelled_witnesses": witnesses,
            "scope": "Exact precision-one provider with cost linear in the explicit symmetric difference"}
