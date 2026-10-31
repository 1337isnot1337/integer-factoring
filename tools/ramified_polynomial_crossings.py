#!/usr/bin/env python3
"""Polynomial endpoint-preserving exterior perturbations over F4.

This constructor starts from an already polynomial tensor. It does not
justify applying the same operation to a numerator with an uncanceled
denominator. Every correction coefficient vanishes at the crossing point.
"""
from __future__ import annotations

import hashlib
import json

from check_rectangular_f4_repair import ZERO, add, scale, times_binary, binary_power


def ramify_binary(polynomial, factor):
    if not isinstance(factor, int) or factor < 1:
        raise ValueError("The ramification must be a positive integer.")
    result = 0
    while polynomial:
        bit = polynomial & -polynomial
        polynomial ^= bit
        result |= 1 << ((bit.bit_length()-1)*factor)
    return result


def endpoint_perturbation(seed, monomials, *, ramification=1):
    """Multiply C(w^d) by 1 + sum a*w^j*(1+w)^k*z^mask.

    Each term is (mask, a, j, k), where a is an F4 coefficient, j is
    integral, and k>=1. One common monomial clears negative j. Exterior
    variables are ordered x, y, u, v; the positions
    within each leg use the ordinary bit-mask order. This is an affine
    sum of exterior monomials, not a product of their factors.
    """
    dimensions = list(seed["dimensions"])
    if len(dimensions) != 4 or len(set(dimensions)) != 1:
        raise ValueError("Use four equal exterior alphabets.")
    q = dimensions[0]
    if q < 2 or q & (q-1):
        raise ValueError("Each alphabet must be a positive power of two.")
    if seed.get("denominator_degree", 0):
        raise ValueError("The source must be polynomial, without a denominator.")
    count = q**4
    if len(seed["values"]) != count:
        raise ValueError("The source coefficient count is inconsistent.")
    monomials = list(monomials)
    for mask, coefficient, shift, order in monomials:
        if not isinstance(mask, int) or not 0 < mask < count:
            raise ValueError("Use a nonempty exterior monomial in the four alphabets.")
        if coefficient not in (1, 2, 3) or not isinstance(shift, int) or not isinstance(order, int) or order < 1:
            raise ValueError("Use a nonzero F4 coefficient, integral shift, and positive endpoint order.")
    common_power = max(0, -min((term[2] for term in monomials), default=0))
    values = [tuple(ramify_binary(int(plane), ramification) for plane in entry)
              for entry in seed["values"]]
    result = [tuple(plane << common_power for plane in entry) for entry in values]
    specifications = []
    for mask, coefficient, shift, order in monomials:
        polynomial = binary_power(3, order) << (shift+common_power)
        for target in range(count):
            if target & mask == mask:
                correction = scale(times_binary(values[target ^ mask], polynomial), coefficient)
                result[target] = add(result[target], correction)
        specifications.append([mask, coefficient, shift, order])
    degree = max((max(entry).bit_length()-1 for entry in result), default=-1)
    rendered = json.dumps(result, separators=(",", ":"))
    return {"field": 4, "dimensions": dimensions, "q_x": q, "q_u": q,
            "denominator_degree": 0, "values": result, "degree": degree,
            "global_w_power": common_power,
            "ramification": ramification, "perturbations": specifications,
            "coefficient_sha256": hashlib.sha256(rendered.encode()).hexdigest(),
            "strict_endpoint": bool(seed.get("strict_endpoint", False)),
            "typed_endpoint": seed.get("typed_endpoint", "inherited from the polynomial source")}


def occupation_charge(mask, modes):
    leg_mask = (1 << modes)-1
    counts = [((mask >> (leg*modes)) & leg_mask).bit_count() for leg in range(4)]
    return counts[1]+counts[3]-counts[0]-counts[2]
