"""Exact source tensors with arbitrary operator parameter of weight at most h.

The physical states have common weight h. The source angle factor for a
parameter monomial of weight w is z^-w (1+z)^(h+w), which is already
included in the neutral kernel. This module makes no rate assertion.
"""

from functools import lru_cache
from math import gcd

from check_modular_crossing import (
    determinant_operator_two, homogeneous_family_two, regular_partitions,
)
from check_modular_parameter_descent import reduced_dual
from neutral_reflected_crossing import kernel_entry, laurent_add, laurent_multiply


def source_monomial_entry(period, parameter, y, v, u):
    h = sum(v)
    if sum(y) != h or sum(u) != h or sum(parameter) > h:
        raise ValueError("physical states must have weight h and parameter weight at most h")
    if period < 2 or period & (period - 1):
        raise ValueError("the auxiliary period must be a power of two")
    for parts in (parameter, y, v, u):
        if tuple(sorted(set(parts))) != tuple(parts) or any(k < 1 or k >= period or k % 2 == 0 for k in parts):
            raise ValueError("use distinct positive odd parts below the auxiliary period")
    fields = (tuple((0, period - m) for m in parameter)
              + tuple((1, m) for m in v) + tuple((2, m) for m in y)
              + tuple((3, m) for m in u))
    if len(fields) % 2:
        return 0, 0
    entries = {(i, j): kernel_entry(period, *fields[i], *fields[j])
               for i in range(len(fields)) for j in range(i + 1, len(fields))}

    @lru_cache(None)
    def pfaffian(mask):
        if not mask:
            return 0, 1
        bit = mask & -mask
        i = bit.bit_length() - 1
        rest = mask ^ bit
        remaining, result = rest, (0, 0)
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            j = bit.bit_length() - 1
            if entries[i, j][1]:
                result = laurent_add(result, laurent_multiply(
                    entries[i, j], pfaffian(rest ^ bit)))
        return result

    return pfaffian((1 << len(fields)) - 1)


def full_weight_completion(weight):
    """All parameter directions below a nonsingular full homogeneous weight."""
    labels, duals = homogeneous_family_two(weight)
    q = len(labels)
    parameters = tuple(m for w in range(weight + 1) for m in regular_partitions(2, w))
    period = 1 << weight.bit_length()
    raw = {}
    for a in parameters:
        raw[a] = [source_monomial_entry(period, a, y, v, u)
                  for v in labels for u in labels for y in labels]

    def combine(terms):
        result = [(0, 0)] * q ** 3
        for a in terms:
            result = [laurent_add(first, second) for first, second in zip(result, raw[a])]
        return result

    baseline = [combine(terms) for terms in duals]
    kernel_parameters = []
    kernel_slabs = []
    for parameter in parameters:
        if sum(parameter) == weight:
            continue
        terms = {parameter}
        for j, b in enumerate(labels):
            if determinant_operator_two(parameter, (), b):
                terms.symmetric_difference_update(duals[j])
        terms = tuple(sorted(terms))
        assert all(sum(determinant_operator_two(a, (), b) for a in terms) % 2 == 0
                   for b in labels)
        kernel_parameters.append(terms)
        kernel_slabs.append(combine(terms))
    integral_baseline = [combine(tuple(a for a, coefficient in reduced_dual(2, b) if coefficient))
                         for b in labels]
    all_slabs = baseline + kernel_slabs + integral_baseline
    least = min(shift for slab in all_slabs for shift, value in slab if value)
    step = 0
    for slab in all_slabs:
        for shift, value in slab:
            while value:
                bit = value & -value
                value ^= bit
                step = gcd(step, shift + bit.bit_length() - 1 - least)
    step = max(step, 1)

    def pack(slab):
        result = []
        for shift, value in slab:
            packed = 0
            while value:
                bit = value & -value
                value ^= bit
                power = shift + bit.bit_length() - 1 - least
                assert power % step == 0
                packed ^= 1 << (power // step)
            result.append(packed)
        return result

    def dense(slabs):
        result = [0] * q ** 4
        for x, slab in enumerate(slabs):
            for j, value in enumerate(pack(slab)):
                result[x + q * j] = value
        return result

    return {'constructor': 'bounded_source_dual_completion', 'weight': weight,
            'dimensions': [q] * 4, 'labels': labels, 'homogeneous_duals': duals,
            'kernel_parameters': kernel_parameters, 'kernel_slabs': [pack(slab) for slab in kernel_slabs],
            'baseline': dense(baseline), 'integral_baseline': dense(integral_baseline),
            'common_laurent_shift': least, 'verified_exponent_step': step,
            'full_parameter_dimension': len(parameters),
            'dual_parameter_dimension': q * len(kernel_parameters), 'strict_endpoint': True}
