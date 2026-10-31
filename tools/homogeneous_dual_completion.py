"""Homogeneous dual-completion seeds from the modular neutral Wick compiler.

This is an exact constructor. It does not certify a strict crossing pressure.
Every coefficient alteration remains in the kernel of the vacuum pairing.
"""

from functools import lru_cache
from math import gcd

from check_modular_crossing import determinant_operator_two, regular_partitions
from neutral_reflected_crossing import (
    kernel_entry, laurent_add, laurent_multiply, sector_labels,
)


def homogeneous_entry(period, parameter, y, v, u):
    """Complete Laurent source entry for an arbitrary homogeneous parameter."""
    weight = sum(parameter)
    labels = (parameter, y, v, u)
    if any(sum(label) != weight for label in labels):
        raise ValueError("all four arguments must have equal weight")
    if period < 2 or period & (period - 1):
        raise ValueError("period must be a power of two")
    if any(tuple(sorted(set(label))) != tuple(label)
           or any(k <= 0 or k >= period or not k % 2 for k in label)
           for label in labels):
        raise ValueError("arguments must be distinct odd modes below the period")
    first = tuple(sorted(period - k for k in parameter))
    fields = (tuple((0, k) for k in first) + tuple((1, k) for k in v)
              + tuple((2, k) for k in y) + tuple((3, k) for k in u))
    size = len(fields)
    entries = {(i, j): kernel_entry(period, *fields[i], *fields[j])
               for i in range(size) for j in range(i + 1, size)}

    @lru_cache(None)
    def pfaffian(mask):
        if not mask:
            return (0, 1)
        first = (mask & -mask).bit_length() - 1
        rest = mask ^ (1 << first)
        remaining, result = rest, (0, 0)
        while remaining:
            second = (remaining & -remaining).bit_length() - 1
            remaining ^= 1 << second
            entry = entries[first, second]
            if entry[1]:
                result = laurent_add(result, laurent_multiply(
                    entry, pfaffian(rest ^ (1 << second))))
        return result

    return pfaffian((1 << size) - 1)


def completion_basis(label_period=16, length=4, weight=32):
    """A right inverse of the selected pairing and its full kernel basis."""
    labels = sector_labels(label_period, length, weight)
    if 2 * weight != label_period * length:
        raise ValueError("this baseline uses the balanced reflected dual")
    reflected = tuple(tuple(sorted(label_period - k for k in mu)) for mu in labels)
    basis = regular_partitions(2, weight)
    complements = []
    for parameter in basis:
        if parameter in reflected:
            continue
        terms = (parameter,) + tuple(reflected[j] for j, nu in enumerate(labels)
                                    if determinant_operator_two(parameter, (), nu))
        for nu in labels:
            assert sum(determinant_operator_two(term, (), nu) for term in terms) % 2 == 0
        complements.append(terms)
    return labels, reflected, tuple(complements)


def build_completion_seed(label_period=16, length=4, weight=32):
    """Return baseline and kernel slabs with one common exact normalization.

    Each slab has dense order y + d*u + d**2*v. A kernel slab can be
    inserted at any logical x without changing the crossing endpoint.
    """
    labels, reflected, complements = completion_basis(label_period, length, weight)
    dimension = len(labels)
    period = 1 << max(weight, label_period - 1).bit_length()
    parameters = tuple(sorted({term for terms in complements for term in terms}
                              | set(reflected)))
    slabs = {}
    for parameter in parameters:
        values = [(0, 0)] * dimension ** 3
        for yi, y in enumerate(labels):
            for vi, v in enumerate(labels):
                for ui, u in enumerate(labels):
                    values[yi + dimension * ui + dimension ** 2 * vi] = (
                        homogeneous_entry(period, parameter, y, v, u))
        slabs[parameter] = values
    complement_slabs = []
    for terms in complements:
        values = [(0, 0)] * dimension ** 3
        for term in terms:
            values = [laurent_add(a, b) for a, b in zip(values, slabs[term])]
        complement_slabs.append(values)
    all_slabs = [slabs[term] for term in reflected] + complement_slabs
    least = min(shift for slab in all_slabs for shift, bits in slab if bits)
    step = 0
    for slab in all_slabs:
        for shift, bits in slab:
            while bits:
                bit = bits & -bits
                step = gcd(step, shift + bit.bit_length() - 1 - least)
                bits ^= bit
    step = max(step, 1)

    def pack(slab):
        values = []
        for shift, bits in slab:
            packed = 0
            while bits:
                bit = bits & -bits
                exponent = shift + bit.bit_length() - 1 - least
                assert exponent % step == 0
                packed ^= 1 << (exponent // step)
                bits ^= bit
            values.append(packed)
        return values

    baseline_slabs = [pack(slabs[term]) for term in reflected]
    kernel_slabs = [pack(slab) for slab in complement_slabs]
    baseline = [0] * dimension ** 4
    for xi, slab in enumerate(baseline_slabs):
        for j, value in enumerate(slab):
            baseline[xi + dimension * j] = value
    for slab in kernel_slabs:
        assert all(value.bit_count() % 2 == 0 for value in slab)
    return {
        "dimensions": [dimension] * 4,
        "label_period": label_period, "operator_period": period,
        "length": length, "weight": weight,
        "labels": labels, "reflected_parameters": reflected,
        "kernel_parameters": complements,
        "baseline": baseline, "kernel_slabs": kernel_slabs,
        "common_laurent_shift": least, "verified_exponent_step": step,
        "dual_parameter_dimension": dimension * len(kernel_slabs),
        "strict_endpoint": True,
    }


def specialize(seed, changes):
    """Apply (logical_x, kernel_index, F4 coefficient) updates exactly."""
    dimension = seed["dimensions"][0]
    values = [[value, 0] for value in seed["baseline"]]
    for logical, kernel_index, scalar in changes:
        if not 0 <= logical < dimension or not 0 <= kernel_index < len(seed["kernel_slabs"]):
            raise ValueError("invalid completion coordinate")
        if scalar not in (0, 1, 2, 3):
            raise ValueError("F4 coefficient must be encoded by two bits")
        for j, value in enumerate(seed["kernel_slabs"][kernel_index]):
            index = logical + dimension * j
            if scalar & 1:
                values[index][0] ^= value
            if scalar & 2:
                values[index][1] ^= value
    return values
