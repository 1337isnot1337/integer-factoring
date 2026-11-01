#!/usr/bin/env python3
"""Exact paired full-letter Laurent transvections of reflection crossings.

Every candidate is expanded before poles and its actual degree are checked.
The two paired legs carry inverse-transpose frames. A frame may have a pole
while its action on the entire tensor is polynomial.
"""
from __future__ import annotations

from hashlib import sha256
import json

from reflection_dual_rail_gaussian import tensor


def scale(value, scalar):
    a, b = value
    if scalar == 0:
        return 0, 0
    if scalar == 1:
        return a, b
    if scalar == 2:
        return b, a ^ b
    if scalar == 3:
        return a ^ b, a
    raise ValueError('Use a coefficient in F4.')


def add(first, second):
    return first[0] ^ second[0], first[1] ^ second[1]


def valuation(value):
    union = value[0] | value[1]
    return (union & -union).bit_length() - 1 if union else None


def degree(value):
    return (value[0] | value[1]).bit_length() - 1


def reverse_polynomial(value, bound):
    result = 0
    while value:
        bit = value & -value
        result ^= 1 << (bound - bit.bit_length() + 1)
        value ^= bit
    return result


def reverse_values(values, bound):
    return [(reverse_polynomial(a, bound), reverse_polynomial(b, bound))
            for a, b in values]


def frame_terms(values, q, leg, target, source):
    """Return linear and quadratic terms for E_target,source on x or y.

    Little-endian physical order is x,y,u,v. On the paired leg the
    inverse-transpose is E_source,target, since E_target,source squared
    is zero. Characteristic two makes the inverse sign positive.
    """
    if leg not in (0, 1) or not (0 <= target < q and 0 <= source < q) or target == source:
        raise ValueError('Use x/y and two distinct physical letters.')
    paired = 3 - leg
    stride, paired_stride = q ** leg, q ** paired
    linear, quadratic = {}, {}
    for index, value in enumerate(values):
        symbol = (index // stride) % q
        paired_symbol = (index // paired_stride) % q
        if symbol == target:
            original = index + (source - target) * stride
            if values[original] != (0, 0):
                linear[index] = values[original]
        if paired_symbol == source:
            original = index + (target - source) * paired_stride
            entry = add(linear.get(index, (0, 0)), values[original])
            if entry != (0, 0):
                linear[index] = entry
            elif index in linear:
                del linear[index]
        if symbol == target and paired_symbol == source:
            original = index + (source - target) * stride + (target - source) * paired_stride
            if values[original] != (0, 0):
                quadratic[index] = values[original]
    return linear, quadratic


def admissible_poles(linear, quadratic, bound, scalar=1):
    """All k>=1 for which w^-k L + w^-2k Q has no pole.

    The intersection expression is tested after exact cancellation; separate
    sufficient valuation inequalities would incorrectly discard candidates.
    """
    if not linear and not quadratic:
        return []
    square = (0, 1, 3, 2)[scalar]
    keys = set(linear) | set(quadratic)
    entries = [(scale(linear.get(i, (0, 0)), scalar),
                scale(quadratic.get(i, (0, 0)), square)) for i in keys]
    upper = bound
    for index, value in linear.items():
        if index not in quadratic:
            upper = min(upper, valuation(value))
    result = []
    for k in range(1, upper + 1):
        mask = (1 << (2 * k)) - 1
        if all(not (((left[0] << k) ^ right[0]) & mask)
               and not (((left[1] << k) ^ right[1]) & mask)
               for left, right in entries):
            result.append(k)
    return result


def transformed(seed, leg, target, source, exponent, scalar=1, constant=False):
    """Apply f=c*w^exponent, optionally plus c, on the two paired legs."""
    q = seed['q_x']
    if seed['q_u'] != q:
        raise ValueError('This frame adapter currently uses a square alphabet.')
    linear, quadratic = frame_terms(seed['values'], q, leg, target, source)
    square = (0, 1, 3, 2)[scalar]
    offset = max(0, -2 * exponent)
    divisor_mask = (1 << offset) - 1
    values = list(seed['values'])
    for index in set(linear) | set(quadratic):
        left = scale(linear.get(index, (0, 0)), scalar)
        right = scale(quadratic.get(index, (0, 0)), square)
        old = values[index]
        result = []
        for component in (0, 1):
            value = ((old[component] << offset)
                     ^ (left[component] << (offset + exponent))
                     ^ (right[component] << (offset + 2 * exponent)))
            if constant:
                value ^= (left[component] ^ right[component]) << offset
            if value & divisor_mask:
                raise ArithmeticError('The paired Laurent frame leaves an actual tensor pole.')
            result.append(value >> offset)
        values[index] = tuple(result)
    for index, (a, b) in enumerate(values):
        x, y, u, v = ((index // q ** leg_number) % q for leg_number in range(4))
        if (a.bit_count() % 2, b.bit_count() % 2) != (int(x == v and y == u), 0):
            raise ArithmeticError('The paired physical frame did not preserve the strict endpoint.')
    result = dict(seed)
    result['values'] = values
    result['degree'] = max(degree(value) for value in values)
    result['monomial_content'] = min(valuation(value) for value in values if value != (0, 0))
    result['frame'] = {'leg': leg, 'target': target, 'source': source,
                       'exponent': exponent, 'scalar': scalar, 'constant': constant}
    return result


def primitive(seed):
    """Remove complete possible new scalar content (a power of w).

    Paired Laurent frames are invertible at every nonzero finite place.
    Thus a primitive starting seed cannot acquire any other common factor.
    """
    content = min(valuation(value) for value in seed['values'] if value != (0, 0))
    result = dict(seed)
    result['values'] = [(a >> content, b >> content) for a, b in seed['values']]
    result['degree'] = seed['degree'] - content
    result['removed_monomial_content'] = content
    return result


def bell_normalize(seed):
    """Normalize an exactly factored constant tensor Q_xy (Q^-1)_uv."""
    from f16_identity_operator_profiles import inverse_matrix
    from physical_split_crossings import identity_split, retract

    q = seed['q_x']
    if seed['q_u'] != q:
        raise ValueError('Use a square physical alphabet.')
    constants = [(a & 1) + 2 * (b & 1) for a, b in seed['values']]
    pivot = next(i for i, value in enumerate(constants) if value)
    x0, y0, u0, v0 = ((pivot // q ** leg) % q for leg in range(4))
    if constants[pivot] != 1:
        raise ValueError('This exact Bell normalizer currently expects a unit pivot equal to one.')
    first = [[constants[x + q * (y + q * (u0 + q * v0))] for y in range(q)] for x in range(q)]
    second = [[constants[x0 + q * (y0 + q * (u + q * v))] for v in range(q)] for u in range(q)]
    for index, value in enumerate(constants):
        x, y, u, v = ((index // q ** leg) % q for leg in range(4))
        if value != (first[x][y] & second[u][v]):
            # The current reflection constants lie in F2. Refuse a wider
            # field input rather than using Boolean multiplication on it.
            raise ArithmeticError('The constant tensor does not have the asserted F2 product factorization.')
    if any(value not in (0, 1) for row in first + second for value in row):
        raise ValueError('The reflected constant factors must lie in F2.')
    first_inverse = inverse_matrix(first)
    if [list(row) for row in first_inverse] != second:
        raise ArithmeticError('The two constant Bell factors are not inverse matrices.')
    result = retract(seed, [list(row) for row in first_inverse], [list(row) for row in zip(*first)],
                     identity_split(q), identity_split(q))
    for index, (a, b) in enumerate(result['values']):
        x, y, u, v = ((index // q ** leg) % q for leg in range(4))
        if (a & 1, b & 1) != (int(x == y and u == v), 0):
            raise ArithmeticError('The normalized constant tensor is not the exact Bell product.')
    return result, {'first_factor': first, 'inverse_factor': second}


def bell_restore(seed, factors):
    from physical_split_crossings import identity_split, retract

    q = seed['q_x']
    return retract(seed, factors['first_factor'],
                   [list(row) for row in zip(*factors['inverse_factor'])],
                   identity_split(q), identity_split(q))


def coupled_terms(seed, target, source, scalar=1):
    """Exact coefficients in f of A_x A^-T_y A_u A^-T_v, A=I+fE.

    The result retains all four-leg interactions. No intermediate
    polynomiality assertion is made when f will be a negative monomial.
    """
    q = seed['q_x']
    if seed['q_u'] != q or not (0 <= target < q and 0 <= source < q) or target == source:
        raise ValueError('Use a square alphabet and two distinct letters.')
    output = [{i: value for i, value in enumerate(seed['values']) if value != (0, 0)}]
    output += [{} for _ in range(4)]
    for leg, row, column in ((0, target, source), (1, source, target),
                             (2, target, source), (3, source, target)):
        stride = q ** leg
        for order in range(3, -1, -1):
            for index, value in list(output[order].items()):
                if (index // stride) % q != column:
                    continue
                destination = index + (row - column) * stride
                changed = add(output[order + 1].get(destination, (0, 0)), scale(value, scalar))
                if changed != (0, 0):
                    output[order + 1][destination] = changed
                elif destination in output[order + 1]:
                    del output[order + 1][destination]
    return output


def coupled_admissible_poles(terms, maximum):
    """Test the complete four-leg polynomial before removing its shift."""
    keys = set().union(*(set(row) for row in terms[1:]))
    if not keys:
        return []
    rows = [[terms[order].get(index, (0, 0)) for order in range(1, 5)] for index in keys]
    allowed = []
    for k in range(1, maximum + 1):
        mask = (1 << (4 * k)) - 1
        valid = True
        for row in rows:
            low = high = 0
            for order, (a, b) in enumerate(row, 1):
                low ^= a << ((4 - order) * k)
                high ^= b << ((4 - order) * k)
            if (low | high) & mask:
                valid = False
                break
        if valid:
            allowed.append(k)
    return allowed


def coupled_transformed(seed, target, source, exponent, scalar=1, terms=None):
    """Apply all four legs with a common global Laurent shift, then test."""
    terms = coupled_terms(seed, target, source, scalar) if terms is None else terms
    shift = max(0, -4 * exponent)
    mask = (1 << shift) - 1
    values = [(0, 0)] * len(seed['values'])
    for order, coefficients in enumerate(terms):
        offset = shift + order * exponent
        for index, (a, b) in coefficients.items():
            before = values[index]
            values[index] = before[0] ^ (a << offset), before[1] ^ (b << offset)
    if any((a | b) & mask for a, b in values):
        raise ArithmeticError('The complete coupled frame retains an actual pole.')
    values = [(a >> shift, b >> shift) for a, b in values]
    q = seed['q_x']
    for index, (a, b) in enumerate(values):
        x, y, u, v = ((index // q ** leg) % q for leg in range(4))
        if (a.bit_count() % 2, b.bit_count() % 2) != (int(x == v and y == u), 0):
            raise ArithmeticError('The coupled physical frame changed the strict endpoint.')
    result = dict(seed, values=values, degree=max(degree(value) for value in values))
    result['coupled_frame'] = {'target': target, 'source': source, 'exponent': exponent, 'scalar': scalar}
    result['monomial_content'] = min(valuation(value) for value in values if value != (0, 0))
    return result


def digest(values):
    return sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest()


def catalogue(interval):
    """Exact catalogue of nontrivial degree-preserving single paired frames."""
    seed = tensor(interval)
    q, bound = seed['q_x'], seed['degree']
    reversals = reverse_values(seed['values'], bound)
    output = []
    for leg in (0, 1):
        for target in range(q):
            for source in range(q):
                if target == source:
                    continue
                for place, values in (('zero', seed['values']), ('infinity', reversals)):
                    linear, quadratic = frame_terms(values, q, leg, target, source)
                    for scalar in (1, 2, 3):
                        allowed = admissible_poles(linear, quadratic, bound, scalar)
                        if allowed:
                            output.append({'leg': leg, 'target': target, 'source': source,
                                           'place': place, 'scalar': scalar, 'pole_orders': allowed})
    return seed, output
