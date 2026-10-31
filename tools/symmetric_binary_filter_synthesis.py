#!/usr/bin/env python3
"""Exact occupation-symmetric polynomial filters in only r+1 weight states.

The output f_w is alpha_{1-k} for any occupation state k of weight w.
Thus degrees are already in the order used by the physical expectation.
Coefficient integers represent additive binary coordinates in any fixed
extension field large enough to contain every nonzero endpoint color.
"""

from functools import lru_cache

from polynomial_state_filter_synthesis import residue_patterns
from weyl_lattice_filter_synthesis import endpoint_image


def submasks(mask):
    current = mask
    while True:
        yield current
        if not current:
            return
        current = (current - 1) & mask


@lru_cache(maxsize=16)
def difference_columns(modes):
    """Column h contributes C^(w-h) when binomial(w,h) is odd."""
    return tuple(tuple((weight, tuple(submasks(weight ^ source)))
                       for weight in range(source, modes + 1)
                       if weight & source == source)
                 for source in range(modes + 1))


@lru_cache(maxsize=16384)
def local_endpoint_image(modes, degrees, orders):
    variables = tuple((weight, degree) for weight, bound in enumerate(degrees)
                      for degree in range(max(0, bound + 1)))
    masks = [0] * (modes + 1)
    equations = {}
    columns = difference_columns(modes)
    for column, (source, degree) in enumerate(variables):
        bit = 1 << column
        masks[source] |= bit
        for weight, shifts in columns[source]:
            for shift in shifts:
                at = degree + shift
                if at < orders[weight]:
                    key = weight, at
                    equations[key] = equations.get(key, 0) ^ bit
    pivots = {}
    for equation in equations.values():
        while equation:
            pivot = equation.bit_length() - 1
            if pivot in pivots:
                equation ^= pivots[pivot]
            else:
                pivots[pivot] = equation
                break
    return variables, endpoint_image(variables, pivots, masks)


def transform(polynomials, m):
    result = []
    for weight in range(len(polynomials)):
        value = {}
        for source in submasks(weight):
            for shift in submasks(weight ^ source):
                for degree, coefficient in polynomials[source].items():
                    at = degree + m * shift
                    value[at] = value.get(at, 0) ^ coefficient
        result.append({degree: coefficient for degree, coefficient in value.items() if coefficient})
    return tuple(result)


def synthesize(modes, m, degrees, orders, construct=True):
    size = modes + 1
    degrees, orders = tuple(degrees), tuple(orders)
    if (modes < 1 or m < 1 or len(degrees) != size or len(orders) != size
            or any(not isinstance(value, int) or value < 0 for value in degrees + orders)):
        raise ValueError('Use positive modes and m, and one nonnegative integer bound per weight.')
    basis = {}
    patterns = residue_patterns(m, degrees, orders)
    for residue, local_degrees, local_orders in patterns:
        variables, image = local_endpoint_image(modes, local_degrees, local_orders)
        for endpoint, coefficients in image:
            lift = {(weight, residue + m * degree)
                    for column, (weight, degree) in enumerate(variables)
                    if (coefficients >> column) & 1} if construct else set()
            while endpoint:
                pivot = endpoint.bit_length() - 1
                if pivot in basis:
                    old, old_lift = basis[pivot]
                    endpoint ^= old
                    if construct:
                        lift.symmetric_difference_update(old_lift)
                else:
                    basis[pivot] = endpoint, lift
                    break
    support = 0
    for endpoint, lift in basis.values():
        support |= endpoint
    if support != (1 << size) - 1:
        return None
    if not construct:
        return True

    field_size = 1 << size.bit_length()
    endpoint = [0] * size
    polynomial = {}
    for vector, lift in basis.values():
        union = tuple(index for index in range(size) if endpoint[index] or (vector >> index) & 1)
        if all(endpoint[index] for index in union):
            continue
        scalar = next(value for value in range(field_size)
                      if all(endpoint[index] ^ (value if (vector >> index) & 1 else 0)
                             for index in union))
        for index in range(size):
            endpoint[index] ^= scalar if (vector >> index) & 1 else 0
        for key in lift:
            polynomial[key] = polynomial.get(key, 0) ^ scalar
            if not polynomial[key]:
                del polynomial[key]
        if all(endpoint):
            break
    assert all(endpoint)
    values = tuple({degree: coefficient for (weight, degree), coefficient in polynomial.items()
                    if weight == index} for index in range(size))
    beta = transform(values, m)
    actual_degrees = tuple(max(value) for value in values)
    actual_orders = tuple(min(value) for value in beta)
    assert all(actual <= bound for actual, bound in zip(actual_degrees, degrees))
    assert all(actual >= bound for actual, bound in zip(actual_orders, orders))
    assert all(value and 0 < color < field_size for value, color in zip(beta, endpoint))
    return {
        'modes': modes, 'm': m, 'coefficient_field_order': field_size,
        'degree_bounds': degrees, 'order_bounds': orders,
        'actual_degrees': actual_degrees, 'actual_orders': actual_orders,
        'endpoint_values': tuple(endpoint), 'endpoint_image_dimension': len(basis),
        'residue_pattern_count': len(patterns), 'f': values, 'beta': beta,
    }
