#!/usr/bin/env python3
"""Largest regular left physical projection preserving a scalar terminal line.

This computes a lower certificate for full-word regular-channel pressure.
It does not compute an upper rank bound or a type-restricted rank lower bound.
All rational arithmetic is over F_2(tau), with exact DVR regularity checks.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product

import adaptive_smith_lattices as adaptive
import check_cancellation_module_automata as arithmetic
import check_exact_module_orbits as lattices
import check_minimal_boundary_modules as linear


ZERO, ONE = arithmetic.ZERO, arithmetic.ONE


def saturated_kernel(matrix, width):
    """Return an O-basis of ker(matrix: O^width -> O^height).

    The regular Schur right factor has a unit maximal minor.  Solving on
    that minor leaves arbitrary regular coordinates in the free columns,
    so the returned vectors generate the whole saturated kernel.
    """
    assert all(len(row) == width for row in matrix)
    assert all(lattices.regular(value) for row in matrix for value in row)
    factors = adaptive.rational_schur_factors(matrix)
    right = tuple(factor[2] for factor in factors)
    rank = len(right)
    if not rank:
        return tuple(tuple(ONE if i == j else ZERO for j in range(width))
                     for i in range(width))

    pivots, independent = [], {}
    for column in range(width):
        residue = sum((arithmetic.jet(right[i][column], 1) & 1) << i
                      for i in range(rank))
        while residue:
            leading = residue.bit_length() - 1
            if leading in independent:
                residue ^= independent[leading]
            else:
                independent[leading] = residue
                pivots.append(column)
                break
    assert len(pivots) == rank
    minor = tuple(tuple(row[column] for column in pivots) for row in right)
    inverse = linear.inverse(minor)
    assert all(lattices.regular(value) for row in inverse for value in row)

    result = []
    for free in range(width):
        if free in pivots:
            continue
        coefficients = linear.column_product(inverse, tuple(row[free] for row in right))
        vector = [ZERO] * width
        vector[free] = ONE
        for column, value in zip(pivots, coefficients):
            vector[column] = value
        assert all(lattices.regular(value) for value in vector)
        assert all(value == ZERO for value in linear.column_product(matrix, vector))
        result.append(tuple(vector))
    assert len(result) == width - rank
    return tuple(result)


def terminal_scalar_projection(local, degree, *, alphabet=2, cut=0, side=0,
                               infinity=False, initial=None, terminal=None):
    """Compute the maximal regular left scalar quotient at one endpoint.

    The input follows the ordered-crossing convention of the adaptive
    engines.  Physical gauges may be applied before this call; in that
    case use infinity=False for the already transformed endpoint tensor
    and restore its external scalar shift when interpreting the pressure.

    For fixed boundaries f,e with f.e=1, set S=fKe.  Stack the physical
    residual matrices Ke-eS, then annihilate them on the left.  The rows
    of P span every regular left projection satisfying PKe=(PS)e.
    """
    dimension = len(next(iter(local.values())))
    assert alphabet >= 2 and degree >= 0 and cut in (0, 1) and side in (0, 1)
    assert set(local) == set(product(range(alphabet), repeat=4))
    initial = tuple(initial or (ONE,) + (ZERO,) * (dimension - 1))
    terminal = tuple(terminal or (ONE,) + (ZERO,) * (dimension - 1))
    assert len(initial) == len(terminal) == dimension
    assert all(value in (ZERO, ONE) for value in initial + terminal)
    assert linear.scalar_product(initial, terminal) == ONE
    width = alphabet * alphabet
    scalar = [[ZERO] * width for _ in range(width)]
    terminal_action = [[[ZERO] * width for _ in range(width)]
                       for _ in range(dimension)]
    residual_transpose = [[ZERO] * width for _ in range(dimension * width)]

    for (x, y, u, v), matrix in local.items():
        assert len(matrix) == dimension and all(len(row) == dimension for row in matrix)
        assert all(isinstance(value, int) and 0 <= value < 1 << (degree + 1)
                   for row in matrix for value in row)
        left, right = ((alphabet * x + y, alphabet * u + v) if cut == 0
                       else (alphabet * x + u, alphabet * y + v))
        physical, other = (left, right) if side == 0 else (right, left)
        if infinity:
            transformed = tuple(tuple(arithmetic.rat(sum(
                1 << (degree - k) for k in range(value.bit_length()) if value >> k & 1
            )) for value in row) for row in matrix)
        else:
            transformed = tuple(tuple(arithmetic.rat(value) for value in row)
                                for row in matrix)
        column = linear.column_product(transformed, terminal)
        entry = linear.scalar_product(initial, column)
        scalar[physical][other] = entry
        for i, value in enumerate(column):
            terminal_action[i][physical][other] = value
            residual_transpose[i * width + other][physical] = arithmetic.radd(
                value, arithmetic.rmul(terminal[i], entry))

    projection = saturated_kernel(tuple(tuple(row) for row in residual_transpose), width)
    projected = linear.matrix_product(projection, scalar) if projection else ()
    identity_checks = 0
    for i in range(dimension):
        actual = linear.matrix_product(projection, terminal_action[i]) if projection else ()
        expected = tuple(tuple(arithmetic.rmul(terminal[i], value) for value in row)
                         for row in projected)
        assert actual == expected
        identity_checks += len(projection) * width
    factors = adaptive.rational_schur_factors(projected) if projection else ()
    return {
        "projection": projection,
        "scalar_matrix": projected,
        "orders": tuple(factor[0] for factor in factors),
        "physical_width": width,
        "projection_rows": len(projection),
        "scalar_rank": len(factors),
        "regular_terminal_identity_entries": identity_checks,
        "infinity": infinity,
        "degree": degree,
    }


def scalar_pressure(certificate, point):
    """The exact lower rate, in the endpoint convention of the input call."""
    assert isinstance(point, Fraction) and point > 1
    shift = certificate["degree"] if certificate["infinity"] else 0
    return sum((point ** (shift - order) for order in certificate["orders"]), Fraction(0))


def boundary_scalar_projections(local, degree, *, alphabet=2, cut=0, side=0,
                                infinity=False, initial=None, terminal=None,
                                maximum_constant_lines=255):
    """Screen every eligible constant auxiliary line, using bounded boundaries.

    A line must be reachable from the actual terminal and visible from the
    actual initial row.  Its scalar projection lower-bounds every uniform
    exponential full-word rate, not necessarily the count at the same length.
    Rational linear combinations of bounded suffix words can have different
    lengths; those bounded shifts do not change a uniform exponential rate.

    Only constant binary lines are enumerated.  This is not a search over
    every F_2(tau)-line.  Exceeding the declared enumeration budget is
    inconclusive, rather than a negative certificate.
    """
    dimension = len(next(iter(local.values())))
    if (1 << dimension) - 1 > maximum_constant_lines:
        raise adaptive.InconclusiveSearch("constant auxiliary line budget exceeded")
    initial = tuple(initial or (ONE,) + (ZERO,) * (dimension - 1))
    terminal = tuple(terminal or (ONE,) + (ZERO,) * (dimension - 1))
    assert len(initial) == len(terminal) == dimension
    assert all(value in (ZERO, ONE) for value in initial + terminal)
    assert linear.scalar_product(initial, terminal) == ONE
    matrices = tuple(tuple(tuple(arithmetic.rat(value) for value in row)
                           for row in matrix) for matrix in local.values())
    reachable, pivots = linear.closure(terminal, matrices, column=True)
    observable, _ = linear.closure(initial, matrices)
    certificates, excluded = [], []
    for mask in range(1, 1 << dimension):
        vector = tuple(ONE if mask >> i & 1 else ZERO for i in range(dimension))
        coordinates = tuple(vector[pivot] for pivot in pivots)
        is_reachable = linear.linear_combination(coordinates, reachable) == vector
        is_visible = any(linear.scalar_product(row, vector) != ZERO for row in observable)
        if not is_reachable or not is_visible:
            excluded.append({"line": mask, "reachable": is_reachable, "visible": is_visible})
            continue
        first = next(i for i, value in enumerate(vector) if value != ZERO)
        dual = tuple(ONE if i == first else ZERO for i in range(dimension))
        certificate = terminal_scalar_projection(
            local, degree, alphabet=alphabet, cut=cut, side=side, infinity=infinity,
            initial=dual, terminal=vector)
        certificates.append({"line": mask, "certificate": certificate})
    return {
        "scope": "uniform full-word exponential rates; constant auxiliary lines only",
        "reachable_dimension": len(reachable),
        "observable_dimension": len(observable),
        "constant_lines_examined": (1 << dimension) - 1,
        "certificates": certificates,
        "excluded_lines": excluded,
    }


def boundary_rate_lower_bound(screen, point):
    """Largest exact rate obstruction returned by boundary_scalar_projections."""
    return max((scalar_pressure(item["certificate"], point)
                for item in screen["certificates"]), default=Fraction(0))
