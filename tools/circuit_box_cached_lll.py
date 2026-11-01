"""Exact cached Gram--Schmidt updates for the existing small-root engine.

The reduction decisions, nearest-integer tie rule, stopping test, and lattice
certificate agree with the reference engine.  The adapters clone function
namespaces; importing this module does not alter either reference module.
"""

from __future__ import annotations

from fractions import Fraction
from types import FunctionType

import circuit_box_centered_windows as _windows
import circuit_box_ratio_small_roots as _roots


def reduce_lattice_cached(rows, threshold=None, *, audit_cache=False):
    """Return exact LLL data, optionally stopping at the reference threshold.

    ``threshold=None`` asks for a fully reduced basis.  Otherwise the first
    row must have strict l1 norm below the supplied threshold on return.
    The optional audit independently recomputes every updated cache.
    """
    if not rows or not rows[0]:
        raise ValueError("a nonempty lattice basis is required")
    dimension, width = len(rows), len(rows[0])
    if any(len(row) != width or any(type(value) is not int for value in row)
           for row in rows):
        raise ValueError("rectangular integer rows are required")
    if threshold is not None and (type(threshold) is not int or threshold <= 0):
        raise ValueError("the optional strict norm threshold must be positive")

    original = [list(row) for row in rows]
    basis = [row[:] for row in original]
    transform = [[int(i == j) for j in range(dimension)] for i in range(dimension)]
    coefficients, norms = _roots._gram_schmidt(basis)
    index, swaps, reductions, audits, updates = 1, 0, 0, 0, 0

    def audit():
        nonlocal audits
        if audit_cache:
            check_coefficients, check_norms = _roots._gram_schmidt(basis)
            if check_coefficients != coefficients or check_norms != norms:
                raise AssertionError("an exact Gram--Schmidt update changed its value")
            audits += 1

    while index < dimension:
        for j in range(index - 1, -1, -1):
            value = coefficients[index][j]
            quotient = (2 * value.numerator + value.denominator) // (2 * value.denominator)
            if quotient:
                basis[index] = [x - quotient * y for x, y in zip(basis[index], basis[j])]
                transform[index] = [x - quotient * y for x, y in zip(transform[index], transform[j])]
                for column in range(j):
                    coefficients[index][column] -= quotient * coefficients[j][column]
                coefficients[index][j] -= quotient
                updates += j + 1
                reductions += 1
                audit()

        mu = coefficients[index][index - 1]
        if norms[index] >= (Fraction(3, 4) - mu * mu) * norms[index - 1]:
            index += 1
        else:
            old_left, old_right = norms[index - 1], norms[index]
            new_left = old_right + mu * mu * old_left
            new_mu = mu * old_left / new_left
            norms[index - 1] = new_left
            norms[index] = old_left * old_right / new_left
            coefficients[index][index - 1] = new_mu

            for column in range(index - 1):
                coefficients[index - 1][column], coefficients[index][column] = (
                    coefficients[index][column], coefficients[index - 1][column]
                )
            for row in range(index + 1, dimension):
                old_mu_left = coefficients[row][index - 1]
                old_mu_right = coefficients[row][index]
                new_mu_right = old_mu_left - mu * old_mu_right
                coefficients[row][index] = new_mu_right
                coefficients[row][index - 1] = old_mu_right + new_mu * new_mu_right
            updates += 3 + 2 * (dimension - 2)

            basis[index], basis[index - 1] = basis[index - 1], basis[index]
            transform[index], transform[index - 1] = transform[index - 1], transform[index]
            swaps += 1
            audit()
            index = max(index - 1, 1)

        # Keep the location of this test identical to the reference engine.
        if threshold is not None and sum(map(abs, basis[0])) < threshold:
            break

    if not any(basis[0]):
        raise AssertionError("an independent lattice produced a zero first row")
    if threshold is not None and sum(map(abs, basis[0])) >= threshold:
        raise AssertionError("the supplied LLL threshold was not attained")
    if any(sum(transform[i][j] * original[j][column] for j in range(dimension))
           != basis[i][column]
           for i in range(dimension) for column in range(width)):
        raise AssertionError("incorrect integer lattice transformation")
    return {
        "basis": basis,
        "transform": transform,
        "lattice_swaps": swaps,
        "lattice_reductions": reductions,
        "full_gram_evaluations": 1,
        "cache_audits": audits,
        "cache_entries_updated": updates,
    }


def short_lattice_vector_cached(basis, threshold):
    result = reduce_lattice_cached(basis, threshold)
    return (result["basis"][0], result["transform"][0],
            result["lattice_swaps"], result["lattice_reductions"])


def _clone(function, namespace):
    result = FunctionType(function.__code__, namespace, function.__name__,
                          function.__defaults__, function.__closure__)
    result.__kwdefaults__ = function.__kwdefaults__
    result.__doc__ = function.__doc__
    return result


_root_namespace = dict(_roots.__dict__)
_root_namespace["_short_lattice_vector"] = short_lattice_vector_cached
residue_short_polynomial = _clone(_roots.residue_short_polynomial, _root_namespace)
_root_namespace["residue_short_polynomial"] = residue_short_polynomial
dyadic_residue_factors = _clone(_roots.dyadic_residue_factors, _root_namespace)
_root_namespace["dyadic_residue_factors"] = dyadic_residue_factors
smallest_factor_by_small_roots = _clone(_roots.smallest_factor_by_small_roots, _root_namespace)

_window_namespace = dict(_windows.__dict__)
_window_namespace["residue_short_polynomial"] = residue_short_polynomial
centered_window_factors = _clone(_windows.centered_window_factors, _window_namespace)
_window_namespace["centered_window_factors"] = centered_window_factors
centered_factor_rectangle = _clone(_windows.centered_factor_rectangle, _window_namespace)
smallest_factor_by_centered_windows = _clone(_windows.smallest_factor_by_centered_windows, _window_namespace)
