#!/usr/bin/env python3
"""Strict Gaussian crossings and affine sums of boundary relabelings.

All constructors use exact F4 polynomial arithmetic. Floating pressure
queries are exploratory; a strict improvement requires exact replay.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product
import math

import numpy as np
from scipy.special import logsumexp

from check_rectangular_f4_repair import ZERO, add, scale, evaluate_one
from check_schur_downset_crossing import sparse_kernel
from check_schur_gaussian_reduction import polynomial_multiply
from native_matrix_smith_f4 import MatrixOracle


def binary_inverse(matrix):
    count = len(matrix)
    work = [sum(value << j for j, value in enumerate(row)) | (1 << (count + i))
            for i, row in enumerate(matrix)]
    for column in range(count):
        pivot = next((row for row in range(column, count) if work[row] >> column & 1), None)
        if pivot is None:
            raise ArithmeticError('The endpoint Pascal matrix is singular.')
        work[column], work[pivot] = work[pivot], work[column]
        for row in range(count):
            if row != column and work[row] >> column & 1:
                work[row] ^= work[column]
    return [[work[row] >> (count + column) & 1 for column in range(count)]
            for row in range(count)]


@lru_cache(None)
def gaussian_tensor(indices, ambient, strict=True):
    """Return coefficients indexed by x + q*y + q^2*u + q^3*v.

    The raw endpoint is delta(y,u) times the exterior Pascal matrix on
    (v,x). Multiplying the one-particle v rows by its inverse repairs all
    exterior degrees simultaneously and introduces no denominator.
    """
    indices = tuple(indices)
    count = len(indices)
    kernel = sparse_kernel(indices, ambient)
    one_particle = [row[:2 * count] for row in kernel[2 * count:]]
    pascal = [[int(b & ~a == 0) for b in indices] for a in indices]
    if strict:
        inverse = binary_inverse(pascal)
        repaired = []
        for row in inverse:
            repaired.append([0] * (2 * count))
            for position, coefficient in enumerate(row):
                if coefficient:
                    repaired[-1] = [a ^ b for a, b in zip(repaired[-1], one_particle[position])]
        one_particle[:count] = repaired

    @lru_cache(None)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        bit = rows & -rows
        row = bit.bit_length() - 1
        value = 0
        remaining = columns
        while remaining:
            chosen = remaining & -remaining
            remaining ^= chosen
            column = chosen.bit_length() - 1
            value ^= polynomial_multiply(one_particle[row][column], minor(rows ^ bit, columns ^ chosen))
        return value

    q = 1 << count
    answer = tuple(minor(v | (y << count), u | (x << count))
                   for v, u, y, x in product(range(q), repeat=4))
    if strict:
        for v, u, y, x in product(range(q), repeat=4):
            value = answer[x | (y << count) | (u << (2 * count)) | (v << (3 * count))]
            assert (value.bit_count() & 1) == int(x == v and y == u)
    return answer


def relabel(values, permutation):
    q = len(permutation)
    if sorted(permutation) != list(range(q)) or q & (q - 1):
        raise ValueError('Supply a permutation of a power-of-two alphabet.')
    count = q.bit_length() - 1
    return tuple(values[permutation[x] | (permutation[y] << count)
                        | (permutation[u] << (2 * count)) | (permutation[v] << (3 * count))]
                 for v, u, y, x in product(range(q), repeat=4))


def substitute_power(polynomial, degree):
    if degree < 1:
        raise ValueError('Substitution degree must be positive.')
    result = 0
    while polynomial:
        bit = polynomial & -polynomial
        polynomial ^= bit
        result |= 1 << (degree * (bit.bit_length() - 1))
    return result


def affine_relabeling(base, permutation, scalar=2, degree=1, vanishing=False, monomial=0):
    """C(w^d)+gamma*w^k*[1+w if requested]*(C^pi(w^d)+C(w^d)).

    A negative k is cleared by a common monomial. Its zero and infinity
    shifts cancel in the product of the two profile partition functions.
    """
    transformed = relabel(base, permutation)
    common = max(0, -monomial)
    factor = (3 if vanishing else 1) << (monomial + common)
    result = []
    for original, conjugate in zip(base, transformed):
        original = substitute_power(original, degree)
        # Form the difference before substitution: ``original`` is already substituted.
        difference = original ^ substitute_power(conjugate, degree)
        correction = scale((polynomial_multiply(difference, factor), 0), scalar)
        result.append(add((original << common, 0), correction))
    return tuple(result)


def affine_sum(base, permutations, coefficients):
    if len(permutations) != len(coefficients) or not permutations:
        raise ValueError('One coefficient is required for each relabeling.')
    total = 0
    for coefficient in coefficients:
        total ^= coefficient
    if total != 1:
        raise ValueError('Affine coefficients must sum to one in F4.')
    result = [ZERO] * len(base)
    for permutation, coefficient in zip(permutations, coefficients):
        source = relabel(base, permutation)
        result = [add(value, scale((other, 0), coefficient)) for value, other in zip(result, source)]
    return tuple(result)


def charge_projection(base, charges, modulus=3):
    """Project onto conservation of an arbitrary physical-letter charge.

    For modulus three this is the sum of the three F4 diagonal-phase
    conjugates. The returned coefficients are already in F2. The exact
    crossing endpoint conserves every supplied charge automatically.
    """
    q = len(charges)
    if q & (q - 1) or modulus < 2:
        raise ValueError('Use a power-of-two alphabet and modulus at least two.')
    count = q.bit_length() - 1
    result = []
    for mask, value in enumerate(base):
        x, y, u, v = ((mask >> (count * leg)) & (q - 1) for leg in range(4))
        keep = (charges[x] + charges[u] - charges[y] - charges[v]) % modulus == 0
        result.append((value if keep else 0, 0))
    return tuple(result)


def check_endpoint(values, count):
    q = 1 << count
    for v, u, y, x in product(range(q), repeat=4):
        value = values[x | (y << count) | (u << (2 * count)) | (v << (3 * count))]
        assert evaluate_one(value) == int(x == v and y == u)
    return q ** 4


def flatten(values, count, infinity):
    q = 1 << count
    matrix = [[ZERO] * (q * q) for _ in range(q * q)]
    for mask, value in enumerate(values):
        x, y, u, v = ((mask >> (count * leg)) & (q - 1) for leg in range(4))
        row, column = (q * x + u, q * y + v) if infinity else (q * x + y, q * u + v)
        matrix[row][column] = value
    return matrix


def natural_heights(indices, ambient):
    indices = tuple(indices)
    count = len(indices)
    q = 1 << count
    zero = ([2 * a for a in indices], [2 * a + 2 - 2 * ambient for a in indices],
            [2 * ambient - 1] * count, [1 - 2 * ambient] * count)
    infinity = ([2 * ambient - 1] * count, [1 - 2 * ambient] * count,
                [2 * a for a in indices], [2 * a + 2 - 2 * ambient for a in indices])
    return np.array([[[sum(leg[a] for a in range(count) if letter >> a & 1)
                       for letter in range(q)] for leg in place]
                     for place in (zero, infinity)], dtype=float)


class RelabelingPressure:
    def __init__(self, values, indices, ambient, degree=1, entry_shifts=None):
        self.indices = tuple(indices)
        self.count = len(self.indices)
        self.q = 1 << self.count
        self.ambient = ambient
        self.degree = degree
        self.endpoint_entries = check_endpoint(values, self.count)
        if entry_shifts is None:
            entry_shifts = [0] * len(values)
        if len(entry_shifts) != len(values) or any(not isinstance(shift, int) for shift in entry_shifts):
            raise ValueError('Supply an integer Laurent exponent for every tensor entry.')
        self.entry_shifts = tuple(entry_shifts)
        self.shift_matrices = [flatten(self.entry_shifts, self.count, infinity) for infinity in (False, True)]
        support = [(value, shift) for value, shift in zip(values, entry_shifts) if value != ZERO]
        minimum = min(min(valuation for plane in value if plane
                          for valuation in [(plane & -plane).bit_length() - 1]) + shift
                      for value, shift in support)
        self.clearing_order = max(0, -minimum)
        maximum_degree = max(max(value).bit_length() - 1 for value in values)
        spread = max(entry_shifts) - min(entry_shifts)
        precision = max(128, 1 << max(0, maximum_degree + spread + 64).bit_length())
        if precision > 8192:
            raise ValueError('The polynomial degree exceeds the selected finite profile engine.')
        self.matrices = [flatten(values, self.count, infinity) for infinity in (False, True)]
        self.oracles = [MatrixOracle(matrix, self.q ** 2, shifts, precision=precision)
                        for matrix, shifts in zip(self.matrices, self.shift_matrices)]
        self.natural = natural_heights(self.indices, ambient)
        self.queries = 0
        self.maximum_degree = max(max(value).bit_length() - 1 + shift + self.clearing_order
                                  for value, shift in support)

    def evaluate(self, log_s, offsets=None, gradient=False):
        if not math.isfinite(log_s) or log_s <= 0:
            raise ValueError('Use positive finite log(s).')
        offsets = np.zeros_like(self.natural) if offsets is None else np.asarray(offsets, dtype=float)
        if offsets.shape != self.natural.shape:
            raise ValueError('Offsets have shape (two places, four legs, alphabet).')
        heights = (self.natural + offsets) * log_s
        log_b = 0.0
        profiles = []
        profile_gradient = np.zeros_like(heights)
        scale_derivative = 0.0
        for place, oracle in enumerate(self.oracles):
            h = heights[place]
            first, second = ((0, 1), (2, 3)) if place == 0 else ((0, 2), (1, 3))
            rows = (h[first[0]][:, None] + h[first[1]][None, :]).reshape(-1)
            columns = (h[second[0]][:, None] + h[second[1]][None, :]).reshape(-1)
            while True:
                try:
                    profile = oracle.profile(rows, columns, scale=log_s / self.degree, infinity=bool(place))
                    break
                except ArithmeticError as error:
                    if 'status 3' not in str(error) or oracle.precision >= 8192:
                        raise
                    oracle = MatrixOracle(self.matrices[place], self.q ** 2, self.shift_matrices[place],
                                          precision=2 * oracle.precision)
                    self.oracles[place] = oracle
            orders = np.asarray(profile['orders'])
            log_partition = float(logsumexp(-orders))
            log_b += log_partition
            probabilities = np.exp(-orders - log_partition)
            chosen_rows = np.asarray(profile['rows'])
            chosen_columns = np.asarray(profile['cols'])
            for leg, labels in ((first[0], chosen_rows // self.q), (first[1], chosen_rows % self.q),
                                (second[0], chosen_columns // self.q), (second[1], chosen_columns % self.q)):
                profile_gradient[place, leg] -= np.bincount(labels, weights=probabilities, minlength=self.q)
            scale_derivative -= float(probabilities @ orders) / log_s
            profiles.append(profile)
        gx = heights[:, 0, :].sum(axis=0) + heights[:, 3, :].sum(axis=0)
        gy = heights[:, 1, :].sum(axis=0) + heights[:, 2, :].sum(axis=0)

        def equation(r):
            return log_b - r * (float(logsumexp(-gx / r)) + float(logsumexp(-gy / r)))

        lower, upper = 1e-10, 4.0
        if equation(lower) <= 0:
            raise ArithmeticError('Nonpositive pressure root requires a separate audit.')
        while equation(upper) > 0:
            upper *= 2
            if upper > 1024:
                raise ArithmeticError('Could not bracket the pressure root.')
        for _ in range(60):
            middle = (lower + upper) / 2
            if equation(middle) > 0:
                lower = middle
            else:
                upper = middle
        r = (lower + upper) / 2
        self.queries += 1
        px = np.exp(-gx / r - logsumexp(-gx / r))
        py = np.exp(-gy / r - logsumexp(-gy / r))
        details = {'log_s_before_substitution': log_s,
                   'log_s_actual': log_s / self.degree,
                   'log_profile_product': log_b,
                   'orders': [profile['orders'] for profile in profiles],
                   'type_x': px.tolist(), 'type_y': py.tolist(),
                   'maximum_polynomial_degree': self.maximum_degree,
                   'common_monomial_clearing_order': self.clearing_order,
                   'profile_precisions': [oracle.precision for oracle in self.oracles]}
        if not gradient:
            return r / 2, details
        entropy = -sum(float(p[p > 0] @ np.log(p[p > 0])) for p in (px, py))
        if entropy < 1e-12:
            raise ArithmeticError('The optimized type entropy is too small for a stable derivative.')
        physical = np.stack([np.stack([px, py, py, px])] * 2)
        offsets_gradient = log_s * (profile_gradient + physical) / (2 * entropy)
        scale_derivative += (float(px @ gx) + float(py @ gy)) / log_s
        log_s_gradient = scale_derivative / (2 * entropy)
        details['entropy'] = entropy
        return r / 2, (offsets_gradient, log_s_gradient), details
