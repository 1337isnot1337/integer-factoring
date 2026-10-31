#!/usr/bin/env python3
"""Floating full pressures from exact joint beta valuations and degrees.

Selected indices may be real for exploratory optimization.  Actual
crossings must replay at distinct integer indices in the proved range.
"""
from itertools import product
from math import fsum, sqrt

import numpy as np
from scipy.special import logsumexp

from prime_state_filter_pressure import gaussian_pressure


def relative_log_cost(prime, n, ramification, charge, scale, selected, beta_orders, degrees):
    states = np.asarray(tuple(product(range(prime), repeat=len(selected))), dtype=float)
    selected = np.asarray(selected, dtype=float)
    a = 2 * selected + 1 - charge
    kappa = 2 * selected + 1 - 2 * n
    zero = 2 * np.asarray(beta_orders) / ramification - states @ a
    infinity = states @ kappa - 2 * np.asarray(degrees) / ramification
    gaussian_zero = states @ (2 * n - a)
    gaussian_infinity = states @ kappa
    minima = [float(np.min(value)) for value in (zero, infinity, gaussian_zero, gaussian_infinity)]
    ground = fsum((-minima[0], -minima[1], minima[2], minima[3]))
    thermal = fsum(sign * float(logsumexp(-(value - minimum) * scale))
                   for sign, value, minimum in zip((1, 1, -1, -1),
                                                  (zero, infinity, gaussian_zero, gaussian_infinity), minima))
    return ground * scale + thermal


def pressure(prime, n, ramification, charge, scale, selected, beta_orders, degrees,
             rho=1 / sqrt(2)):
    return (gaussian_pressure(prime, n, charge, scale, rho)
            + relative_log_cost(prime, n, ramification, charge, scale,
                                selected, beta_orders, degrees))


def index_gradient(prime, n, ramification, charge, scale, selected, beta_orders, degrees):
    """Relative log cost and exact floating gradient in the selected indices."""
    states = np.asarray(tuple(product(range(prime), repeat=len(selected))), dtype=float)
    selected = np.asarray(selected, dtype=float)
    a = 2 * selected + 1 - charge
    kappa = 2 * selected + 1 - 2 * n
    arrays = (2 * np.asarray(beta_orders) / ramification - states @ a,
              states @ kappa - 2 * np.asarray(degrees) / ramification,
              states @ (2 * n - a), states @ kappa)
    means = []
    for array in arrays:
        weights = np.exp(-scale * array - logsumexp(-scale * array))
        means.append(weights @ states)
    gradient = 2 * scale * (means[0] - means[1] - means[2] + means[3])
    return relative_log_cost(prime, n, ramification, charge, scale, selected, beta_orders, degrees), gradient


def compiled_binary_beta_orders(n, ramification, modes, exponents):
    """Compile complete beta polynomials for alpha_d=color_d*x^(-exponent_d)."""
    from prime_power_identity_filters import finite_difference_terms
    states = tuple(product(range(2), repeat=modes))
    lookup = {state: index for index, state in enumerate(states)}
    if len(exponents) != len(states):
        raise ValueError('Use one exponent per joint state.')
    offset = max(exponents)
    basis = []
    for k, terms in finite_difference_terms(2, modes).items():
        row = []
        for state, power, scalar in terms:
            index = lookup[state]
            polynomial = 0
            for degree in range(power + 1):
                if degree & power == degree:
                    polynomial ^= 1 << (n * ramification * degree - exponents[index] + offset)
            row.append((index, polynomial))
        basis.append(tuple(row))

    def evaluate(colors):
        if len(colors) != len(states) or any(value <= 0 for value in colors):
            raise ValueError('Use one nonzero binary-extension color per joint state.')
        bits = max(colors).bit_length()
        result = []
        for row in basis:
            planes = [0] * bits
            for index, polynomial in row:
                for bit in range(bits):
                    if (colors[index] >> bit) & 1:
                        planes[bit] ^= polynomial
            union = 0
            for plane in planes:
                union |= plane
            assert union
            result.append((union & -union).bit_length() - 1 - offset)
        return tuple(result)

    return evaluate
