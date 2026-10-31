#!/usr/bin/env python3
"""Floating pressure evaluations for the exact prime-state filter profiles.

Ground energies are separated from the bounded thermal terms. This
avoids cancellation between large log partition functions when a common
Laurent monomial clears the filter. These values guide searches; a
negative value still needs an independent rigorous pressure certificate.

The scale is t=P*log(s). beta_orders and infinity_active are exact z
orders from prime_state_filter_profiles or polynomial synthesis.
"""
from __future__ import annotations

from math import fsum, sqrt
import numpy as np
from scipy.special import logsumexp


def digit_thermal(prime, energies, scale):
    energies = np.asarray(energies, dtype=float)
    return logsumexp(-np.abs(energies)[..., None] * scale * np.arange(prime), axis=-1)


def gaussian_pressure(prime, n, charge, scale, rho=1 / sqrt(2)):
    """Full paired-type Gaussian pressure at the given charge and scale."""
    if n < 1 or scale <= 0 or rho <= 0:
        raise ValueError('Use a positive interval length, scale, and exponent.')
    ratio = 2 * rho
    indices = np.arange(n)
    zero = np.concatenate((2 * indices + 1 - charge, 2 * indices + 1 + charge))
    infinity = 2 * np.arange(2 * n) + 1 - 2 * n
    physical = np.concatenate((4 * indices - 2 * (n - 1) + charge,
                               4 * indices - 2 * (n - 1) - charge))
    ground = (prime - 1) * fsum(
        [min(float(value), 0) for value in physical]
        + [-min(float(value), 0) for value in zero]
        + [-min(float(value), 0) for value in infinity])
    thermal = fsum(digit_thermal(prime, zero, scale))
    thermal += fsum(digit_thermal(prime, infinity, scale))
    thermal -= ratio * fsum(digit_thermal(prime, physical, scale / ratio))
    return ground * scale + thermal


def relative_log_cost(prime, n, ramification, charge, scale,
                      beta_orders, infinity_active):
    """Exact-profile log-cost difference from the same Gaussian background."""
    if len(beta_orders) != prime or len(infinity_active) != prime:
        raise ValueError('One finite-difference order and infinity order per digit is required.')
    digits = np.arange(prime)
    zero = np.asarray(beta_orders, dtype=float) / ramification - (1 - charge) * digits
    infinity = np.asarray(infinity_active, dtype=float) / ramification
    gaussian_zero = (2 * n - 1 + charge) * digits
    gaussian_infinity = (1 - 2 * n) * digits
    minima = [float(np.min(values)) for values in
              (zero, infinity, gaussian_zero, gaussian_infinity)]
    ground = fsum((-minima[0], -minima[1], minima[2], minima[3]))
    thermal = 0.0
    for sign, values, minimum in zip((1, 1, -1, -1),
                                    (zero, infinity, gaussian_zero, gaussian_infinity), minima):
        thermal += sign * float(logsumexp(-(values - minimum) * scale))
    return ground * scale + thermal


def filter_pressure(prime, n, ramification, charge, scale,
                    beta_orders, infinity_active, rho=1 / sqrt(2)):
    return (gaussian_pressure(prime, n, charge, scale, rho)
            + relative_log_cost(prime, n, ramification, charge, scale,
                                beta_orders, infinity_active))
