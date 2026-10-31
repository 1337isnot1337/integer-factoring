#!/usr/bin/env python3
"""Exact homogeneous Boolean joint filters and their analytic cost bounds.

All floating helpers evaluate the proved formulas; they are not interval
certificates.  The literal alpha/beta construction is over F2.
"""
from itertools import product
from math import exp, expm1, log, log1p, sqrt

import numpy as np
from scipy.special import logsumexp


def truth_mask(modes, monomials):
    value = 0
    for monomial in monomials:
        value ^= sum(1 << state for state in range(1 << modes)
                     if state & monomial == monomial)
    return value


def filter_data(modes, n, ramification, monomials, final_degree):
    monomials = tuple(monomials)
    if (modes < 1 or n < modes or n & (n - 1)
            or ramification < 1 or ramification % 2 == 0
            or not monomials or len(set(monomials)) != len(monomials)
            or any(value <= 0 or value >= 1 << modes for value in monomials)):
        raise ValueError('Use distinct nonconstant Boolean monomials and admissible binary parameters.')
    degree = monomials[0].bit_count()
    if any(value.bit_count() != degree for value in monomials):
        raise ValueError('The Boolean coefficient pattern must be homogeneous.')
    m = n * ramification
    if final_degree <= m * degree:
        raise ValueError('The final polynomial degree must exceed m times the Boolean degree.')
    mask = truth_mask(modes, monomials)
    states = tuple(product(range(2), repeat=modes))
    size = len(states)
    reversed_alpha = tuple(1 ^ (((1 << (m * degree)) ^ (1 << final_degree))
                                  if (mask >> state) & 1 else 0)
                           for state in range(size))
    alpha = tuple(reversed(reversed_alpha))
    orders = tuple(m * state.bit_count() + (final_degree - m * degree
                                          if (mask >> state) & 1 else 0)
                   for state in range(size))
    degrees = tuple(value.bit_length() - 1 for value in alpha)
    return {'states': states, 'degree': degree, 'truth_mask': mask,
            'reversed_alpha': reversed_alpha, 'alpha': alpha,
            'beta_orders': orders, 'degrees': degrees}


def log_support_probability(energies, scale, mask):
    states = np.asarray(tuple(product(range(2), repeat=len(energies))), dtype=float)
    weights = -scale * (states @ np.asarray(energies, dtype=float))
    support = np.asarray([(mask >> i) & 1 for i in range(len(states))], dtype=bool)
    return float(logsumexp(weights[support]) - logsumexp(weights))


def relative_formula(n, ramification, selected, charge, scale, degree, final_degree, mask):
    energies = 2 * n - 1 - 2 * np.asarray(selected, dtype=float)
    log_p = log_support_probability(energies, scale, mask)
    log_q = log_support_probability(energies + charge, scale, mask)
    log_a = -2 * scale * (final_degree - n * ramification * degree) / ramification
    log_h = 2 * n * scale * degree
    def complement(value):
        return log(-expm1(value)) if value < 0 else -np.inf
    return float(np.logaddexp(complement(log_q), log_q + log_a)
                 + np.logaddexp(complement(log_p), log_p + log_h - log_a))


def scale_threshold(n):
    """Unique positive t satisfying exp(2nt)(1-exp(-t))=1."""
    lower, upper = 0.0, log(2)
    for _ in range(80):
        middle = (lower + upper) / 2
        if middle < log1p(exp(-(2 * n - 1) * middle)):
            lower = middle
        else:
            upper = middle
    return (lower + upper) / 2


def distinct_index_log_lower_bound(n, degree, scale):
    return (scale * degree * degree
            - sum(log1p(exp(-scale * (2 * n - 1 - 2 * i))) for i in range(degree)))


def bounded_charge_small_scale_lower_bound(n, degree, charge_bound):
    if n < 3 or degree < 1 or charge_bound < 0:
        raise ValueError('Use n>=3, positive Boolean degree, and a nonnegative charge bound.')
    upper = log(2 * n) / (2 * n - 1)
    first = sqrt(2) * (n + 1)
    tail = 2 * log(2) * exp(-first * upper) * (-expm1(-(2 * n - first) * upper)) / upper
    return tail - degree * log(2) - n * upper * upper * charge_bound * charge_bound / (4 * sqrt(2))
