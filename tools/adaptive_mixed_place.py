#!/usr/bin/env python3
"""Exact adaptive Hasse-jet flags and exploratory profile-rate calculations.

The exact routines use only Python's standard library. Numerical rate curves
are guidance for producing finite rational certificates, never certificates
by themselves.
"""
from __future__ import annotations

from collections import Counter
from functools import lru_cache
from heapq import heappop, heappush
from math import log


@lru_cache(maxsize=256)
def _counts(orders, copies, cutoff):
    if not orders or copies < 1 or cutoff < 0 or min(orders) < 0:
        raise ValueError("nonempty nonnegative profile, positive copies, nonnegative cutoff required")
    counts = Counter(orders)
    # Kronecker packing is exact: every coefficient of every intermediate
    # power is at most len(orders)**copies, strictly below the chosen radix.
    bits = (len(orders) ** copies).bit_length()
    radix_mask = (1 << bits) - 1
    truncation_mask = (1 << (bits * (cutoff + 1))) - 1
    packed = sum(count << (bits * order) for order, count in counts.items()
                 if order <= cutoff)
    result, exponent = 1, copies
    while exponent:
        if exponent & 1:
            result = (result * packed) & truncation_mask
        exponent //= 2
        if exponent:
            packed = (packed * packed) & truncation_mask
    cumulative, output = 0, []
    for _ in range(cutoff + 1):
        cumulative += result & radix_mask
        output.append(cumulative)
        result >>= bits
    return tuple(output)


def channel_counts(orders, copies, cutoff):
    """Counts of channel words with total order <= j, for every 0<=j<=cutoff."""
    return list(_counts(tuple(sorted(orders)), copies, cutoff))


def greedy_flag(profiles, copies, total_degree):
    """The cheapest total_degree+1 admissible initial jets.

    ``profiles`` is a nonempty list of equal-size complete local order lists.
    ``total_degree`` is the degree of the tensor power (usually D*copies).
    Every prefix preserves the initial segment at each physical place.
    The prefix cost is capped by the full cut dimension.
    """
    if not profiles or any(len(p) != len(profiles[0]) for p in profiles):
        raise ValueError("nonempty profiles of equal size required")
    dimension = len(profiles[0]) ** copies
    counts = [channel_counts(p, copies, total_degree) for p in profiles]
    heap = [(c[0], i, 0) for i, c in enumerate(counts)]
    # Entries begin with j=0; a heap built by pushes avoids any assumptions
    # about the order of the individual profile zero counts.
    queue = []
    for entry in heap:
        heappush(queue, entry)
    prefix, result = 0, []
    for _ in range(total_degree + 1):
        cost, place, order = heappop(queue)
        prefix += cost
        result.append({'place_index': place, 'order': order,
                       'cost': cost, 'prefix_cost': min(prefix, dimension)})
        if order < total_degree:
            next_order = order + 1
            heappush(queue, ((next_order + 1) * counts[place][next_order],
                            place, next_order))
    return result


def finite_rank_cost(row_profiles, column_profiles, copies, degree):
    """Exact worst product of the adaptive flag rank upper bounds."""
    full_degree = copies * degree
    row = greedy_flag(row_profiles, copies, full_degree)
    column = greedy_flag(column_profiles, copies, full_degree)
    costs = [row[k]['prefix_cost'] * column[full_degree - k]['prefix_cost']
             for k in range(full_degree + 1)]
    index = max(range(len(costs)), key=costs.__getitem__)
    return {'product': costs[index], 'index': index, 'row': row, 'column': column}


@lru_cache(maxsize=256)
def capacity_curve(orders, entropy_steps=8193, temperature_steps=16385):
    """Numerical inverse lower-tail entropy on a uniform entropy grid.

    Interpolation has no outward rounding guarantee. This routine is used
    only to locate promising allocations for an independent exact check.
    """
    import numpy as np
    from scipy.special import logsumexp

    values, counts = np.unique(np.array(orders, dtype=float), return_counts=True)
    n = len(orders)
    top_entropy = log(n)
    entropies = np.linspace(0, top_entropy, entropy_steps)
    if len(values) == 1:
        capacity = np.full(entropy_steps, float(values[0]))
        return entropies, capacity
    # Log-spaced temperatures resolve both small perturbations from the full
    # channel distribution and the exponentially small low-tail regime.
    temperatures = np.concatenate(([0.0], np.geomspace(1e-7, 80, temperature_steps)))
    logits = np.log(counts)[None, :] - temperatures[:, None] * values[None, :]
    logs = logsumexp(logits, axis=1)
    probabilities = np.exp(logits - logs[:, None])
    means = probabilities @ values
    channel_entropy = logs + temperatures * means
    # Reversal puts the entropy coordinate in increasing order. Repeated
    # floating values at the zero-temperature boundary are removed.
    reversed_entropy = channel_entropy[::-1]
    unique_entropy, indices = np.unique(reversed_entropy, return_index=True)
    reversed_means = means[::-1][indices]
    capacity = np.interp(entropies, unique_entropy, reversed_means,
                         left=values[0], right=float(means[0]))
    return entropies, capacity


def rate_curve(profiles, degree, entropy_steps=8193, rate_steps=4097):
    """Exploratory optimal adaptive prefix rate on 0<=x<=degree."""
    import numpy as np
    if not profiles:
        raise ValueError("nonempty place set required")
    curves = [capacity_curve(tuple(sorted(p)), entropy_steps) for p in profiles]
    entropies = curves[0][0]
    capacity = sum((c[1] for c in curves), np.zeros_like(entropies))
    unique, indices = np.unique(capacity, return_index=True)
    x = np.linspace(0, degree, rate_steps)
    rate = np.interp(x, unique, entropies[indices],
                     left=float(entropies[indices[0]]), right=log(len(profiles[0])))
    free = sum(min(p) for p in profiles)
    rate[x < free - 1e-12] = -np.inf
    # At the end of the free prefix one additional jet is required. Its
    # exponential multiplicity need not be one, including after content shifts.
    rate[np.abs(x - free) < 1e-12] = min(log(p.count(min(p))) for p in profiles)
    return x, rate


def adaptive_rate(row_profiles, column_profiles, degree, alphabet,
                  entropy_steps=8193, rate_steps=4097):
    import numpy as np
    x, row = rate_curve(row_profiles, degree, entropy_steps, rate_steps)
    _, column = rate_curve(column_profiles, degree, entropy_steps, rate_steps)
    total = row + column[::-1]
    index = int(np.argmax(total))
    return {'rho': float(total[index] / (4 * log(alphabet))),
            'critical_prefix': float(x[index]),
            'row_rate': float(row[index]),
            'column_rate': float(column[-index - 1]),
            'entropy_grid_points': entropy_steps,
            'prefix_grid_points': rate_steps,
            'scope': 'floating interpolated rate; not an exact certificate'}
