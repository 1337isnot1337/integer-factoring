#!/usr/bin/env python3
"""Exact typed adaptive flags using a catalogue of valid local gauges.

The physical basis and paired physical types are common to all catalogues.
Only the local diagonal gauge used for a rank estimate may change per jet.
Profile/mean-shift validity is supplied by an independent tensor checker.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from heapq import heappop, heappush
from math import comb

from adaptive_mixed_place import channel_counts


@dataclass(frozen=True)
class GaugeProfile:
    orders: tuple[int, ...]
    ramification: int
    mean_shift: Fraction

    def __post_init__(self):
        if not self.orders or any(not isinstance(value, int) for value in self.orders):
            raise ValueError('A complete nonempty list of integral ramified Smith orders is required.')
        if not isinstance(self.ramification, int) or self.ramification < 1:
            raise ValueError('The local ramification is a positive integer.')
        if not isinstance(self.mean_shift, (int, Fraction)):
            raise TypeError('The physical mean shift must be exact.')


def type_class_size(counts):
    if not counts or any(not isinstance(value, int) or value < 0 for value in counts):
        raise ValueError('Nonnegative integer type counts are required.')
    remaining, result = sum(counts), 1
    for count in counts:
        result *= comb(remaining, count)
        remaining -= count
    return result


def jet_costs(profile, copies, total_degree):
    """Sum_w max(L*j + copies*A - sum_w nu + 1,0), for every jet j."""
    if copies < 1 or total_degree < 0:
        raise ValueError('Positive tensor length and nonnegative full degree are required.')
    shift = copies * Fraction(profile.mean_shift)
    if shift.denominator != 1:
        raise ValueError('The tensor length must clear the physical type denominator.')
    minimum = min(profile.orders)
    offset = int(shift) - copies * minimum
    cutoff = profile.ramification * total_degree + offset
    if cutoff < 0:
        return [0] * (total_degree + 1)
    nonnegative = [value - minimum for value in profile.orders]
    cumulative_channels = channel_counts(nonnegative, copies, cutoff)
    integrated, total = [], 0
    for count in cumulative_channels:
        total += count
        integrated.append(total)
    return [integrated[profile.ramification * j + offset]
            if profile.ramification * j + offset >= 0 else 0
            for j in range(total_degree + 1)]


def catalogue_flag(catalogues, copies, total_degree, cut_dimension):
    """Cheapest complete initial-jet flag, with an independent gauge per jet."""
    if not catalogues or any(not gauges for gauges in catalogues) or cut_dimension < 1:
        raise ValueError('Use nonempty local catalogues and a positive exact cut dimension.')
    size = len(catalogues[0][0].orders)
    if any(len(gauge.orders) != size for gauges in catalogues for gauge in gauges):
        raise ValueError('Every gauge must contain the same complete generic cut dimension.')
    costs, selected = [], []
    for gauges in catalogues:
        alternatives = [jet_costs(gauge, copies, total_degree) for gauge in gauges]
        values, indices = [], []
        for j in range(total_degree + 1):
            value, index = min((row[j], i) for i, row in enumerate(alternatives))
            values.append(value)
            indices.append(index)
        if values != sorted(values):
            raise ArithmeticError('The gauge-catalogue jet costs must be nondecreasing.')
        costs.append(values)
        selected.append(indices)
    queue = []
    for place, row in enumerate(costs):
        heappush(queue, (row[0], place, 0))
    prefix, result = 0, []
    for _ in range(total_degree + 1):
        cost, place, order = heappop(queue)
        prefix += cost
        result.append({'place_index': place, 'order': order,
                       'gauge_index': selected[place][order], 'cost': cost,
                       'prefix_cost': min(prefix, cut_dimension)})
        if order < total_degree:
            heappush(queue, (costs[place][order + 1], place, order + 1))
    return result


def finite_rank_cost(row_catalogues, column_catalogues, copies, degree, xcounts, ycounts):
    if sum(xcounts) != copies or sum(ycounts) != copies or len(xcounts) != len(ycounts):
        raise ValueError('One common paired type on each physical alphabet is required.')
    alphabet = len(xcounts)
    if any(len(gauge.orders) != alphabet * alphabet
           for catalogues in (row_catalogues, column_catalogues)
           for gauges in catalogues for gauge in gauges):
        raise ValueError('The profiles must contain every channel of the physical seed cuts.')
    dimension = type_class_size(xcounts) * type_class_size(ycounts)
    total_degree = copies * degree
    row = catalogue_flag(row_catalogues, copies, total_degree, dimension)
    column = catalogue_flag(column_catalogues, copies, total_degree, dimension)
    products = [row[k]['prefix_cost'] * column[total_degree - k]['prefix_cost']
                for k in range(total_degree + 1)]
    index = max(range(total_degree + 1), key=lambda k: products[k])
    return {'product': products[index], 'critical_index': index, 'cut_dimension': dimension,
            'row_flag': row, 'column_flag': column}
