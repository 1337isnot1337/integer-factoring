#!/usr/bin/env python3
"""Exact finite weighted adaptive ranks and a tensor-power amplification gate.

Input complete Smith profiles and physical gauges require their own exact
local certificate. This module checks the remaining finite arithmetic,
including physical type means, every positive jet, full dimension caps,
zero pieces, and the rank-table amplification penalty. It does not infer
valid local profiles from floating optimization output.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from math import isqrt
import json
from pathlib import Path

from weighted_adaptive_mixed_place import GaugeProfile, catalogue_flag, type_class_size


def gauge_from_heights(orders, ramification, heights, xcounts, ycounts):
    """Build the exact mean for integral, already ramified physical heights."""
    if len(heights) != 4:
        raise ValueError('There must be one height list per physical leg.')
    copies = sum(xcounts)
    if copies < 1 or sum(ycounts) != copies:
        raise ValueError('Both physical types must have the same positive tensor length.')
    counts = (xcounts, ycounts, ycounts, xcounts)
    if any(len(h) != len(c) or any(not isinstance(v, int) for v in h)
           for h, c in zip(heights, counts)):
        raise ValueError('Each ramified height must be an integer on the correct alphabet.')
    mean = Fraction(sum(value * count for h, c in zip(heights, counts)
                        for value, count in zip(h, c)), copies)
    return GaugeProfile(tuple(orders), ramification, mean)


def ceiling_square_root(value, precision):
    """Return the numerator of ceil(2^precision sqrt(value))."""
    if not isinstance(value, int) or value < 0 or not isinstance(precision, int) or precision < 0:
        raise ValueError('A nonnegative integer and nonnegative binary precision are required.')
    scaled = value << (2 * precision)
    root = isqrt(scaled)
    return root + int(root * root != scaled)


def amplification_bounds(rank_pairs, external_product, exponent=(140, 99), precision=64):
    """Exact sufficient balanced gates with a rational exponent below sqrt(2)."""
    pairs = [tuple(pair) for pair in rank_pairs]
    if any(len(pair) != 2 or any(not isinstance(value, int) or value < 0 for value in pair)
           for pair in pairs):
        raise ValueError('Every cut rank bound must be a nonnegative integer.')
    if not isinstance(external_product, int) or external_product < 2:
        raise ValueError('The product external label count must be an integer at least two.')
    a, b = exponent
    if (not isinstance(a, int) or not isinstance(b, int) or min(a, b) < 1
            or a * a >= 2 * b * b):
        raise ValueError('Use a positive rational exponent strictly below sqrt(2).')
    active = [r * c for r, c in pairs if r and c]
    if not active:
        raise ValueError('A strict crossing cannot have every piece certified zero.')
    count, largest = len(active), max(active)
    coarse = largest * count * count
    denominator = 1 << precision
    numerator = sum(ceiling_square_root(value, precision) for value in active)
    coarse_gate = pow(coarse, b) < pow(external_product, a)
    sharp_gate = pow(numerator, 2 * b) < pow(denominator, 2 * b) * pow(external_product, a)
    return {'active_piece_upper_bound': count, 'zero_piece_bounds': len(pairs) - count,
            'maximum_rank_product': largest, 'coarse_base': coarse,
            'square_root_sum_upper_numerator': numerator,
            'square_root_sum_upper_denominator': denominator,
            'rational_sqrt2_lower_bound': [a, b], 'exponent_square_test': a * a < 2 * b * b,
            'coarse_gate_passed': coarse_gate, 'rank_sum_gate_passed': sharp_gate,
            'certified_rho_strict_upper_bound_if_passed': [a, 2 * b]}


def finite_certificate(row_catalogues, column_catalogues, copies, degree,
                       xcounts, ycounts, exponent=(140, 99), precision=64):
    """Every finite prefix and rank product for fixed paired physical types."""
    if (not isinstance(copies, int) or copies < 1 or not isinstance(degree, int) or degree < 0
            or sum(xcounts) != copies or sum(ycounts) != copies):
        raise ValueError('Use positive length, nonnegative degree, and exact common-length types.')
    dx, dy = type_class_size(xcounts), type_class_size(ycounts)
    if min(dx, dy) < 2:
        raise ValueError('Both finite external alphabets must contain at least two labels.')
    channels = len(xcounts) * len(ycounts)
    if any(len(gauge.orders) != channels
           for side in (row_catalogues, column_catalogues) for place in side for gauge in place):
        raise ValueError('Each profile must contain all qX*qY generic cut channels.')
    dimension, total_degree = dx * dy, copies * degree
    row = catalogue_flag(row_catalogues, copies, total_degree, dimension)
    column = catalogue_flag(column_catalogues, copies, total_degree, dimension)
    pairs = [(row[k]['prefix_cost'], column[total_degree - k]['prefix_cost'])
             for k in range(total_degree + 1)]
    gate = amplification_bounds(pairs, dimension, exponent, precision)
    critical = max(range(total_degree + 1), key=lambda k: pairs[k][0] * pairs[k][1])
    return {'copies': copies, 'seed_degree': degree, 'tensor_degree': total_degree,
            'type_counts_x': list(xcounts), 'type_counts_y': list(ycounts),
            'external_label_counts': [dx, dy], 'cut_dimension': dimension,
            'critical_index': critical, 'rank_pairs': pairs, 'amplification': gate,
            'row_flag': row, 'column_flag': column}


def read_catalogues(data):
    """Parse externally certified complete local records, with disjoint labels."""
    xcounts, ycounts = data['type_counts_x'], data['type_counts_y']
    catalogues, seen = [], set()
    for side in ('row_places', 'column_places'):
        result = []
        for place in data[side]:
            label = json.dumps(place['place'], sort_keys=True)
            if label in seen:
                raise ValueError('Each labelled geometric place may occur on only one side.')
            seen.add(label)
            result.append([gauge_from_heights(g['orders'], g['ramification'], g['heights'], xcounts, ycounts)
                           for g in place['gauges']])
        catalogues.append(result)
    return catalogues


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--exponent-numerator', type=int, default=140)
    parser.add_argument('--exponent-denominator', type=int, default=99)
    parser.add_argument('--sqrt-bits', type=int, default=64)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    data = json.loads(raw)
    row, column = read_catalogues(data)
    result = finite_certificate(row, column, sum(data['type_counts_x']), data['degree'],
                                data['type_counts_x'], data['type_counts_y'],
                                (args.exponent_numerator, args.exponent_denominator), args.sqrt_bits)
    names = (Path(__file__).name, 'weighted_adaptive_mixed_place.py', 'adaptive_mixed_place.py')
    result.update(scope=__doc__, input_sha256=sha256(raw).hexdigest(),
                  local_profile_provenance=data.get('local_profile_provenance'),
                  source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                 for name in names})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['amplification']))


if __name__ == '__main__':
    main()
