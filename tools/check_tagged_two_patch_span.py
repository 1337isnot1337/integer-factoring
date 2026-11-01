#!/usr/bin/env python3
"""Bounded exact controls for a single weight-two trace coordinate."""
import argparse
import hashlib
import json
import math
import os
import time
from pathlib import Path

from check_trace_bilinear_kernel import Field, irreducibility_by_division, rank
from check_two_patch_span import cut_ranks, dense_cut_ranks, shortest_paths


def check_deadline(deadline):
    if time.monotonic() > deadline:
        raise TimeoutError('bounded tagged-span time limit reached')


def tagged_kernels(t, polynomial, patch, coordinate):
    field = Field(polynomial)
    assert field.degree == t
    irreducibility_by_division(polynomial)
    d = 1 << t
    powers = [field.power(field.alpha, j) for j in range(t)]
    coordinates = [sum(field.trace(field.mul(a, z)) << j
                       for j, a in enumerate(powers)) for z in range(d)]
    assert len(set(coordinates)) == d
    coefficients = [0] * (2*t+1)
    for s in range(d):
        for w in range(d):
            error = s ^ w
            b_values = (field.mul(s, error), field.mul(s ^ 1, error))
            bits = 1
            for side, value in enumerate(b_values):
                for j in range(t):
                    if coordinates[value] >> j & 1:
                        weight = 2 if (side, j) == (patch, coordinate) else 1
                        bits ^= bits << weight
            assert bits.bit_length() <= 2*t+1
            assert bits.bit_count() % 2 == int(s == w)
            for j in range(2*t+1):
                if bits >> j & 1:
                    coefficients[j] |= 1 << (s*d+w)
    assert coefficients[-1] != 0
    basis, basis_indices = [], []
    for index, kernel in enumerate(coefficients):
        if rank(basis+[kernel]) > len(basis):
            basis.append(kernel)
            basis_indices.append(index)
    echelon = {}
    for index, kernel in enumerate(basis):
        value, mask = kernel, 1 << index
        while value:
            pivot = value.bit_length()-1
            if pivot not in echelon:
                echelon[pivot] = (value, mask)
                break
            row, expression = echelon[pivot]
            value ^= row
            mask ^= expression
        assert value

    def representation(kernel):
        value, mask = kernel, 0
        while value:
            row, expression = echelon[value.bit_length()-1]
            value ^= row
            mask ^= expression
        return mask

    coefficient_masks = [representation(kernel) for kernel in coefficients]
    target_kernel = sum(1 << (s*d+s) for s in range(d))
    target = representation(target_kernel)
    summed_mask = 0
    for mask in coefficient_masks:
        summed_mask ^= mask
    assert summed_mask == target
    tensors = [0] * (1 << len(basis))
    for mask in range(1, len(tensors)):
        bit = mask & -mask
        tensors[mask] = tensors[mask ^ bit] ^ basis[bit.bit_length()-1]
    assert len(set(tensors)) == len(tensors)
    assert tensors[target] == target_kernel
    return tensors, target, basis_indices, coefficient_masks


def profile(t, polynomial, patch, coordinate, cache, deadline):
    d = 1 << t
    tensors, target, basis_indices, coefficient_masks = tagged_kernels(
        t, polynomial, patch, coordinate)
    ranks = []
    for kernel in tensors:
        check_deadline(deadline)
        ranks.append(cut_ranks(kernel, d, cache))
    assert ranks[target] == [d*d, d*d]
    if (patch, coordinate) == (0, 0):
        dense_masks = list(range(len(tensors)))
    else:
        dense_masks = sorted({0, target, 3, 5, 6, *coefficient_masks})
    for mask in dense_masks:
        check_deadline(deadline)
        assert dense_cut_ranks(tensors[mask], d) == ranks[mask]
    denominator = 1
    while True:
        check_deadline(deadline)
        costs = [math.isqrt(r*c*denominator**2) for r, c in ranks]
        distances, previous = shortest_paths(costs)
        path, state = [], target
        while state:
            parent, edge = previous[state]
            path.append(edge)
            state = parent
        assert len(path) <= len(basis_indices)
        upper = 0
        path_mask = 0
        for edge in path:
            path_mask ^= edge
            r, c = ranks[edge]
            square = r*c*denominator**2
            floor = math.isqrt(square)
            upper += floor + int(floor*floor != square)
        assert path_mask == target
        if distances[target] == d*d*denominator:
            assert upper == d*d*denominator
            status = 'exact_optimum_d_squared'
            break
        if upper < d*d*denominator:
            status = 'strict_improvement_over_one_piece'
            break
        if denominator >= 1 << 16:
            status = 'unresolved_interval'
            break
        denominator *= 2
    return {
        't': t, 'field_polynomial': polynomial,
        'tagged_patch': patch, 'tagged_coordinate': coordinate,
        'labels': d, 'degree': 2*t,
        'basis_coefficient_indices': basis_indices,
        'coefficient_masks_in_basis': coefficient_masks,
        'span_dimension': len(basis_indices),
        'span_tensors_checked': len(tensors),
        'target_coefficient_mask': target,
        'rank_pairs_by_mask': ranks,
        'dense_masks_checked': dense_masks,
        'dense_cut_comparisons': 2*len(dense_masks),
        'lower_cost_denominator': denominator,
        'lower_cost_numerators_by_mask': costs,
        'lower_potential_by_mask': distances,
        'potential_inequalities_checked': len(tensors)*(len(tensors)-1),
        'lower_optimum_numerator': distances[target],
        'upper_decomposition_masks': path,
        'upper_decomposition_rank_pairs': [ranks[edge] for edge in path],
        'upper_cost_numerator': upper,
        'status': status,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--time-limit', type=float, default=110)
    args = parser.parse_args()
    deadline = time.monotonic()+args.time_limit
    here = Path(__file__).resolve().parent
    sources = [Path(__file__), here/'check_two_patch_span.py',
               here/'check_trace_bilinear_kernel.py']
    evidence = {
        'status': 'partial', 'profiles': [],
        'expected_profiles': 14,
        'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sources},
        'scope': 'Specified F2 spans for one weight-two tagged trace coordinate, t=3,4. No larger width, other weights, other trace bases, or extension-coefficient exclusion.',
    }

    def save():
        temporary = args.output.with_suffix(args.output.suffix+'.tmp')
        with temporary.open('w') as out:
            json.dump(evidence, out, sort_keys=True, separators=(',', ':'))
            out.write('\n')
            out.flush()
            os.fsync(out.fileno())
        temporary.replace(args.output)

    try:
        for t, polynomial in ((3, 11), (4, 19)):
            cache = {0: 0}
            for patch in range(2):
                for coordinate in range(t):
                    result = profile(t, polynomial, patch, coordinate,
                                     cache, deadline)
                    evidence['profiles'].append(result)
                    save()
                    print(json.dumps({k: result[k] for k in (
                        't', 'tagged_patch', 'tagged_coordinate',
                        'span_dimension', 'span_tensors_checked', 'status',
                        'lower_cost_denominator', 'lower_optimum_numerator',
                        'upper_cost_numerator')}), flush=True)
    except TimeoutError:
        evidence['status'] = 'partial_time_limit'
        save()
        print(json.dumps({'status': evidence['status'],
                          'complete_profiles': len(evidence['profiles'])}))
        raise SystemExit(2)
    assert len(evidence['profiles']) == evidence['expected_profiles']
    evidence['status'] = 'passed'
    save()
    assert json.loads(args.output.read_text())['status'] == 'passed'


if __name__ == '__main__':
    main()
