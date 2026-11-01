#!/usr/bin/env python3
"""Bounded exact bivariate trace-span control at t=3, tag (0,0)."""
import argparse
import hashlib
import json
import math
import os
import time
from pathlib import Path

from check_trace_bilinear_kernel import Field, irreducibility_by_division, rank
from check_two_patch_span import cut_ranks, dense_cut_ranks, shortest_paths


def bivariate_span():
    t, polynomial, patch, coordinate = 3, 11, 0, 0
    d = 1 << t
    field = Field(polynomial)
    irreducibility_by_division(polynomial)
    powers = [field.power(field.alpha, j) for j in range(t)]
    coordinates = [sum(field.trace(field.mul(a, z)) << j
                       for j, a in enumerate(powers)) for z in range(d)]
    assert len(set(coordinates)) == d
    coefficients = [0] * (4*t)
    for s in range(d):
        for w in range(d):
            error = s ^ w
            b_values = (field.mul(s, error), field.mul(s ^ 1, error))
            rest_polynomial, tagged_bit = 1, None
            for side, value in enumerate(b_values):
                for j in range(t):
                    bit = coordinates[value] >> j & 1
                    if (side, j) == (patch, coordinate):
                        tagged_bit = bit
                    elif bit:
                        rest_polynomial ^= rest_polynomial << 1
            assert tagged_bit in (0, 1)
            assert rest_polynomial.bit_length() <= 2*t
            endpoint = (rest_polynomial.bit_count() & 1) * (1 ^ tagged_bit)
            assert endpoint == int(s == w)
            for j in range(2*t):
                if rest_polynomial >> j & 1:
                    coefficients[j] |= 1 << (s*d+w)
                    if tagged_bit:
                        coefficients[2*t+j] |= 1 << (s*d+w)
    assert coefficients[-1] == 0
    basis, basis_indices = [], []
    for index, kernel in enumerate(coefficients):
        if rank(basis+[kernel]) > len(basis):
            basis.append(kernel)
            basis_indices.append(index)
    assert len(basis) <= 4*t-1
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
    assert len(tensors) <= 2048 and len(set(tensors)) == len(tensors)
    assert tensors[target] == target_kernel
    return tensors, target, basis_indices, coefficient_masks


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
        'status': 'partial', 't': 3, 'field_polynomial': 11,
        'tagged_patch': 0, 'tagged_coordinate': 0, 'labels': 8,
        'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sources},
        'progress': {'block_rank_tensors': 0, 'dense_rank_tensors': 0},
        'scope': 'Bivariate F2 coefficient span at t=3 and tag (0,0) only. No other tags, widths, trace bases, or extension-field recombinations.',
    }

    def save():
        temporary = args.output.with_suffix(args.output.suffix+'.tmp')
        with temporary.open('w') as out:
            json.dump(evidence, out, sort_keys=True, separators=(',', ':'))
            out.write('\n')
            out.flush()
            os.fsync(out.fileno())
        temporary.replace(args.output)

    def check_deadline():
        if time.monotonic() > deadline:
            raise TimeoutError('bounded bivariate-span time limit reached')

    save()
    try:
        check_deadline()
        tensors, target, basis_indices, coefficient_masks = bivariate_span()
        evidence.update({
            'span_dimension': len(basis_indices),
            'span_tensors_expected': len(tensors),
            'basis_coefficient_indices': basis_indices,
            'coefficient_index_convention': 'A_0,...,A_5,B_0,...,B_5',
            'coefficient_masks_in_basis': coefficient_masks,
            'target_coefficient_mask': target,
        })
        ranks, cache = [], {0: 0}
        for index, kernel in enumerate(tensors):
            check_deadline()
            ranks.append(cut_ranks(kernel, 8, cache))
            evidence['progress']['block_rank_tensors'] = index+1
        evidence['rank_pairs_by_mask'] = ranks
        assert ranks[target] == [64, 64]
        save()
        for index, kernel in enumerate(tensors):
            check_deadline()
            assert dense_cut_ranks(kernel, 8) == ranks[index]
            evidence['progress']['dense_rank_tensors'] = index+1
        evidence['dense_cut_comparisons'] = 2*len(tensors)
        save()
        denominator = 1
        while True:
            check_deadline()
            costs = [math.isqrt(r*c*denominator**2) for r, c in ranks]
            distances, previous = shortest_paths(costs)
            check_deadline()
            path, state = [], target
            while state:
                parent, edge = previous[state]
                path.append(edge)
                state = parent
            assert len(path) <= len(basis_indices)
            upper, path_mask = 0, 0
            for edge in path:
                path_mask ^= edge
                r, c = ranks[edge]
                square = r*c*denominator**2
                floor = math.isqrt(square)
                upper += floor+int(floor*floor != square)
            assert path_mask == target
            evidence.update({
                'lower_cost_denominator': denominator,
                'lower_cost_numerators_by_mask': costs,
                'lower_potential_by_mask': distances,
                'potential_inequalities_checked': len(tensors)*(len(tensors)-1),
                'lower_optimum_numerator': distances[target],
                'upper_decomposition_masks': path,
                'upper_decomposition_rank_pairs': [ranks[edge] for edge in path],
                'upper_cost_numerator': upper,
            })
            save()
            if distances[target] == 64*denominator:
                assert upper == 64*denominator
                evidence['optimization_status'] = 'exact_optimum_d_squared'
                break
            if upper < 64*denominator:
                evidence['optimization_status'] = 'strict_improvement_over_one_piece'
                break
            if denominator >= 1 << 16:
                evidence['optimization_status'] = 'unresolved_interval'
                break
            denominator *= 2
    except TimeoutError:
        evidence['status'] = 'partial_time_limit'
        save()
        print(json.dumps({'status': evidence['status'],
                          'progress': evidence['progress']}))
        raise SystemExit(2)
    evidence['status'] = 'passed'
    save()
    saved = json.loads(args.output.read_text())
    assert saved['progress']['block_rank_tensors'] == saved['span_tensors_expected']
    assert saved['progress']['dense_rank_tensors'] == saved['span_tensors_expected']
    print(json.dumps({k: saved[k] for k in (
        'status', 'span_dimension', 'span_tensors_expected',
        'dense_cut_comparisons', 'potential_inequalities_checked',
        'optimization_status', 'lower_cost_denominator',
        'lower_optimum_numerator', 'upper_cost_numerator',
        'upper_decomposition_masks')}))


if __name__ == '__main__':
    main()
