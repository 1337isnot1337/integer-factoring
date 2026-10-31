#!/usr/bin/env python3
"""Exact complementary Gaussian columns aligned with the degree-18 q5 seed.

Complementary four-particle states are inserted only in the x boundary.
The other three boundaries retain the five-state reflection code.  The
strict endpoint therefore makes every extracted slab vanish at one.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import time

from reflection_dual_rail_gaussian import reflection_letters
from saturated_dual_kernel_directions import divide, polynomial_gcd
from search_f4_recentered_crossings import digest, five_letter_source
from selected_gaussian_coefficients import selected_tensor


def construct(degree=18, progress=None):
    code = tuple(reflection_letters(8)[letter] for letter in (0, 2, 4, 8, 12))
    reference = five_letter_source()[0]
    old_reference = tuple(sum(((value >> k) & 1) << (16+2*k) for k in range(19))
                          for value, high in reference)
    raw_slabs = []
    primitives = []
    zero_slabs = []
    baseline_entries = endpoint_entries = 0
    for positions in combinations(range(8), 4):
        complement = sum(1 << position for position in positions)
        if complement in code:
            continue
        source = selected_tensor(tuple(range(8)), 8, code+(complement,), code)
        values = source['values']
        slab = []
        for v in range(5):
            for u in range(5):
                for y in range(5):
                    for x in range(5):
                        value, high = values[x+6*(y+5*(u+5*v))]
                        assert high == 0 and value == old_reference[x+5*(y+5*(u+5*v))]
                        baseline_entries += 1
                    value, high = values[5+6*(y+5*(u+5*v))]
                    assert high == 0 and value.bit_count() % 2 == 0
                    assert all(not ((value >> k) & 1) for k in range(1, value.bit_length(), 2))
                    compressed = sum(((value >> (2*k)) & 1) << k
                                     for k in range((value.bit_length()+1)//2))
                    slab.append(compressed)
                    endpoint_entries += 1
        if not any(slab):
            zero_slabs.append(complement)
            continue
        common = 0
        for value in slab:
            common = polynomial_gcd(common, value)
        normalized = []
        for value in slab:
            quotient, remainder = divide(value, common)
            assert not remainder
            normalized.append(quotient ^ (quotient << 1))
        item = {'complement_mask': complement, 'complement_modes': positions,
                'raw_slab_in_z': slab, 'raw_degree': max(value.bit_length()-1 for value in slab),
                'entry_gcd': common, 'entry_gcd_degree': common.bit_length()-1,
                'normalized_slab': normalized, 'normalized_degree': max(value.bit_length()-1 for value in normalized),
                'direct_degree18_direction': all(not (value & 255) for value in slab)
                    and max(value.bit_length()-1 for value in slab) <= 26,
                'native_complete_degree_bound_in_w': source['coefficient_degree_bound'],
                'native_exact_precision': source['exact_precision']}
        raw_slabs.append(item)
        if item['normalized_degree'] <= degree:
            primitives.append({key: value for key, value in item.items() if key != 'raw_slab_in_z'})
        if progress:
            progress(len(raw_slabs)+len(zero_slabs), item)
    descriptors = [(primitive, logical, shift)
                   for primitive, item in enumerate(primitives) for logical in range(5)
                   for shift in range(degree-item['normalized_degree']+1)]
    return {'constructor': 'aligned_selected_gaussian_complement_directions', 'dimensions': [5]*4,
            'degree': degree, 'ambient_modes': 8, 'code_masks': code,
            'reference_values': reference, 'reference_digest': digest(reference),
            'reference_identity': 'G_code(w)=w^16*T5(w^2), checked coefficient by coefficient.',
            'baseline_coefficient_comparisons': baseline_entries,
            'complementary_endpoint_entries': endpoint_entries,
            'complementary_four_particle_masks': 65, 'zero_complement_masks': zero_slabs,
            'raw_slabs': raw_slabs, 'primitives': primitives, 'direction_descriptors': descriptors,
            'direction_count_before_rank_reduction': len(descriptors),
            'scope': 'Exact Gaussian complements at the actual q5 reference; optional gcd normalization gives polynomial endpoint-zero directions outside the original split family.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degree', type=int, default=18)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()

    def progress(count, item):
        if count % 10 == 0:
            print(json.dumps({'event': 'complements', 'count': count, 'normalized_degree': item['normalized_degree'],
                              'seconds': time.monotonic()-started}), flush=True)

    result = construct(args.degree, progress)
    names = (Path(__file__).name, 'selected_gaussian_coefficients.py', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp', 'saturated_dual_kernel_directions.py',
             'search_f4_recentered_crossings.py', 'reflection_dual_rail_gaussian.py')
    result['source_sha256'] = {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}
    result['elapsed_seconds'] = time.monotonic()-started
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'event': 'finished', 'primitives': len(result['primitives']),
                      'directions_before_rank_reduction': len(result['direction_descriptors']),
                      'direct_degree18_complements': sum(item['direct_degree18_direction'] for item in result['raw_slabs']),
                      'normalized_degrees': [item['normalized_degree'] for item in result['primitives']],
                      'seconds': time.monotonic()-started}), flush=True)


if __name__ == '__main__':
    main()
