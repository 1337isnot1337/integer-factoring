#!/usr/bin/env python3
"""Compare the exact hot degree specialization with the generic constructor."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from check_coherent_q32_pressure import displayed_maps
from coherent_hot_compounds import profiles as generic_profiles
from coherent_hot_compounds_fast import profiles as fast_profiles
from coherent_inverse_systems import inverse_system


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    code = inverse_system([(1 << 6) ^ (1 << 9)], 4)
    cases = [('inverse_six', code['embedding'], code['retraction'], 4),
             ('inverse_six_transpose', code['retraction'], code['embedding'], 4)]
    for coherent in (False, True):
        embedding, retraction = displayed_maps(coherent)
        cases.append((f'q32_coherent_{coherent}', embedding, retraction, 8))
    counts, records = Counter(), []
    for name, embedding, retraction, modes in cases:
        generic = generic_profiles(embedding, retraction, modes, keep_blocks=True)
        fast = fast_profiles(embedding, retraction, modes, keep_blocks=True)
        assert generic['zero_orders_t'] == fast['zero_orders_t']
        assert generic['infinity_orders_t'] == fast['infinity_orders_t']
        counts['complete_profile_comparisons'] += 2
        for old_cut, new_cut in zip(generic['cuts'], fast['cuts']):
            for old_block, new_block in zip(old_cut['blocks'], new_cut['blocks']):
                assert old_block['pairs'] == new_block['pairs']
                assert old_block['orders_t'] == new_block['orders_t']
                assert old_block.pop('matrix') == new_block.pop('matrix')
                counts['complete_polynomial_comparisons'] += old_block['size'] ** 2
                counts['complete_smith_order_comparisons'] += old_block['size']
        records.append({'name': name, 'modes': modes, 'q': len(embedding),
                        'generic_seconds': generic['seconds'], 'fast_seconds': fast['seconds'],
                        'generic_precisions': [[b['precision_t'] for b in c['blocks']]
                                                for c in generic['cuts']],
                        'fast_precisions': [[b['precision_t'] for b in c['blocks']]
                                             for c in fast['cuts']],
                        'fast_result': fast})
    names = (Path(__file__).name, 'coherent_hot_compounds_fast.py', 'coherent_hot_compounds_fast.cpp',
             'coherent_hot_compounds.py', 'coherent_hot_compounds.cpp', 'coherent_inverse_systems.py',
             'check_coherent_q32_pressure.py', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp')
    result = {'status': 'pass', 'scope': __doc__, 'counts': dict(counts), 'records': records,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}, 'factoring_exponent_improved': False}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': 'pass', 'counts': dict(counts),
                      'timings': [{key: r[key] for key in ('name', 'generic_seconds', 'fast_seconds')}
                                  for r in records]}, indent=2))


if __name__ == '__main__':
    main()
