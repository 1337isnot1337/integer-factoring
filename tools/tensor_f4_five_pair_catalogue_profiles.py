#!/usr/bin/env python3
"""Complete five-column F4 orthogonal profile classes at a chosen cutoff.

The exact coefficient catalogue has twelve ordered kernels modulo row
operations, paired into six Frobenius classes. This driver computes every
zero-cut profile. The infinity profile is supplied by the proved orthogonal
equal-profile identity and is explicitly marked as derived, not recomputed.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from coherent_differential_kernels_f4 import differential_kernel
from coherent_hot_compounds_f4 import compound_block, tensor_vector, validate_maps
from orthogonal_differential_pairs_f4 import orthogonal_retraction
from search_coherent_hot_retractions import pressure


def orthogonal_zero_profile(code, left):
    right, modes = code['embedding'], code['modes']
    hs, charges = validate_maps(right, left, modes)
    q = len(hs)
    orders, blocks = [], []
    start = time.monotonic()
    for charge in sorted({a + b for a in charges for b in charges}):
        pairs = [(x, y) for x in range(q) for y in range(q) if charges[x] + charges[y] == charge]
        rows = [tensor_vector(left[x], left[y], modes) for x, y in pairs]
        columns = [tensor_vector(right[y], right[x], modes) for x, y in pairs]
        block = compound_block(modes, rows, columns)
        blocks.append(block)
        orders.extend(block['orders_t'])
    assert len(orders) == q * q
    orders.sort()
    return {'q': q, 'modes': modes, 'heights': hs, 'charges': charges,
            'zero_orders_t': orders, 'infinity_orders_t': orders[:],
            'computed_cuts': ['zero'], 'computed_smith_orders': q * q,
            'infinity_profile_source': 'orthogonal equal-profile theorem; not separately recomputed',
            'zero_blocks': blocks, 'seconds': time.monotonic() - start}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--height', type=int, default=20)
    ap.add_argument('--catalogue', type=Path, default=Path('evidence/tensor-f4-rank-two-catalogue.json'))
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    catalogue = json.loads(args.catalogue.read_text())
    assert catalogue['passed'] and catalogue['ordered_kernels_modulo_row_operations'] == 12
    matrices = catalogue['frobenius_profile_representatives']
    assert len(matrices) == 6
    pairs = [(i, 9 - i) for i in range(5)]
    names = [Path(__file__).name, 'coherent_differential_kernels_f4.py',
             'coherent_hot_compounds_f4.py', 'coherent_hot_compounds_f4.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp',
             'orthogonal_differential_pairs_f4.py', 'search_coherent_hot_retractions.py']
    report = {'scope': __doc__, 'height_cap': args.height, 'modes': 14, 'pairs': pairs,
              'catalogue_path': args.catalogue.as_posix(),
              'catalogue_sha256': sha256(args.catalogue.read_bytes()).hexdigest(),
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'records': [], 'completed': False, 'in_progress': None}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()

    def save():
        report['elapsed_seconds'] = time.monotonic() - start
        temporary = args.output.with_suffix(args.output.suffix + '.tmp')
        temporary.write_text(json.dumps(report, separators=(',', ':')) + '\n')
        temporary.replace(args.output)

    for index, C in enumerate(matrices):
        operators = [{(1 << a) | (1 << b): row[i] for i, (a, b) in enumerate(pairs) if row[i]}
                     for row in C]
        code = differential_kernel(operators, 14, args.height)
        left = orthogonal_retraction(code)
        assert left is not None
        report['in_progress'] = {'index': index, 'coefficients': C, 'q': code['q']}
        save()
        print(json.dumps({'starting': index, 'q': code['q']}), flush=True)
        profile = orthogonal_zero_profile(code, left)
        fit = pressure(profile)
        report['records'].append({'index': index, 'coefficients': C, 'code': code,
                                  'orthogonal_retraction': left, 'profile': profile, 'pressure': fit})
        report['in_progress'] = None
        save()
        print(json.dumps({'completed': index, 'rho': fit['rho'], 'seconds': profile['seconds']}), flush=True)
    report['completed'] = True
    save()


if __name__ == '__main__':
    main()
