#!/usr/bin/env python3
"""Bounded complete-profile searches using the auxiliary Clifford reduction."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
import json
from pathlib import Path
import time

from tensor_affine_clifford_profiles import affine, nilpotent_masks, profile
from tensor_literal_full_tail_affine_frames import pressure


def run_case(parameters):
    M, degree, a_mask, b_mask, low_height = parameters
    started = time.monotonic()
    L, P = affine(a_mask, degree), affine(b_mask, degree)
    cuts = [profile(M, L, P, cut, low_height) for cut in ('zero', 'infinity')]
    fit = pressure(M, cuts[0]['orders'], cuts[1]['orders'], low_height)
    compact = [{k: v for k, v in cut.items() if k not in
                ('orders', 'active_full_mode_indices', 'active_ambient_order', 'active_generic_rank',
                 'auxiliary_smith_multiplicity', 'clifford_matrix_order')} for cut in cuts]
    return {'M': M, 'degree': degree, 'a_mask': a_mask, 'b_mask': b_mask,
            'low_heights': [0, low_height], 'profiles': compact, 'pressure_fit': fit,
            'physical_alphabet_dimension': 1 << M, 'seconds': time.monotonic()-started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes', type=int, nargs='+', default=[16])
    parser.add_argument('--degrees', type=int, nargs='+', default=[1])
    parser.add_argument('--a-masks', type=int, nargs='+')
    parser.add_argument('--b-masks', type=int, nargs='+')
    parser.add_argument('--all-matrices', action='store_true')
    parser.add_argument('--frobenius-representatives', action='store_true')
    parser.add_argument('--low-height', type=int, default=1)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.workers <= 4 or args.low_height < 0:
        raise ValueError('Use one through four workers and a nonnegative low height.')
    defaults = list(range(256)) if args.all_matrices else nilpotent_masks()
    a_masks = defaults if args.a_masks is None else args.a_masks
    b_masks = defaults if args.b_masks is None else args.b_masks
    tasks = [(M, degree, a, b, args.low_height) for M in args.modes for degree in args.degrees
             for a in a_masks for b in b_masks]
    assert len(set(tasks)) == len(tasks)
    normalized_cases = len(tasks)
    if args.frobenius_representatives:
        def frobenius(mask):
            return sum((c if c < 2 else c ^ 1) << (2*i)
                       for i in range(4) for c in [(mask >> (2*i)) & 3])
        if set(a_masks) != {frobenius(a) for a in a_masks} or set(b_masks) != {frobenius(b) for b in b_masks}:
            raise ValueError('Frobenius reduction requires invariant matrix catalogues.')
        tasks = [task for task in tasks if task[2:4] <= (frobenius(task[2]), frobenius(task[3]))]
    names = [Path(__file__).name, 'tensor_affine_clifford_profiles.py',
             'tensor_literal_full_tail_core.py', 'tensor_literal_full_tail_affine_frames.py',
             'check_tensor_hot_compounds_f16.py', 'tensor_hot_compounds_f16.py',
             'tensor_orthogonal_scaling.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'native_weighted_smith_f4.py',
             'native_rectangular_smith_f4.cpp']
    report = {'scope': __doc__, 'parameters': {k: v for k, v in vars(args).items() if k != 'output'},
              'a_masks': a_masks, 'b_masks': b_masks, 'declared_cases': len(tasks),
              'normalized_cases_covered': normalized_cases,
              'complete_profile_encoding': 'minimum_order plus every subset sum of binomial_gaps',
              'orders_per_cut': 256, 'capped_or_omitted_orders': 0,
              'completed': False, 'records': [],
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    started, last_save = time.monotonic(), time.monotonic()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        report['elapsed_seconds'] = time.monotonic()-started
        temporary = args.output.with_suffix(args.output.suffix+'.tmp')
        temporary.write_text(json.dumps(report, separators=(',', ':'))+'\n')
        temporary.replace(args.output)

    def retain(record):
        nonlocal last_save
        report['records'].append(record)
        if len(report['records']) % 1024 == 0 or time.monotonic()-last_save > 10:
            save()
            last_save = time.monotonic()
            best = min(report['records'], key=lambda r: r['pressure_fit']['rho'])
            print(json.dumps({'completed_cases': len(report['records']), 'declared_cases': len(tasks),
                              'best_fit_so_far': best['pressure_fit']['rho'],
                              'best_parameters': [best[k] for k in ('M', 'degree', 'a_mask', 'b_mask')],
                              'elapsed_seconds': report['elapsed_seconds']}), flush=True)

    save()
    if args.workers == 1:
        for task in tasks:
            retain(run_case(task))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(run_case, task) for task in tasks]
            for future in as_completed(futures):
                retain(future.result())
    assert len(report['records']) == len(tasks)
    report['records'].sort(key=lambda r: (r['M'], r['degree'], r['a_mask'], r['b_mask']))
    report['completed'] = True
    save()
    print(json.dumps({'completed': True, 'cases': len(tasks), 'seconds': report['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
