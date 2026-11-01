#!/usr/bin/env python3
"""Complete natural profiles for normalized orthogonal F4 paired kernels.

Every retained trial computes both cuts separately. Floating pressures are
discovery diagnostics, not exact bounds or a new factoring exponent.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

from coherent_differential_kernels_f4 import differential_kernel
from coherent_hot_compounds_f4 import mul, profiles
from orthogonal_differential_pairs_f4 import subset_gram_certificate, orthogonal_retraction
from search_coherent_hot_retractions import pressure


def normalized_scalings():
    """Global nonzero scaling preserves the common kernel; fix its first entry."""
    for tail in product((1, 2, 3), repeat=3):
        scales = (1,) + tail
        coefficients = [[mul(scales[i], 1) for i in range(4)],
                        [mul(scales[i], i) for i in range(4)]]
        certificate = subset_gram_certificate(coefficients)
        if certificate['all_subsets_nondegenerate']:
            yield scales, coefficients, certificate


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--caps', type=int, nargs='+', default=[16, 20])
    ap.add_argument('--modes', type=int, default=12)
    ap.add_argument('--only-scale', default=None)
    ap.add_argument('--max-seconds', type=float, default=None)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    pairs = [(0, 7), (1, 6), (2, 5), (3, 4)]
    names = [Path(__file__).name, 'coherent_differential_kernels_f4.py',
             'coherent_hot_compounds_f4.py', 'coherent_hot_compounds_f4.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp',
             'orthogonal_differential_pairs_f4.py', 'search_coherent_hot_retractions.py']
    scales = list(normalized_scalings())
    assert len(scales) == 6
    if args.only_scale:
        wanted = tuple(map(int, args.only_scale.split(',')))
        scales = [item for item in scales if item[0] == wanted]
        if not scales:
            raise ValueError('The requested normalized scale is not admissible.')
    report = {'scope': __doc__, 'base_commit': '941cc1d8f56e6e80ec45938f1d660532de2f5046',
              'pairs': pairs, 'modes': args.modes, 'caps': args.caps,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'normalized_scales': [list(s) for s, _, _ in scales],
              'records': [], 'completed': False, 'in_progress': None}
    start = time.monotonic()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        report['elapsed_seconds'] = time.monotonic() - start
        temp = args.output.with_suffix(args.output.suffix + '.tmp')
        temp.write_text(json.dumps(report, indent=2) + '\n')
        temp.replace(args.output)

    for cap in args.caps:
        for scaling, C, certificate in scales:
            if args.max_seconds is not None and time.monotonic() - start >= args.max_seconds:
                report['stop_reason'] = 'time budget between complete trials'
                save()
                return
            operators = [{(1 << a) | (1 << b): row[i] for i, (a, b) in enumerate(pairs) if row[i]}
                         for row in C]
            code = differential_kernel(operators, args.modes, cap)
            left = orthogonal_retraction(code)
            assert left is not None
            report['in_progress'] = {'height_cap': cap, 'scaling': scaling, 'q': code['q']}
            save()
            print(json.dumps({'starting': report['in_progress']}), flush=True)
            result = profiles(code['embedding'], left, args.modes)
            assert result['zero_orders_t'] == result['infinity_orders_t']
            fit = pressure(result)
            report['records'].append({'height_cap': cap, 'scaling': scaling, 'coefficients': C,
                                      'subset_gram_certificate': certificate, 'code': code,
                                      'orthogonal_retraction': left, 'profile': result,
                                      'pressure': fit, 'independent_cut_profiles_agree': True})
            report['in_progress'] = None
            save()
            print(json.dumps({'completed': {'cap': cap, 'scaling': scaling, 'q': code['q'],
                                            'rho': fit['rho'], 'seconds': result['seconds']}}), flush=True)
    report['completed'] = True
    save()


if __name__ == '__main__':
    main()
