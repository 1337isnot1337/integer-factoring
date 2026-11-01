#!/usr/bin/env python3
"""Complete profiles for the maximal five-column orthogonal F4 pair code.

Floating pressure fits are exploratory. Both cuts are computed separately;
the orthogonal identity is then checked on their complete order lists.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from coherent_differential_kernels_f4 import differential_kernel
from coherent_hot_compounds_f4 import profiles
from orthogonal_differential_pairs_f4 import orthogonal_retraction, subset_gram_certificate
from search_coherent_hot_retractions import pressure


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--caps', type=int, nargs='+', default=[20, 24])
    ap.add_argument('--coordinate', action='store_true')
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    C = [[1, 1, 2, 2, 0], [0, 1, 3, 1, 3]]
    pairs = [(i, 9 - i) for i in range(5)]
    operators = [{(1 << a) | (1 << b): row[i] for i, (a, b) in enumerate(pairs) if row[i]}
                 for row in C]
    certificate = subset_gram_certificate(C)
    assert certificate['all_subsets_nondegenerate']
    names = [Path(__file__).name, 'coherent_differential_kernels_f4.py',
             'coherent_hot_compounds_f4.py', 'coherent_hot_compounds_f4.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp',
             'orthogonal_differential_pairs_f4.py', 'search_coherent_hot_retractions.py']
    report = {'scope': __doc__, 'coefficients': C, 'pairs': pairs,
              'subset_gram_certificate': certificate, 'modes': 14,
              'retraction_kind': 'coordinate' if args.coordinate else 'orthogonal',
              'caps': args.caps, 'source_sha256': {
                  n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'records': [], 'completed': False, 'in_progress': None}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()

    def save():
        report['elapsed_seconds'] = time.monotonic() - start
        temp = args.output.with_suffix(args.output.suffix + '.tmp')
        temp.write_text(json.dumps(report, indent=2) + '\n')
        temp.replace(args.output)

    for cap in args.caps:
        code = differential_kernel(operators, 14, cap)
        left = code['retraction'] if args.coordinate else orthogonal_retraction(code)
        assert left is not None
        report['in_progress'] = {'height_cap': cap, 'q': code['q']}
        save()
        print(json.dumps({'starting': report['in_progress'], 'retraction': report['retraction_kind']}), flush=True)
        result = profiles(code['embedding'], left, 14)
        if not args.coordinate:
            assert result['zero_orders_t'] == result['infinity_orders_t']
        fit = pressure(result)
        report['records'].append({'height_cap': cap, 'code': code, 'retraction': left,
                                  'profile': result, 'pressure': fit})
        report['in_progress'] = None
        save()
        print(json.dumps({'completed': cap, 'q': code['q'], 'pressure': fit,
                          'seconds': result['seconds']}), flush=True)
    report['completed'] = True
    save()


if __name__ == '__main__':
    main()
