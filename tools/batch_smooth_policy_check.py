#!/usr/bin/env python3
"""Matched batch-filter timings for unpromoted and selective cofactor policies."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
from time import perf_counter, process_time

from batch_smooth_check import FIXTURES, ROOT, SOURCES, compare, payload_digest
from check_prime_cofactor_saturation import verify_result
from prime_cofactor_saturation import factor_with_prime_cofactor_saturation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bits', type=int, choices=tuple(FIXTURES), default=144)
    parser.add_argument('--batch-size', type=int, default=32)
    parser.add_argument('--promotion-multiplier', type=int, default=16)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.batch_size < 2 or args.promotion_multiplier < 1:
        parser.error('Require batch size at least two and a positive promotion multiplier.')
    integer, cutoff, window = FIXTURES[args.bits]
    sources = (*SOURCES, 'tools/batch_smooth_policy_check.py',
               'tools/check_prime_cofactor_saturation.py')
    report = {'claim': 'Matched cofactor-policy batch comparisons',
              'scope': 'Finite completed traversals; no uniform complexity or timing claim.',
              'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest()
                                for name in sources},
              'environment': {'python': platform.python_version(),
                              'machine': platform.machine(), 'system': platform.system()},
              'policies': []}
    for policy, promote in (('unpromoted', False), ('selective', True)):
        options = {'promote_residuals': promote,
                   'promotion_cutoff': args.promotion_multiplier * cutoff}
        reference, runs = None, []
        for batch in (1, args.batch_size):
            start, cpu_start = perf_counter(), process_time()
            result = factor_with_prime_cofactor_saturation(
                integer, cutoff, window, **options, smooth_batch_size=batch,
                verify_candidates=True, include_relations=False)
            elapsed, cpu_elapsed = perf_counter() - start, process_time() - cpu_start
            if reference is None:
                reference = result
                verification = verify_result(reference)
            else:
                compare(reference, result, 'saturation', batch)
            row = {'input_bits': args.bits, 'policy': policy, 'batch_size': batch,
                   'seconds': elapsed, 'cpu_seconds': cpu_elapsed,
                   'factor': result['factor'], 'norm_queries': result['norm_queries'],
                   'scanned_norm_queries': result.get('scanned_norm_queries', result['norm_queries']),
                   'generated_ahead_norms': result.get('generated_ahead_norms', 0),
                   'relation_rank': result['relation_rank'],
                   'prime_product_bits': result['maximum_saturation_product_bits'],
                   'full_mathematical_payload_sha256': payload_digest(result, 'saturation'),
                   'matches_reference_witnesses_and_relations': True}
            runs.append(row)
            print(json.dumps(row), flush=True)
        report['policies'].append({'policy': policy, 'options': options,
                                   'runs': runs, 'reference_verification': verification})
    args.output.write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
