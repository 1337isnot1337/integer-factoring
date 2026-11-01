#!/usr/bin/env python3
"""Batched search and full replay on an already certified fresh input."""

import argparse
import gc
from hashlib import sha256
import json
import lzma
from pathlib import Path
import resource
from time import perf_counter, process_time

from check_fresh_prime_cofactor_control import verify_input
from check_prime_cofactor_saturation import summary, verify_result
from prime_cofactor_saturation import factor_with_prime_cofactor_saturation

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('batch_smooth_fresh_run.py', 'batch_smooth_parts.py',
           'prime_cofactor_saturation.py', 'leading_power_smooth_relations.py',
           'check_fresh_prime_cofactor_control.py', 'check_prime_cofactor_saturation.py',
           'certified_small_factor.py', 'cofactor_collision_relations.py',
           'relation_two_saturation.py')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    defaults = dict(bits=192, cutoff=65521, window=1600000,
                    max_queries=120000000, max_generators=250000,
                    promotion_multiplier=16, batch_size=32, memory_limit_mib=3584)
    for name, value in defaults.items():
        parser.add_argument('--' + name.replace('_', '-'), type=int, default=value)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.suffix != '.xz' or args.memory_limit_mib < 512:
        parser.error('Require an .xz output and an address-space limit of at least 512 MiB.')
    limit = args.memory_limit_mib * 2**20
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    input_path = ROOT / 'evidence' / f'prime-cofactor-fresh-{args.bits}-input.json'
    input_bytes = input_path.read_bytes()
    data = json.loads(input_bytes)
    assert data['bits'] == args.bits
    input_verification = verify_input(data)
    hashes = {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
              for name in SOURCES}
    started, cpu_started, next_report = perf_counter(), process_time(), 100000

    def progress(state):
        nonlocal next_report
        if state['norm_queries'] < next_report:
            return
        row = {key: state[key] for key in ('norm_queries', 'generator_count',
               'relation_rank', 'current_generator', 'maximum_saturation_product_bits')}
        row.update(scanned_norm_queries=state.get('scanned_norm_queries', state['norm_queries']),
                   elapsed_seconds=perf_counter() - started,
                   cpu_seconds=process_time() - cpu_started,
                   maximum_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
        print(json.dumps(row), flush=True)
        next_report = (state['norm_queries'] // 100000 + 1) * 100000

    print(json.dumps({'event': 'started', 'policy': {name: getattr(args, name) for name in defaults},
                      'input_verification': input_verification, 'source_sha256': hashes}), flush=True)
    # Only the input integer, not its certified factors, is passed to the search.
    result = factor_with_prime_cofactor_saturation(
        data['integer'], args.cutoff, args.window, max_queries=args.max_queries,
        max_generators=args.max_generators, promotion_cutoff=args.promotion_multiplier * args.cutoff,
        verify_candidates=True, include_relations=False,
        smooth_batch_size=args.batch_size, progress=progress)
    elapsed, cpu_elapsed = perf_counter() - started, process_time() - cpu_started
    gc.collect()
    verification = verify_result(result)
    assert result['factoring_oracle_calls'] == result['residual_primality_tests'] == 0
    assert result['factor'] is None or result['factor'] in data['input_primes']
    report = {'claim': 'Batched search on a precommitted certified input', 'mode': 'selective',
              'input_file': input_path.name, 'input_sha256': sha256(input_bytes).hexdigest(),
              'input_verification': input_verification, 'elapsed_seconds': elapsed,
              'cpu_seconds': cpu_elapsed, 'result': result, 'verification': verification,
              'source_sha256': hashes,
              'scope': 'Finite deterministic traversal; no uniform complexity improvement.'}
    # Stream compression avoids allocating another entire uncompressed certificate.
    digest = sha256()
    temporary = args.output.with_name(args.output.name + '.partial')
    with lzma.open(temporary, 'wb', preset=6) as stream:
        for piece in json.JSONEncoder(sort_keys=True, separators=(',', ':')).iterencode(report):
            block = piece.encode()
            digest.update(block)
            stream.write(block)
        digest.update(b'\n')
        stream.write(b'\n')
    temporary.replace(args.output)
    print(json.dumps({'event': 'completed', 'elapsed_seconds': elapsed,
                      'cpu_seconds': cpu_elapsed, 'summary': summary(result),
                      'scanned_norm_queries': result.get('scanned_norm_queries', result['norm_queries']),
                      'generated_ahead_norms': result.get('generated_ahead_norms', 0),
                      'verification': verification,
                      'certificate_uncompressed_sha256': digest.hexdigest()}), flush=True)


if __name__ == '__main__':
    main()
