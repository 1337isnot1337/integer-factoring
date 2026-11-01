#!/usr/bin/env python3
"""Generate certified fresh inputs and replay prime-cofactor controls.

Input construction is separate from the factoring algorithm. Recursive
Lucas certificates prove both input primes; the search receives only N.
"""

import argparse
from hashlib import sha256
import json
from math import prod
from pathlib import Path
from time import perf_counter

from certified_small_factor import CertifiedFactorizer, verify_prime_certificates
from check_prime_cofactor_saturation import normalize_record, summary, verify_result
from prime_cofactor_saturation import factor_with_prime_cofactor_saturation


ROOT = Path(__file__).resolve().parents[1]
SEED = 'integer-factoring-fresh-prime-cofactor-20261101'
SOURCES = ('tools/check_fresh_prime_cofactor_control.py',
           'tools/prime_cofactor_saturation.py', 'tools/check_prime_cofactor_saturation.py',
           'tools/certified_small_factor.py')


def starts(bits):
    assert bits >= 32 and bits % 2 == 0
    half = bits//2
    return [(3 << (half-2)) | (int.from_bytes(sha256(
        f'{SEED}:{bits}:{side}'.encode()).digest(), 'big') & ((1 << (half-2))-1)) | 1
        for side in ('left', 'right')]


def generate_input(bits):
    certifier, primes = CertifiedFactorizer(), []
    for first in starts(bits):
        candidate = first
        while not (certifier._probable(candidate) and certifier._certificate(candidate)):
            candidate += 2
        primes.append(candidate)
    return {'seed': SEED, 'bits': bits, 'starts': starts(bits), 'input_primes': primes,
            'integer': prod(primes), 'prime_certificates': certifier.certificates,
            'input_generation_statistics': certifier.stats,
            'input_generation_is_not_part_of_the_factoring_algorithm': True}


def verify_input(data):
    bits = data['bits']
    assert data['seed'] == SEED and data['starts'] == starts(bits)
    primes = data['input_primes']
    assert len(primes) == 2 and primes[0] != primes[1]
    assert prod(primes) == data['integer'] and data['integer'].bit_length() == bits
    verified = verify_prime_certificates(data['prime_certificates'])
    assert set(primes) <= verified
    rejected = 0
    for first, prime in zip(starts(bits), primes):
        assert prime >= first and (prime-first) % 2 == 0
        for candidate in range(first, prime, 2):
            assert not CertifiedFactorizer._probable(candidate)
            rejected += 1
    return {'prime_certificates': len(verified), 'rejected_candidates': rejected}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bits', type=int, default=160)
    parser.add_argument('--input', type=Path)
    parser.add_argument('--generate-input', action='store_true')
    parser.add_argument('--mode', choices=('fixed', 'selective', 'all'))
    parser.add_argument('--cutoff', type=int, default=16381)
    parser.add_argument('--window', type=int, default=400000)
    parser.add_argument('--max-queries', type=int, default=20000000)
    parser.add_argument('--max-generators', type=int, default=100000)
    parser.add_argument('--promotion-multiplier', type=int, default=16)
    parser.add_argument('--replay', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    input_path = args.input or ROOT/'evidence'/f'prime-cofactor-fresh-{args.bits}-input.json'
    if args.generate_input:
        data = generate_input(args.bits)
        verification = verify_input(data)
        input_path.write_text(json.dumps(data, sort_keys=True, separators=(',', ':'))+'\n')
        print(json.dumps({'input': str(input_path), 'integer': data['integer'],
                          'primes': data['input_primes'], 'verification': verification}), flush=True)
        if not args.mode:
            return
    data = json.loads(input_path.read_text())
    input_verification = verify_input(data)
    hashes = {name: sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES}
    if args.replay:
        report = json.loads(args.replay.read_text())
        assert report['input_file'] == input_path.name
        assert report['input_sha256'] == sha256(input_path.read_bytes()).hexdigest()
        result = normalize_record(report['result'])
        assert result['integer'] == data['integer']
        verification = verify_result(result)
        if result['factor'] is not None:
            assert result['factor'] in data['input_primes']
        print(json.dumps({'mode': report['mode'], 'summary': summary(result),
                          'input_verification': input_verification,
                          'verification': verification}, indent=2))
        return
    if not args.mode or not args.output:
        parser.error('Supply --mode and --output to run a search, or --replay.')
    started = perf_counter()
    progress_path = Path('/tmp')/('integer-factoring-'+args.output.stem+'-progress.json')
    def progress(state):
        if state['norm_queries'] % 100000 == 0:
            row = {key: state[key] for key in ('norm_queries', 'generator_count',
                                              'relation_rank', 'current_generator')}
            row['elapsed_seconds'] = perf_counter()-started
            progress_path.write_text(json.dumps(row)+'\n')
            print(json.dumps(row), flush=True)
    cutoff = args.cutoff
    promotion = cutoff*cutoff if args.mode == 'all' else args.promotion_multiplier*cutoff
    result = factor_with_prime_cofactor_saturation(
        data['integer'], cutoff, args.window, max_queries=args.max_queries,
        max_generators=args.max_generators, promote_residuals=args.mode != 'fixed',
        promotion_cutoff=promotion, verify_candidates=True, include_relations=False,
        progress=progress)
    elapsed = perf_counter()-started
    verification = verify_result(result)
    assert result['factoring_oracle_calls'] == result['residual_primality_tests'] == 0
    if result['factor'] is not None:
        assert result['factor'] in data['input_primes']
    report = {'claim': 'Fresh certified-input prime-cofactor control', 'mode': args.mode,
              'input_file': input_path.name, 'input_sha256': sha256(input_path.read_bytes()).hexdigest(),
              'input_verification': input_verification, 'elapsed_seconds': elapsed,
              'result': result, 'verification': verification, 'source_sha256': hashes,
              'scope': 'Finite deterministic traversal; no uniform complexity improvement.'}
    args.output.write_text(json.dumps(report, sort_keys=True, separators=(',', ':'))+'\n')
    print(json.dumps({'mode': args.mode, 'elapsed_seconds': elapsed,
                      'summary': summary(result), 'verification': verification}), flush=True)


if __name__ == '__main__':
    main()
