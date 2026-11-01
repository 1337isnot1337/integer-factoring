"""Generate and replay an unrelated fixed-seed 124-bit factoring control."""

import argparse
from hashlib import sha256
import json
from math import prod
from pathlib import Path

from certified_small_factor import CertifiedFactorizer, verify_prime_certificates
from check_leading_power_smooth_relations import canonical, verify_witnesses
from leading_power_smooth_relations import factor_with_smooth_leading_powers


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / 'evidence/fresh-leading-power-control.json'
SEED = 'integer-factoring-leading-power-fresh-2026-11-01'
SOURCES = ('tools/check_fresh_leading_power_control.py',
           'tools/leading_power_smooth_relations.py',
           'tools/certified_small_factor.py')


def starts():
    return [(1 << 61) | (int.from_bytes(sha256((SEED + ':' + side).encode()).digest(),
                                       'big') & ((1 << 61) - 1)) | 1
            for side in ('left', 'right')]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--regenerate', action='store_true')
    args = parser.parse_args()
    hashes = {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in SOURCES}
    if args.regenerate:
        certifier = CertifiedFactorizer()
        primes = []
        for start in starts():
            candidate = start
            while not (certifier._probable(candidate) and certifier._certificate(candidate)):
                candidate += 2
            primes.append(candidate)
        evidence = canonical({'seed': SEED, 'starts': starts(), 'input_primes': primes,
                              'prime_certificates': certifier.certificates,
                              'source_sha256': hashes})
    else:
        evidence = json.loads(EVIDENCE.read_text())
        assert evidence['source_sha256'] == hashes
    assert evidence['seed'] == SEED and evidence['starts'] == starts()
    primes = evidence['input_primes']
    assert primes == [4509576403369997671, 3701143448763469987]
    verified = verify_prime_certificates(evidence['prime_certificates'])
    assert set(primes).issubset(verified)
    rejected = 0
    for start, prime in zip(starts(), primes):
        for candidate in range(start, prime, 2):
            assert not CertifiedFactorizer._probable(candidate)
            rejected += 1
    assert rejected == 39
    integer = prod(primes)
    assert integer == 16690589162031198237853308486068400277
    result = canonical(factor_with_smooth_leading_powers(integer, 4093, 100000))
    if args.regenerate:
        evidence['search'] = result
        EVIDENCE.write_text(json.dumps(evidence, separators=(',', ':')) + '\n')
    else:
        assert evidence['search'] == result
    verify_witnesses(result)
    assert (result['factor'], result['base_size'], result['relation_rank'],
            result['norm_queries'], result['trial_divisions'],
            result['maximum_exponent_tested'], result['maximum_leading_bits']) == (
                3701143448763469987, 262, 248, 2301916, 48870, 51101, 63)
    print(json.dumps({'status': 'passed', 'verified_input_prime_certificates': len(verified),
                      'rejected_input_candidates': rejected,
                      'search': {key: result[key] for key in
                                 ('integer', 'factor', 'base_size', 'relation_rank',
                                  'norm_queries', 'modular_smoothness_tests',
                                  'trial_divisions', 'maximum_leading_bits')}}, indent=2))


if __name__ == '__main__':
    main()
