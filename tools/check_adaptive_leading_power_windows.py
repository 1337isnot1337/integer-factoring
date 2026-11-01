"""Check resumable power windows and the fixed-seed 144-bit factorization.

Default: complete small searches, exact large-prefix witnesses, and a fresh
154884-query continuation. --full also repeats the 8077846-query first phase.
--regenerate reconstructs both phases and the input-prime certificates.
"""

import argparse
from copy import deepcopy
from hashlib import sha256
import json
from math import prod
from pathlib import Path

from adaptive_leading_power_windows import (factor_with_expanding_power_windows,
                                            resume_leading_power_windows)
from arithmetic_two_group import power
from certified_small_factor import CertifiedFactorizer, verify_prime_certificates
from check_leading_power_smooth_relations import canonical, verify_witnesses
from leading_power_smooth_relations import factor_with_smooth_leading_powers
from quadratic_ideal_hnf import compose_forms_hnf
from weighted_norm_relations import least_power_above

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / 'evidence/adaptive-leading-power-windows.json'
SEED = 'integer-factoring-144-bit-leading-power-2026-11-01'
SOURCES = ('tools/adaptive_leading_power_windows.py',
           'tools/check_adaptive_leading_power_windows.py',
           'tools/leading_power_smooth_relations.py',
           'tools/certified_small_factor.py',
           'tools/dynamic_leading_norm_relations.py',
           'tools/check_leading_power_smooth_relations.py')


def starts():
    return [(3 << 70) | (int.from_bytes(sha256((SEED + ':' + side).encode()).digest(),
                                      'big') & ((1 << 70) - 1)) | 1
            for side in ('left', 'right')]


def small_controls():
    records = []
    independently_enumerated_candidates = resumed_phases = 0
    for integer, cutoff in ((77, 3), (97, 7), (103, 7), (211, 7), (307, 13),
                            (1009, 13), (10007, 13), (10807, 13),
                            (1022117, 7), (1000036000099, 31)):
        result = canonical(factor_with_expanding_power_windows(integer, cutoff, (1, 4, 16)))
        verify_witnesses(result)
        resumed_phases += len(result.get('resume_phases', []))
        delta = result['discriminant']
        identity = (1, delta % 2, (delta % 2 - delta) // 4)
        prime_records = {int(p): (r[0], tuple(r[1]))
                         for p, r in result['prime_records'].items()}
        # Independent HNF steps and literal divisions find the first accepted
        # exponent for every retained target, across the phase boundaries.
        for witness in result['witnesses']:
            prime = witness['prime']
            first = least_power_above(prime, result['height_bound'])
            form = power(prime_records[prime][1], first, identity, compose_forms_hnf)
            for exponent in range(first, witness['exponent'] + 1):
                remaining = form[0]
                for p in result['active_primes']:
                    while remaining % p == 0:
                        remaining //= p
                independently_enumerated_candidates += 1
                if remaining == 1:
                    assert exponent == witness['exponent']
                    break
                form = compose_forms_hnf(form, prime_records[prime][1])
            else:
                raise AssertionError('The independent prefix did not find its accepted norm.')
        records.append({k: result[k] for k in ('integer', 'factor', 'base_size',
                                               'relation_rank', 'norm_queries')})
    assert resumed_phases > 0
    control = canonical(factor_with_smooth_leading_powers(4951760154835678088235319297, 31, 4))
    assert control['factor'] is None
    rejections = 0
    for field in ('norm_queries', 'trial_divisions', 'maximum_exponent_tested'):
        bad = deepcopy(control)
        bad[field] += 1
        try:
            resume_leading_power_windows(bad, 8)
        except ValueError:
            rejections += 1
        else:
            raise AssertionError('Corrupted retained work counter was accepted.')
    for windows in ((), (0,), (4, 4), (4, 2)):
        try:
            factor_with_expanding_power_windows(103, 7, windows)
        except ValueError:
            rejections += 1
        else:
            raise AssertionError('Malformed window schedule was accepted.')
    return {'records': records, 'resumed_phases': resumed_phases,
            'independent_hnf_prefix_candidates': independently_enumerated_candidates,
            'invalid_input_rejections': rejections}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true')
    parser.add_argument('--regenerate', action='store_true')
    args = parser.parse_args()
    hashes = {p: sha256((ROOT / p).read_bytes()).hexdigest() for p in SOURCES}
    if args.regenerate:
        certifier = CertifiedFactorizer()
        primes = []
        for start in starts():
            candidate = start
            while not (certifier._probable(candidate) and certifier._certificate(candidate)):
                candidate += 2
            primes.append(candidate)
        data = canonical({'seed': SEED, 'starts': starts(), 'input_primes': primes,
                          'prime_certificates': certifier.certificates,
                          'source_sha256': hashes})
    else:
        data = json.loads(EVIDENCE.read_text())
        assert data['source_sha256'] == hashes
    assert data['seed'] == SEED and data['starts'] == starts()
    assert data['input_primes'] == [4108225714090758818561, 4224614633883335075341]
    verified_primes = verify_prime_certificates(data['prime_certificates'])
    assert set(data['input_primes']).issubset(verified_primes)
    for start, prime in zip(starts(), data['input_primes']):
        for candidate in range(start, prime, 2):
            assert not CertifiedFactorizer._probable(candidate)
    integer = prod(data['input_primes'])
    assert integer == 17355670471043633865640786569418441084204301
    assert integer.bit_length() == 144
    if args.full or args.regenerate:
        initial = canonical(factor_with_smooth_leading_powers(integer, 8191, 50000))
        if args.regenerate:
            data['first_phase'] = initial
        else:
            assert initial == data['first_phase']
    initial = data['first_phase']
    verify_witnesses(initial)
    assert initial['factor'] is None and initial['relation_rank'] == 494
    result = canonical(resume_leading_power_windows(initial, 200000))
    verify_witnesses(result)
    if args.regenerate:
        data['search'] = result
    else:
        assert result == data['search']
    assert (result['factor'], result['base_size'], result['relation_rank'],
            result['norm_queries'], result['trial_divisions'],
            result['maximum_leading_bits']) == (
                4108225714090758818561, 513, 499, 8232730, 182539, 73)
    phase = result['resume_phases'][0]
    assert phase['additional_norm_queries'] == 154884
    assert phase['additional_relations'] == 5
    assert result['norm_queries'] == initial['norm_queries'] + 154884
    controls = small_controls()
    if args.regenerate:
        data['small_controls'] = controls
        EVIDENCE.write_text(json.dumps(data, separators=(',', ':')) + '\n')
    elif 'small_controls' in data:
        assert data['small_controls'] == controls
    print(json.dumps({'status': 'pass', 'input_bits': 144,
                      'factor': result['factor'], 'total_norm_queries': result['norm_queries'],
                      'resumed_norm_queries': 154884, 'accepted_relations': 499,
                      'base_size': 513, 'input_prime_certificates': len(verified_primes),
                      'first_phase_search_repeated': bool(args.full or args.regenerate),
                      'small_controls': controls}, indent=2))


if __name__ == '__main__':
    main()
