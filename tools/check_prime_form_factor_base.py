#!/usr/bin/env python3
"""Independent arithmetic controls for prime-form root construction."""
import json
from math import isqrt
from pathlib import Path

from binary_quadratic_forms import inverse_form, reduce_form
from factored_norm_representation import represent_factored_norm
from prime_form_factor_base import build_prime_form_factor_base
from prime_form_factor_base import least_quadratic_nonresidue
from prime_form_factor_base import modular_square_root, nonresidue_scan_bound
from prime_form_factor_base import primes_through
from quadratic_ideal_hnf import represented_norm_from_word


def elementary_prime(value):
    return value >= 2 and all(value % divisor for divisor in range(2, isqrt(value) + 1))


def square_root_checks():
    primes = primes_through(1000)
    assert primes == [value for value in range(2, 1001) if elementary_prime(value)]
    decisions = 0
    stats = {}
    largest_small_nonresidue = 0
    for prime in primes:
        residues = {value * value % prime for value in range(prime)}
        for value in range(prime):
            root = modular_square_root(value, prime, stats)
            assert (root is not None) == (value in residues)
            if root is not None:
                assert 0 <= root < prime and root * root % prime == value
            decisions += 1
        if prime > 2:
            expected = next(value for value in range(2, prime) if value not in residues)
            actual = least_quadratic_nonresidue(prime)
            assert actual == expected <= nonresidue_scan_bound(prime)
            largest_small_nonresidue = max(largest_small_nonresidue, actual)

    large_controls = []
    for prime in [65537, 998244353, 2013265921]:
        assert elementary_prime(prime)
        local_stats = {}
        for value in [0, 1, 2, 3, 5, 17, 257, 123456789 % prime]:
            root = modular_square_root(value, prime, local_stats)
            is_square = value == 0 or pow(value, (prime - 1) // 2, prime) == 1
            assert (root is not None) == is_square
            if root is not None:
                assert root * root % prime == value
        large_controls.append({'prime': prime, 'stats': local_stats})
    return {'small_prime_count': len(primes), 'complete_residue_decisions': decisions,
            'largest_small_nonresidue': largest_small_nonresidue,
            'stats': stats, 'large_two_power_controls': large_controls}


def factor_base_checks():
    discriminants = [-3, -4, -7, -8, -11, -15, -16, -20, -23, -24, -27,
                    -31, -36, -40, -44, -48, -63, -84, -92, -120,
                    -41612, -4145955690243]
    cutoff = 499
    records_tested = 0
    examples = []
    for discriminant in discriminants:
        result = build_prime_form_factor_base(discriminant, cutoff)
        by_prime = {record.prime: record for record in result['records']}
        for prime in primes_through(cutoff):
            expected_roots = [root for root in range(2 * prime)
                              if root % 2 == discriminant % 2
                              and (root * root - discriminant) % (4 * prime) == 0]
            eligible = discriminant % prime != 0 and bool(expected_roots)
            assert (prime in by_prime) == eligible
            assert (prime in result['excluded_primes']) == (discriminant % prime == 0)
            if not eligible:
                continue
            record = by_prime[prime]
            assert record.root in expected_roots
            expected_forms = {reduce_form((prime, root,
                                           (root * root - discriminant) // (4 * prime)))
                              for root in expected_roots}
            assert {record.positive, record.negative} == expected_forms
            assert record.negative == inverse_form(record.positive)
            witness = represented_norm_from_word(discriminant, [(prime, 1, record.root)])
            assert witness.norm == prime
            assert witness.form == record.positive
            records_tested += 1
        examples.append({'discriminant': discriminant, 'stats': result['stats']})
    return {'discriminant_count': len(discriminants), 'prime_cutoff': cutoff,
            'records_compared': records_tested, 'examples': examples}


def integration_checks():
    result = build_prime_form_factor_base(-23, 19)
    by_prime = {record.prime: record for record in result['records']}
    positive = represent_factored_norm((2, -1, 3), [(2, 2, by_prime[2].root)])
    assert positive['representation'] is not None
    assert positive['representation'].norm == 4
    negative = represent_factored_norm((2, -1, 3), [])
    assert negative['representation'] is None

    discriminant = -4145955690243
    base = build_prime_form_factor_base(discriminant, 101)['records']
    terms = [(record.prime, (1 if index % 2 == 0 else -1) * (1 + index % 3), record.root)
             for index, record in enumerate(base[:8])]
    witness = represented_norm_from_word(discriminant, terms)
    factors = [(prime, abs(exponent), root) for prime, exponent, root in terms]
    answer = represent_factored_norm(witness.form, factors)
    assert answer['representation'] is not None
    assert answer['representation'].norm == witness.norm
    return {'small_positive_and_negative_controls': 2,
            'large_input_norm_bits': witness.norm.bit_length(),
            'large_input_distinct_primes': len(terms),
            'large_decoder_stats': answer['stats']}


def main():
    evidence = {'claim_scope': 'Root and factor-base arithmetic; no norm selection algorithm.',
                'square_roots': square_root_checks(),
                'factor_bases': factor_base_checks(),
                'orientation_decoder': integration_checks()}
    destination = Path(__file__).resolve().parents[1] / 'evidence' / 'prime-form-factor-base.json'
    destination.write_text(json.dumps(evidence, indent=2) + '\n')
    print(json.dumps({'status': 'pass',
                      'residue_decisions': evidence['square_roots']['complete_residue_decisions'],
                      'factor_base_records': evidence['factor_bases']['records_compared'],
                      'evidence': str(destination)}, indent=2))


if __name__ == '__main__':
    main()
