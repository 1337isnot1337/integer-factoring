#!/usr/bin/env python3
"""Exact finite controls for a unary factoring predicate, not a fast count."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_carry_fourier_elimination import irreducible
from one_trace_inverse_words import (
    inverse_factor_prefix_value, inverse_query_parameters,
    odd_inverse_polynomial, smallest_factor_from_inverse_counts,
)
from teichmuller_cyclic_carries import UnramifiedRing
from witt_trace_word_circuits import trace_value

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/one-trace-inverse-word-controls.json'


def prefixes(bits):
    return [''.join(p) for k in range(bits+1) for p in product('01', repeat=k)]


def literal(target, prefix, x):
    a = (target.bit_length()+1)//2
    return int(x >= 2 and target % x == 0 and target//x >= 2
               and format(x, f'0{a}b').startswith(prefix))


def run():
    inverse_controls = 0
    for h in range(1, 10):
        for x in range(1, 1 << h, 2):
            assert odd_inverse_polynomial(x, h) == pow(x, -1, 1 << h)
            inverse_controls += 1

    digest, pointwise, queries = sha256(), 0, 0
    scalar_records, recovery_controls = [], 0
    for target in range(3, 256, 2):
        plan = inverse_query_parameters(target)
        table = {}
        for prefix in prefixes(plan['candidate_bits']):
            total = 0
            for x in range(1 << plan['candidate_bits']):
                # Nonzero upper trace bits must not alter the canonical word.
                supplied = x+((3+target % 7) << plan['candidate_bits'])
                answer = inverse_factor_prefix_value(target, prefix, supplied)
                assert answer == literal(target, prefix, x)
                total += answer
                pointwise += 1
                digest.update(json.dumps([target, prefix, x, answer]).encode())
            table[prefix] = total
            queries += 1
        result, used = smallest_factor_from_inverse_counts(target, lambda n, p: table[p])
        expected = next(x for x in range(2, target+1) if target % x == 0)
        assert result == expected and used <= plan['candidate_bits']+1
        recovery_controls += 1
        scalar_records.append(dict(target=target, candidate_bits=plan['candidate_bits'],
                                   accepted_candidates=table[''],
                                   prefix_table_sha256=sha256(json.dumps(table).encode()).hexdigest(),
                                   recovered_factor=result, recovery_queries=used))

    # This is an independently measured complete trace domain, not a uniform
    # fibre assumption. The weighted oracle is deliberately finite/exhaustive.
    degree, polynomial, precision = 9, 0x211, 18
    assert irreducible(polynomial)
    ring = UnramifiedRing(precision, [(polynomial >> j) & 1 for j in range(degree+1)])
    trace_values = [trace_value(ring, ring.teichmuller(x)) for x in range(1 << degree)]
    multiplicities = Counter(x % 4 for x in trace_values)
    assert set(multiplicities) == set(range(4))
    weighted_records = []
    for target in (5, 7, 9, 11, 13, 15):
        plan = inverse_query_parameters(target, degree)
        assert plan['starting_precision'] <= precision and plan['candidate_bits'] == 2
        table = {}
        for prefix in prefixes(2):
            actual = sum(inverse_factor_prefix_value(target, prefix, t, degree) for t in trace_values)
            expected = sum(multiplicities[x]*literal(target, prefix, x) for x in range(4))
            assert actual == expected and 0 <= actual < (1 << plan['count_precision'])
            table[prefix] = actual
            queries += 1
        result, used = smallest_factor_from_inverse_counts(target, lambda n, p: table[p])
        assert result == next(x for x in range(2, target+1) if target % x == 0)
        weighted_records.append(dict(target=target, parameters=plan,
                                     complete_weighted_counts=table,
                                     recovered_factor=result, recovery_queries=used))

    # Large scalar points test the same actual polynomial/halving implementation;
    # no larger field is enumerated and no complete finite-factorization claim
    # is made for these supplied candidate divisors.
    large_records = []
    for p, q in ((1000003, 1000033), (1000000007, 1000000009),
                 (18446744073709551557, 18446744073709551533)):
        target = p*q
        plan = inverse_query_parameters(target)
        for x in (1, 2, min(p, q), min(p, q)+2):
            assert x < (1 << plan['candidate_bits'])
            for prefix in ('', format(x, f"0{plan['candidate_bits']}b")[:3]):
                value = inverse_factor_prefix_value(target, prefix, x+(5 << plan['candidate_bits']))
                assert value == literal(target, prefix, x)
                pointwise += 1
        large_records.append(dict(target=target, supplied_candidate=min(p, q), parameters=plan))

    # With only b inverse bits, x=3 falsely passes the range test for prime N=7.
    wrapped = dict(target=7, candidate=3, insufficient_modulus=8,
                   modular_cofactor=7*pow(3, -1, 8) % 8)
    assert wrapped['modular_cofactor'] == 5 and 3*5 != 7
    assert inverse_factor_prefix_value(7, '', 3) == 0
    even_records = []
    for target in (2, 4, 6, 100, 1 << 128):
        def forbidden(n, prefix):
            raise AssertionError('Even input must not call the count oracle')
        assert smallest_factor_from_inverse_counts(target, forbidden) == (2, 0)
        even_records.append(target)
        recovery_controls += 1

    rejected = []
    for target, degree in ((2, None), (9, 1), (9, 4)):
        try:
            inverse_query_parameters(target, degree)
        except ValueError:
            rejected.append([target, degree])
        else:
            raise AssertionError('Invalid query accepted')

    names = ['tools/one_trace_inverse_words.py', 'tools/check_one_trace_inverse_words.py',
             'tools/witt_trace_word_circuits.py', 'tools/teichmuller_cyclic_carries.py',
             'tools/check_carry_fourier_elimination.py']
    return dict(scope='Unary half-width pointwise inverse predicates; exhaustive small scalar and one complete measured trace model. The finite weighted count oracle is exhaustive; no fast contractor or stronger complete factoring bound.',
                scalar_controls=scalar_records,
                trace_model=dict(degree=9, polynomial=polynomial,
                                 precision=precision, complete_points=len(trace_values),
                                 word_multiplicities=[multiplicities[x] for x in range(4)],
                                 trace_sha256=sha256(json.dumps(trace_values).encode()).hexdigest()),
                weighted_controls=weighted_records, large_scalar_controls=large_records,
                insufficient_modulus_negative=wrapped, even_inputs=even_records,
                rejected_parameters=rejected,
                counts=dict(inverse_polynomial_controls=inverse_controls,
                            scalar_pointwise_queries=pointwise, complete_trace_points=len(trace_values),
                            prefix_count_queries=queries, scalar_minimum_recoveries=recovery_controls),
                pointwise_sha256=digest.hexdigest(),
                source_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in names},
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    args = parser.parse_args()
    receipt = run()
    if args.generate:
        RECEIPT.write_text(json.dumps(receipt, indent=2)+'\n')
    else:
        assert receipt == json.loads(RECEIPT.read_text())
    print(json.dumps(dict(scope=receipt['scope'], counts=receipt['counts'])))


if __name__ == '__main__':
    main()
