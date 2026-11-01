#!/usr/bin/env python3
"""Finite trace images and pointwise native-word circuit controls, not a fast count."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_carry_fourier_elimination import irreducible
from teichmuller_cyclic_carries import UnramifiedRing
from witt_trace_word_circuits import word_query_parameters, trace_value, extract_bits, factor_prefix_value

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/witt-trace-word-controls.json'


def prefixes(bits):
    return [''.join(x) for size in range(bits+1) for x in product('01', repeat=size)]


def native_predicate(target, prefix, x, y):
    b = target.bit_length()
    return int(x >= 2 and y >= 2 and x*y == target and format(x, f'0{b}b').startswith(prefix))


def recover_from_counts(target, table):
    queries = 1
    if not table['']:
        return target, queries
    prefix = ''
    for _ in range(target.bit_length()):
        queries += 1
        prefix += '0' if table[prefix+'0'] else '1'
        assert table[prefix]
    return int(prefix, 2), queries


def run():
    models = [(5, 0x25, 15), (9, 0x211, 29)]
    traces, model_records = {}, []
    for degree, polynomial, precision in models:
        assert irreducible(polynomial)
        coefficients = [(polynomial >> j) & 1 for j in range(degree+1)]
        ring = UnramifiedRing(precision, coefficients)
        values, digest = [], sha256()
        for residue in range(1 << degree):
            point = ring.teichmuller(residue)
            value = trace_value(ring, point)
            assert ring.frobenius(value) == value
            values.append(value)
            digest.update(json.dumps([residue, point, value]).encode())
        traces[degree] = values
        model_records.append(dict(degree=degree, polynomial=polynomial, precision=precision,
                                  complete_field_points=len(values), point_trace_sha256=digest.hexdigest()))
    images, decoded = [], 0
    banks = {}
    for bits, degree in ((1, 5), (2, 9), (3, 9)):
        plan = word_query_parameters(bits, degree)
        multiplicities, representatives = Counter(), {}
        for value in traces[degree]:
            digits, remaining = extract_bits(value, plan['starting_precision'], bits)
            word = sum(bit << j for j, bit in enumerate(digits))
            assert word == value % (1 << bits)
            assert remaining == plan['count_precision']+plan['zero_bits']
            multiplicities[word] += 1
            representatives.setdefault(word, value)
            decoded += 1
        assert set(multiplicities) == set(range(1 << bits))
        banks[bits] = (degree, multiplicities, representatives)
        images.append(dict(parameters=plan, complete_trace_inputs=len(traces[degree]),
                           word_multiplicities=[multiplicities[x] for x in range(1 << bits)]))
    queries, pointwise, recovery_queries, digest = [], 0, 0, sha256()
    for target in range(2, 8):
        b = target.bit_length()
        degree, multiplicities, representatives = banks[b]
        plan = word_query_parameters(b, degree)
        table = {}
        for prefix in prefixes(b):
            actual, expected = 0, 0
            for x, y in product(range(1 << b), repeat=2):
                value = factor_prefix_value(target, prefix, representatives[x], representatives[y], degree)
                independent = native_predicate(target, prefix, x, y)
                assert value == independent
                weight = multiplicities[x]*multiplicities[y]
                actual += weight*value
                expected += weight*independent
                pointwise += 1
                digest.update(json.dumps([target, prefix, x, y, value]).encode())
            assert actual == expected and actual < (1 << plan['count_precision'])
            assert bool(actual) == any(native_predicate(target, prefix, x, y)
                                       for x, y in product(range(1 << b), repeat=2))
            table[prefix] = actual
        recovered, used = recover_from_counts(target, table)
        assert recovered == next(x for x in range(2, target+1) if target % x == 0)
        recovery_queries += used
        queries.append(dict(target=target, degree=degree, count_precision=plan['count_precision'],
                            weighted_counts=table, recovered_factor=recovered, recovery_queries=used))
    # These are scalar word-pair controls, without enumerating the larger field.
    scalar_controls = []
    for target in (9, 15):
        b = target.bit_length()
        degree = 4*b+1
        table = {}
        for prefix in prefixes(b):
            total = 0
            for x, y in product(range(1 << b), repeat=2):
                value = factor_prefix_value(target, prefix, x+(3 << b), y+(5 << b), degree)
                assert value == native_predicate(target, prefix, x, y)
                total += value
                pointwise += 1
                digest.update(json.dumps([target, prefix, x, y, value]).encode())
            table[prefix] = total
        recovered, used = recover_from_counts(target, table)
        assert recovered == 3
        recovery_queries += used
        scalar_controls.append(dict(target=target, degree=degree, counts=table, recovered_factor=recovered,
                                     scope='Scalar word-pair circuit controls; the field domain is not enumerated.'))
    # An undersized trace domain really can omit native words.
    binary = UnramifiedRing(3, (0, 1))
    missed = sorted(set(range(4))-{trace_value(binary, binary.teichmuller(x)) % 4 for x in range(2)})
    assert missed == [2, 3]
    rejected = []
    for bits, degree in ((2, 1), (2, 4)):
        try:
            word_query_parameters(bits, degree)
        except ValueError:
            rejected.append([bits, degree])
        else:
            raise AssertionError('Uncertified trace model accepted')
    names = ['tools/witt_trace_word_circuits.py', 'tools/check_witt_trace_word_circuits.py',
             'tools/teichmuller_cyclic_carries.py', 'tools/check_carry_fourier_elimination.py']
    return dict(scope='Pointwise known-exact-division word circuits and finite trace images. Weighted small-query controls use complete measured fibres; no fast contraction or larger-field enumeration is claimed.',
                models=model_records, trace_images=images, weighted_native_queries=queries,
                scalar_word_controls=scalar_controls, rejected_parameter_pairs=rejected,
                negative_missing_words=missed,
                counts=dict(complete_field_points=sum(x['complete_field_points'] for x in model_records),
                            decoded_trace_words=decoded, pointwise_word_pair_circuits=pointwise,
                            weighted_prefix_queries=sum(len(x['weighted_counts']) for x in queries),
                            scalar_prefix_queries=sum(len(x['counts']) for x in scalar_controls),
                            recovery_queries=recovery_queries), query_sha256=digest.hexdigest(),
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
    print(json.dumps(dict(scope=receipt['scope'], counts=receipt['counts'], negative_missing_words=receipt['negative_missing_words'])))


if __name__ == '__main__':
    main()
