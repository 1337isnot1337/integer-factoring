#!/usr/bin/env python3
"""Independent literal controls for bijective field-coordinate candidates."""
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import argparse
import json

from check_carry_fourier_elimination import irreducible
from coordinate_inverse_words import (
    cofactor_bit_circuit, coordinate_bits, coordinate_factor_prefix_value, coordinate_query_parameters,
    smallest_factor_from_coordinate_counts, trace_dual_basis,
)
from one_trace_inverse_words import inverse_query_parameters
from teichmuller_cyclic_carries import UnramifiedRing
from witt_trace_word_circuits import trace_value

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/coordinate-inverse-word-controls.json'
MODELS = {1:0x3, 2:0x7, 3:0xB, 4:0x13, 5:0x25, 6:0x43}


def model(degree, precision):
    polynomial = MODELS[degree]
    assert polynomial == 3 or irreducible(polynomial)
    coefficients = [(polynomial >> i) & 1 for i in range(degree+1)]
    dual = trace_dual_basis(coefficients)
    field = UnramifiedRing(1, coefficients)
    assert [[trace_value(field, field.mul(c, 1 << j)) for j in range(degree)]
            for c in dual] == [[int(i == j) for j in range(degree)]
                              for i in range(degree)]
    ring = UnramifiedRing(precision, coefficients)
    return ring, tuple(ring.teichmuller(c) for c in dual)


def literal(target, prefix, word, degree):
    x = 1+2*word
    return int(x > 1 and target % x == 0 and target//x > 1
               and format(word, f'0{degree}b').startswith(prefix))


def run():
    bijections, checks, recoveries = [], 0, []
    inverse_controls = 0
    for width in range(2,8):
        for target in range(1 << width):
            for x in range(1,1 << min(width,4),2):
                digits = tuple((x >> j) & 1 for j in range(min(width,4)))
                bits = cofactor_bit_circuit(target,digits,width,1+target % 3)
                assert set(bits) <= {0,1}
                assert sum(bit << j for j,bit in enumerate(bits)) == target*pow(x,-1,1 << width) % (1 << width)
                inverse_controls += 1
    for degree in MODELS:
        ring, dual = model(degree, 9)
        histogram = Counter()
        upper_traces = 0
        for word in range(1 << degree):
            point = ring.teichmuller(word)
            digits = coordinate_bits(ring, dual, point)
            decoded = sum(bit << i for i, bit in enumerate(digits))
            assert decoded == word
            histogram[decoded] += 1
            upper_traces += sum(trace_value(ring, ring.mul(c, point)) > 1 for c in dual)
        assert set(histogram.values()) == {1}
        bijections.append(dict(degree=degree, polynomial=MODELS[degree],
                               complete_points=len(histogram), all_fibres_one=True,
                               lifted_traces_with_upper_digits=upper_traces))

    cases = [(n,1) for n in range(5,128,2)]
    cases += [(143,3), (221,3), (289,4), (323,4), (437,4)]
    digest = sha256()
    all_prefix_targets = {9,15,35,49,105,323}
    for target,t in cases:
        assert all(target % p for p in range(2,(1 << t)+1))
        plan = coordinate_query_parameters(target,t)
        d = plan['degree']
        ring, dual = model(d,plan['starting_precision'])
        points = [ring.teichmuller(w) for w in range(1 << d)]
        table = {}

        def query(n,prefix,trial_bits):
            nonlocal checks
            assert n == target and trial_bits == t
            if prefix not in table:
                values = [coordinate_factor_prefix_value(n,prefix,ring,dual,u,t)
                          for u in points]
                expected = [literal(n,prefix,w,d) for w in range(1 << d)]
                assert values == expected
                count = sum(values)
                assert count < (1 << plan['count_precision'])
                table[prefix] = count
                checks += len(points)
                digest.update(json.dumps([n,t,prefix,values]).encode())
            return table[prefix]

        if target in all_prefix_targets:
            for k in range(d+1):
                for p in product('01',repeat=k):
                    query(target,''.join(p),t)
        found,used = smallest_factor_from_coordinate_counts(target,query,t)
        expected = next(x for x in range(2,target+1) if target % x == 0)
        assert found == expected and used <= d+1
        recoveries.append(dict(target=target,trial_bits=t,parameters=plan,
                               exact_count=table[''],recovered_factor=found,
                               queries_used=used))

    # An unjustified roughness promise can wrap a positive count to zero.
    # The complete pretest prevents that count call on this actual input.
    n,t = 105,4
    plan = coordinate_query_parameters(n,t)
    count = sum(literal(n,'',w,plan['degree']) for w in range(1 << plan['degree']))
    assert count == 4 and count % (1 << plan['count_precision']) == 0
    def forbidden(*args):
        raise AssertionError('Pretested input must not call the count oracle')
    assert smallest_factor_from_coordinate_counts(n,forbidden,t) == (3,0)
    for n in (2,3,4,6,1 << 128):
        expected = 3 if n == 3 else 2
        assert smallest_factor_from_coordinate_counts(n,forbidden) == (expected,0)

    rejected = 0
    for n,t in ((3,1),(4,1),(9,0)):
        try:
            coordinate_query_parameters(n,t)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError('Invalid parameters accepted')
    ring,dual = model(2,coordinate_query_parameters(35)['starting_precision'])
    for prefix in ('2','000'):
        try:
            coordinate_factor_prefix_value(35,prefix,ring,dual,0)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError('Invalid query accepted')

    names = ['tools/coordinate_inverse_words.py','tools/check_coordinate_inverse_words.py',
             'tools/one_trace_inverse_words.py','tools/witt_trace_word_circuits.py',
             'tools/teichmuller_cyclic_carries.py','tools/check_carry_fourier_elimination.py']
    return dict(scope='Bijective trace-dual odd words, actual lifted pointwise predicates, complete finite prefix counts and charged pretest; exhaustive test oracle only, no fast contractor or stronger factoring exponent.',
                bijections=bijections,complete_recoveries=recoveries,
                division_free_inverse_controls=inverse_controls,
                pointwise_controls=checks,query_digest=digest.hexdigest(),
                false_roughness_guard=dict(target=105,trial_bits=4,actual_count=count,
                                           unjustified_precision=plan['count_precision'],
                                           wrapped_count=0,pretest_factor=3),
                comparison=dict(target=323,trial_bits=4,
                                old_weighted_trace_plan=inverse_query_parameters(323),
                                bijective_plan=coordinate_query_parameters(323,4)),
                rejected_parameters=rejected,
                source_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in names})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate',action='store_true')
    args = parser.parse_args()
    receipt = run()
    if args.generate:
        RECEIPT.write_text(json.dumps(receipt,indent=2)+'\n')
    else:
        assert receipt == json.loads(RECEIPT.read_text())
    print(json.dumps(dict(scope=receipt['scope'],pointwise_controls=receipt['pointwise_controls'],
                         complete_recoveries=len(receipt['complete_recoveries']),
                         bijection_points=sum(x['complete_points'] for x in receipt['bijections']))))


if __name__ == '__main__':
    main()
