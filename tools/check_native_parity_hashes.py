#!/usr/bin/env python3
"""Actual pure-field native states and independent sparse parity recovery."""
from hashlib import sha256
from pathlib import Path
import argparse
import json

from check_bit_tagged_witness_recovery import field_for_bound, hash_horner, hash_word
from check_carry_fourier_elimination import irreducible
from coordinate_inverse_words import trace_dual_basis
from native_parity_hashes import (
    masked_parity_value, native_parity_state, smallest_factor_from_native_parities,
)
from teichmuller_cyclic_carries import UnramifiedRing

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/native-parity-hash-controls.json'
MODELS = {1:3,2:7,3:11,4:19,10:0x409}


def run():
    records,pointwise,masks,queries_total = [],0,0,0
    digest = sha256()
    cases = [(n,1) for n in range(5,64,2)]
    cases += [(105,1),(143,3),(221,3),(323,4),(67*131*257,6)]
    for target,t in cases:
        b = target.bit_length()
        d = (b+1)//2-1
        polynomial = MODELS[d]
        assert polynomial == 3 or irreducible(polynomial)
        assert all(target % p for p in range(2,(1 << t)+1))
        coefficients = [(polynomial >> i) & 1 for i in range(d+1)]
        ring = UnramifiedRing(1,coefficients)
        dual = tuple(ring.teichmuller(c) for c in trace_dual_basis(coefficients))
        active,expected_words = [],[]
        for word in range(1 << d):
            state = native_parity_state(target,ring,dual,ring.teichmuller(word))
            digits,flag = state
            assert sum(bit << i for i,bit in enumerate(digits)) == word
            x = 1+2*word
            expected = int(x > 1 and target % x == 0 and target//x > 1)
            assert flag == expected
            if flag:
                active.append(state)
                expected_words.append(word)
            pointwise += 1
        bound = 1 << min(d,(b+t-1)//t)
        control = field_for_bound(d,bound)
        isolated = set()
        for point in range(control.size):
            basis = control.hash_basis(point,d)
            buckets = {}
            for word in expected_words:
                h = hash_horner(word,point,control,d)
                assert h == hash_word(word,basis)
                buckets.setdefault(h,[]).append(word)
                digits = tuple((word >> i) & 1 for i in range(d))
                assert masked_parity_value(digits,1,basis,control.degree,h) == 1
                assert masked_parity_value(digits,1,basis,control.degree,h ^ 1) == 0
                masks += 2
            isolated.update(values[0] for values in buckets.values() if len(values) == 1)
        assert isolated == set(expected_words)
        calls = 0
        def query(n,basis,h,bucket,tag):
            nonlocal calls
            assert n == target and h == control.degree
            actual = sum(masked_parity_value(*state,basis,h,bucket,tag)
                         for state in active) % 2
            expected = sum(int(hash_word(word,basis) == bucket)
                           *(1 if tag is None else (word >> tag & 1))
                           for word in expected_words) % 2
            assert actual == expected
            calls += 1
            digest.update(json.dumps([n,basis,h,bucket,tag,actual]).encode())
            return actual
        factor,used = smallest_factor_from_native_parities(target,query,t)
        assert factor == next(x for x in range(2,target+1) if target % x == 0)
        assert calls == used <= control.size**2*(d+1)
        queries_total += calls
        record = dict(target=target,trial_bits=t,degree=d,polynomial=polynomial,
                      coefficient_precision=1,complete_native_points=1 << d,
                      witnesses=len(active),bound=bound,control_degree=control.degree,
                      control_size=control.size,parity_queries=calls,recovered_factor=factor)
        if target == 67*131*257:
            basis = control.hash_basis(0,d)
            collided = [w for w in expected_words if hash_word(w,basis) == 1]
            assert [1+2*w for w in collided] == [67,131]
            assert sum(masked_parity_value(*state,basis,control.degree,1)
                       for state in active) % 2 == 0
            record['point_zero_even_bucket'] = dict(bucket=1,hidden_factors=[67,131],parity=0)
        records.append(record)

    # A multi-witness bucket can have scalar one and misleading XOR tags.
    # This deliberately tests a zero mask, not a point of the acquired family.
    words = [1,3,5]
    digits = [tuple(w >> i & 1 for i in range(3)) for w in words]
    scalar = sum(masked_parity_value(bits,1,(0,0,0),1,0) for bits in digits) % 2
    tags = [sum(masked_parity_value(bits,1,(0,0,0),1,0,i) for bits in digits) % 2
            for i in range(3)]
    fake = 1+2*sum(bit << i for i,bit in enumerate(tags))
    assert scalar == 1 and fake == 15 and 231 % fake != 0
    def forbidden(*args):
        raise AssertionError('Pretested input must not request a parity')
    for target,t,factor in ((2,1,2),(3,1,3),(105,4,3),(143,4,11)):
        assert smallest_factor_from_native_parities(target,forbidden,t) == (factor,0)

    names = ['tools/native_parity_hashes.py','tools/check_native_parity_hashes.py',
             'tools/coordinate_inverse_words.py','tools/check_bit_tagged_witness_recovery.py',
             'tools/teichmuller_cyclic_carries.py','tools/witt_trace_word_circuits.py',
             'tools/check_carry_fourier_elimination.py']
    return dict(scope='Complete finite pure-field native predicates and directly masked parity recovery; the cached exhaustive state oracle is not a fast contractor.',
                records=records,pointwise_controls=pointwise,hash_mask_controls=masks,
                parity_queries=queries_total,query_digest=digest.hexdigest(),
                false_singleton_guard=dict(target=231,actual_candidates=[3,7,11],
                                           zero_mask_outside_acquired_family=True,
                                           scalar=scalar,decoded_candidate=fake,verified=False),
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
    print(json.dumps({k:receipt[k] for k in ('scope','pointwise_controls','hash_mask_controls','parity_queries')}))


if __name__ == '__main__':
    main()
