#!/usr/bin/env python3
"""Independent literal-divisor controls for the supplied semiprime decoder."""
import argparse
from hashlib import sha256
from itertools import combinations
import json
from math import gcd, isqrt, prod
from pathlib import Path

from semiprime_single_field import component_data, squared_factor_residues, recover_semiprime


def supplied_component(integer, prime):
    # Finite fixture construction only. The decoder never calls this function.
    _, _, _, logs, powers = component_data(prime)
    result = 0
    for divisor in range(1, isqrt(integer)+1):
        if integer % divisor == 0:
            result ^= powers[logs[divisor % prime]]
            other = integer//divisor
            if other != divisor:
                result ^= powers[logs[other % prime]]
    return result


def run():
    primes = [2,3,5,7,11,13,17,19,23,29,31,37,43]
    records, digest, count, doubles = [], sha256(), 0, 0
    for level in [7,11,19,23,31,47,503]:
        local, repeated = 0, 0
        for left, right in combinations(primes, 2):
            integer = left*right
            if integer % level == 0:
                continue
            observed = supplied_component(integer, level)
            residues = squared_factor_residues(integer, level, observed)
            assert residues == tuple(sorted([left*left % level, right*right % level]))
            repeated += int(residues[0] == residues[1])
            digest.update(json.dumps([integer,level,observed,residues]).encode())
            local += 1
        order, field, _, _, _ = component_data(level)
        records.append(dict(prime=level,odd_order=order,field_polynomial=field,
                            field_degree=field.bit_length()-1,semiprimes=local,
                            equal_character_roots=repeated))
        count += local
        doubles += repeated
    assert doubles > 0
    crt = []
    for left, right in [(101,113),(1009,1013),(2,1009)]:
        integer = left*right
        supplied = [(q,supplied_component(integer,q))
                    for q in [7,11,19,23,31,43,47,59,67,71,79,83] if gcd(q,integer)==1]
        assert prod(q for q,_ in supplied) > integer**2+1
        assert recover_semiprime(integer,supplied) == (left,right)
        crt.append(dict(integer=integer,factors=[left,right],levels=[q for q,_ in supplied]))
    # Independent input checks, rather than silently accepting a short lift.
    rejected = []
    for label, call in [
            ('wrong_level',lambda:squared_factor_residues(15,13,0)),
            ('noncoprime',lambda:squared_factor_residues(21,7,0)),
            ('noncanonical',lambda:squared_factor_residues(15,7,4)),
            ('short_crt',lambda:recover_semiprime(15,[(7,supplied_component(15,7))])),
            ('duplicate_level',lambda:recover_semiprime(15,[(7,0),(7,0)]))]:
        try:
            call()
        except ValueError:
            rejected.append(label)
        else:
            raise AssertionError(label)
    return dict(scope='Supplied-coefficient decoder under the distinct-semiprime promise; literal divisors construct finite fixtures only. No fast coefficient engine or arbitrary-integer bound.',
                component_controls=records,semiprimes_checked=count,
                equal_character_roots=doubles,digest=digest.hexdigest(),
                crt_controls=crt,rejected=rejected,
                source_sha256={name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                               for name in ('semiprime_single_field.py','check_semiprime_single_field.py',
                                            'sparse_faithful_character.py','cyclotomic_character_decoder.py',
                                            'odd_character_residues.py')})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='passed',semiprimes_checked=result['semiprimes_checked'],
                          equal_character_roots=result['equal_character_roots'],
                          crt_controls=len(result['crt_controls']),factoring_exponent='1/2; unchanged')))


if __name__ == '__main__':
    main()
