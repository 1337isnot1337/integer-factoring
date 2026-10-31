#!/usr/bin/env python3
"""Exact controls for the universal multipoint local-jet Smith floor."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations_with_replacement,product
import json
from pathlib import Path
from random import Random

from check_multipoint_hermite_crossings import (
    check_power,local_profiles,multiply,
)
from search_multipoint_rank_one_seeds import KEYS,SWAP,pure_tensors,representatives

ROOT = Path(__file__).resolve().parents[1]


def mean_budget(first,second):
    assert len(first) == len(second) and min(first) == min(second)
    return Fraction(sum(first)+sum(second),2*len(first))


def distribution_control(first,second):
    n = len(first)
    d = mean_budget(first,second)
    minimum = min(first)
    small,large = sorted((first,second),key=sum)
    assert Fraction(sum(small),n) <= d
    assert minimum + Fraction(sum(x-minimum for x in large),2*n) <= d
    return d


def exact_pressure(first,second,s,r):
    n = len(first)
    d = distribution_control(first,second)
    u = sum(Fraction(s)**(-x) for x in first)
    v = sum(Fraction(r)**(-x) for x in second)
    left = (u*v)**(2*n)*Fraction(s*r)**int(2*n*d)
    right = Fraction((2*n*n)**(2*n),(n+1)**(n+1))
    assert left >= right


def abstract_controls():
    profiles = [(0,)+tail for tail in combinations_with_replacement(range(6),3)]
    count = 0
    for first,second in product(profiles,repeat=2):
        for s,r in ((1,1),(2,3),(7,2)):
            exact_pressure(first,second,s,r)
            count += 1
    random = Random(1337)
    for n in (2,3,9,16,25):
        for _ in range(100):
            minimum = random.randrange(4)
            first = [minimum]+[minimum+random.randrange(20) for _ in range(n-1)]
            second = [minimum]+[minimum+random.randrange(20) for _ in range(n-1)]
            exact_pressure(first,second,2,3)
            count += 1
    # A different common minimum on the two cuts destroys the statement.
    first,second = (1,)*4,(0,)*4
    d = Fraction(sum(first)+sum(second),8)
    left = (4*Fraction(100)**-1*4)**8*100**int(8*d)
    right = Fraction(32**8,5**5)
    assert left < right
    try:
        mean_budget(first,second)
    except AssertionError:
        pass
    else:
        raise AssertionError('unequal-minimum negative control was accepted')
    # Sharp abstract profile corner, with s^(2n)=n+1 and r=1.
    for n in range(2,129):
        d = distribution_control([0]+[2*n]*(n-1),[0]*n)
        assert d == n-1
        attained = Fraction(2*n*n,n+1)**(2*n)*(n+1)**(n-1)
        assert attained == Fraction((2*n*n)**(2*n),(n+1)**(n+1))
        assert 4**n*n**n > (n+1)**(n+1)
    return count


def literal_controls():
    random = Random(314159)
    cases = profiles_checked = inequalities = 0
    for prime,alphabet in ((3,2),(5,2),(7,2),(3,3)):
        keys = tuple(product(range(alphabet),repeat=4))
        endpoint = {key:int(key[0] == key[3] and key[1] == key[2]) for key in keys}
        for degree in range(1,5 if alphabet == 2 else 3):
            for _ in range(4):
                seed = {}
                for key in keys:
                    value = [random.randrange(prime) for _ in range(degree)]
                    value.append((endpoint[key]-sum(value)) % prime)
                    seed[key] = value
                if cases % 4 == 0:
                    # A literal scalar factor checks positive common minima.
                    scalar = [0,1]
                    seed = {key:multiply(value,scalar,prime) for key,value in seed.items()}
                    actual_degree = degree+1
                else:
                    actual_degree = degree
                assert all(sum(value) % prime == endpoint[key] for key,value in seed.items())
                points = [x for x in range(prime) if x != 1]+[None]
                local = [local_profiles(seed,alphabet,actual_degree,point,prime) for point in points]
                n = alphabet**2
                for cut in (0,1):
                    assert sum(sum(pair[cut]) for pair in local) <= n*actual_degree
                budgets = [mean_budget(*pair) for pair in local]
                assert sum(budgets) <= actual_degree
                for pair in local:
                    exact_pressure(*pair,2,3)
                    profiles_checked += 2
                    inequalities += 1
                # Nonuniform rational degree allocations, including omitted points.
                for _ in range(12):
                    weights = [random.randrange(4) for _ in points]
                    if not sum(weights):
                        weights[0] = 1
                    allocation = [Fraction(actual_degree*w,sum(weights)) for w in weights]
                    assert any(d > 0 and d >= b for d,b in zip(allocation,budgets))
                cases += 1
    return cases,profiles_checked,inequalities


def cubic_controls(path):
    saved = json.loads(path.read_text())
    assert saved['complete'] and saved['processed'] == saved['total'] == 1579014
    assert saved['signatures'] == len(saved['catalogue']) == 200
    assert sum(record['count'] for record in saved['catalogue']) == saved['total']
    assert saved['prefix_certificate_failures'] == 0
    tensors,reps = pure_tensors(),representatives()
    for record in saved['catalogue']:
        ai,bi,ci = record['indices']
        a,b,c = reps[ai],tensors[bi],tensors[ci]
        seed = {key:(x,(y-e-z)%3,(-x-y-e)%3,z)
                for key,x,y,z,e in zip(KEYS,a,b,c,SWAP)}
        budgets = []
        for j,point in enumerate((0,2,None)):
            pair = local_profiles(seed,2,3,point,3)
            assert pair == record['profiles'][2*j:2*j+2]
            budgets.append(mean_budget(*pair))
            exact_pressure(*pair,2,3)
            numerator,denominator,k,l = record['prefix_certificates'][j]
            bound = Fraction(numerator,denominator)
            assert k*l >= 8
            assert Fraction(sum(pair[0][:k]),k) <= bound
            assert Fraction(sum(pair[1][:l]),l) <= bound
        assert sum(budgets) <= 3
        assert sum(Fraction(a,b) for a,b,_,_ in record['prefix_certificates']) <= 3
        check_power(seed,2,3,2,[0,2,None],[3,2,2],3)
    return len(saved['catalogue'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    abstract = abstract_controls()
    cases,profiles_checked,inequalities = literal_controls()
    cubic = cubic_controls(ROOT/'evidence/multipoint-pure-cubics.json')
    result = {'status':'pass','abstract_exact_pressure_comparisons':abstract,
              'sharp_abstract_profile_corners':127,'literal_polynomial_seeds':cases,
              'literal_cut_profiles':profiles_checked,
              'literal_exact_pressure_comparisons':inequalities,
              'complete_native_cubic_search_seeds':1579014,
              'independent_cubic_toeplitz_profiles':6*cubic,
              'reconstructed_two_copy_cubic_crossings':cubic,
              'complete_two_copy_cubic_entries':cubic*2**8,
              'universal_rho_floor':'strictly greater than 3/4 for every fixed finite alphabet',
              'scope':'The unweighted local-jet multipoint Smith pressure criterion; no factoring lower bound or exclusion of sharper actual-jet bounds.',
              'source_sha256':{name:sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in ('search_multipoint_pure_cubics.cpp',
                                            'check_multipoint_smith_floor.py',
                                            'check_multipoint_hermite_crossings.py')},
              'cubic_catalogue_sha256':sha256((ROOT/'evidence/multipoint-pure-cubics.json').read_bytes()).hexdigest()}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
