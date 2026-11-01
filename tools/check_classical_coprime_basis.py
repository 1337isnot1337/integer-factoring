#!/usr/bin/env python3
"""Check classical coprime refinement against independent valuation oracles."""
import argparse
from functools import reduce
from hashlib import sha256
import json
import lzma
from math import gcd, prod
from pathlib import Path
from random import Random
from time import perf_counter, process_time

from classical_coprime_basis import CoprimeArithmetic, factor_into_coprimes, classical_oriented_coprime_basis

ROOT = Path(__file__).resolve().parents[1]


def known_factorization(value):
    factors, divisor = {}, 2
    while divisor*divisor<=value:
        while value%divisor==0:
            factors[divisor] = factors.get(divisor,0)+1
            value //= divisor
        divisor += 1
    if value>1:
        factors[value] = factors.get(value,0)+1
    return factors


def valuation_oracle(primes, rows):
    groups = {}
    for column,prime in enumerate(primes):
        vector = [row[column] for row in rows]
        content = reduce(gcd,vector,0)
        if content:
            direction = tuple(value//content for value in vector)
            groups[direction] = groups.get(direction,1)*prime**content
    return sorted(groups.values())


def pair_controls():
    factors = {n:known_factorization(n) for n in range(1,101)}
    checks, greatest_depth = 0, 0
    for a in range(1,101):
        for b in range(1,101):
            primes = sorted(set(factors[a])|set(factors[b]))
            rows = [[factors[value].get(p,0) for p in primes] for value in (a,b)]
            arithmetic = CoprimeArithmetic()
            actual = sorted(arithmetic.pair_basis(a,b))
            assert actual==valuation_oracle(primes,rows),(a,b,actual)
            greatest_depth = max(greatest_depth,arithmetic.statistics['maximum_pair_depth'])
            checks += 1
    return {'two_integer_pairs':checks,'maximum_pair_depth':greatest_depth}


def family_controls():
    rng, checks, words, max_bits = Random(202611011414), 0, 0, 0
    primes = [2,3,5,7,11,13,17,19]
    for trial in range(240):
        width = rng.randrange(1,len(primes)+1)
        rows = [[rng.randrange(8) if rng.randrange(3) else 0 for _ in range(width)]
                for _ in range(rng.randrange(1,24))]
        if trial%3==0 and width>=2:
            for row in rows:
                row[1] = 2*row[0]
        values = [prod(prime**e for prime,e in zip(primes,row)) for row in rows]
        values.extend(values[:2])
        rows.extend(rows[:2])
        result = factor_into_coprimes(values)
        assert result['basis']==valuation_oracle(primes[:width],rows)
        assert sum(map(len,result['words'])) <= sum(value.bit_length()-1 for value in values)
        max_bits = max(max_bits,max(map(int.bit_length,values)))
        checks += 1
        words += len(values)
    for values in ([],[1],[1,1,1],[12],[12,18,12,1],[36,216]):
        factor_into_coprimes(values)
        checks += 1
        words += len(values)
    rejected = 0
    for values in ([0],[-1],[1,2.5]):
        try:
            factor_into_coprimes(values)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError(values)
    for value,basis in ((3,[9]),(3,[15]),(15,[3]),(25,[5,7,11])):
        if value==25:
            assert CoprimeArithmetic().factor_family([25],basis)=={25:[(5,2)]}
            continue
        try:
            CoprimeArithmetic().factor_family([value],basis)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError((value,basis))
    return {'integer_families':checks,'complete_sparse_words':words,
            'maximum_input_bits':max_bits,'malformed_or_insufficient_inputs_rejected':rejected}


def power_controls():
    cases = []
    for scale in (8,32,128,512,2048,8192):
        rows = [(scale,scale-1),(scale-1,scale//2)]
        values = [2**a*3**b for a,b in rows]
        result = factor_into_coprimes(values)
        assert result['basis']==valuation_oracle([2,3],rows)
        cases.append({'maximum_input_bits':max(map(int.bit_length,values)),
                      'pair_depth':result['statistics']['maximum_pair_depth'],
                      'gcd_calls':result['statistics']['gcd_calls'],
                      'gcd_operand_bits':result['statistics']['gcd_operand_bits'],
                      'norm_words':result['words']})
    return cases


def ideal_controls():
    cases = []
    for delta in (-4,-23,-84,-260,-575):
        keys = [(a,b) for a in range(1,120,2) if gcd(a,delta)==1
                for b in range(2*a) if (b*b-delta)%(4*a)==0]
        result = classical_oriented_coprime_basis(delta,keys)
        cases.append({'discriminant':delta,'keys':len(keys),'atoms':len(result['atoms']),
                      'orientation_splits':result['orientation_splits'],
                      'verification':result['verification']})
    return cases


def fixture(bits, prefix=None):
    name = f'evidence/sparse-oriented-power-{bits}.json.xz' if bits in (64,80) else f'evidence/optimized-oriented-power-{bits}.json.xz'
    saved = json.loads(lzma.decompress((ROOT/name).read_bytes()))['result']['refinement']
    keys = saved['keys'] if prefix is None else saved['keys'][:prefix]
    started, cpu = perf_counter(),process_time()
    result = classical_oriented_coprime_basis(saved['discriminant'],keys,verify=False)
    wall, elapsed_cpu = perf_counter()-started,process_time()-cpu
    from oriented_norm_basis_lift import verify_lifted_oriented_basis
    verification = verify_lifted_oriented_basis(result)
    matched = None
    if prefix is None:
        matched = all(json.loads(json.dumps(result[key]))==saved[key] for key in ('keys','atoms','sparse_words'))
    return {'bits':bits,'input_path':name,'input_sha256':sha256((ROOT/name).read_bytes()).hexdigest(),
            'prefix':prefix,'keys':len(keys),'atoms':len(result['atoms']),
            'same_retained_complete_basis':matched,'wall_seconds':wall,'cpu_seconds':elapsed_cpu,
            'norm_construction':result['norm_construction'],
            'orientation_gcd_calls':result['statistics']['orientation_gcd_calls'],
            'verification':verification}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture-bits',type=int,nargs='*',default=[64])
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    started = perf_counter()
    report = {'scope':'Classical supplied-integer coprime construction and oriented ideal output; no new factoring exponent.',
              'pairs':pair_controls(),'families':family_controls(),'power_ladders':power_controls(),
              'ideals':ideal_controls(),'fixtures':[fixture(bits) for bits in args.fixture_bits]}
    sources = ['classical_coprime_basis.py','check_classical_coprime_basis.py','oriented_norm_basis_lift.py',
               'oriented_coprime_basis.py','binary_quadratic_forms.py','quadratic_ideal_hnf.py']
    report['source_sha256'] = {'tools/'+name:sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in sources}
    report['seconds'] = perf_counter()-started
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
