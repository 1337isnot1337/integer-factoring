#!/usr/bin/env python3
"""Check orientation lifting, literal ideals, and amortized sparse writes."""
import argparse
from copy import deepcopy
from hashlib import sha256
import json
import lzma
from math import gcd, prod
from pathlib import Path
from random import Random
from time import perf_counter
from incremental_oriented_coprime_basis import IncrementalOrientedBasis
from oriented_norm_basis_lift import lift_coprime_norm_basis, verify_coprime_norms, verify_lifted_oriented_basis

ROOT = Path(__file__).resolve().parents[1]


def merge_equal_incidence(norms, words):
    columns = [[] for _ in norms]
    for i, word in enumerate(words):
        for j, exponent in word:
            columns[j].append((i, abs(exponent)))
    groups = {}
    for norm, column in zip(norms, columns):
        assert column
        groups.setdefault(tuple(column), []).append(norm)
    merged, output = [], [[] for _ in words]
    for column, factors in sorted(groups.items(), key=lambda pair: prod(pair[1])):
        index = len(merged)
        merged.append(prod(factors))
        for row, exponent in column:
            output[row].append((index, exponent))
    return merged, output


def small_norm_basis(keys):
    # Trial division supplies only these small test fixtures.
    factored = []
    for norm, _ in keys:
        factors, divisor = {}, 2
        while divisor <= norm:
            while norm % divisor == 0:
                factors[divisor] = factors.get(divisor, 0)+1
                norm //= divisor
            divisor += 1
        factored.append(factors)
    norms = sorted({p for factors in factored for p in factors})
    indices = {p:i for i,p in enumerate(norms)}
    words = [[(indices[p],e) for p,e in sorted(factors.items())] for factors in factored]
    return merge_equal_incidence(norms, words)


def write_bound(delta, keys):
    state, logarithms, nonunits = IncrementalOrientedBasis(delta), 0, 0
    for norm, root in keys:
        state.add(norm, root)
        logarithms += norm.bit_length()-1
        nonunits += norm > 1
        assert state.usage_writes <= 2*logarithms-nonunits
    return {'prefixes':len(keys), 'usage_writes':state.usage_writes,
            'proved_integer_upper_bound':2*logarithms-nonunits}


def small_controls():
    rng, coprime = Random(202611011313), 0
    for _ in range(800):
        values = [rng.randrange(2,1000) for _ in range(rng.randrange(15))]
        expected = all(gcd(a,b)==1 for i,a in enumerate(values) for b in values[i+1:])
        try:
            verify_coprime_norms(values)
            actual = True
        except ValueError:
            actual = False
        assert actual == expected
        coprime += 1
    batches = identities = prefixes = splits = 0
    for delta in (-4,-23,-84,-260,-575):
        keys = [(a,b) for a in range(1,140,2) if gcd(a,delta)==1
                for b in range(2*a) if (b*b-delta) % (4*a)==0]
        shuffled = list(keys)
        rng.shuffle(shuffled)
        for ordered in (keys,list(reversed(keys)),shuffled+shuffled[:20]):
            norms, words = small_norm_basis(ordered)
            result = lift_coprime_norm_basis(delta,ordered,norms,words)
            assert result['statistics']['orientation_gcd_calls'] <= sum(a.bit_length()-1 for a,_ in ordered)
            prefixes += write_bound(delta,ordered)['prefixes']
            identities += len(ordered)
            splits += result['orientation_splits']
            batches += 1
    assert not lift_coprime_norm_basis(-23,[],[],[])['atoms']
    return {'coprime_batches':coprime,'ordered_ideal_batches':batches,
            'literal_ideal_identities':identities,'orientation_splits':splits,
            'prefix_sparse_write_bounds':prefixes}


def mixed_powers():
    primes = [5,13,17,29,37,41,53,61]
    keys, words = [], []
    for mask in range(256):
        exponent, root, norm = 1+mask%5, 0, 1
        for i,prime in enumerate(primes):
            local = next(r for r in range(prime) if (r*r+4)%prime==0)
            modulus = prime
            for _ in range(1,exponent):
                local += (((-4-local*local)//modulus)*pow(2*local,-1,prime)%prime)*modulus
                modulus *= prime
            if (mask>>i)&1:
                local = -local%modulus
            root += ((local-root)*pow(norm,-1,modulus)%modulus)*norm
            norm *= modulus
        root += norm if root%2 else 0
        keys.append((norm,root))
        words.append([(0,exponent)])
    result = lift_coprime_norm_basis(-4,keys,[prod(primes)],words)
    assert [a['norm'] for a in result['atoms']] == primes
    return {'keys':len(keys),'maximum_norm_bits':max(a.bit_length() for a,_ in keys),
            'coarse_norms':1,'final_atoms':8,'orientation_splits':result['orientation_splits'],
            'orientation_gcd_calls':result['statistics']['orientation_gcd_calls'],
            'literal_ideal_identities':len(keys),'sparse_write_control':write_bound(-4,keys)}


def corruptions():
    keys, checks = [(65,16),(65,36)], 0
    good = lift_coprime_norm_basis(-4,keys,[65],[[(0,1)],[(0,1)]])
    for norms,words in [([65],[[(0,1)],[]]),([65],[[(0,2)],[(0,1)]]),
                        ([65,17],[[(0,1)],[(0,1)]]),([5,65],[[(1,1)],[(1,1)]])]:
        try:
            lift_coprime_norm_basis(-4,keys,norms,words)
        except (ValueError,AssertionError):
            checks += 1
        else:
            raise AssertionError((norms,words))
    for kind in ('exponent','orientation','duplicate'):
        bad = deepcopy(good)
        if kind=='duplicate':
            bad['atoms'].append(deepcopy(bad['atoms'][0]))
        else:
            j,e = bad['sparse_words'][0][0]
            bad['sparse_words'][0][0] = (j,2*e if kind=='exponent' else -e)
        try:
            verify_lifted_oriented_basis(bad)
        except (ValueError,AssertionError):
            checks += 1
        else:
            raise AssertionError(kind)
    return checks


def fixture(bits):
    name = f'evidence/sparse-oriented-power-{bits}.json.xz' if bits==80 else f'evidence/optimized-oriented-power-{bits}.json.xz'
    path = ROOT/name
    saved = json.loads(lzma.decompress(path.read_bytes()))['result']['refinement']
    norms = [a['norm'] for a in saved['atoms']]
    words = [[(j,abs(e)) for j,e in word] for word in saved['sparse_words']]
    result = lift_coprime_norm_basis(saved['discriminant'],saved['keys'],norms,words)
    for key in ('keys','atoms','sparse_words'):
        assert json.loads(json.dumps(result[key])) == saved[key],key
    merged, merged_words = merge_equal_incidence(norms,words)
    coarse = lift_coprime_norm_basis(saved['discriminant'],saved['keys'],merged,merged_words)
    bound = 2*sum(a.bit_length()-1 for a,_ in saved['keys'])-sum(a>1 for a,_ in saved['keys'])
    assert saved['statistics']['usage_writes'] <= bound
    return {'bits':bits,'input_path':name,'input_sha256':sha256(path.read_bytes()).hexdigest(),
            'keys':len(saved['keys']),'original_atoms':len(norms),
            'original_atoms_and_words_identical':True,
            'original_basis_orientation_gcds':result['statistics']['orientation_gcd_calls'],
            'merged_norm_basis':len(merged),'refined_merged_basis':len(coarse['atoms']),
            'merged_orientation_splits':coarse['orientation_splits'],
            'literal_ideal_identities':len(saved['keys'])*2,
            'retained_incremental_usage_writes':saved['statistics']['usage_writes'],
            'proved_integer_usage_write_upper_bound':bound}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    started = perf_counter()
    report = {'scope':'Supplied norm-basis orientations and linear sparse write charges; no relation-supply bound.',
              'small_controls':small_controls(),'mixed_prime_powers':mixed_powers(),
              'corruptions_rejected':corruptions(),'fixtures':[fixture(bits) for bits in (80,96)]}
    sources = ['oriented_norm_basis_lift.py','check_oriented_norm_basis_lift.py',
               'incremental_oriented_coprime_basis.py','oriented_coprime_basis.py',
               'binary_quadratic_forms.py','quadratic_ideal_hnf.py']
    report['source_sha256'] = {'tools/'+name:sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in sources}
    report['seconds'] = perf_counter()-started
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
