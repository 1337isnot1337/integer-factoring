#!/usr/bin/env python3
"""Independent HNF controls and large-discriminant compact-basis checks."""
import argparse
from copy import deepcopy
from hashlib import sha256
import json,lzma
from math import gcd,prod
from pathlib import Path
from time import perf_counter,process_time
import sys
from binary_quadratic_forms import reduce_form
from circuit_box_oriented_batch import two_pass_oriented_basis
from check_circuit_box_oriented_batch import _mixed_keys,_prime_words
from compact_oriented_basis import compact_oriented_coprime_basis,validate_compact_keys,verify_compact_oriented_basis
from oriented_coprime_basis import _key
from oriented_norm_basis_lift import verify_lifted_oriented_basis

ROOT = Path(__file__).resolve().parents[1]


def compare_case(discriminant,keys,*,raw=False,prime_check=True):
    compact = compact_oriented_coprime_basis(discriminant,keys,unitary_preconditioning=not raw)
    original = two_pass_oriented_basis(discriminant,keys,unitary_preconditioning=not raw)
    assert compact['keys']==original['keys'] and compact['sparse_words']==original['sparse_words']
    assert compact['atoms']==[{k:a[k] for k in ('norm','root')} for a in original['atoms']]
    expanded = dict(compact,atoms=[dict(a,form=reduce_form((a['norm'],a['root'],
        (a['root']**2-discriminant)//(4*a['norm'])))) for a in compact['atoms']])
    hnf = verify_lifted_oriented_basis(expanded)
    if prime_check: _prime_words(compact['keys'],compact)
    check = verify_compact_oriented_basis(json.loads(json.dumps(compact)),
        expected_discriminant=discriminant,expected_keys=keys)
    assert check['arithmetic_literal_ideal_identities']==hnf['literal_ideal_identities']
    return {'keys':len(keys),'atoms':len(compact['atoms']),'raw_differences':raw,
            'input_bits':compact['total_input_bits'],'literal_hnf_identities':hnf['literal_ideal_identities'],
            'independent_prime_expansions':len(keys) if prime_check else 0,
            'same_complete_compact_payload':True,'validation':compact['input_validation']}


def small_controls():
    cases = []
    for delta in (-4,-23,-84,-260,-575):
        keys = [(r,b) for r in range(1,140,2) if gcd(r,delta)==1
                for b in range(2*r) if (b*b-delta)%(4*r)==0]
        cases.append(compare_case(delta,keys))
        shifted = [(r,b+2*r*((-1 if i%2 else 1)*(1<<(i%101)))) for i,(r,b) in enumerate(keys)]
        cases.append(compare_case(delta,shifted,raw=True))
    for keys in ([],[(1,0)]*7,[(65,16),(65,36)],[(325,114),(325,614)],[(169,140),(13,10),(169,198)]):
        cases.append(compare_case(-4,keys))
    _,mixed = _mixed_keys()
    cases.append(compare_case(-4,mixed,prime_check=False))
    cases.append(compare_case(-4,mixed,raw=True,prime_check=False))
    return cases


def large_discriminants():
    primes = [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]
    keys = [(p,1 if i%2 else 2*p-1) for i in range(2048) for p in [primes[i%len(primes)]]]
    observations = []
    for shift in (4096,65536,1048576):
        delta = 1-4*prod(primes)*(1<<shift)
        start,cpu = perf_counter(),process_time()
        expected = [_key(delta,*key) for key in keys]
        old_wall,old_cpu = perf_counter()-start,process_time()-cpu
        start,cpu = perf_counter(),process_time()
        actual,validation = validate_compact_keys(delta,keys)
        new_wall,new_cpu = perf_counter()-start,process_time()-cpu
        assert actual==expected
        full = compact_oriented_coprime_basis(delta,keys)
        assert [a['norm'] for a in full['atoms']]==primes
        assert all(a['root']==1 for a in full['atoms'])
        for (norm,root),word in zip(keys,full['sparse_words']):
            assert word==[(primes.index(norm),1 if root==1 else -1)]
        observations.append({'shift':shift,'discriminant_bits':delta.bit_length(),'keys':len(keys),
            'definition':'Delta=1-4*product(primes)*2^shift','primes':primes,
            'old_validation_wall_seconds':old_wall,'old_validation_cpu_seconds':old_cpu,
            'batch_validation_wall_seconds':new_wall,'batch_validation_cpu_seconds':new_cpu,
            'same_normalized_keys':True,'complete_compact_basis_verified':True,
            'validation':validation,'compact_verification':full['verification']})
    return observations


def rejections():
    rejected = 0
    for delta,keys in [(1,[]),(-3,[(2,1)]),(-4,[(5,1)]),(-20,[(5,0)]),(-4,[(0,0)])]:
        try: compact_oriented_coprime_basis(delta,keys)
        except ValueError: rejected += 1
        else: raise AssertionError('Invalid input accepted.')
    original = compact_oriented_coprime_basis(-4,[(65,16),(65,36)])
    changes = [
        lambda r:r['sparse_words'][0].__setitem__(0,(0,1<<100000)),
        lambda r:r['sparse_words'][0].append(r['sparse_words'][0][0]),
        lambda r:r['sparse_words'][0].__setitem__(0,(0,0)),
        lambda r:r['sparse_words'][0].__setitem__(0,(99,1)),
        lambda r:r['sparse_words'][0].__setitem__(0,(0,-r['sparse_words'][0][0][1])),
        lambda r:r['atoms'][0].__setitem__('root',0),
        lambda r:r['atoms'][0].__setitem__('norm',1),
        lambda r:r['atoms'][0].__setitem__('form',(1,0,1)),
        lambda r:r['atoms'].append(dict(r['atoms'][0])),
    ]
    corruptions = 0
    for change in changes:
        damaged = deepcopy(original);change(damaged)
        try: verify_compact_oriented_basis(damaged,expected_discriminant=-4,expected_keys=[(65,16),(65,36)])
        except (ValueError,AssertionError,IndexError,KeyError): corruptions += 1
        else: raise AssertionError('Corrupted certificate accepted.')
    return {'invalid_inputs':rejected,'corruptions':corruptions,'huge_binary_exponent_rejected_before_power':True}


def fixtures():
    cases = []
    for bits,limit in [(64,None),(80,128),(96,128)]:
        relative = f'evidence/sparse-oriented-power-{bits}.json.xz' if bits in (64,80) else 'evidence/optimized-oriented-power-96.json.xz'
        raw = (ROOT/relative).read_bytes();saved = json.loads(lzma.decompress(raw))['result']['refinement']
        keys = saved['keys'][:limit]
        result = compare_case(saved['discriminant'],keys,prime_check=False)
        cases.append({'input_path':relative,'input_sha256':sha256(raw).hexdigest(),
            'full_key_count':len(saved['keys']),'checked_key_count':len(keys),**result})
    return cases


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    sources=sorted({Path(m.__file__).resolve() for m in list(sys.modules.values()) if getattr(m,'__file__',None)
        and Path(m.__file__).resolve().parent==ROOT/'tools' and Path(m.__file__).suffix=='.py'})
    guards={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sources}
    start=perf_counter()
    report={'scope':'Compact supplied-ideal construction and verification; no factoring-exponent claim.',
        'small_cases':small_controls(),'large_discriminants':large_discriminants(),
        'rejections':rejections(),'fixtures':fixtures(),'source_sha256':guards}
    for p,h in guards.items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    report['seconds']=perf_counter()-start
    if args.output:args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'small_cases':len(report['small_cases']),
        'hnf_identities':sum(x['literal_hnf_identities'] for x in report['small_cases']+report['fixtures']),
        'large_discriminant_bits':[x['discriminant_bits'] for x in report['large_discriminants']],
        'rejections':report['rejections'],'seconds':report['seconds']}),flush=True)


if __name__=='__main__':main()
