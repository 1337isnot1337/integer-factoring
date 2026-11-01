#!/usr/bin/env python3
"""Finite exact controls and charged equivalent-norm stencil comparisons."""
import argparse
from copy import deepcopy
from hashlib import sha256
import json
import lzma
from math import gcd,prod
from pathlib import Path
import sys
from time import perf_counter,process_time
from binary_quadratic_forms import reduced_forms,reduce_form
from cofactor_collision_relations import partial_leading_word,verify_partial_identity
from equivalent_norm_stencils import STENCILS,proper_transform,partial_represented_word,analyze_power_prefix,verify_prefix_analysis
from prime_form_factor_base import build_prime_form_factor_base

ROOT = Path(__file__).resolve().parents[1]


def small_controls():
    forms = transforms = splits = leading = unusable = 0
    for integer in range(35,206,2):
        delta = -integer if integer%4==3 else -4*integer
        base = build_prime_form_factor_base(delta,7)
        records = {r.prime:(r.root,r.positive) for r in base['records']}
        if integer%4==1: records[2] = (2,reduce_form((2,2,(integer+1)//2)))
        primes = sorted(records)
        for form in reduced_forms(delta):
            forms += 1
            for j,matrix in enumerate(STENCILS['four_norms']):
                transformed = proper_transform(form,matrix)
                assert reduce_form(transformed)==form
                (x,u),(y,v) = matrix
                for s,t in [(1,0),(0,1),(1,1),(2,-3)]:
                    X,Y = x*s+u*t,y*s+v*t
                    assert transformed[0]*s*s+transformed[1]*s*t+transformed[2]*t*t==form[0]*X*X+form[1]*X*Y+form[2]*Y*Y
                transforms += 1
                if gcd(transformed[0],integer)!=1:
                    unusable += 1
                    continue
                split = partial_represented_word(integer,transformed,primes,records)
                assert verify_partial_identity(integer,form,split,records)
                assert split['smooth_part']*split['residual']==transformed[0]
                assert gcd(split['residual'],prod(primes))==1
                splits += 1
                if j==0:
                    assert split==partial_leading_word(integer,form,primes,records)
                    leading += 1
    invalid = 0
    for matrix in [((1,0),(0,-1)),((2,0),(0,1)),((0,0),(0,0)),((1,0,0),(0,1,0))]:
        try: proper_transform((1,1,6),matrix)
        except ValueError: invalid += 1
        else: raise AssertionError('Invalid matrix accepted.')
    return {'reduced_forms':forms,'proper_transformations':transforms,
        'independent_coordinate_evaluations':4*transforms,'class_identities':splits,
        'complete_leading_decoder_equalities':leading,'charged_unusable_norms':unusable,
        'invalid_matrices_rejected':invalid}


def read_reference(bits):
    name = f'evidence/sparse-oriented-power-{bits}.json.xz' if bits in (64,80) else f'evidence/optimized-oriented-power-{bits}.json.xz'
    path = ROOT/name
    with lzma.open(path,'rt') as f: reference = json.load(f)['result']
    return name,sha256(path.read_bytes()).hexdigest(),reference


def corruption_controls(reference):
    original = analyze_power_prefix(reference,32,stencil='four_norms')
    assert verify_prefix_analysis(reference,original)['complete_declared_stencil_replayed']
    changes = [
        lambda r:r.__setitem__('norm_queries',r['norm_queries']-1),
        lambda r:r['witnesses'].pop(),
        lambda r:r['witnesses'][0].__setitem__('represented_form',(1,0,1)),
        lambda r:r['witnesses'][0].__setitem__('exponent',r['witnesses'][0]['exponent']+1),
        lambda r:r.__setitem__('old_base_rank',r['old_base_rank']+1),
        lambda r:r.__setitem__('independent_rank',r['independent_rank']+1),
        lambda r:r.__setitem__('factor',2),
        lambda r:r['transported_rows'][0].__setitem__(0,999),
    ]
    count = 0
    for change in changes:
        damaged = deepcopy(original); change(damaged)
        try: verify_prefix_analysis(reference,damaged)
        except (AssertionError,KeyError,ValueError,IndexError): count += 1
        else: raise AssertionError('Corrupted certificate accepted.')
    return count


def fixture_comparisons(bits,budgets):
    name,digest,reference = read_reference(bits)
    cases = []
    for budget in budgets:
        if budget>len(reference['witnesses']): raise ValueError('Budget exceeds the retained leading prefix.')
        for stencil in STENCILS:
            length = len(STENCILS[stencil])
            if budget%length: raise ValueError('A norm budget must contain complete stencils.')
            result = analyze_power_prefix(reference,budget//length,stencil=stencil)
            wall,cpu = perf_counter(),process_time()
            verification = verify_prefix_analysis(reference,json.loads(json.dumps(result)))
            replay_wall,replay_cpu = perf_counter()-wall,process_time()-cpu
            case = {k:result[k] for k in ('stencil','power_prefix','norm_queries','trial_divisions',
                'maximum_norm_bits','old_base_rank','core_columns','independent_rank',
                'factor','relation_factor','analysis_wall_seconds','analysis_cpu_seconds')}
            case.update(keys=len(result['keys']),atoms=len(result['refinement']['atoms']),
                retained_norms=len(result['witnesses']),unusable_norms=len(result['rejected_norms']),
                replay_wall_seconds=replay_wall,replay_cpu_seconds=replay_cpu,verification=verification)
            cases.append(case)
            print(json.dumps({'bits':bits,**{k:case[k] for k in ('stencil','power_prefix','norm_queries','independent_rank','factor')}}),flush=True)
    return {'bits':bits,'reference_file':name,'reference_sha256':digest,'norm_budgets':budgets,'cases':cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bits',type=int,nargs='*',default=[64,80])
    parser.add_argument('--budgets',type=int,nargs='*')
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    source_paths = sorted({Path(m.__file__).resolve() for m in list(sys.modules.values())
        if getattr(m,'__file__',None) and Path(m.__file__).resolve().parent==ROOT/'tools' and Path(m.__file__).suffix=='.py'})
    hashes = {str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in source_paths}
    started = perf_counter()
    controls = small_controls()
    _,_,reference = read_reference(64)
    corruptions = corruption_controls(reference)
    defaults = {64:[256,512,1024,2048],80:[1024,2048,4096,8192],96:[2048,4096,8192,16384]}
    fixtures = [fixture_comparisons(bits,args.budgets or defaults[bits]) for bits in args.bits]
    for p,digest in hashes.items(): assert sha256((ROOT/p).read_bytes()).hexdigest()==digest,p
    report = {'claim':'Proper variable changes give exact extra norm candidates on supplied powers.',
        'controls':controls,'corrupted_records_rejected':corruptions,'fixtures':fixtures,'source_sha256':hashes,
        'seconds':perf_counter()-started,'scope':'Finite supplied-prefix analysis. All norms are charged; no relation-supply or factoring-exponent claim.',
        'timing_scope':'Analysis includes every declared norm and one final presentation. Source-power generation, loading, and independent replay are separate; this is not an end-to-end search comparison.'}
    if args.output: args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'controls':controls,'corruptions':corruptions,'seconds':report['seconds']}),flush=True)


if __name__=='__main__': main()
