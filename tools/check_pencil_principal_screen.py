#!/usr/bin/env python3
"""Check every principal-screen input/curve key and sample independent profiles."""
import argparse
from hashlib import sha256
from itertools import product
import json
import lzma
from pathlib import Path

from check_joint_energy_conjugations import cut_matrix
from check_scaled_asymmetric_pressure import inputs
from pencil_principal_alphabets import pressure_key, restrict
from weighted_smith_oracle_f4 import WeightedSmithOracle


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    path = root/'evidence/pencil-principal-screen-complete.json.xz'
    packet = json.loads(lzma.decompress(path.read_bytes()))
    summary = packet['summary']
    assert summary['status']=='complete_registered_screen' and not summary['unresolved_profiles']
    assert sha256(packet['profiles_jsonl'].encode()).hexdigest()==summary['profile_ledger_sha256']
    assert sha256(packet['scalar_fits_json'].encode()).hexdigest()==summary['scalar_fits_sha256']
    for n,h in summary['source_sha256'].items():
        assert sha256(summary['as_run_sources'][n].encode()).hexdigest()==h
    manifest = summary['registration']
    assert sha256((root/manifest['source_manifest']).read_bytes()).hexdigest()==manifest['source_manifest_sha256']
    original,cases = inputs(root)
    subsets = [[i for i in range(5) if mask >> i & 1] for mask in range(32) if mask.bit_count()>=2]
    assert manifest['subset_masks']==[sum(1<<i for i in s) for s in subsets]
    assert len(cases)==67 and len(subsets)==26
    records = [json.loads(line) for line in packet['profiles_jsonl'].splitlines()]
    fits = json.loads(packet['scalar_fits_json'])
    assert len(records)==summary['completed_cases']==manifest['registered_cases']==45292
    seen,order_count,position,proper = {},0,0,0
    for (case,values),ax,ay in product(cases,subsets,subsets):
        record = records[position]
        assert record['source_case_index']==case['index'] and record['alphabet_x']==ax and record['alphabet_y']==ay
        seed,parameters,gx,gy = restrict(values,original['natural_heights_and_scale'],ax,ay)
        assert record['parameters']==parameters and record['physical_energies']==[gx,gy]
        assert sha256(json.dumps(seed['values'],separators=(',',':')).encode()).hexdigest()==record['literal_values_sha256']
        profiles = record['orders']
        assert len(profiles)==2 and all(len(p)==len(ax)*len(ay) and p==sorted(p) for p in profiles)
        assert all(type(v) is int for p in profiles for v in p)
        key,index = pressure_key(profiles,gx,gy),record['scalar_fit_index']
        assert 0<=index<len(fits)
        if key in seen:
            assert seen[key]==index
        else:
            assert index==len(seen)
            seen[key]=index
        if ax==ay==list(range(5)):
            assert sum(profiles,[])==case['natural_orders']
        order_count += sum(map(len,profiles))
        proper += int(len(ax)<5 or len(ay)<5)
        position += 1
    assert order_count==summary['resolved_orders']==manifest['expected_complete_orders']==753750
    assert len(seen)==len(fits)==summary['scalar_curve_classes']==6490
    assert proper==45225
    # Extract from the separate complete five-letter matrices; do not build
    # the independent matrices through the selected-entry constructor.
    independent,independent_orders = [],0
    for index in [i*len(records)//17 for i in range(1,17)]:
        r = records[index]
        ax,ay = r['alphabet_x'],r['alphabet_y']
        a,b = len(ax),len(ay)
        source = cases[r['source_case_index']][1]
        stride,offsets = 2*(a+b-2),[0,a-1,a+b-2,a+2*b-3]
        exact = []
        for place in range(2):
            h = [[0]+r['parameters'][place*stride+o:place*stride+o+q-1]
                 for o,q in zip(offsets,(a,b,b,a))]
            full = cut_matrix(source,place==0)
            ri = [x+5*u for u in ay for x in ax]
            ci = [y+5*v for v in ax for y in ay]
            if place:
                ri = [x+5*y for y in ay for x in ax]
                ci = [v+5*u for u in ay for v in ax]
            matrix = [[full[i][j] for j in ci] for i in ri]
            row = [h[0][i%a]+h[2 if place==0 else 1][i//a] for i in range(a*b)]
            col = ([h[1][i%b]+h[3][i//b] for i in range(a*b)] if place==0 else
                   [h[3][i%a]+h[2][i//a] for i in range(a*b)])
            result,_ = WeightedSmithOracle(matrix,infinity=place==0).profile(row,col)
            assert list(result)==r['orders'][place]
            exact.append(list(result))
            independent_orders += len(result)
        independent.append({'record_index':index,'source_case_index':r['source_case_index'],
                            'alphabet_x':ax,'alphabet_y':ay,'independent_profiles':exact})
    names = ['tools/check_pencil_principal_screen.py','tools/pencil_principal_alphabets.py',
             'tools/check_joint_energy_conjugations.py','tools/weighted_smith_oracle_f4.py']
    report = {'status':'pass','registered_cases':45292,'proper_restrictions':proper,
              'resolved_native_orders':order_count,'exact_inherited_scalar_curve_classes':len(seen),
              'independent_profile_samples':16,'independent_untruncated_orders':independent_orders,
              'independent_controls':independent,'archive_sha256':sha256(path.read_bytes()).hexdigest(),
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'scope':'Checks all coverage identities, literal-input hashes, restored parameters '
                      'and exact scalar-curve cache keys. Sixteen profiles in both cuts are '
                      'independently reconstructed; the remaining native orders are not. '
                      'No global numerical optimality or improved factoring bound.'}
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('status','registered_cases','resolved_native_orders',
        'exact_inherited_scalar_curve_classes','independent_untruncated_orders')}))


if __name__ == '__main__':
    main()
