#!/usr/bin/env python3
"""Validate the complete low-order census and its higher-precision selected checks."""
import itertools,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
from probe_f8_hyperoval_retraction import maps

def strip_times(x):
    if isinstance(x,dict):return {k:strip_times(v)for k,v in x.items()if k!='seconds'}
    if isinstance(x,list):return list(map(strip_times,x))
    return x

def main():
    data=json.loads(Path('evidence/f8-hyperoval-charge3-retraction-screen.json').read_text());manifest=json.loads(Path('evidence/f8-hyperoval-retraction-census-manifest.json').read_text())
    assert data['completed']and data['completed_trials']==120 and manifest['completed']
    assert sha256(json.dumps(strip_times(data),sort_keys=True,separators=(',',':')).encode()).hexdigest()==manifest['semantic_data_sha256']
    for n,h in manifest['source_sha256'].items():assert sha256(Path('tools',n).read_bytes()).hexdigest()==h,n
    triples=list(itertools.combinations(range(10),3));triples.remove((0,8,9));triples=[(0,8,9)]+triples
    assert data['registered_trials']==list(map(list,triples))and [r['pivot_columns']for r in data['records']]==list(map(list,triples))
    classes=Counter();best=[];reference=data['records'][0]
    for r in data['records']:
        maps(r['pivot_columns'])
        key=[]
        for c in r['cuts']:
            assert all(0<=int(k)<=63 and type(v)is int and v>0 for k,v in c['counts'].items())
            assert sum(c['counts'].values())+c['tail_count']==932
            key.append(tuple(sorted((int(k),v)for k,v in c['counts'].items())))
        classes[tuple(key)]+=1
        if r['orders63_total']==89:best.append(r)
    assert sorted(classes.values())==[4,8,8,16,32,52]
    assert [r['pivot_columns']for r in best]==[[6,7,8],[6,7,9],[6,8,9],[7,8,9]]
    for r in best:
        for i,c in enumerate(r['cuts']):
            expected={int(k):v for k,v in reference['cuts'][i]['counts'].items()}
            if i==0:expected[61]-=1;expected[63]-=1
            else:expected[63]-=2
            assert {int(k):v for k,v in c['counts'].items()}==expected
            for other in data['records']:
                assert all(v<=other['cuts'][i]['counts'].get(str(k),0)for k,v in expected.items())
    for tag,S in [('678',[6,7,8]),('789',[7,8,9])]:
        full=json.loads(Path('evidence/f8-hyperoval-retraction-'+tag+'.json').read_text());assert full['completed']and full['pivot_columns']==S
        rec=next(r for r in data['records']if r['pivot_columns']==S)
        charges=sorted({a+b for a in full['profile']['charges']for b in full['profile']['charges']});index=charges.index(3)
        for i,cut in enumerate(full['construction']['cuts']):
            block=cut['blocks'][index];assert block['size']==932 and block['tail_count']==0
            counts=Counter(v for v in block['retained_orders_t']if v<=63)
            assert counts==Counter({int(k):v for k,v in rec['cuts'][i]['counts'].items()})
    print(json.dumps({'completed':True,'coordinate_retractions':120,'low_profile_classes':6,'jointly_best_truncated_profiles':4,'higher_precision_complete_charge3_checks':4,'full_census_recomputed_this_run':False,'global_full_profile_optimum_claimed':False}))
if __name__=='__main__':main()
