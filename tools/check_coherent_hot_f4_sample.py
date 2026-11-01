#!/usr/bin/env python3
"""Exact all-scale pressure comparisons for the saved 100-map F4 sample.

This checks integer polynomial identities on saved complete profiles. It does
not rerun the native profiles or exclude the full 5,440-map F4 catalogue.
"""
import json
from pathlib import Path
from hashlib import sha256
from collections import Counter
from check_coherent_q32_catalogue import certificate


def main():
    source=Path('evidence/coherent-hot-f4-q32.json');data=json.loads(source.read_text())
    assert data['completed'] and not data['catalogue_exhausted'] and len(data['records'])==100
    for name,expected in data['source_sha256'].items():assert sha256(Path('tools',name).read_bytes()).hexdigest()==expected,name
    old=Path('evidence/coherent-hot-q32-retractions.json');binary=json.loads(old.read_text());winner=binary['records'][55]
    assert winner['description']['kernel']==7 and winner['description']['covector']==15
    counts=Counter();proofs=[]
    for i,row in enumerate(data['records']):
        assert sorted(row['heights'])==sorted(winner['heights'])
        assert len(row['zero_orders_t'])==len(row['infinity_orders_t'])==1024
        cert=certificate(row,winner);counts[cert['method']]+=1;proofs.append({'index':i,**cert})
    report={'scope':__doc__,'completed':True,'sample_sha256':sha256(source.read_bytes()).hexdigest(),'binary_reference_sha256':sha256(old.read_bytes()).hexdigest(),'counts':dict(counts),'comparisons':proofs,'all_scale_interval':'0 < r < 1','native_profiles_recomputed':False,'full_F4_catalogue_excluded':False,'factoring_exponent_improved':False}
    Path('evidence/coherent-hot-f4-sample-comparisons.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'counts':dict(counts),'status':'passed'}))

if __name__=='__main__':main()
