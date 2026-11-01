#!/usr/bin/env python3
"""Reconstruct the F8 q96 hyperoval profile or validate its source-pinned replay."""
import argparse,json,itertools,subprocess,sys,time
from pathlib import Path
from hashlib import sha256
from coherent_f8_maps import height,mul,basis,validate_maps
from coherent_hot_compounds_f8 import partial_profiles
SOURCES=['replay_f8_hyperoval_q96.py','check_f8_hyperoval_entries.py','coherent_f8_maps.py','coherent_hot_compounds_f8.py','coherent_hot_compounds_f8.cpp','selected_gaussian_coefficients.cpp','native_rectangular_smith_f4.cpp','coherent_hot_compounds.py']
def hashes():return {n:sha256(Path('tools',n).read_bytes()).hexdigest()for n in SOURCES}
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--full',action='store_true');args=ap.parse_args()
    path=Path('evidence/f8-hyperoval-q96-full-profile.json');saved=json.loads(path.read_text());assert saved['completed']
    lower=[w for w in range(4096)if height(w)<25];boundary=sorted(w for w in range(4096)if height(w)==25 and w.bit_count()==3)
    assert len(lower)==92 and len(boundary)==10
    columns=[(1,t,mul(t,t))for t in range(8)]+[(0,1,0),(0,0,1)]
    for three in itertools.combinations(columns,3):assert len(basis([{i:c for i,c in enumerate(v)if c}for v in three]))==3
    p=[{w:1}for w in lower+[31]]+[{boundary[j]:columns[j][i]for j in range(10)if columns[j][i]}for i in range(3)]
    l=[{w:1}for w in lower+[31]]+[{boundary[j]:1}for j in (0,8,9)]
    assert json.loads(json.dumps(p))==saved['embedding'] and json.loads(json.dumps(l))==saved['retraction']
    hs,ch=validate_maps(p,l,12);assert hs==saved['heights']and ch==saved['charges']
    assert all(len(saved[k])==9216 for k in ['zero_orders_t','infinity_orders_t'])
    initial=hashes();receipt=Path('evidence/f8-hyperoval-q96-clean-replay.json')
    if args.full:
        start=time.monotonic();result=partial_profiles(p,l,12,127,progress=lambda v:print(json.dumps(v),flush=True))
        assert result['all_orders_resolved']
        assert all(c['tail_count']==0 for c in result['cuts'])
        for key,c in zip(['zero_orders_t','infinity_orders_t'],result['cuts']):assert c['retained_orders_t']==saved[key]
        assert hashes()==initial,'Sources changed during replay'
        record={'completed':True,'all_full_orders_replayed':18432,'all_120_hyperoval_minors_checked':True,'all_endpoints_checked_by_constructor':True,'cutoff':127,'tail_counts':[0,0],'profile_sha256':sha256(path.read_bytes()).hexdigest(),'source_sha256':initial,'seconds':time.monotonic()-start,'scope':'Both complete finite profiles were reconstructed by the exact F8 engine. Not a factoring-exponent improvement or independent mathematical review.'}
        receipt.write_text(json.dumps(record,indent=2)+'\n')
    else:
        record=json.loads(receipt.read_text());assert record['completed']and record['profile_sha256']==sha256(path.read_bytes()).hexdigest()and record['source_sha256']==initial
    subprocess.run([sys.executable,'tools/check_f8_hyperoval_entries.py'],check=True)
    print(json.dumps({'completed':True,'full_native_replay_this_run':args.full,'source_pinned_full_replay_valid':True,'channels_per_cut':9216,'hyperoval_minors':120}))
if __name__=='__main__':main()
