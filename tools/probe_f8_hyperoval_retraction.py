#!/usr/bin/env python3
"""Complete source-pinned q96 F8 hyperoval profiles for a chosen coordinate retract."""
import argparse,json,itertools,time
from hashlib import sha256
from pathlib import Path
from coherent_f8_maps import mul,height,validate_maps
from coherent_hot_compounds_f8 import partial_profiles
from search_coherent_hot_retractions import pressure

def inverse(matrix):
    n=len(matrix);a=[list(row)+[int(i==j)for j in range(n)]for i,row in enumerate(matrix)]
    for j in range(n):
        i=next(i for i in range(j,n)if a[i][j]);a[i],a[j]=a[j],a[i]
        inv=(0,1,5,6,7,2,3,4)[a[j][j]];a[j]=[mul(inv,v)for v in a[j]]
        for i in range(n):
            if i!=j:
                c=a[i][j];a[i]=[v^mul(c,w)for v,w in zip(a[i],a[j])]
    return [row[n:]for row in a]

def maps(pivots):
    if len(pivots)!=3 or len(set(pivots))!=3 or any(type(i)is not int or not 0<=i<10 for i in pivots):raise ValueError('Choose three distinct columns from0 through9.')
    columns=[(1,t,mul(t,t))for t in range(8)]+[(0,1,0),(0,0,1)]
    inv=inverse([[columns[j][i]for j in pivots]for i in range(3)]);rows=[]
    for row in inv:
        coefficients=[]
        for c in columns:
            v=0
            for a,b in zip(row,c):v^=mul(a,b)
            coefficients.append(v)
        rows.append(coefficients)
    lower=[w for w in range(4096)if height(w)<25];boundary=sorted(w for w in range(4096)if height(w)==25 and w.bit_count()==3)
    p=[{w:1}for w in lower+[31]]+[{boundary[j]:rows[i][j]for j in range(10)if rows[i][j]}for i in range(3)]
    l=[{w:1}for w in lower+[31]]+[{boundary[j]:1}for j in pivots]
    validate_maps(p,l,12)
    return p,l,rows,boundary

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pivots',required=True,help='Comma-separated triple of column numbers.')
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    pivots=list(map(int,args.pivots.split(',')));p,l,rows,boundary=maps(pivots)
    names=[Path(__file__).name,'coherent_f8_maps.py','coherent_hot_compounds_f8.py','coherent_hot_compounds_f8.cpp','selected_gaussian_coefficients.cpp','native_rectangular_smith_f4.cpp','coherent_hot_compounds.py','search_coherent_hot_retractions.py']
    fingerprints=lambda:{n:sha256(Path('tools',n).read_bytes()).hexdigest()for n in names}
    initial=fingerprints();data={'completed':False,'pivot_columns':pivots,'normalized_embedding_rows':rows,'boundary_words':boundary,'embedding':p,'retraction':l,'source_sha256':initial,'scope':'Complete finite q96 profile attempt; numerical pressure alone is not an exact bound certificate.'}
    args.output.write_text(json.dumps(data,indent=2)+'\n');start=time.monotonic()
    result=partial_profiles(p,l,12,127,progress=lambda v:print(json.dumps({'pivots':pivots,**v}),flush=True))
    assert fingerprints()==initial,'Sources changed during construction'
    assert result['all_orders_resolved']and all(c['tail_count']==0 for c in result['cuts'])
    profile={'q':96,'modes':12,'heights':result['heights'],'charges':result['charges'],'zero_orders_t':result['cuts'][0]['retained_orders_t'],'infinity_orders_t':result['cuts'][1]['retained_orders_t']}
    fit=pressure(profile);data.update(completed=True,profile=profile,construction=result,numerical_pressure=fit,seconds=time.monotonic()-start)
    args.output.write_text(json.dumps(data,indent=2)+'\n');print(json.dumps({'pivots':pivots,'completed':True,'pressure':fit,'seconds':data['seconds']}),flush=True)
if __name__=='__main__':main()
