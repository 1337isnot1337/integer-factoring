#!/usr/bin/env python3
"""All-field 17/24 floor for a protected low core inside height36."""
import argparse,json
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from coherent_hot_compounds import height
from coherent_hot_f4_partial import partial_profiles
from check_hot_boundary_core_barrier import integer_at,integer_hash
from certify_natural_pressure_interval import polynomial_interval,log_interval

def strip_times(x):
    if isinstance(x,dict):return {k:strip_times(v) for k,v in x.items() if k!='seconds'}
    if isinstance(x,list):return list(map(strip_times,x))
    return x

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay-core',action='store_true');args=ap.parse_args()
    source=Path('evidence/height24-core-low-orders.json');d=json.loads(source.read_text())
    assert d['completed']
    for n,h in d['source_sha256'].items():assert sha256(Path('tools',n).read_bytes()).hexdigest()==h,n
    words=[w for w in range(4096)if height(w)<24];assert d['words']==words and len(words)==81
    p=d['profile'];assert p['q']==81 and p['modes']==12 and p['cutoff']==64
    if args.replay_core:
        P=[{w:1}for w in words]
        assert strip_times(partial_profiles(P,P,12,64))==p
    cuts=p['cuts'];orders=cuts[0]['retained_orders_t'];assert orders==cuts[1]['retained_orders_t']
    assert len(orders)==3359 and all(0<=v<=64 for v in orders)
    for c in cuts:
        assert c['tail_count']==3202 and len(c['retained_orders_t'])+c['tail_count']==81**2
        assert sorted(v for b in c['blocks'] for v in b['retained_orders_t'])==orders
    heights=[height(w)for w in range(4096)if height(w)<=36];assert len(heights)==294
    Z,Q=Counter(orders),Counter(heights);a,b=17,24;G={}
    for i,ci in Z.items():
        for j,cj in Q.items():
            if i or j:
                k=a*i+b*j-a;assert k>=0;G[k]=G.get(k,0)+(i-2*j)*ci*cj
    G={k:v for k,v in G.items()if v};degree=max(G);transformed=[0]*(degree+1)
    for i,c in G.items():
        choose=1;n=degree-i
        for j in range(n+1):
            transformed[i+j]+=c*choose
            if j<n:choose=choose*(n-j)//(j+1)
    signs=[(x>0)-(x<0)for x in transformed];nz=[x for x in signs if x]
    variations=sum(x!=y for x,y in zip(nz,nz[1:]));assert degree==1935 and variations==3 and G[0]>0 and sum(G.values())<0
    boxes=[(9081026427,9081026428),(9876550127,9876550128),(9975603842,9975603843)];den=10**10;root_controls=[]
    for i,(ln,hn) in enumerate(boxes):
        lo,hi=Fraction(ln,den),Fraction(hn,den);gl=integer_at(G,lo.numerator,lo.denominator);gh=integer_at(G,hi.numerator,hi.denominator)
        assert gl*gh<0 and (gl>0)==(i!=1)
        root_controls.append({'box':[[lo.numerator,lo.denominator],[hi.numerator,hi.denominator]],'signs':[1 if gl>0 else -1,1 if gh>0 else -1],'value_hashes':[integer_hash(gl),integer_hash(gh)]})
    lo,hi=map(lambda n:Fraction(n,den),boxes[1]);bits=192;S=1<<bits
    zl=polynomial_interval([a*x for x in orders],lo,bits)[0];qu=polynomial_interval([b*x for x in heights],hi,bits)[1]
    gap=b*log_interval(Fraction(zl,S),bits,96)[0]-2*a*log_interval(Fraction(qu,S),bits,96)[1]
    assert gap>0 and 3359**24>294**34 and 2*a*a>b*b
    out={'status':'passed','core_q':81,'maximum_physical_height':36,'maximum_height_inventory_q':294,'known_core_orders_per_cut':3359,'unknown_core_orders_per_cut':3202,'rho_floor':[17,24],'derivative_degree':degree,'mobius_sign_variations':variations,'root_controls':root_controls,'precision_bits':bits,'minimum_log_gap_lower':gap,'endpoint_integer_power_gap_positive':True,'core_receipt_sha256':sha256(source.read_bytes()).hexdigest(),'scope':'Balanced natural pressure of homogeneous coherent retractions on twelve modes retaining every coordinate state of height<24 and supported in height<=36, over any characteristic-two field. Not an arbitrary-weight or asymmetric obstruction.'}
    Path('evidence/height36-core-floor.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({**out,'core_native_replayed_this_run':args.replay_core}))
if __name__=='__main__':main()
