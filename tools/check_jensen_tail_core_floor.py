#!/usr/bin/env python3
"""Global protected-core floors from an aggregate bound on unresolved orders."""
import json
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from coherent_hot_compounds import height,one_particle
from check_hot_boundary_core_barrier import integer_at,integer_hash
from certify_natural_pressure_interval import polynomial_interval,log_interval

def main():
    path=Path('evidence/height24-core-low-orders.json');d=json.loads(path.read_text());assert d['completed']
    for n,h in d['source_sha256'].items():assert sha256(Path('tools',n).read_bytes()).hexdigest()==h,n
    p=d['profile'];words=[w for w in range(4096)if height(w)<24];assert d['words']==words and len(words)==81
    assert p['heights']==list(map(height,words))and p['q']==81 and p['cutoff']==64
    orders=p['cuts'][0]['retained_orders_t'];assert orders==p['cuts'][1]['retained_orders_t']and len(orders)==3359
    tail=3202;assert all(c['tail_count']==tail for c in p['cuts'])
    B=one_particle(12)
    for i in range(24):
        for j in range(24):
            assert B[i][j].bit_length()-1 <= (2*(i%12)+1)+(2*(j%12)+1)
            assert B[i][j].bit_count()%2==int(i//12!=j//12 and i%12==j%12)
    determinant_degree_bound=2*sum(a+b for a in p['heights']for b in p['heights'])
    assert determinant_degree_bound==4*81*sum(p['heights'])==433836
    known_sum=sum(orders);remaining_sum_bound=determinant_degree_bound-known_sum
    assert known_sum==173663 and remaining_sum_bound==260173 and remaining_sum_bound>=65*tail
    exponent=(remaining_sum_bound+tail-1)//tail;assert exponent==82
    # Aggregate Jensen lower bound only: these are not individual tail orders.
    lower_orders=orders+[exponent]*tail;Z=Counter(lower_orders);assert sum(Z.values())==6561
    cases=[(37,17,24,[(9081026427,9081026428),(9902443571,9902443572)],2265),
           (38,29,41,[(9454699347,9454699348),(9946866891,9946866892)],3907)]
    results=[]
    for cap,a,b,boxes,expected_degree in cases:
        Q=Counter({0:1})
        for odd in range(1,cap+1,2):
            for h,c in list(Q.items()):
                if h+odd<=cap:Q[h+odd]+=c
        heights=sorted(h for h,c in Q.items()for _ in range(c));assert len(heights)==(372 if cap==37 else 409)
        G={}
        for i,ci in Z.items():
            for j,cj in Q.items():
                if i or j:
                    k=a*i+b*j-a;assert k>=0;G[k]=G.get(k,0)+(i-2*j)*ci*cj
        G={k:v for k,v in G.items()if v};degree=max(G);assert degree==expected_degree
        transformed=[0]*(degree+1)
        for i,c in G.items():
            choose=1;n=degree-i
            for j in range(n+1):
                transformed[i+j]+=c*choose
                if j<n:choose=choose*(n-j)//(j+1)
        signs=[(v>0)-(v<0)for v in transformed];nz=[v for v in signs if v];variations=sum(x!=y for x,y in zip(nz,nz[1:]));assert variations==2 and G[0]>0 and sum(G.values())>0
        controls=[]
        for i,(ln,hn)in enumerate(boxes):
            lo,hi=Fraction(ln,10**10),Fraction(hn,10**10);gl=integer_at(G,lo.numerator,lo.denominator);gh=integer_at(G,hi.numerator,hi.denominator);assert gl*gh<0 and (gl>0)==(i==0)
            controls.append({'box':[[lo.numerator,lo.denominator],[hi.numerator,hi.denominator]],'signs':[1 if gl>0 else-1,1 if gh>0 else-1],'value_hashes':[integer_hash(gl),integer_hash(gh)]})
        lo,hi=[Fraction(v,10**10)for v in boxes[1]];bits=192;S=1<<bits
        zl=polynomial_interval([a*v for v in lower_orders],lo,bits)[0];qu=polynomial_interval([b*h for h in heights],hi,bits)[1]
        gap=b*log_interval(Fraction(zl,S),bits,96)[0]-2*a*log_interval(Fraction(qu,S),bits,96)[1]
        assert gap>0 and 6561**b>len(heights)**(2*a) and 2*a*a>b*b
        results.append({'maximum_height':cap,'ambient_height_inventory_q':len(heights),'rho_floor':[a,b],'derivative_degree':degree,'sign_variations':variations,'root_controls':controls,'minimum_gap_lower':gap,'precision_bits':bits,'endpoint_integer_power_gap_positive':True})
    names=['check_jensen_tail_core_floor.py','coherent_hot_compounds.py','check_hot_boundary_core_barrier.py','certify_natural_pressure_interval.py']
    out={'completed':True,'core_q':81,'known_orders_per_cut':3359,'unknown_orders_per_cut':3202,'known_order_sum':known_sum,'core_determinant_degree_bound':determinant_degree_bound,'unknown_sum_upper_bound':remaining_sum_bound,'jensen_tail_exponent':exponent,'tail_orders_claimed_individually':False,'cases':results,'core_receipt_sha256':sha256(path.read_bytes()).hexdigest(),'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest()for n in names},'scope':'Constant height-homogeneous paired splits retaining the full height<24 coordinate core, over any characteristic-two field and with all hot modes allowed. Separate alphabets/types allowed; common balanced natural cut scale only. No deletion closure required above the core.'}
    Path('evidence/jensen-tail-core-floors.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'tail_mean_upper_bound':[remaining_sum_bound,tail],'rounded_tail_exponent':exponent,'cases':results}))
if __name__=='__main__':main()
