#!/usr/bin/env python3
"""An all-field MDS-family pressure floor using only the binary 28-label core."""
from collections import Counter
from fractions import Fraction
from pathlib import Path
from hashlib import sha256
import json
from coherent_hot_compounds_f4 import profiles,height
from check_hot_boundary_core_barrier import integer_at,integer_hash
from certify_natural_pressure_interval import polynomial_interval,log_interval


def physical_heights():
    Q=Counter()
    for mask in range(16):
        term=Counter({0:1})
        for _ in range(max(0,mask.bit_count()-2)):
            new=Counter(term)
            for h,c in term.items():new[h+16]+=c
            term=new
        for i in range(4):
            if not mask>>i&1:
                new=Counter()
                for h,c in term.items():
                    new[h+2*i+1]+=c;new[h+15-2*i]+=c
                term=new
        Q.update(term)
    assert sum(Q.values())==92
    return sorted(h for h,c in Q.items() for _ in range(c))


def main():
    words=[w for w in range(256) if height(w)<16];P=[{w:1} for w in words];assert len(P)==28
    result=profiles(P,P,8);assert result['zero_orders_t']==result['infinity_orders_t']
    orders=result['zero_orders_t'];heights=physical_heights();Z=Counter(orders);Q=Counter(heights);a,b=32,45
    assert Counter(h for h in heights if h<16)==Counter(height(w) for w in words)
    G={}
    for i,ci in Z.items():
        for j,cj in Q.items():
            if i or j:
                k=a*i+b*j-a;assert k>=0;G[k]=G.get(k,0)+(i-2*j)*ci*cj
    G={k:v for k,v in G.items() if v};degree=max(G);B=[0]*(degree+1)
    for i,c in G.items():
        choose=1;n=degree-i
        for j in range(n+1):
            B[i+j]+=c*choose
            if j<n:choose=choose*(n-j)//(j+1)
    signs=[(v>0)-(v<0) for v in B];nz=[x for x in signs if x];variations=sum(x!=y for x,y in zip(nz,nz[1:]))
    assert variations==3 and G[0]>0 and sum(G.values())<0
    lo=Fraction(990386623,10**9);hi=Fraction(990386624,10**9)
    gl=integer_at(G,lo.numerator,lo.denominator);gh=integer_at(G,hi.numerator,hi.denominator);assert gl<0<gh
    endpoint=len(orders)**b-len(heights)**(2*a);assert endpoint>0
    bits=192;S=1<<bits;zl=polynomial_interval([a*v for v in orders],lo,bits)[0];qu=polynomial_interval([b*v for v in heights],hi,bits)[1]
    direct=pow(zl,b)*pow(S,2*a-b)-pow(qu,2*a);assert direct>0
    gap=b*log_interval(Fraction(zl,S),bits,96)[0]-2*a*log_interval(Fraction(qu,S),bits,96)[1];assert gap>0
    names=[Path(__file__).name,'coherent_hot_compounds_f4.py','coherent_hot_compounds_f4.cpp','selected_gaussian_coefficients.cpp','native_rectangular_smith_f4.cpp','check_hot_boundary_core_barrier.py','certify_natural_pressure_interval.py']
    report={'status':'passed','rho_floor':[a,b],'core_q':28,'full_q':92,'binary_core_orders':orders,'uniform_matroid_heights':heights,'derivative_degree':degree,'derivative_coefficients':G,'mobius_signs':signs,'mobius_variations':variations,'derivative_at_zero_positive':True,'derivative_at_one_negative':True,'minimum_box':[[lo.numerator,lo.denominator],[hi.numerator,hi.denominator]],'box_signs':[-1,1],'box_sign_hashes':[integer_hash(gl),integer_hash(gh)],'endpoint_power_gap_positive':True,'endpoint_power_gap_sha256':integer_hash(endpoint),'direct_minimum_power_gap_sha256':integer_hash(direct),'precision_bits':bits,'minimum_log_gap_lower':gap,'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest() for n in names},'scope':'Any homogeneous coherent retraction retaining the 28 low coordinate labels with this 92-label height inventory; balanced natural pressure only.'}
    Path('evidence/paired-mds-low-core-floor.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','rho_floor':[a,b],'core_q':28,'derivative_degree':degree,'variations':variations,'minimum_log_gap_lower':gap}))

if __name__=='__main__':main()
