#!/usr/bin/env python3
"""All-scale balanced natural-pressure floor for a fixed orthogonal code core."""
from collections import Counter
from fractions import Fraction
from pathlib import Path
from hashlib import sha256
import json
from coherent_differential_kernels_f4 import differential_kernel
from orthogonal_differential_pairs_f4 import orthogonal_retraction
from coherent_hot_compounds_f4 import profiles
from check_hot_boundary_core_barrier import integer_at,integer_hash
from certify_natural_pressure_interval import polynomial_interval,log_interval


def main():
    C=[[1,1,2,2],[0,1,3,1]]
    ops=[{(1<<i)|(1<<(7-i)):row[i] for i in range(4) if row[i]} for row in C]
    low=differential_kernel(ops,8,20);full=differential_kernel(ops,8,64)
    L=orthogonal_retraction(low);LF=orthogonal_retraction(full)
    assert low['q']==43 and full['q']==92 and L is not None and LF is not None
    assert full['embedding'][:43]==low['embedding'] and LF[:43]==L
    profile=profiles(low['embedding'],L,8)
    assert profile['zero_orders_t']==profile['infinity_orders_t']
    orders=profile['zero_orders_t'];heights=full['heights'];Z=Counter(orders);Q=Counter(heights);a,b=32,45
    G={}
    for i,ci in Z.items():
        for j,cj in Q.items():
            if i:
                k=a*(i-1)+b*j;G[k]=G.get(k,0)+i*ci*cj
            if j:
                k=a*i+b*(j-1)+(b-a);G[k]=G.get(k,0)-2*j*ci*cj
    G={k:v for k,v in G.items() if v};degree=max(G);transformed=[0]*(degree+1)
    for i,c in G.items():
        choose=1;n=degree-i
        for j in range(n+1):
            transformed[i+j]+=c*choose
            if j<n:choose=choose*(n-j)//(j+1)
    signs=[(x>0)-(x<0) for x in transformed];nonzero=[x for x in signs if x]
    variations=sum(x!=y for x,y in zip(nonzero,nonzero[1:]));assert variations==2 and G[0]>0
    lo=Fraction(990294031,10**9);hi=Fraction(990294032,10**9)
    gl=integer_at(G,lo.numerator,lo.denominator);gh=integer_at(G,hi.numerator,hi.denominator)
    assert gl<0<gh
    bits=192;S=1<<bits
    zl=polynomial_interval([a*v for v in orders],lo,bits)[0]
    qu=polynomial_interval([b*v for v in heights],hi,bits)[1]
    direct=pow(zl,b)*pow(S,2*a-b)-pow(qu,2*a);assert direct>0
    gap=b*log_interval(Fraction(zl,S),bits,96)[0]-2*a*log_interval(Fraction(qu,S),bits,96)[1];assert gap>0
    assert 2*a*a>b*b and a*10**9>b*710005861
    names=[Path(__file__).name,'coherent_differential_kernels_f4.py','orthogonal_differential_pairs_f4.py','coherent_hot_compounds_f4.py','coherent_hot_compounds_f4.cpp','selected_gaussian_coefficients.cpp','native_rectangular_smith_f4.cpp','check_hot_boundary_core_barrier.py','certify_natural_pressure_interval.py']
    report={'status':'passed','coefficient_matrix':C,'core_q':43,'full_q':92,'core_maps_agree':True,'core_both_complete_profiles_recomputed':True,'full_92_profile_recomputed':False,'core_orders':orders,'full_heights':heights,'rho_floor':[a,b],'derivative_degree':degree,'derivative_coefficients':G,'mobius_coefficient_signs':signs,'mobius_sign_variations':variations,'minimum_box':[[lo.numerator,lo.denominator],[hi.numerator,hi.denominator]],'derivative_box_signs':[-1,1],'derivative_hashes':[integer_hash(gl),integer_hash(gh)],'precision_bits':bits,'minimum_box_log_gap_lower':gap,'direct_power_gap_positive':True,'direct_power_gap_sha256':integer_hash(direct),'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest() for n in names},'scope':'Balanced natural profiles for the specified retained core and height inventory. General gauges, asymmetric supports and other low cores are outside the claim.'}
    Path('evidence/orthogonal-pair-core-floor.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','rho_floor':[a,b],'core_q':43,'full_q':92,'derivative_degree':degree,'sign_variations':variations,'minimum_box_log_gap_lower':gap,'full_profile_recomputed':False}))

if __name__=='__main__':main()
