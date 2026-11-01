#!/usr/bin/env python3
"""Exact global natural-pressure floor for every height-25 boundary retraction.

The common 92-label core supplies a partition lower bound. Descartes' rule
and an exact minimum-box inequality prove the floor 17/24 at every scale.
This does not exclude interior-code changes, general gauges or other places.
"""
import json
from pathlib import Path
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from certify_natural_pressure_interval import polynomial_interval,log_interval


def integer_at(coefficients,n,d):
    degree=max(coefficients);value=coefficients.get(degree,0);power=1
    for i in range(degree-1,-1,-1):
        power*=d;value=value*n+coefficients.get(i,0)*power
    return value


def integer_hash(value):
    return sha256(bytes([int(value<0)])+abs(value).to_bytes(max(1,(abs(value).bit_length()+7)//8),'big')).hexdigest()


def main():
    path=Path('evidence/hot-height25-common-core.json');core=json.loads(path.read_text());assert core['completed'] and core['q']==92 and len(core['orders_t'])==92**2
    for name,expected in core['source_sha256'].items():assert sha256(Path('tools',name).read_bytes()).hexdigest()==expected,name
    Z=Counter(core['orders_t']);Q=Counter(core['heights']+[25]*4);a,b=17,24
    # F=24 log Z(z^17)-34 log Q(z^24); its derivative has the
    # sign of G=Z'(z^17)Q(z^24)-2 z^7 Q'(z^24)Z(z^17).
    G={}
    for i,ci in Z.items():
        for j,cj in Q.items():
            if i:
                k=a*(i-1)+b*j;G[k]=G.get(k,0)+i*ci*cj
            if j:
                k=a*i+b*(j-1)+(b-a);G[k]=G.get(k,0)-2*j*ci*cj
    G={k:v for k,v in G.items() if v};degree=max(G);bernstein=[0]*(degree+1)
    for i,c in G.items():
        choose=1;n=degree-i
        for j in range(n+1):
            bernstein[i+j]+=c*choose
            if j<n:choose=choose*(n-j)//(j+1)
    signs=[(x>0)-(x<0) for x in bernstein];nonzero=[x for x in signs if x];variations=sum(x!=y for x,y in zip(nonzero,nonzero[1:]));assert variations==2
    assert G.get(0,0)>0
    left=Fraction(98447622,10**8);right=Fraction(98447623,10**8)
    gl=integer_at(G,left.numerator,left.denominator);gr=integer_at(G,right.numerator,right.denominator)
    assert gl<0<gr
    bits=192;S=1<<bits
    zl=polynomial_interval([a*x for x in core['orders_t']],left,bits)[0]
    qu=polynomial_interval([b*x for x in core['heights']+[25]*4],right,bits)[1]
    direct=pow(zl,b)*pow(S,2*a-b)-pow(qu,2*a);assert direct>0
    lz=log_interval(Fraction(zl,S),bits,96)[0];lq=log_interval(Fraction(qu,S),bits,96)[1]
    loggap=b*lz-2*a*lq;assert loggap>0
    assert 2*a*a>b*b # 17/24 is strictly larger than 1/sqrt(2).
    names=[Path(__file__).name,'certify_natural_pressure_interval.py']
    report={'scope':__doc__,'status':'passed','core_sha256':sha256(path.read_bytes()).hexdigest(),'rho_uniform_lower_bound':[a,b],'comparison_above_post_half':2*a*a-b*b,'derivative_degree':degree,'derivative_power_coefficients':G,'mobius_coefficient_signs':signs,'mobius_sign_variations':variations,'minimum_box':[[left.numerator,left.denominator],[right.numerator,right.denominator]],'derivative_box_signs':[-1,1],'derivative_integer_hashes':[integer_hash(gl),integer_hash(gr)],'precision_bits':bits,'minimum_box_Z_lower':zl,'minimum_box_Q_upper':qu,'direct_power_gap_positive':True,'direct_power_gap_bits':direct.bit_length(),'direct_power_gap_sha256':integer_hash(direct),'minimum_box_log_gap_lower':loggap,'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest() for n in names},'scope_excludes_interior_changes_general_gauges_other_places':True}
    Path('evidence/hot-boundary-core-barrier.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'status':'passed','rho_floor':[a,b],'descartes_variations':variations,'minimum_box_log_gap_lower':loggap,'degree':degree}))

if __name__=='__main__':main()
