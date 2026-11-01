#!/usr/bin/env python3
"""Exact q96 F8 pressure bound and all-scale comparison with the F4 code."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import accumulate
import json
from pathlib import Path
from certify_coherent_hot_f4_hyperoval import certify,certify_small_type

def main():
    path=Path('evidence/f8-hyperoval-q96-full-profile.json');p=json.loads(path.read_text())
    replay_path=Path('evidence/f8-hyperoval-q96-clean-replay.json');replay=json.loads(replay_path.read_text())
    assert p['completed']and replay['completed']and replay['tail_counts']==[0,0]
    assert replay['profile_sha256']==sha256(path.read_bytes()).hexdigest()
    for n,h in replay['source_sha256'].items():assert sha256(Path('tools',n).read_bytes()).hexdigest()==h,n
    rho=Fraction(142001167,200000000);base=Fraction(557123173,557123174)
    c=certify(p,rho,base);high=certify(p,rho,base,224,112);assert c['passed']and high['passed']
    negative=certify(p,Fraction(355002917,500000000),Fraction(1392807933,1392807934));assert negative['gap_4a_logQ_minus_b_logZ0Zinf'][1]<0
    typ=certify_small_type(p,c)
    old_path=Path('evidence/coherent-hot-f4-hyperoval-q96.json');old=json.loads(old_path.read_text());assert old['completed'];old=old['profile']
    assert Counter(p['heights'])==Counter(old['heights'])and len(p['heights'])==96
    comparisons={}
    for key in ['zero_orders_t','infinity_orders_t']:
        assert len(p[key])==len(old[key])==9216
        differences=[a-b for a,b in zip(accumulate(p[key]),accumulate(old[key]))]
        assert min(differences)>=0 and max(differences)>0
        comparisons[key]={'every_prefix_sum_at_least_old':True,'first_strict_prefix':next(i+1 for i,v in enumerate(differences)if v),'strict_prefix_count':sum(v>0 for v in differences),'total_order_gain':differences[-1]}
    assert rho<Fraction(710005861,10**9)and 2*rho*rho>1
    names=['certify_f8_hyperoval_q96.py','certify_coherent_hot_f4_hyperoval.py','certify_natural_pressure_interval.py']
    out={'completed':True,'certificate':c,'higher_precision':high,'negative_nearby_parameter_control':negative,'small_rational_type':typ,'both_natural_profiles_dominate_F4_at_all_positive_scales':comparisons,'profile_sha256':sha256(path.read_bytes()).hexdigest(),'clean_replay_sha256':sha256(replay_path.read_bytes()).hexdigest(),'old_F4_profile_sha256':sha256(old_path.read_bytes()).hexdigest(),'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest()for n in names},'large_native_profiles_recomputed_this_run':False,'factoring_exponent':'1/2; unchanged'}
    Path('evidence/f8-hyperoval-q96-pressure.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'strict_rho_bound':[rho.numerator,rho.denominator],'decimal_bound':'0.710005835','gap_lower':c['gap_4a_logQ_minus_b_logZ0Zinf'][0],'all_scale_comparison':True,'type_denominator':typ['denominator'],'factoring_exponent':'1/2; unchanged'}))
if __name__=='__main__':main()
