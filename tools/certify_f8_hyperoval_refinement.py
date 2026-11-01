#!/usr/bin/env python3
"""Exact pressure and all-scale comparison for the refined F8 coordinate retract."""
import json
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import accumulate
from pathlib import Path
from certify_coherent_hot_f4_hyperoval import certify,certify_small_type
from coherent_f8_maps import mul

def main():
    paths={tag:Path('evidence/f8-hyperoval-retraction-'+tag+'.json')for tag in('678','789')};data={k:json.loads(p.read_text())for k,p in paths.items()}
    for tag,d in data.items():
        assert d['completed']and d['construction']['all_orders_resolved']and all(c['tail_count']==0 for c in d['construction']['cuts'])
        for n,h in d['source_sha256'].items():assert sha256(Path('tools',n).read_bytes()).hexdigest()==h,n
        for key,c in zip(['zero_orders_t','infinity_orders_t'],d['construction']['cuts']):assert d['profile'][key]==c['retained_orders_t']==sorted(c['retained_orders_t'])and len(c['retained_orders_t'])==9216
        entry=json.loads(Path('evidence/f8-hyperoval-retraction-'+tag+'-entry-controls.json').read_text());assert entry['completed']and entry['full_polynomial_entries']==166 and entry['profile_sha256']==sha256(paths[tag].read_bytes()).hexdigest()
        for n,h in entry['source_sha256'].items():assert sha256(Path('tools',n).read_bytes()).hexdigest()==h,n
    p=data['789']['profile'];rho=Fraction(355002893599,500000000000);base=Fraction(1392830724074,1392830724075)
    certificate=certify(p,rho,base);higher=certify(p,rho,base,224,112);assert certificate['passed']and higher['passed'];typed=certify_small_type(p,certificate)
    negative=certify(p,Fraction(710005787197,10**12),Fraction(2785661448145,2785661448146));assert negative['gap_4a_logQ_minus_b_logZ0Zinf'][1]<0
    old_path=Path('evidence/f8-hyperoval-q96-full-profile.json');old=json.loads(old_path.read_text());comparisons={}
    for name,reference in [('previous_F8',old),('alternative_678',data['678']['profile'])]:
        assert Counter(reference['heights'])==Counter(p['heights']);cuts={}
        for key in ['zero_orders_t','infinity_orders_t']:
            assert reference[key]==sorted(reference[key])
            diff=[a-b for a,b in zip(accumulate(p[key]),accumulate(reference[key]))];assert min(diff)>=0 and max(diff)>0
            cuts[key]={'first_strict_prefix':next(i+1 for i,v in enumerate(diff)if v),'strict_prefix_count':sum(v>0 for v in diff),'total_order_gain':diff[-1],'every_prefix_at_least_reference':True}
        comparisons[name]=cuts
    assert rho<Fraction(142001167,200000000)and 2*rho*rho>1
    assert data['789']['embedding'][:93]==old['embedding'][:93]and data['789']['retraction'][:93]==old['retraction'][:93]
    G0=[[c[i]for c in old['columns']]for i in range(3)];G1=data['789']['normalized_embedding_rows'];eta=7;eta2=mul(eta,eta)
    assert G1[0]==G0[0]
    assert G1[1]==[b^mul(eta,a)for a,b in zip(G0[0],G0[1])]
    assert G1[2]==[b^mul(eta2,a)for a,b in zip(G0[0],G0[2])]
    phi={0:1,7:1,8:eta,9:eta2}
    for row in G0:
        value=0
        for j,c in phi.items():value^=mul(c,row[j])
        assert value==0
    for r in range(10):
        for c in range(10):
            before=0;after=0
            for i,j in enumerate([0,8,9]):
                if c==j:before^=G0[i][r]
            for i,j in enumerate([7,8,9]):
                if c==j:after^=G1[i][r]
            assert before^after==mul(G0[0][r],phi.get(c,0))
    rank_one={'entry_checks':100,'field_shift':eta,'shift_square':eta2,'annihilator_covector_on_boundary':phi,'same_embedding_image':True}

    names=['certify_f8_hyperoval_refinement.py','certify_coherent_hot_f4_hyperoval.py','certify_natural_pressure_interval.py','coherent_f8_maps.py']
    out={'completed':True,'pivot_columns':[7,8,9],'certificate':certificate,'higher_precision_certificate':higher,'negative_nearby_parameter_control':negative,'small_rational_type':typed,'all_scale_natural_profile_comparisons':comparisons,'rank_one_projector_update':rank_one,'profile_sha256':{k:sha256(p.read_bytes()).hexdigest()for k,p in paths.items()},'previous_F8_profile_sha256':sha256(old_path.read_bytes()).hexdigest(),'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest()for n in names},'full_native_profiles_recomputed_this_run':False,'global_retraction_optimum_claimed':False,'factoring_exponent':'1/2; unchanged'}
    assert out['certificate']['passed']
    Path('evidence/f8-hyperoval-retraction-refinement-pressure.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'completed':True,'strict_rho_bound':[rho.numerator,rho.denominator],'decimal_bound':'0.710005787198','gap_lower':certificate['gap_4a_logQ_minus_b_logZ0Zinf'][0],'type_denominator':typed['denominator'],'minimum_type_count':typed['minimum_count'],'factoring_exponent':'1/2; unchanged'}))
if __name__=='__main__':main()
