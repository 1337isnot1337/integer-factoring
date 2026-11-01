#!/usr/bin/env python3
"""Exact two-cut pressure certificate on the complete F4 q96 hyperoval profile.

This arithmetic replay does not reconstruct the large polynomial matrices.
"""
from pathlib import Path
from fractions import Fraction
from hashlib import sha256
import json
from certify_natural_pressure_interval import polynomial_interval,log_interval,power_interval


def certify(profile,rho,base,bits=192,terms=96):
    a,b=rho.numerator,rho.denominator;S=1<<bits
    zero=profile['zero_orders_t'];inf=profile['infinity_orders_t'];heights=profile['heights']
    assert len(zero)==len(inf)==len(heights)**2
    assert all(type(x)is int and x>=0 for x in zero+inf+heights) and 0 in zero and 0 in inf and 0 in heights
    ps=[polynomial_interval([a*x for x in orders],base,bits) for orders in (zero,inf)]
    q=polynomial_interval([b*h for h in heights],base,bits)
    logs=[(log_interval(Fraction(lo,S),bits,terms)[0],log_interval(Fraction(hi,S),bits,terms)[1]) for lo,hi in ps]
    lq=(log_interval(Fraction(q[0],S),bits,terms)[0],log_interval(Fraction(q[1],S),bits,terms)[1])
    gap=(4*a*lq[0]-b*sum(x[1] for x in logs),4*a*lq[1]-b*sum(x[0] for x in logs))
    return {'rho_strict_bound':[a,b],'rational_base':[base.numerator,base.denominator],'precision_bits':bits,'logarithm_terms':terms,'two_channel_sum_intervals':ps,'physical_sum_interval':q,'two_log_sum_intervals':logs,'physical_log_interval':lq,'gap_4a_logQ_minus_b_logZ0Zinf':gap,'passed':gap[0]>0,'channel_counts':[len(zero),len(inf)],'physical_labels':len(heights)}


def certify_small_type(profile,certificate,total=1<<36):
    a,b=certificate['rho_strict_bound'];base=Fraction(*certificate['rational_base']);bits=certificate['precision_bits'];terms=certificate['logarithm_terms'];S=1<<bits
    Qhi=certificate['physical_sum_interval'][1]
    counts=[total*power_interval(base,b*h,bits)[0]//Qhi for h in profile['heights']]
    deficit=total-sum(counts);assert 0<=deficit<=len(counts)
    for i in range(deficit):counts[i]+=1
    assert sum(counts)==total and min(counts)>0
    logtotal=log_interval(Fraction(total),bits,terms)
    logcounts={k:log_interval(Fraction(k),bits,terms) for k in set(counts)}
    entropy=(total*logtotal[0]-sum(k*logcounts[k][1] for k in counts),total*logtotal[1]-sum(k*logcounts[k][0] for k in counts))
    lb=log_interval(base,bits,terms);scale=(-a*lb[1],-a*lb[0]);mean=sum(k*h for k,h in zip(counts,profile['heights']))
    logs=certificate['two_log_sum_intervals']
    gap=(4*a*entropy[0]-4*b*scale[1]*mean-b*total*sum(v[1] for v in logs),4*a*entropy[1]-4*b*scale[0]*mean-b*total*sum(v[0] for v in logs))
    assert gap[0]>0
    return {'counts':counts,'denominator':total,'minimum_count':min(counts),'entropy_numerator_interval':entropy,'scale_t_dyadic_interval':scale,'mean_height_numerator':mean,'strict_gap_integer_interval':gap,'gap_denominator':b*total*S,'passed':True,'floating_point_used':False}


def main():
    path=Path('evidence/coherent-hot-f4-hyperoval-q96.json');data=json.loads(path.read_text());assert data['completed']
    for name,expected in data['source_sha256'].items():assert sha256(Path('tools',name).read_bytes()).hexdigest()==expected,name
    result=certify(data['profile'],Fraction(710005861,10**9),Fraction(2785578335,2785578336))
    assert result['passed'],result
    # Independent precision and truncation settings must also establish the sign.
    repeated=certify(data['profile'],Fraction(710005861,10**9),Fraction(2785578335,2785578336),bits=224,terms=112)
    assert repeated['passed']
    # A bound below the numerical candidate is a negative control, not evidence
    # of a lower bound. It must fail this specific rational-base certificate.
    rejected=certify(data['profile'],Fraction(710005860,10**9),Fraction(2785578335,2785578336))
    assert not rejected['passed']
    small_type=certify_small_type(data['profile'],result)
    from check_coherent_q32_pressure import counting_polynomial,multiply
    coordinate_path=Path('evidence/hot96-charge-profile.json');coordinate=json.loads(coordinate_path.read_text())
    assert all(x%2==0 for x in coordinate['orders_w'])
    assert sorted(x//4 for x in coordinate['energies_x'])==sorted(data['profile']['heights'])
    cp=counting_polynomial([x//2 for x in coordinate['orders_w']]);old_product=multiply(cp,cp)
    new_product=multiply(counting_polynomial(data['profile']['zero_orders_t']),counting_polynomial(data['profile']['infinity_orders_t']))
    diff=[(old_product[i] if i<len(old_product) else 0)-(new_product[i] if i<len(new_product) else 0) for i in range(max(len(old_product),len(new_product)))]
    first=next(i for i,v in enumerate(diff) if v);assert first==63 and diff[first]==1
    signs=[]
    for point in (Fraction(1,2),Fraction(9,10)):
        value=Fraction(0)
        for coefficient in reversed(diff):value=value*point+coefficient
        signs.append((value>0)-(value<0))
    assert signs==[1,-1]
    comparison={'coordinate_profile_sha256':sha256(coordinate_path.read_bytes()).hexdigest(),'first_difference_order':first,'first_difference_coefficient':diff[first],'difference_coefficients':diff,'exact_sign_at_one_half':signs[0],'exact_sign_at_nine_tenths':signs[1],'all_scale_dominance_claimed':False}
    report={'coordinate_comparison':comparison,'small_rational_type':small_type,'scope':__doc__,'profile_sha256':sha256(path.read_bytes()).hexdigest(),'certificate':result,'higher_precision_certificate':repeated,'negative_control':rejected,'profile_source_pins_match':True,'large_profile_recomputed':False,'factoring_exponent_improved':False,'source_sha256':{name:sha256(Path('tools',name).read_bytes()).hexdigest() for name in [Path(__file__).name,'certify_natural_pressure_interval.py']}}
    Path('evidence/coherent-hot-f4-hyperoval-pressure.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'strict_rho_bound':result['rho_strict_bound'],'gap_lower_integer':result['gap_4a_logQ_minus_b_logZ0Zinf'][0],'precision_bits':192,'type_denominator':small_type['denominator'],'minimum_type_count':small_type['minimum_count'],'passed':True,'factoring_exponent_improved':False}))

if __name__=='__main__':main()
