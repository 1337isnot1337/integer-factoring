#!/usr/bin/env python3
"""Replace a Gibbs type by explicit small-denominator counts using a KL bound."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from certify_natural_pressure_interval import certify, load_profile, power_interval


def encode(value):
    return [value.numerator,value.denominator]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile',type=Path,default=Path('evidence/hot96-charge-profile.json'))
    parser.add_argument('--certificate',type=Path,default=Path('evidence/hot96-interval-pressure.json'))
    parser.add_argument('--denominator-bits',type=int,default=36)
    parser.add_argument('--output',type=Path,default=Path('evidence/hot96-small-denominator-type.json'))
    args=parser.parse_args()
    profile=load_profile(args.profile)
    previous=json.loads(args.certificate.read_text())
    assert previous['input_sha256']==sha256(args.profile.read_bytes()).hexdigest()
    engine=Path(__file__).with_name('certify_natural_pressure_interval.py')
    assert previous['source_sha256']==sha256(engine.read_bytes()).hexdigest()
    rho,base=Fraction(*previous['rho_strict_bound']),Fraction(*previous['rational_base'])
    certificate=certify(profile['orders_w'],profile['energies_x'],rho,base,
                        bits=previous['precision_bits'],terms=previous['logarithm_terms'])
    assert certificate['passed']
    for key,value in certificate.items():
        assert previous[key]==value
    assert args.denominator_bits>=1 and profile['energies_x'][0]==0
    denominator=1<<args.denominator_bits
    a,b=certificate['support_parameter']
    qlo,qhi=certificate['polynomial_q_dyadic_interval']
    intervals=[]
    for energy in profile['energies_x']:
        lo,hi=power_interval(base,b*energy,certificate['precision_bits'])
        assert lo>0
        intervals.append((Fraction(lo,qhi),Fraction(hi,qlo)))
    counts=[0]+[(denominator*lo).numerator//(denominator*lo).denominator for lo,hi in intervals[1:]]
    counts[0]=denominator-sum(counts)
    assert min(counts)>0 and sum(counts)==denominator
    bound=sum(max(abs(Fraction(count,denominator)-lo),abs(Fraction(count,denominator)-hi))**2/lo
              for count,(lo,hi) in zip(counts,intervals))
    gap=Fraction(certificate['a_log_q_minus_b_log_p_dyadic_interval'][0],1<<certificate['precision_bits'])
    remaining=gap-a*bound
    result={'scope':'Exact type rounding: KL(p||g)<=sum((p-g)^2/g). Both paired alphabets use these positive counts.',
            'profile_sha256':sha256(args.profile.read_bytes()).hexdigest(),
            'certificate_sha256':sha256(args.certificate.read_bytes()).hexdigest(),
            'source_sha256':{Path(__file__).name:sha256(Path(__file__).read_bytes()).hexdigest(),engine.name:sha256(engine.read_bytes()).hexdigest()},
            'rho_strict_bound':encode(rho),'type_denominator':denominator,'counts_in_profile_order':counts,
            'gibbs_intervals':[[encode(lo),encode(hi)] for lo,hi in intervals],
            'chi_square_upper':encode(bound),'original_log_gap_lower':encode(gap),
            'remaining_log_gap_lower':encode(remaining),'passed':remaining>0}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':result['passed'],'type_denominator':denominator,'minimum_count':min(counts)}))
    if not result['passed']:
        raise SystemExit('The rounded type consumes the available strict margin.')


if __name__=='__main__':
    main()
