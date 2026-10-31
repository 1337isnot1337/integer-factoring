#!/usr/bin/env python3
"""Literal all-minor controls for the universal binary Smith pressure floor."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import random

from charged_gaussian_floor import Dyadic
from check_binary_fixed_ray_floor import add, encode_polynomials, log_one_plus, log_partition
from check_binary_two_cut_lattices import literal_profiles, pmul, endpoint
from search_binary_two_cut_lattices import labels


F=Fraction


def point(arithmetic,value):
    value=F(value)
    return arithmetic.rational(value.numerator,value.denominator)


def general_floor(arithmetic,A,X,Y):
    a,x,y=map(lambda value:point(arithmetic,value),(A,X,Y))
    minimum=point(arithmetic,min(X,Y))
    first=arithmetic.exp_negative(arithmetic.add(a,minimum))
    last_exponent=arithmetic.divide_integer(arithmetic.subtract(arithmetic.add(
        arithmetic.multiply(point(arithmetic,3),a),
        arithmetic.multiply(point(arithmetic,4),arithmetic.add(x,y))),minimum),2)
    last=arithmetic.multiply(point(arithmetic,2),arithmetic.exp_negative(last_exponent))
    positive=arithmetic.add(a,arithmetic.multiply(point(arithmetic,2),
        log_one_plus(arithmetic,arithmetic.add(first,last))))
    physical=arithmetic.multiply(arithmetic.root_two,arithmetic.add(
        arithmetic.softplus(arithmetic.multiply(arithmetic.root_two,x)),
        arithmetic.softplus(arithmetic.multiply(arithmetic.root_two,y))))
    return arithmetic.subtract(positive,physical)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    arithmetic=Dyadic(80);rng=random.Random(311721)
    records=[];profiles=endpoints=pressure_comparisons=flips=0
    for degree in range(7):
        for trial in range(12):
            values=[]
            for index in range(16):
                entry={power:value for power in range(1,degree+1) if (value:=rng.randrange(4))}
                x,y,u,v=labels(index)
                constant=int(x==v and y==u)
                for value in entry.values():constant^=value
                if constant:entry[0]=constant
                values.append(entry)
            heights=[rng.randrange(-4,5) for _ in range(8)]
            packed=encode_polynomials(values);endpoint(packed);endpoints+=16
            original=literal_profiles(packed,heights=heights);profiles+=2
            total=[heights[i]+heights[i+4] for i in range(4)]
            flip_mask=0
            for i,j in ((0,3),(1,2)):
                if total[i]+total[j]>0:flip_mask|=(1<<i)|(1<<j)
            constants=[sum(heights[4*side+i] for i in range(4) if (flip_mask>>i)&1)
                       for side in range(2)]
            canonical_heights=[-value if (flip_mask>>(i%4))&1 else value
                               for i,value in enumerate(heights)]
            canonical=[values[index^flip_mask] for index in range(16)]
            infinity,zero=literal_profiles(encode_polynomials(canonical),heights=canonical_heights)
            profiles+=2
            assert infinity==tuple(value-constants[0] for value in original[0])
            assert zero==tuple(value-constants[1] for value in original[1])
            flips+=bool(flip_mask)
            total=[canonical_heights[i]+canonical_heights[i+4] for i in range(4)]
            a=-F(total[0]+total[3],2);b=-F(total[1]+total[2],2)
            assert a>=0 and b>=0
            clearing=rng.randrange(degree+4)
            infinity=tuple(value+clearing for value in infinity)
            zero=tuple(value-clearing for value in zero)
            laurent=[{power-clearing:value for power,value in entry.items()} for entry in canonical]
            minimum0,minimum_infinity=zero[0],infinity[0]
            deficit=sum(total)-minimum0-minimum_infinity
            assert deficit>=0
            normalized=[]
            for index in range(16):
                bits=labels(index)
                shift=-minimum0+sum(canonical_heights[4+i]*(1-bits[i]) for i in range(4))
                entry={power+shift:value for power,value in laurent[index^15].items()}
                bound=deficit-sum(total[i]*bits[i] for i in range(4))
                assert all(0<=power<=bound for power in entry)
                normalized.append(entry)
            packed=encode_polynomials(normalized);endpoint(packed);endpoints+=16
            other,first=literal_profiles(packed,heights=tuple(-value for value in total)+(0,)*4)
            second=tuple(value+deficit for value in other);profiles+=2
            assert first==tuple(value-minimum0 for value in zero)
            assert second==tuple(value-minimum_infinity for value in infinity)
            assert first[0]==second[0]==0
            assert sum(first)<=4*deficit+4*(a+b)
            assert sum(second)<=4*deficit+4*(a+b)
            assert first[1]+second[1]<=2*deficit+2*min(a,b)
            supports=[]
            for indices,pair_energy in (((0,1,8,9),a),((0,2,4,6),b)):
                i,j,k,l=indices
                delta=add(pmul(normalized[i],normalized[l]),pmul(normalized[j],normalized[k]))
                assert delta and max(delta)<=2*deficit+2*pair_energy
                color=0
                for value in delta.values():color^=value
                assert color==1
                assert first[1]+second[1]<=2*deficit+2*pair_energy+min(delta)-max(delta)
                supports.append([min(delta),max(delta)])
            for scale in (F(1,10),F(1,2),F(1)):
                t=point(arithmetic,scale)
                actual=arithmetic.add(log_partition(arithmetic,zero,t),log_partition(arithmetic,infinity,t))
                physical=arithmetic.add(point(arithmetic,2*(a+b)*scale),
                    arithmetic.multiply(arithmetic.root_two,arithmetic.add(
                        arithmetic.softplus(arithmetic.multiply(arithmetic.root_two,point(arithmetic,a*scale))),
                        arithmetic.softplus(arithmetic.multiply(arithmetic.root_two,point(arithmetic,b*scale))))))
                actual=arithmetic.subtract(actual,physical)
                bound=general_floor(arithmetic,deficit*scale,a*scale,b*scale)
                assert actual[0]>=bound[1]
                assert actual[0]>0
                pressure_comparisons+=1
            records.append({'degree':degree,'original_heights_infinity_then_zero':heights,
                'paired_label_flip_mask':flip_mask,'Laurent_clearing':clearing,
                'canonical_paired_energy_halves':[str(a),str(b)],'deficit':deficit,
                'normalized_zero_profile':first,'normalized_infinity_profile':second,
                'shared_minor_supports':supports})
    result={'status':'Every arbitrary-height normalization and shared-minor pressure control passed.',
        'scope':'All binary Laurent-polynomial strict crossing seeds, every coefficient field, and all independent physical-letter heights; no improved factoring exponent.',
        'literal_endpoint_coefficients':endpoints,'complete_profiles_by_all_minors':profiles,
        'arbitrary_integer_height_cases':len(records),'nontrivial_paired_label_flips':flips,
        'outward_actual_pressure_comparisons':pressure_comparisons,
        'compact_analytic_certificate':'binary-fixed-ray-floor.json',
        'critical_lower_bound':'exp(-max(a,b)*t)/28, with paired energies normalized to -2a and -2b',
        'literal_controls':records,
        'source_sha256':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in (Path(__file__).name,'check_binary_fixed_ray_floor.py',
                         'check_binary_two_cut_lattices.py','charged_gaussian_floor.py')}}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key not in ('literal_controls','source_sha256')},indent=2))


if __name__=='__main__':main()
