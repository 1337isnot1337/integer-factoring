#!/usr/bin/env python3
"""All-seed binary fixed-ray pressure floor: exact compact and tensor controls."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import factorial
from pathlib import Path
import random

from charged_gaussian_floor import Dyadic
from check_global_homogeneous_boolean_filter_floor import construct_cover,audit_cover
from check_binary_two_cut_lattices import literal_profiles,polynomial,pmul,endpoint
from search_binary_two_cut_lattices import labels


F=Fraction


def interval(arithmetic,lo,hi):
    return (arithmetic.rational(lo.numerator,lo.denominator)[0],
            arithmetic.rational(hi.numerator,hi.denominator)[1])


def log_one_plus(arithmetic,value):
    # Monotonic endpoint evaluation avoids dependency inflation in the
    # atanh-series transform when the argument interval is initially wide.
    return (arithmetic.log_one_plus((value[0],value[0]))[0],
            arithmetic.log_one_plus((value[1],value[1]))[1])


def floor_interval(arithmetic,rectangle,keep_last_pair=True):
    lo_a,hi_a,lo_t,hi_t=rectangle
    a=interval(arithmetic,lo_a,hi_a);t=interval(arithmetic,lo_t,hi_t)
    first=arithmetic.exp_negative(arithmetic.add(a,t))
    argument=first
    if keep_last_pair:
        second=arithmetic.exp_negative(arithmetic.divide_integer(arithmetic.add(
            arithmetic.multiply(arithmetic.rational(3),a),
            arithmetic.multiply(arithmetic.rational(7),t)),2))
        argument=arithmetic.add(first,arithmetic.multiply(arithmetic.rational(2),second))
    positive=arithmetic.add(a,arithmetic.multiply(arithmetic.rational(2),log_one_plus(arithmetic,argument)))
    negative=arithmetic.multiply(arithmetic.multiply(arithmetic.rational(2),arithmetic.root_two),
                               arithmetic.softplus(arithmetic.multiply(arithmetic.root_two,t)))
    return arithmetic.subtract(positive,negative)


def split_box(rectangle):
    la,ha,lt,ht=rectangle
    if ha-la>=ht-lt:
        middle=(la+ha)/2
        return ('0',(la,middle,lt,ht)),('1',(middle,ha,lt,ht))
    middle=(lt+ht)/2
    return ('2',(la,ha,lt,middle)),('3',(la,ha,middle,ht))


def add(a,b):
    out=dict(a)
    for exponent,value in b.items():
        out[exponent]=out.get(exponent,0)^value
        if not out[exponent]:del out[exponent]
    return out


def encode_polynomials(values):
    return [(sum((value&1)<<power for power,value in entry.items()),
             sum((value>>1)<<power for power,value in entry.items())) for entry in values]


def log_partition(arithmetic,orders,t):
    minimum=min(orders)
    total=arithmetic.zero
    for value in orders:
        total=arithmetic.add(total,arithmetic.exp_negative(
            arithmetic.multiply(arithmetic.rational(value-minimum),t)))
    return arithmetic.subtract(log_one_plus(arithmetic,arithmetic.subtract(total,arithmetic.one)),
                               arithmetic.multiply(arithmetic.rational(minimum),t))


def tensor_controls(arithmetic):
    rng=random.Random(311709)
    profiles=endpoints=normalizations=pressure_checks=0
    records=[]
    for degree in range(9):
        for trial in range(8):
            values=[]
            for index in range(16):
                entry={power:value for power in range(1,degree+1) if (value:=rng.randrange(4))}
                x,y,u,v=labels(index)
                constant=int(x==v and y==u)
                for value in entry.values():constant^=value
                if constant:entry[0]=constant
                values.append(entry)
            packed=encode_polynomials(values);endpoint(packed);endpoints+=16
            infinity,zero=literal_profiles(packed)
            clearing=rng.randrange(degree+3)
            infinity=tuple(value+clearing for value in infinity)
            zero=tuple(value-clearing for value in zero)
            laurent=[{power-clearing:value for power,value in entry.items()} for entry in values]
            delta=add(pmul(laurent[6],laurent[15]),pmul(laurent[7],laurent[14]))
            assert delta and min(delta)<=max(delta)
            color=0
            for value in delta.values():color^=value
            assert color==1
            assert sum(infinity[:2])+sum(zero[:2])<=-6+min(delta)-max(delta)<=-6
            a,b=zero[0],infinity[0];d=-4-a-b
            assert d>=0
            normalized=[]
            for index in range(16):
                x,y,u,v=labels(index)
                entry={power-a-2+x+y:value for power,value in laurent[index^15].items()}
                assert all(0<=power<=d+x+y+u+v for power in entry)
                normalized.append(entry)
            packed_normalized=encode_polynomials(normalized);endpoint(packed_normalized);endpoints+=16
            other,first=literal_profiles(packed_normalized,heights=(1,1,1,1,0,0,0,0))
            second=tuple(value+d for value in other)
            assert first==tuple(value-a for value in zero)
            assert second==tuple(value-b for value in infinity)
            assert first[0]==second[0]==0
            assert sum(first)<=4*d+8 and sum(second)<=4*d+8
            assert first[1]+second[1]<=2*d+2
            for scale in (F(1,10),F(1,2),F(1)):
                t=arithmetic.rational(scale.numerator,scale.denominator)
                actual=arithmetic.add(log_partition(arithmetic,zero,t),log_partition(arithmetic,infinity,t))
                physical=arithmetic.add(arithmetic.multiply(arithmetic.rational(4),t),
                    arithmetic.multiply(arithmetic.multiply(arithmetic.rational(2),arithmetic.root_two),
                        arithmetic.softplus(arithmetic.multiply(arithmetic.root_two,t))))
                actual=arithmetic.subtract(actual,physical)
                lower=floor_interval(arithmetic,(d*scale,d*scale,scale,scale))
                assert actual[0]>=lower[1]
                pressure_checks+=1
            profiles+=4;normalizations+=1
            records.append({'degree':degree,'common_clearing':clearing,'deficit':d,
                            'infinity_profile':infinity,'zero_profile':zero,
                            'shared_minor_support':[min(delta),max(delta)]})
    return {'literal_endpoint_coefficients':endpoints,'all_minor_profiles':profiles,
            'complete_normalizations':normalizations,'outward_pressure_comparisons':pressure_checks,
            'literal_tensor_controls':records}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();arithmetic=Dyadic(80)
    function=lambda rectangle:arithmetic.subtract(floor_interval(arithmetic,rectangle),arithmetic.rational(1,16))
    compact=construct_cover(function,(F(0),F(3,2),F(0),F(3,2)),split_box)
    nodes=audit_cover(function,compact,split_box)
    for mutation in (dict(compact,cover=compact['cover'][1:]),
                     dict(compact,cover=[[path,str(int(value)+(index==0))]
                                        for index,(path,value) in enumerate(compact['cover'])])):
        try:audit_cover(function,mutation,split_box)
        except AssertionError:pass
        else:raise AssertionError('Malformed cover accepted.')
    assert F(7,5)**2<2<F(10,7)**2
    assert F(5,3)**5<16
    assert sum(F(3)**k/factorial(k) for k in range(9))>16
    assert 2*(F(7,8)-F(6,7))==F(1,28)
    assert 1-F(1,4)-F(1,2)==F(1,4)
    omitted=floor_interval(arithmetic,(F(6,5),F(6,5),F(0),F(0)),keep_last_pair=False)
    assert omitted[1]<0
    controls=tensor_controls(arithmetic)
    result={'status':'Universal binary fixed-ray floor certificate and complete tensor controls passed.',
        'scope':'Every Laurent-polynomial crossing seed over any field at the specified fixed physical height ray; no improved factoring exponent.',
        'infinity_heights':[0,0,-1,-1],'zero_heights':[-1,-1,0,0],
        'critical_pressure_lower_bound':'exp(-t)/28 for every t>0',
        'dyadic_precision':80,'common_denominator':str(arithmetic.scale),
        'compact_margin':'1/16','compact_certificate':compact,'audited_nodes':nodes,
        'malformed_certificate_negative_controls':2,
        'missing_rank_channels_negative_control':[str(x) for x in omitted],**controls,
        'source_sha256':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in [Path(__file__).name,'charged_gaussian_floor.py',
                                       'check_global_homogeneous_boolean_filter_floor.py',
                                       'check_binary_two_cut_lattices.py']}}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in
        ['compact_certificate','literal_tensor_controls','source_sha256']},indent=2))
    print('compact_leaves',compact['positive_leaves'],flush=True)


if __name__=='__main__':main()
