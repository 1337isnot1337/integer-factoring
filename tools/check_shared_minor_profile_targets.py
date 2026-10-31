#!/usr/bin/env python3
"""Exact feasible targets for necessary shared-minor profile inequalities.

These are points of a relaxation. No polynomial tensor realizing either
negative-pressure target is supplied or inferred.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from charged_gaussian_floor import Dyadic
from check_binary_fixed_ray_floor import log_one_plus


F=Fraction


def prefix_caps(deficit,energy_x,energy_y,shared=True):
    q=len(energy_x)
    if len(energy_y)!=q or sorted(energy_x)!=list(energy_x) or sorted(energy_y)!=list(energy_y):
        raise ValueError('Use equal-size sorted paired energy lists.')
    if deficit<0 or energy_x[0]!=0 or energy_y[0]!=0:
        raise ValueError('Normalize each paired minimum to zero and keep nonnegative deficit.')
    caps={0:F(0),1:F(0)}
    for a in range(1,q+1):
        for b in range(1,q+1):
            k=a*b;bound=k*deficit+b*sum(energy_x[:a])+a*sum(energy_y[:b])
            caps[k]=min(caps.get(k,bound),bound)
    if shared:
        for k in range(2,q+1):
            bound=(k*deficit+min(sum(energy_x[:k]),sum(energy_y[:k])))/2
            caps[k]=min(caps.get(k,bound),bound)
    return caps


def envelope(caps):
    hull=[]
    for point in sorted(caps.items()):
        while len(hull)>1 and ((hull[-1][1]-hull[-2][1])*(point[0]-hull[-1][0])
                >=(point[1]-hull[-1][1])*(hull[-1][0]-hull[-2][0])):
            hull.pop()
        hull.append(point)
    orders=[]
    for a,b in zip(hull,hull[1:]):orders.extend([(b[1]-a[1])/(b[0]-a[0])]*(b[0]-a[0]))
    assert orders==sorted(orders)
    assert all(sum(orders[:k])<=bound for k,bound in caps.items())
    return tuple(orders),hull


def pressure(arithmetic,deficit,orders,energy_x,energy_y):
    def point(value):return arithmetic.rational(value.numerator,value.denominator)
    def partition(values,divide_root_two=False):
        total=arithmetic.zero
        for value in values:
            exponent=point(value)
            if divide_root_two:exponent=arithmetic.divide_integer(arithmetic.multiply(exponent,arithmetic.root_two),2)
            total=arithmetic.add(total,arithmetic.exp_negative(exponent))
        return log_one_plus(arithmetic,arithmetic.subtract(total,arithmetic.one))
    raw=arithmetic.add(point(deficit),arithmetic.multiply(point(F(2)),partition(orders)))
    physical=arithmetic.multiply(arithmetic.root_two,
        arithmetic.add(partition(energy_x,True),partition(energy_y,True)))
    return arithmetic.subtract(raw,physical)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();arithmetic=Dyadic(90);rows=[]
    for energies,margin in [((F(0),F(0),F(2),F(8,3)),F(1,50)),
                            ((F(0),F(0),F(2),F(8,3),F(3)),F(9,100))]:
        caps=prefix_caps(F(1),energies,energies);orders,hull=envelope(caps)
        interval=pressure(arithmetic,F(1),orders,energies,energies)
        assert interval[1]<arithmetic.rational(-margin.numerator,margin.denominator)[0]
        rows.append({'alphabet':len(energies),'deficit':'1','paired_energies':list(map(str,energies)),
            'average_normalized_orders':list(map(str,orders)),
            'prefix_caps':[[k,str(value)] for k,value in sorted(caps.items())],
            'convex_hull':[[k,str(value)] for k,value in hull],
            'pressure_dyadic_interval':list(map(str,interval)),
            'pressure_upper_bound':'-'+str(margin),
            'status':'Feasible relaxation point only; no realizing tensor.'})
    energies=(F(0),F(4));without=prefix_caps(F(0),energies,energies,shared=False)
    false_orders,_=envelope(without);correct=prefix_caps(F(0),energies,energies)
    assert sum(false_orders[:2])>correct[2]
    false_pressure=pressure(arithmetic,F(0),false_orders,energies,energies)
    assert false_pressure[1]<0
    result={'status':'Both rational larger-alphabet relaxation targets and the missing-minor negative control passed.',
        'scope':'Necessary prefix inequalities only. Negative relaxed pressure is not a crossing certificate or factoring improvement.',
        'dyadic_precision':90,'denominator':str(arithmetic.scale),'targets':rows,
        'missing_shared_minor_negative_control':{
            'alphabet':2,'average_orders':list(map(str,false_orders)),
            'violated_prefix':2,'correct_cap':str(correct[2]),
            'false_pressure_interval':list(map(str,false_pressure))}}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
