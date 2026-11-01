#!/usr/bin/env python3
"""Exact crossing-prefix and rational cold-tail certificate controls."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
from random import Random
import time

from check_polynomial_crossing_axis import brauer, profiles
from cold_gap_bound import certify


def f(value):
    return math.log1p(math.exp(-value))


def tail(values, scale):
    return math.log1p(sum(math.exp(-scale*(v-min(values))) for v in values[1:]))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    begun,rng=time.monotonic(),Random(0xC01D20261101)
    cases=[]
    for ax,ay in ((2,2),(2,3),(3,3),(3,5)):
        for kind in ('constant_crossing','nonbinary_brauer','random_strict'):
            values=[]
            if kind=='nonbinary_brauer':
                values=brauer(ax,ay,2)
            else:
                for v,u,y,x in product(range(ax),range(ay),range(ay),range(ax)):
                    low,high=(rng.getrandbits(4),rng.getrandbits(4)) if kind=='random_strict' else (0,0)
                    values.append((int(x==v and y==u)^low^(low<<1),high^(high<<1)))
            for mode in ('random_both_places','cold_moving'):
                height=([[[0]+[rng.randrange(-4,5) for _ in range(q-1)]
                           for q in (ax,ay,ay,ax)] for _ in range(2)]
                        if mode=='random_both_places' else
                        [[[0]+[(1000 if p==0 else 7)*i for i in range(1,q)]
                          for q in (ax,ay,ay,ax)] for p in range(2)])
                cases.append((ax,ay,values,height,kind,mode))
    records,certificates=[],[]
    count,checks,shared,saturated=0,0,0,0
    angles=[]
    for z in (Fraction(1,2),Fraction(1,3),Fraction(1,10),Fraction(1,100)):
        angles.append(((1-z*z)/(1+z*z),2*z/(1+z*z)))
    for ax,ay,values,height,kind,mode in cases:
        cuts=profiles(values,ax,ay,height)
        count+=2*ax*ay
        gx=sorted(sum(height[p][leg][i] for p in range(2) for leg in (0,3)) for i in range(ax))
        gy=sorted(sum(height[p][leg][i] for p in range(2) for leg in (1,2)) for i in range(ay))
        e=sorted(x+y for x in gx for y in gy)
        gap=e[1]-e[0]
        for strong in (0,1):
            r,s=cuts[strong],cuts[1-strong]
            slack=e[0]-r[0]-s[0]
            u,v=r[1]-r[0],s[1]-s[0]
            assert slack>=0 and u+v<=gap+2*slack
            shared+=1
            saturated+=int(slack==0)
            for a,b in angles:
                c=a+b
                for t in (.05,.2,1.,3.):
                    physical=tail(gx,t)+tail(gy,t)
                    margin=t*slack+float(a)*tail(r,t/float(a))+float(b)*tail(s,t/float(b))-physical
                    lower=float(c)*f(t*gap/float(c))-physical
                    assert margin+1e-9>=lower,(margin,lower)
                    checks+=1
                if gap>0:
                    # The linear exponential lower bound supplies an exact
                    # test point beyond the logarithmic threshold, at ANY
                    # height magnitude or supporting angle.
                    scaled_gap=(2*(ax+ay-2)-c)/(c-1)
                    certificate=certify(ax,ay,a,b,scaled_gap)
                    assert certificate['status']=='certified_positive_margin'
                    certificates.append(certificate)
        records.append({'kind':kind,'height_mode':mode,'alphabets':[ax,ay],
                        'heights':height,'profiles':cuts,'physical_energies':[gx,gy],
                        'first_physical_pair_gap':gap})
    # Invalid abstract spectra cannot use the theorem's crossing prefixes.
    r,s=[0,10,10,10],[0,10,10,10]
    assert (r[1]-r[0])+(s[1]-s[0])>1
    a,b=Fraction(3,5),Fraction(4,5)
    bad_margin=float(a)*tail(r,10/float(a))+float(b)*tail(s,10/float(b))-2*f(10)
    assert bad_margin<0
    inconclusive=certify(2,2,a,b,Fraction(1,100),maximum_terms=2)
    assert inconclusive['status']=='not_certified'
    names=['cold_gap_bound.py','check_cold_gap_bound.py','check_polynomial_crossing_axis.py',
           'weighted_smith_oracle_f4.py']
    report={'status':'pass','independent_untruncated_orders':count,'shared_gap_controls':shared,
            'saturated_ground_cases':saturated,'floating_lower_bound_controls':checks,
            'exact_rational_cold_certificates':len(certificates),'certificates':certificates,
            'distinct_rational_certificate_parameters':len({(r['alphabet_x'],r['alphabet_y'],
                r['a'],r['b'],r['scaled_physical_gap']) for r in certificates}),
            'cases':records,'invalid_prefix_negative_margin':bad_margin,
            'inconclusive_certificate_control':inconclusive,'seconds':time.monotonic()-begun,
            'source_sha256':{name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
            'scope':'Finite controls of the ordinary all-field proof and rational threshold '
                    'arithmetic. Numerical comparisons are not a uniform proof. No interior '
                    'or complete moving-gauge exclusion, new crossing, or score improvement.'}
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('status','independent_untruncated_orders',
                     'shared_gap_controls','floating_lower_bound_controls',
                     'exact_rational_cold_certificates','seconds')}))


if __name__=='__main__':
    main()
