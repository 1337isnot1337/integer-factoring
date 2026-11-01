#!/usr/bin/env python3
"""Exact orbit normalization, all-scale mass and integration scope controls."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json

from binary_quadratic_forms import compose_forms,reduce_form
from window_orbit_duality import analyze,encode,elements,subtract

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/window-orbit-duality-controls.json'


def report():
    cases=[]
    specified=[((4,4),[(0,0),(1,0),(3,0)],[(1,0)]),
        ((4,4),[(0,0),(1,0),(3,0),(0,1),(0,3)],[(1,0),(1,1),(2,0)]),
        ((6,2),[(0,0),(1,0),(5,0),(0,1),(1,1),(5,1)],[(1,0),(2,0),(1,1)]),
        ((3,3),[(0,0),(1,0),(2,0),(0,1),(0,2)],[(1,0),(1,1)]),
        ((2,2,2),[(0,0,0),(1,0,0),(0,1,0)],[(0,0,1)])]
    for shape,A,targets in specified:
        for target in targets:cases.append(analyze(shape,A,target))
    special=cases[0]
    assert special['orbit_mean_multiplier']==4
    assert special['windows'][0]['error']==F(4,3) and special['windows'][0]['certified']
    # This coset is never in A-A, whereas every one-step target-cycle window is certified.
    bad_coset=[(x,1) for x in range(4)]
    differences={subtract(a,b,(4,4)) for a in special['A'] for b in special['A']}
    assert all(x not in differences for x in bad_coset)
    q=reduce_form((2,2,3));principal=reduce_form((1,0,5))
    assert compose_forms(q,q)==principal and q!=principal
    # 2log2<log5 and log2<log(sqrt5), both certified by the same integer inequality.
    assert 2**2<5
    ramified=dict(discriminant=-20,form=q,square=principal,
        integer_relation=[2],relation_norm_comparison='log4<log5=2R',
        integer_ball_collision=[1,-1],quotient_domain='Z/2: 1=-1 and 2=0')
    pigeonholes=[]
    for shape in ((4,4),(6,2),(3,3),(2,2,2)):
        group=elements(shape);m=len(shape)
        for L in (2,3,4):
            if len(group)<=L**m:continue
            witness=next(k for k in group if any(k) and all(
                min(F(k[i],d),1-F(k[i],d))<F(1,L) for i,d in enumerate(shape)))
            pigeonholes.append(dict(shape=shape,L=L,witness=witness))
    assert min(F(1,2),1-F(1,2))==F(1,2)
    paths=['tools/window_orbit_duality.py','tools/check_window_orbit_duality.py',
           'tools/binary_quadratic_forms.py']
    return encode(dict(scope='General ordinary orbit-restriction proof, with rational finite-group controls. No fresh PARI scan or quantitative arithmetic window bound.',
        factoring_exponent_improved=False,cases=cases,ramified_counterexample=ramified,
        conditional_pigeonhole_controls=pigeonholes,
        counts=dict(group_target_cases=len(cases),fourier_restriction_values=sum(len(c['cycle']) for c in cases),
            window_mass_integrals=sum(len(c['windows']) for c in cases),
            certified_windows=sum(w['certified'] for c in cases for w in c['windows'])),
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();data=report()
    if args.write:RECEIPT.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    else:assert json.loads(RECEIPT.read_text())==data
    print(json.dumps(dict(status='PASS',counts=data['counts'],factoring_exponent_improved=False)))
    print('Exact ramified-relation counterexample and target-cycle normalization; arithmetic W bound remains open.')


if __name__=='__main__':main()
