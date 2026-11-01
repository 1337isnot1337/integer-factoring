#!/usr/bin/env python3
"""Literal independent-height controls for graded original/adjugate box caps."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from random import Random
import time

from box_adjugate_word_caps import inequalities as bare_adjugate_caps
from check_all_field_pascal_crossings import add, multiply, smith_orders, weighted_orders
from check_bare_seven_to_twelve_letters import box_weights
from check_box_adjugate_word_caps import crossed, reverse
from check_higher_jet_word_caps import detpoly
from graded_box_minor_caps import original, adjugate
from higher_jet_word_caps import word_cost_profile

ROOT = Path(__file__).resolve().parents[1]


def shift(value,amount):
    value = list(value)
    if amount < 0:
        assert not any(value[:-amount])
        value = value[-amount:]
    else:
        value = [0]*amount+value
    while value and not value[-1]:
        value.pop()
    return tuple(value)


def source_cases(records):
    result = []
    for record in records:
        result.append(record)
        result.extend(source_cases(record['restricted_controls']))
    return result


def curve(case,h0,hi):
    A,q,D,p = case['source_polynomials'],case['alphabet'],case['degree_bound'],case['prime']
    n = q*q
    letters = [(i//q,i % q,j % q,j//q) for i in range(n) for j in range(n)]
    def charge(h,word):
        return sum(h[leg][symbol] for leg,symbol in enumerate(word))
    minimum0 = min(e+charge(h0,word) for word in letters
                   for e,v in enumerate(A[word[0]*q+word[1]][word[3]*q+word[2]]) if v)
    minimumi = min(-e+charge(hi,word) for word in letters
                   for e,v in enumerate(A[word[0]*q+word[1]][word[3]*q+word[2]]) if v)
    zero,infinity = [],[]
    Ar = reverse(A,D)
    for i in range(n):
        zr,ir = [],[]
        for j in range(n):
            word = (i//q,i % q,j % q,j//q)
            zr.append(shift(A[i][j],charge(h0,word)-minimum0))
            ir.append(shift(Ar[i][j],charge(hi,word)-D-minimumi))
        zero.append(zr)
        infinity.append(ir)
    r,_ = smith_orders(zero,p)
    s,_ = smith_orders(crossed(infinity,q),p)
    assert r[0] == s[0] == 0
    H = [[a+b for a,b in zip(x,y)] for x,y in zip(h0,hi)]
    GX = [H[0][i]+H[3][i] for i in range(q)]
    GY = [H[1][i]+H[2][i] for i in range(q)]
    bx,by = min(range(q),key=GX.__getitem__),min(range(q),key=GY.__getitem__)
    X,Y = [v-GX[bx] for v in GX],[v-GY[by] for v in GY]
    d = GX[bx]+GY[by]-minimum0-minimumi
    assert d >= 0
    plain = original(d,r,s,X,Y)
    composed = adjugate(d,r,s,X,Y)
    assert min(c['slack'] for c in plain) >= 0
    assert min(c['slack'] for c in composed) >= 0
    if not any(X+Y):
        assert composed == bare_adjugate_caps(d,r,s)
    # Literal endpoint box determinants are computed before gauges; their
    # monomial charges at each place are restored separately.
    boxes = []
    for alpha in range(1,q+1):
        for beta in range(1,q+1):
            xs = sorted(range(q),key=lambda i:X[i])[:alpha]
            ys = sorted(range(q),key=lambda i:Y[i])[:beta]
            N = alpha*beta
            weights = box_weights(n,N,min(alpha,beta)**2)
            for kind in ('ordinary','crossed'):
                if kind == 'ordinary':
                    matrix = [[A[x*q+y][v*q+u] for v in xs for u in ys] for x in xs for y in ys]
                    lower0,loweri = sum(r[:N]),sum(a*b for a,b in zip(weights,s))
                else:
                    matrix = [[A[x*q+y][v*q+u] for v in xs for y in ys] for x in xs for u in ys]
                    lower0,loweri = sum(a*b for a,b in zip(weights,r)),sum(s[:N])
                determinant = detpoly(matrix,N*D+1,p)
                nonzero = [i for i,v in enumerate(determinant) if v]
                assert nonzero and sum(determinant) % p == 1
                zcharge = beta*sum(h0[0][i]+h0[3][i] for i in xs)+alpha*sum(h0[1][i]+h0[2][i] for i in ys)
                icharge = beta*sum(hi[0][i]+hi[3][i] for i in xs)+alpha*sum(hi[1][i]+hi[2][i] for i in ys)
                order0,orderi = min(nonzero)+zcharge,-max(nonzero)+icharge
                assert order0 >= lower0+N*minimum0 and orderi >= loweri+N*minimumi
                physical = beta*sum(GX[i] for i in xs)+alpha*sum(GY[i] for i in ys)
                assert order0+orderi <= physical
                boxes.append(dict(alpha=alpha,beta=beta,cut=kind,x_letters=xs,y_letters=ys,
                                  determinant=determinant,raw_orders=[order0,orderi],physical_budget=physical))
    # Acquire the actual graded adjugate from the saved complete ungauged
    # cofactors, using exact row/column monomial charges; verify its identity.
    rowhi = [hi[0][x]+hi[1][y] for x in range(q) for y in range(q)]
    colhi = [hi[3][v]+hi[2][u] for v in range(q) for u in range(q)]
    coefficient = (n-1)*(-D-minimumi)+sum(rowhi)+sum(colhi)
    B0 = case['complete_adjugate_polynomials']
    B = [[shift(B0[i][j],coefficient-colhi[i]-rowhi[j]) for j in range(n)] for i in range(n)]
    determinant = shift(detpoly(Ar,n*D+1,p),n*(-D-minimumi)+sum(rowhi)+sum(colhi))
    for i in range(n):
        for j in range(n):
            value = ()
            for k in range(n):
                value = add(value,multiply(infinity[i][k],B[k][j],p),p)
            assert not add(value,determinant if i == j else (),p,-1)
    eta,_ = smith_orders(crossed(B,q),p)
    c = [[H[leg][i]-H[leg][bx if leg in (0,3) else by] for i in range(q)] for leg in range(4)]
    row = [-c[3][v]-c[2][u] for v in range(q) for u in range(q)]
    column = [-c[0][x]-c[1][y] for x in range(q) for y in range(q)]
    theta,_ = weighted_orders(B,row,column,p,infinity=True)
    C = (n-1)*d+q*(sum(X)+sum(Y))
    assert list(theta) == [sum(r)-v-C for v in r[::-1]]
    assert all(a >= b for a,b in zip(eta,word_cost_profile(q,s)))
    dB = -max(X)-max(Y)-eta[0]-theta[0]
    actual_adj = original(dB,[v-eta[0] for v in eta],[v-theta[0] for v in theta],
                          [max(X)-v for v in X],[max(Y)-v for v in Y])
    assert min(c['slack'] for c in actual_adj) >= 0
    return dict(prime=p,alphabet=q,base_label=case['label'],zero_heights=h0,infinity_heights=hi,
                raw_profile_minima=[minimum0,minimumi],normalized_deficit=d,energy_x=X,energy_y=Y,
                source_profiles=[r,s],original_checks=len(plain),original_minimum_slack=min(c['slack'] for c in plain),
                adjugate_checks=len(composed),adjugate_minimum_slack=min(c['slack'] for c in composed),
                complete_adjugate_polynomials=B,adjugate_profiles=[eta,theta],
                actual_adjugate_box_checks=len(actual_adj),literal_box_determinants=boxes)


def run():
    path = ROOT/'evidence/box-adjugate-word-controls.json'
    cases = source_cases(json.loads(path.read_text())['literal_controls'])
    randomizer,records = Random(2026110118),[]
    for case in cases:
        q = case['alphabet']
        height_sets = [([[0]*q for _ in range(4)],[[0]*q for _ in range(4)]),
                       ([[randomizer.randrange(-2,3) for _ in range(q)] for _ in range(4)],
                        [[randomizer.randrange(-2,3) for _ in range(q)] for _ in range(4)]),
                       ([[((3*leg+2*i) % 5)-2 for i in range(q)] for leg in range(4)],
                        [[((leg+3*i) % 5)-2 for i in range(q)] for leg in range(4)])]
        for h0,hi in height_sets:
            records.append(curve(case,h0,hi))
    names = ['tools/check_graded_box_minor_caps.py','tools/graded_box_minor_caps.py',
             'tools/box_adjugate_word_caps.py','tools/check_box_adjugate_word_caps.py',
             'tools/check_bare_seven_to_twelve_letters.py','tools/check_all_field_pascal_crossings.py',
             'tools/check_higher_jet_word_caps.py','tools/higher_jet_word_caps.py']
    return dict(scope='All-field scalar graded box/adjugate compatibility; exact finite Laurent controls, no global rate floor.',
                curves=records,curve_count=len(records),
                literal_box_determinants=sum(len(r['literal_box_determinants']) for r in records),
                complete_orders=sum(4*r['alphabet']**2 for r in records),
                original_checks=sum(r['original_checks'] for r in records),
                composed_checks=sum(r['adjugate_checks'] for r in records),
                input_sha256=sha256(path.read_bytes()).hexdigest(),
                source_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in names})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    begun = time.monotonic()
    result = json.loads(json.dumps(run(),default=str))
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    else:
        assert result == json.loads((ROOT/'evidence/graded-box-minor-controls.json').read_text())
    print(json.dumps({k:v for k,v in result.items() if k not in ('curves','source_sha256')}))
    print(json.dumps({'seconds':time.monotonic()-begun}))


if __name__ == '__main__':
    main()
