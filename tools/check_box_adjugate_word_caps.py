#!/usr/bin/env python3
"""Complete polynomial controls and exact target rejection for adjugate boxes."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from random import Random
import time

from box_adjugate_word_caps import inequalities, first_jet_linear_caps
from check_all_field_pascal_crossings import add, multiply, smith_orders
from check_bare_seven_to_twelve_letters import constraints, satisfies
from check_higher_jet_word_caps import adjpoly, detpoly
from check_negation_odd_qutrit import cell
from higher_jet_word_caps import word_cost_profile

ROOT = Path(__file__).resolve().parents[1]


def crossed(A,q):
    return [[A[x*q+y][v*q+u] for v in range(q) for y in range(q)]
            for x in range(q) for u in range(q)]


def reverse(A,degree):
    return [[tuple((value[i] if i < len(value) else 0) for i in range(degree,-1,-1))
             for value in row] for row in A]


def profiles(A,q,degree,prime):
    first,_ = smith_orders(A,prime)
    second,_ = smith_orders(crossed(reverse(A,degree),q),prime)
    return first,second


def check_curve(A,q,degree,prime,label):
    n = q*q
    assert all(sum(A[i][j]) % prime == int(i == j) for i in range(n) for j in range(n))
    r,s = profiles(A,q,degree,prime)
    caps = inequalities(degree,r,s)
    assert min(cap['slack'] for cap in caps) >= 0
    C = (n-1)*degree
    T = reverse(A,degree)
    B = adjpoly(T,C+1,prime)  # Full degree bound, without truncating any cofactor.
    determinant = tuple(detpoly(T,n*degree+1,prime))
    for i in range(n):
        for j in range(n):
            value = ()
            for k in range(n):
                value = add(value,multiply(tuple(T[i][k]),tuple(B[k][j]),prime),prime)
            expected = determinant if i == j else ()
            assert not add(value,expected,prime,-1)
    assert all(sum(B[i][j]) % prime == int(i == j) for i in range(n) for j in range(n))
    eta,_ = smith_orders(crossed(B,q),prime)
    infinity,_ = smith_orders(reverse(B,C),prime)
    expected = tuple(sum(r)-v for v in r[::-1])
    assert infinity == expected
    lower = word_cost_profile(q,s)
    assert all(a >= b for a,b in zip(eta,lower))
    assert not satisfies(q,list(infinity),list(eta),C)
    restricted = []
    for m in range(2,q):
        indices = [x*q+y for x in range(m) for y in range(m)]
        small = [[A[i][j] for j in indices] for i in indices]
        sr,ss = profiles(small,m,degree,prime)
        assert all(a >= b for a,b in zip(sr,r)) and all(a >= b for a,b in zip(ss,s))
        restricted.append(check_curve(small,m,degree,prime,label+' / restricted '+str(m)))
    return dict(label=label,prime=prime,alphabet=q,degree_bound=degree,source_polynomials=A,
                source_profiles=[r,s],complete_adjugate_polynomials=B,
                adjugate_profiles_zero_crossed_and_infinity=[eta,infinity],word_lower=lower,
                composed_checks=len(caps),minimum_slack=min(cap['slack'] for cap in caps),
                literal_adjugate_box_checks=len(constraints(q)),restricted_controls=restricted)


def controls():
    randomizer = Random(2026110117)
    records = []
    for prime in (2,3,5):
        for degree in (1,2):
            n = 4
            A = []
            for i in range(n):
                row = []
                for j in range(n):
                    tail = [randomizer.randrange(prime) for _ in range(degree)]
                    row.append(tuple([(int(i == j)-sum(tail)) % prime]+tail))
                A.append(row)
            records.append(check_curve(A,2,degree,prime,'random strict polynomial'))
    q = 2
    A = [[(int(x == y and v == u),int(x == v and y == u)^int(x == y and v == u))
          for v in range(q) for u in range(q)] for x in range(q) for y in range(q)]
    records.append(check_curve(A,2,1,2,'two rank-one leading cuts'))
    identity_times_z = [[(0,int(i == j)) for j in range(4)] for i in range(4)]
    records.append(check_curve(identity_times_z,2,2,3,'both profile minima positive'))
    tensor,_ = cell()
    A = [[tuple(map(int,tensor[:,i//3,i % 3,j//3,j % 3])) for j in range(9)] for i in range(9)]
    records.append(check_curve(A,3,3,7,'qutrit cell with nontrivial word jets'))
    dual = crossed(reverse(A,3),3)
    records.append(check_curve(dual,3,3,7,'reversed opposite qutrit orientation'))
    return records


def run():
    records = controls()
    path = ROOT/'evidence/q13-joint-word-target.json'
    data = json.loads(path.read_text())
    r,degree = data['mu'],data['M']
    caps = inequalities(degree,r,r)
    bad = [cap for cap in caps if cap['slack'] < 0]
    assert len(caps) == 1816 and len(bad) == 52
    worst = min(bad,key=lambda cap:cap['slack'])
    assert (worst['physical_size'],worst['alpha'],worst['beta'],worst['lhs'],worst['rhs']) == (13,7,13,30588956,30576000)
    linear = first_jet_linear_caps(13)
    for side in (0,1):
        for cap in linear:
            lhs = sum(a*b for a,b in zip(cap['first'],r))+sum(a*b for a,b in zip(cap['second'],r))
            matching = next(c for c in caps if c['source_orientation'] == side and c['physical_size'] == cap['physical_size']
                            and (c['box_orientation'],c['alpha'],c['beta']) == tuple(cap['box']))
            assert lhs <= matching['lhs']
    names = ['tools/check_box_adjugate_word_caps.py','tools/box_adjugate_word_caps.py',
             'tools/check_bare_seven_to_twelve_letters.py','tools/higher_jet_word_caps.py',
             'tools/check_higher_jet_word_caps.py','tools/check_all_field_pascal_crossings.py',
             'tools/check_negation_odd_qutrit.py']
    return dict(scope='All-field bare scalar box-adjugate necessary conditions, literal finite controls, and rejection of one q13 point; no universal q13 rate floor.',
                literal_controls=records,target_rejection=dict(input_sha256=sha256(path.read_bytes()).hexdigest(),
                    checks=len(caps),violations=bad,violated_checks=len(bad),worst=worst),
                first_jet_majorization_checks=2*len(linear),
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
        assert result == json.loads((ROOT/'evidence/box-adjugate-word-controls.json').read_text())
    print(json.dumps({'scope':result['scope'],'top_level_controls':len(result['literal_controls']),
                      'target_violations':result['target_rejection']['violated_checks'],
                      'seconds':time.monotonic()-begun}))


if __name__ == '__main__':
    main()
