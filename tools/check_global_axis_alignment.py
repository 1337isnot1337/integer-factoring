#!/usr/bin/env python3
"""Exact finite polynomial witnesses for the homogeneous alignment obstruction."""
import argparse
from hashlib import sha256
from itertools import combinations
import json
import lzma
from pathlib import Path
from random import Random

from high_axis_frobenius_residues import residue, crossblock
from high_axis_plucker import power, MUL
from weighted_smith_oracle_f4 import scalar_multiply, polynomial_multiply
from windowed_f4_smith import pack
from weighted_smith_oracle_f4 import WeightedSmithOracle
from native_matrix_smith_f4 import MatrixOracle

ROOT = Path(__file__).resolve().parents[1]


def add_term(poly,monomial,value):
    new = poly.get(monomial,0) ^ value
    if new:
        poly[monomial] = new
    else:
        poly.pop(monomial,None)


def endpoint(poly):
    value = 0
    for c in poly.values():
        value ^= c
    return value


def P_and_Q(coefficients,h,u):
    P = {}
    for b,poly in enumerate(coefficients):
        for monomial,c in poly.items():
            assert monomial[b] == 0
            value = list(monomial)
            value[b] += u
            add_term(P,tuple(value),c)
    Q = []
    for a in range(h):
        value = P.copy()
        for monomial,c in coefficients[a].items():
            exponents = list(monomial)
            exponents[a] += u
            add_term(value,tuple(exponents),c)
        Q.append(value)
    return P,Q


def check(coefficients,D,u,label):
    h = len(coefficients)
    assert h % 2 and D//(h-1) < u
    for b,poly in enumerate(coefficients):
        assert endpoint(poly) == 1
        assert all(sum(m)==D and m[b]==0 for m in poly)
    P,Q = P_and_Q(coefficients,h,u)
    assert endpoint(P) == 1 and all(sum(m)==D+u for m in P)
    monomial = min(P)
    b = next(i for i,e in enumerate(monomial) if e == u)
    a = min((i for i in range(h) if i != b),key=lambda i:monomial[i])
    assert monomial[a] <= D//(h-1) < u and Q[a].get(monomial) == P[monomial]
    orders = [min((m[i] for m in poly),default=None) for i,poly in enumerate(Q)]
    assert orders[a] <= D//(h-1)
    return dict(label=label,degree=D,power=u,set_size=h,coefficients=[
        [[list(m),c] for m,c in sorted(poly.items())] for poly in coefficients],
        P=[[list(m),c] for m,c in sorted(P.items())],
        Q_orders=orders,witness=dict(monomial=list(monomial),coefficient=P[monomial],
                                    distinguished_index=b,offending_axis=a,order=orders[a],bound=D//(h-1)))


def sharp_case(h):
    D,u = 64,32
    e,j,f = residue(h-1)
    coefficients = []
    for b in range(h):
        indices = [i for i in range(h) if i != b]
        poly = {}
        for chosen in combinations(indices,j):
            monomial = tuple(0 if i==b else e+(f if i in chosen else 0) for i in range(h))
            add_term(poly,monomial,1)
        coefficients.append(poly)
    result = check(coefficients,D,u,'degree-specific sharp Frobenius fixture')
    assert result['Q_orders'] == [D//(h-1)]*h
    return result


def quotient_minor():
    values = [2]+[1]*7
    matrix = crossblock(9,values)
    form = [pack(matrix[mask][0]) for mask in range(256)]
    U = [power(v,32) for v in values]
    product = [0]*256
    for mask in range(256):
        for i in range(8):
            if mask >> i & 1:
                product[mask] ^= scalar_multiply(form[mask ^ (1 << i)],U[i])
    assert product[15] == 1 << 32  # z^16, packed with two bits per field coefficient.
    inverse = (0,1,3,2)[U[0]]
    def change(coefficients):
        result = [0]*256
        for mask,value in enumerate(coefficients):
            if mask & 1:
                rest = mask ^ 1
                result[rest | 1] ^= scalar_multiply(value,inverse)
                for i in range(1,8):
                    if not (rest >> i & 1):
                        result[rest | (1 << i)] ^= scalar_multiply(value,MUL[inverse][U[i]])
            else:
                result[mask] ^= value
        return result
    transformed,transformed_product = change(form),change(product)
    assert all(transformed_product[mask] == (transformed[mask ^ 1] if mask & 1 else 0)
               for mask in range(256))
    cols = [j << 4 for j in range(16)]
    rows = [15 | j for j in cols]
    minor = [[product[i ^ j] if j & ~i == 0 else 0 for j in cols] for i in rows]
    assert all(not minor[i][j] for i in range(16) for j in range(i+1,16))
    determinant = 1
    for i in range(16):
        assert minor[i][i] == product[15]
        determinant = polynomial_multiply(determinant,minor[i][i])
    assert determinant == 1 << 512  # z^256.
    return dict(scope='One exact nonzero quotient-prefix fixture, without a full-axis pressure conclusion.',
                values=values,U=U,exterior_coordinate_change='new epsilon_0=U',
                transformed_form=transformed,coordinate_rows=rows,coordinate_columns=cols,
                coordinate_minor=minor,determinant_packed=determinant,
                quotient_dimension=128,prefix_size=16,prefix_sum_upper_bound=256)


def module_counterexample():
    def matrix(variables,coefficients):
        n = 1 << variables
        return [[(coefficients.get(i ^ j,0),0) if j & ~i == 0 else (0,0)
                 for j in range(n)] for i in range(n)]
    small = matrix(4,{0:1 << 16,3:1 << 10,15:1 << 8})
    b = [8,10,10,12]+[16]*8+[20,22,22,24]
    got,_ = WeightedSmithOracle([[pack(v) for v in row] for row in small]).profile([0]*16,[0]*16)
    assert list(got) == b
    quotient = matrix(7,{0:1 << 64,7:1 << 16,31:1 << 10,127:1 << 8})
    quotient_orders = b+[64]*96+[128-v for v in b[::-1]]
    assert MatrixOracle(quotient,128,precision=256).profile()['orders'] == quotient_orders
    assert [sum(quotient_orders[:m]) for m in (1,4,16)] == [8,40,256]
    full = matrix(8,{0:1 << 64,1:1 << 32,14:1 << 16,62:1 << 10,254:1 << 8})
    low = sorted(b+b+[32]*96)
    orders = low+[128-v for v in low[::-1]]
    assert MatrixOracle(full,256,precision=256).profile()['orders'] == orders
    assert sum(orders[:86]) == 2240 and 3*sum(orders[:86]) > 64*86
    return dict(scope='Exterior-module counterexample to lifting quotient bounds alone; its endpoint is not a crossing.',
                small_untruncated_orders=b,quotient_orders=quotient_orders,
                full_orders=orders,first_86_sum=2240,required_upper_bound_times_three=64*86,
                actual_sum_times_three=3*2240)


def run():
    randomizer,controls = Random(2026110119),[]
    for h in (5,7,9):
        for D in (8,16):
            for _ in range(16):
                coefficients = []
                for b in range(h):
                    indices = [i for i in range(h) if i != b]
                    poly = {}
                    for term in range(8):
                        exponents = [0]*h
                        for unit in range(D):
                            exponents[randomizer.choice(indices)] += 1
                        add_term(poly,tuple(exponents),randomizer.randrange(1,4))
                    correction = endpoint(poly) ^ 1
                    monomial = tuple(D if i == indices[0] else 0 for i in range(h))
                    add_term(poly,monomial,correction)
                    coefficients.append(poly)
                controls.append(check(coefficients,D,D//2,'arbitrary homogeneous F4 coefficients'))
    sharp = [sharp_case(h) for h in (5,7,9)]
    names = ['tools/check_global_axis_alignment.py','tools/high_axis_frobenius_residues.py',
             'tools/high_axis_plucker.py','tools/weighted_smith_oracle_f4.py','tools/windowed_f4_smith.py',
             'tools/native_matrix_smith_f4.py','tools/native_matrix_smith_f4.cpp',
             'tools/native_rectangular_smith_f4.cpp']
    return dict(scope='Uniform homogeneous divisibility and quotient-prefix obstructions; no full-axis pressure or factoring bound.',
                characteristic=2,arbitrary_coefficient_controls=controls,sharp_fixtures=sharp,
                quotient_minor_fixture=quotient_minor(),
                quotient_lifting_counterexample=module_counterexample(),
                control_count=len(controls),sharp_orders=[case['Q_orders'][0] for case in sharp],
                source_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in names})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        payload = (json.dumps(result,indent=2)+'\n').encode()
        if args.output.suffix == '.xz':
            args.output.write_bytes(lzma.compress(payload))
        else:
            args.output.write_bytes(payload)
    else:
        with lzma.open(ROOT/'evidence/global-axis-alignment-controls.json.xz','rt') as stream:
            assert result == json.load(stream)
    print(json.dumps({k:result[k] for k in ('scope','control_count','sharp_orders')}))


if __name__ == '__main__':
    main()
