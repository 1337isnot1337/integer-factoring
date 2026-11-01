#!/usr/bin/env python3
"""Exact block and independent magnetic controls for the enlarged spin crossing.

Requires NumPy. The entry constructor itself uses only the standard library.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path
from random import Random

import numpy as np

from bounded_spin_crossing import paths, path_entry, path_norm, recursive_vector
from bounded_spin_six_vertex import monomial_entry, unveil_coefficients, integral_polynomial_entry
from check_bounded_spin_rank_algebra import coefficient_tensor, rank

ROOT = Path(__file__).resolve().parents[1]


def literal_integral_entry(four_paths):
    n = len(four_paths[0])
    vertices = {(0,0,0,0):(0,0), (1,1,1,1):(0,0), (1,0,0,1):(0,0),
                (0,1,0,1):(1,0), (1,0,1,0):(1,0), (0,1,1,0):(0,1)}
    result = defaultdict(int)
    vectors = [tuple(recursive_vector(p).items()) for p in four_paths]
    for choices in product(*vectors):
        words, value = [w for w,_ in choices], 1
        for _,coefficient in choices:
            value *= coefficient
        ell,c = 0,0
        for i in range(n):
            bits = tuple((w >> i)&1 for w in words)
            if bits not in vertices:
                break
            dt,dw = vertices[bits]
            ell,c = ell+dt,c+dw
        else:
            assert ell % 2 == 0
            result[n//2-ell//2-c,ell//2,c] += value
    scale = 12**n//(path_norm(four_paths[2])*path_norm(four_paths[3]))
    return unveil_coefficients({m:value*scale for m,value in result.items()},n,prime=None)


def pivot_columns(a):
    a = a.copy() % 5
    r, columns = 0, []
    for j in range(a.shape[1]):
        hits = np.flatnonzero(a[r:,j])
        if not len(hits):
            continue
        h = r+int(hits[0])
        a[[r,h]] = a[[h,r]]
        a[r,j:] = a[r,j:]*pow(int(a[r,j]),-1,5) % 5
        a[r+1:,j:] = (a[r+1:,j:]-a[r+1:,j,None]*a[r,None,j:]) % 5
        columns.append(j)
        r += 1
        if r == len(a):
            break
    return columns


def block_tensor(n):
    ps = paths(n)
    words = [w for w in range(1 << n) if w.bit_count() == n//2]
    indices = {w:i for i,w in enumerate(words)}
    frame = np.array([[path_entry(p,w) % 5 for w in words] for p in ps],dtype=np.int64)
    groups = defaultdict(list)
    for x,y in product(words,repeat=2):
        groups[x & y,x ^ y].append((x,y))
    r = n//2
    monomials = [(a,b,r-a-b) for a in range(r+1) for b in range(r-a+1)]
    slots = {m:i for i,m in enumerate(monomials)}
    size = len(ps)**2
    result = np.zeros((len(monomials),size,size),dtype=np.int64)
    for (_,difference),pairs in groups.items():
        s = difference.bit_count()
        a = np.stack([np.outer(frame[:,indices[x]],frame[:,indices[y]]).ravel()
                      for x,y in pairs],axis=1) % 5
        local = np.zeros((s//2+1,len(pairs),len(pairs)),dtype=np.int64)
        for i,(x,_) in enumerate(pairs):
            for j,(_,u) in enumerate(pairs):
                local[(x ^ u).bit_count()//2,i,j] = 1
        for k in range(s//2+1):
            slot = slots[r-s//2,s//2-k,k]
            result[slot] = (result[slot]+(a @ local[k] % 5) @ a.T) % 5
    q = len(ps)
    return monomials,result.reshape(len(monomials),q,q,q,q).transpose(0,1,2,4,3)


def operators(tensor,n):
    q = len(paths(n))
    norms = [path_norm(p) % 5 for p in paths(n)]
    inverse = np.array([pow(x*y,-1,5) for x,y in product(norms,repeat=2)])
    return [tensor.transpose(0,1,2,4,3).reshape(len(tensor),q*q,q*q)*inverse[None,None,:] % 5,
            tensor.transpose(0,1,3,4,2).reshape(len(tensor),q*q,q*q)*inverse[None,None,:] % 5]


def check():
    random = Random(1337)
    cases, controls, integer_controls = [], 0, 0
    for n in (4,6,8):
        ps = paths(n)
        monomials,t = block_tensor(n)
        chosen = (list(product(range(len(ps)),repeat=4)) if n == 4 else
                  [tuple(random.randrange(len(ps)) for _ in range(4)) for _ in range(12)])
        norm = [path_norm(p) for p in ps]
        endpoint = sum(t[j] for j,m in enumerate(monomials) if m[1] == 0) % 5
        for ids in product(range(len(ps)),repeat=4):
            expected = int(ids[0] == ids[3] and ids[1] == ids[2])
            assert int(endpoint[ids])*pow(norm[ids[2]]*norm[ids[3]],-1,5) % 5 == expected
        for ids in chosen:
            selected = tuple(ps[i] for i in ids)
            coefficients = monomial_entry(selected)
            for j,m in enumerate(monomials):
                expected = int(t[(j,)+ids])*pow(norm[ids[2]]*norm[ids[3]],-1,5) % 5
                assert coefficients.get(m,0) == expected
            poly = unveil_coefficients(coefficients,n)
            assert sum(poly) % 5 == int(ids[0] == ids[3] and ids[1] == ids[2])
            controls += 1
            if n == 4:
                integral = integral_polynomial_entry(selected)
                assert integral == literal_integral_entry(selected)
                assert [c % 5 for c in integral] == [pow(12,n,5)*c % 5 for c in poly]
                assert sum(integral) == 12**n*int(ids[0] == ids[3] and ids[1] == ids[2])
                integer_controls += 1
        # Full column rank certifies every monomial tensor is recoverable
        # from the coefficient tensors of the univariate construction.
        substitution = np.array([unveil_coefficients({m:1},n) for m in monomials])
        assert rank(substitution) == len(monomials)
        row_cut,column_cut = operators(t,n)
        slots = {m:i for i,m in enumerate(monomials)}
        for i,(a,b,c) in enumerate(monomials):
            assert np.array_equal(row_cut[i],column_cut[slots[c,b,a]])
        cases.append({'length':n,'labels':len(ps),'monomials':len(monomials),
                      'coefficient_span_dimension':rank(t.reshape(len(t),-1)),
                      'univariate_degree_bound':n//2*(n+1),
                      'complete_cut_reciprocity_identities':len(monomials)})

    # Recover the old pencil as u=1,v=(1-z)^2,w=z^2, entry by entry.
    old = np.zeros((9,)+t.shape[1:],dtype=np.int64)
    for (a,b,c),tensor in zip(monomials,t):
        for j in range(2*b+1):
            old[2*c+j] = (old[2*c+j]+comb(2*b,j)*(-1)**j*tensor) % 5
    assert np.array_equal(old,coefficient_tensor(8))
    old_row,old_column = operators(old,8)
    assert np.array_equal(old_row,old_column[::-1])
    basis = np.concatenate((old[[0,2,3,5,4,6]],t[0,None,:]))
    assert rank(basis.reshape(7,-1)) == 7
    modules = []
    for matrices in operators(basis,8):
        assert rank(matrices[-1]) == 34
        for i in range(6):
            for j in range(i):
                assert not np.any((matrices[i] @ matrices[j]-matrices[j] @ matrices[i]) % 5)
        commutators = [(matrices[-1] @ a-a @ matrices[-1]) % 5 for a in matrices[:-1]]
        images = np.concatenate(commutators,axis=1)
        module = images[:,pivot_columns(images)]
        assert module.shape == (169,28)
        assert rank(np.concatenate([module]+[a @ module % 5 for a in matrices],axis=1)) == 28
        modules.append({'dimension':28,'ambient_dimension':169,
                        'commutator_ranks':[rank(a) for a in commutators]})
    assert [c['coefficient_span_dimension'] for c in cases] == [3,4,7]
    return {'status':'pass','prime':5,'cases':cases,
            'independent_magnetic_entry_controls':controls,
            'integral_entry_controls':integer_controls,'old_pencil_recovered_exactly':True,
            'invariant_commutator_modules':modules,
            'scope':'Uniform enlarged coherent crossing and exact finite span/module controls; no favorable uniform rank theorem or improved factoring bound.',
            'source_sha256':{name:sha256((ROOT/name).read_bytes()).hexdigest() for name in
                             ('tools/bounded_spin_crossing.py','tools/bounded_spin_six_vertex.py',
                              'tools/check_bounded_spin_six_vertex.py','tools/check_bounded_spin_rank_algebra.py')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
