#!/usr/bin/env python3
"""Exact ten-site algebra and correction-module growth; optional FLINT replay."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path
from random import Random

import numpy as np

from bounded_spin_crossing import paths,path_norm
from bounded_spin_six_vertex import monomial_entry
from bounded_spin_batched_blocks import coefficient_tensor
from check_bounded_spin_six_vertex import block_tensor,operators
from check_bounded_spin_rank_algebra import rank as reference_rank
from native_gf5_flint import library,multiply,rank,rref,pivot_columns,inverse

ROOT = Path(__file__).resolve().parents[1]


def combination(coefficients,matrices):
    return (np.tensordot(np.asarray(coefficients,dtype=np.int16),
                        matrices.astype(np.int16),axes=(0,0)) % 5).astype(np.uint8)


def original_coefficients(monomials,t,n):
    result = np.zeros((n+1,)+t.shape[1:],dtype=np.uint8)
    for (_,b,c),tensor in zip(monomials,t):
        for j in range(2*b+1):
            result[2*c+j] = (result[2*c+j]+(comb(2*b,j)*(-1)**j % 5)*tensor) % 5
    return result


def canon(rows):
    if not len(rows):
        return ()
    reduced,r = rref(np.array(rows,dtype=np.int64))
    return tuple(tuple(map(int,row)) for row in reduced[:r])


def kernel(rows,dimension):
    pivots = [next(i for i,x in enumerate(row) if x) for row in rows]
    result = []
    for j in range(dimension):
        if j in pivots:
            continue
        v = [int(i == j) for i in range(dimension)]
        for p,row in zip(pivots,rows):
            v[p] = -row[j] % 5
        result.append(v)
    return result


def dot(a,b):
    return sum(int(x)*int(y) for x,y in zip(a,b)) % 5


def check_strata(model,reciprocity,b,raw_basis):
    characters = np.array(model['characters'],dtype=np.int64)
    e = np.array(model['idempotents'],dtype=np.int64)
    radical = np.array(model['radical'],dtype=np.int64)
    p = next(i for i,x in enumerate(radical) if x)
    derivative = (np.eye(8,dtype=np.int64)[p]-e[:,p] @ characters)*pow(int(radical[p]),-1,5) % 5
    row_forms = np.concatenate((characters,derivative[None,:]))
    forms = np.concatenate((row_forms,row_forms @ reciprocity.T % 5))
    hyperplanes = {tuple(map(int,f*pow(int(next(x for x in f if x)),-1,5) % 5)) for f in forms}
    flats = {()}
    for form in sorted(hyperplanes):
        flats |= {canon(list(flat)+[form]) for flat in list(flats)}
    dims,nil,nilrank = model['dimensions'],model['nil_component'],model['nil_rank']

    def ranks(vanishing):
        result = []
        for cut in (0,1):
            value = sum(d for j,d in enumerate(dims) if not vanishing[cut*8+j])
            if vanishing[cut*8+nil] and not vanishing[cut*8+7]:
                value += nilrank
            result.append(value)
        return result

    threshold = model['optimization_threshold']
    separator,unit = model['optimization_separator'],model['unit']
    assert dot(separator,unit) == 1
    cheaper = []
    for flat in flats:
        space = kernel(flat,8)
        vanishing = [not any(dot(f,v) for v in space) for f in forms]
        r,c = ranks(vanishing)
        if r*c < threshold:
            assert all(dot(separator,v) == 0 for v in space)
            cheaper.extend(space)
    assert len(flats) == 9805 and len(canon(cheaper)) == 7
    witness = np.array(model['optimization_basis'],dtype=np.int64)
    assert rank(witness) == 8
    coefficients = multiply(np.array(unit,dtype=np.uint8)[None,:],inverse(witness))[0]
    pieces = [(int(c)*w % 5).tolist() for c,w in zip(coefficients,witness) if c]
    assert np.array_equal(np.array(pieces).sum(axis=0) % 5,unit)
    pairs = []
    for w in pieces:
        expected = ranks([dot(f,w) == 0 for f in forms])
        direct = combination(w,raw_basis)
        actual = [rank(combination(w,b)),rank(operators(direct[None,:],10)[1][0])]
        assert actual == expected and actual[0]*actual[1] <= threshold
        pairs.append(actual)
    assert max(r*c for r,c in pairs) == threshold
    return {'common_kernel_flats':len(flats),'sharp_all_extension_minmax':threshold,
            'separator':separator,'attaining_piece_rank_pairs':pairs}


def check():
    rng = np.random.default_rng(1337)
    native_cases = 0
    for rows,inner,columns in ((1,2,3),(5,8,4),(31,34,32),(68,71,69)):
        a = rng.integers(-5,5,size=(rows,inner),dtype=np.int64)
        b = rng.integers(-5,5,size=(inner,columns),dtype=np.int64)
        assert np.array_equal(multiply(a,b),a @ b % 5)
        pivots = pivot_columns(a)
        assert rank(a) == reference_rank(a) == len(pivots)
        assert rank(a[:,pivots]) == len(pivots)
        native_cases += 1
    monomials8,t8 = coefficient_tensor(8)
    control_monomials,control = block_tensor(8)
    assert monomials8 == control_monomials and np.array_equal(t8,control)
    n = 10
    monomials,t = coefficient_tensor(n)
    ps,q = paths(n),len(paths(n))
    size = q*q
    random = Random(1337)
    for _ in range(16):
        labels = tuple(random.randrange(q) for _ in range(4))
        expected = monomial_entry(tuple(ps[i] for i in labels))
        denominator = path_norm(ps[labels[2]])*path_norm(ps[labels[3]])
        for j,m in enumerate(monomials):
            assert int(t[(j,)+labels])*pow(denominator,-1,5) % 5 == expected.get(m,0)
    model_path = ROOT/'evidence/bounded-spin-ten-algebra-model.json'
    model = json.loads(model_path.read_text())
    old = original_coefficients(monomials,t,n)
    old_operators = operators(old,n)[0].astype(np.uint8)
    b = old_operators[model['basis_indices']]
    extreme = operators(t[:1],n)[0][0].astype(np.uint8)
    full_basis = np.concatenate((b,extreme[None,:]))
    sample = set(map(int,rng.choice(size*size,512,replace=False)))
    for p,r in product(range(q),repeat=2):
        sample.update(((p*q+r)*size+(p*q+r),(p*q+r)*size+(r*q+p),
                       (p*q+p)*size+(r*q+r)))
    sample = np.array(sorted(sample))
    flat = full_basis.reshape(9,-1)
    chosen = pivot_columns(flat[:,sample])
    assert len(chosen) == 9
    chosen = sample[chosen]
    vi = inverse(flat[:,chosen])
    norm = [path_norm(p) % 5 for p in ps]
    inverse_norms = np.array([pow(x*y,-1,5) for x,y in product(norm,repeat=2)],dtype=np.uint8)
    for tensor in t:
        matrix = tensor.transpose(0,1,3,2).reshape(size,size)*inverse_norms[None,:] % 5
        coords = multiply(matrix.ravel()[chosen][None,:],vi)[0]
        assert np.array_equal(combination(coords,full_basis),matrix)
    positions = sample[pivot_columns(b.reshape(8,-1)[:,sample])]
    vi = inverse(b.reshape(8,-1)[:,positions])
    old_coordinates = []
    for matrix in old_operators:
        coords = multiply(matrix.ravel()[positions][None,:],vi)[0]
        assert np.array_equal(combination(coords,b),matrix)
        old_coordinates.append(coords)
    reciprocity = np.array(old_coordinates)[[n-j for j in model['basis_indices']]]
    assert reciprocity.tolist() == model['optimization_reciprocity']
    unit = np.array(model['unit'],dtype=np.uint8)
    regular = np.array(model['regular'],dtype=np.uint8)
    assert np.array_equal(combination(unit,b),np.eye(size,dtype=np.uint8))
    assert np.array_equal(old_operators.sum(axis=0) % 5,np.eye(size,dtype=np.uint8))
    assert np.array_equal(np.array(old_coordinates).sum(axis=0) % 5,unit)
    for i,s in enumerate(regular):
        assert np.array_equal(multiply(s,unit[:,None])[:,0],np.eye(8,dtype=np.uint8)[:,i])
        for other in regular:
            assert np.array_equal(multiply(s,other),multiply(other,s))
    generators = np.array(model['generators'],dtype=np.uint8)
    actions = [combination(g,regular) for g in generators]
    for gamma,action in zip(generators,actions):
        matrix = combination(gamma,b)
        for j,piece in enumerate(b):
            assert np.array_equal(multiply(matrix,piece),combination(action[:,j],b))
    words = [unit]
    for v in words:
        for action in actions:
            w = multiply(action,v[:,None])[:,0]
            if rank(np.array(words+[w])) > len(words):
                words.append(w)
        if len(words) == 8:
            break
    assert len(words) == 8
    idempotents = np.array(model['idempotents'],dtype=np.uint8)
    radical = np.array(model['radical'],dtype=np.uint8)
    characters = np.array(model['characters'],dtype=np.uint8)
    assert np.array_equal(idempotents.sum(axis=0) % 5,unit)
    assert rank(np.concatenate((idempotents,radical[None,:]))) == 8
    assert np.array_equal(multiply(characters,idempotents.T),np.eye(7,dtype=np.uint8))
    assert np.array_equal(multiply(characters,unit[:,None]),np.ones((7,1),dtype=np.uint8))
    assert not np.any(multiply(characters,radical[:,None]))
    for i,e in enumerate(idempotents):
        action = combination(e,regular)
        for j,f in enumerate(idempotents):
            assert np.array_equal(multiply(action,f[:,None])[:,0],e if i == j else np.zeros(8,dtype=np.uint8))
        assert rank(action) == (2 if i == model['nil_component'] else 1)
    assert not np.any(multiply(combination(radical,regular),radical[:,None]))
    for s,c in zip(regular,np.array(model['characters']).T):
        assert np.array_equal(multiply(s,radical[:,None])[:,0],int(c[model['nil_component']])*radical % 5)
    dimensions = [rank(combination(e,b)) for e in idempotents]
    assert dimensions == model['dimensions'] and sum(dimensions) == size
    assert rank(combination(radical,b)) == model['nil_rank']
    assert len(model['characters']) == 7 > 5
    vertex_columns = extreme[:,pivot_columns(extreme)]
    assert vertex_columns.shape[1] == 116
    module = np.concatenate([multiply(s,vertex_columns) for s in b],axis=1)
    module_dimension = rank(module)
    assert module_dimension == 170
    strata = check_strata(model,reciprocity,b,old[model['basis_indices']])
    # A small entry certificate disproves the linear B+vertex ansatz at twelve sites.
    ps12 = paths(12)
    mons12 = [(a,b,6-a-b) for a in range(7) for b in range(7-a)]
    rows,old_rows,entry_labels = [],[],[]
    random = Random(1337)
    for i in range(150):
        ids = [random.randrange(len(ps12)) for _ in range(4)]
        if i % 4 == 0:
            ids[2:] = ids[1::-1]
        elif i % 4 == 1:
            ids[2:] = ids[:2]
        elif i % 4 == 2:
            ids = [ids[0],ids[0],ids[1],ids[1]]
        values = monomial_entry(tuple(ps12[j] for j in ids))
        row = [values.get(m,0) for m in mons12]
        poly = [0]*13
        for (_,bb,cc),value in zip(mons12,row):
            for j in range(2*bb+1):
                poly[2*cc+j] = (poly[2*cc+j]+value*comb(2*bb,j)*(-1)**j) % 5
        rows.append(row)
        old_rows.append(poly)
        entry_labels.append(ids)
        if rank(np.array(rows)) == 14:
            break
    assert rank(np.array(rows)) == 14 and rank(np.array(old_rows)) == 10
    with_vertex = np.column_stack((old_rows,np.array(rows)[:,0]))
    assert rank(with_vertex) == 11
    files = ('tools/native_gf5_flint.c','tools/native_gf5_flint.py',
             'tools/bounded_spin_batched_blocks.py','tools/check_bounded_spin_growth.py',
             'tools/bounded_spin_crossing.py','tools/bounded_spin_six_vertex.py',
             'tools/check_bounded_spin_six_vertex.py',
             'evidence/bounded-spin-ten-algebra-model.json')
    return {'status':'pass','flint_version':library()[1],'native_reference_controls':native_cases,
            'batched_eight_site_tensor_matches_reference':True,'ten_site_magnetic_entries':16,
            'length':10,'labels':q,'original_span_dimension':8,'six_vertex_span_dimension':9,
            'algebra':'F5^6 times F5[e]/(e^2)','F5_characters':7,'single_F5_generator_impossible':True,
            'verified_F5_generators':2,'primary_dimensions':dimensions,'nilpotent_rank':model['nil_rank'],
            'vertex_rank':116,'vertex_generated_module_dimension':module_dimension,
            'original_optimization':strata,
            'twelve_site_entry_certificate':{'labels':89,'entry_labels':entry_labels,
                                             'full_span_dimension':14,'original_span_dimension':10,
                                             'original_plus_vertex_dimension':11},
            'scope':'Exact finite algebra and correction-module structure, with a conditional uniform module bound; no uniform favorable crossing decomposition.',
            'source_sha256':{name:sha256((ROOT/name).read_bytes()).hexdigest() for name in files}}


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
