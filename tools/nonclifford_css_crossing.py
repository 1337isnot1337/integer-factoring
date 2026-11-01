#!/usr/bin/env python3
"""A fixed order-five local twist of the four-qutrit CSS crossing split."""
from hashlib import sha256
from itertools import product

import numpy as np

from collective_css_pauli import code, compress, polynomial_multiply
from cyclic_pauli_crossings import construct as cell_construct, monomial, parity, literal_seed
from physical_split_crossings import F4_PRODUCT
from regularized_transpose_words import scalar

TABLE = np.asarray(F4_PRODUCT,dtype=np.uint8)


def matmul(a,b):
    result = np.zeros((a.shape[0],b.shape[1]),dtype=np.uint8)
    for k in range(a.shape[1]):
        result ^= TABLE[a[:,k,None],b[k,None,:]]
    return result


def rotation():
    identity = np.eye(3,dtype=np.uint8)
    v,w = np.array([1,2,3],dtype=np.uint8),np.array([1,1,0],dtype=np.uint8)
    O = matmul(identity ^ TABLE[v[:,None],v[None,:]], identity ^ TABLE[w[:,None],w[None,:]])
    assert np.array_equal(matmul(O.T,O),identity) and not np.array_equal(O,identity)
    power = identity
    for _ in range(5):
        power = matmul(power,O)
    assert np.array_equal(power,identity)
    weyls = []
    for a,b,c in product(range(3),range(3),range(1,4)):
        permutation,phases = monomial((a,),(b,))
        matrix = np.zeros((3,3),dtype=np.uint8)
        matrix[permutation,np.arange(3)] = phases
        weyls.append(TABLE[c,matrix])
    permutation,phases = monomial((1,),(0,))
    X = np.zeros((3,3),dtype=np.uint8)
    X[permutation,np.arange(3)] = phases
    conjugate = matmul(matmul(O.T,X),O)
    assert not any(np.array_equal(conjugate,W) for W in weyls)
    return O,conjugate


def twist_cell(planes,O):
    tensor = [p.reshape(3,3,3,3).copy() for p in planes]
    for axis in range(4):
        moved = [np.moveaxis(p,axis,0) for p in tensor]
        result = [np.zeros_like(p) for p in moved]
        for i,j in product(range(3),repeat=2):
            low,high = scalar(moved[0][j],moved[1][j],int(O[j,i]))
            result[0][i] ^= low
            result[1][i] ^= high
        tensor = [np.moveaxis(p,0,axis) for p in result]
    return [p.reshape(9,9) for p in tensor]


def twisted_embedding(O):
    _,embedding = code()
    original = embedding.reshape(3,27,9)
    result = np.zeros_like(original)
    for i,j in product(range(3),repeat=2):
        result[i] ^= TABLE[O[i,j],original[j]]
    result = result.reshape(81,9)
    assert np.array_equal(matmul(result.T,result),np.eye(9,dtype=np.uint8))
    detections = 0
    for site,shift,phase in product(range(4),range(3),range(3)):
        if not (shift or phase):
            continue
        a,b = [0]*4,[0]*4
        a[site],b[site] = shift,phase
        permutation,phases = monomial(tuple(a),tuple(b))
        acted = np.zeros_like(result)
        acted[permutation] = TABLE[phases[:,None],result]
        assert not np.any(matmul(result.T,acted))
        detections += 1
    assert detections == 32
    return result,detections


def construct():
    O,conjugate = rotation()
    embedding,detections = twisted_embedding(O)
    words,_ = code()
    original,cell_metadata = cell_construct(1)
    changed = twist_cell(original,O)
    assert np.array_equal(parity(changed[0]),np.eye(9)) and not np.any(parity(changed[1]))
    logical = np.array(list(product(range(9),repeat=4)),dtype=np.int64)
    shifts = np.array(list(product(range(3),repeat=4)),dtype=np.int64)
    physical = [words[logical[:,leg,None],shifts[None,:,leg],:] for leg in range(4)]
    result = [np.ones((9**4,3**4),dtype=np.uint64),np.zeros((9**4,3**4),dtype=np.uint64)]
    for site in range(4):
        local = changed if site == 0 else original
        x,y,v,u = (s[:,:,site] for s in physical)
        result = polynomial_multiply(result,[p[3*x+y,3*v+u] for p in local],3)
    planes = [np.bitwise_xor.reduce(p,axis=1).reshape(81,81) for p in result]
    assert np.array_equal(parity(planes[0]),np.eye(81)) and not np.any(parity(planes[1]))
    assert not any(np.any(p >> np.uint64(13)) for p in planes)
    old,old_metadata = compress()
    odd_counts = {str(j):int(np.count_nonzero(((planes[0] | planes[1]) >> np.uint64(j)) & 1))
                  for j in range(1,13,2)}
    metadata = dict(q=9,degree=12,physical_length=4,logical_length=2,rotation=O.tolist(),
                    conjugated_X=conjugate.tolist(),rotation_order=5,
                    single_site_nonidentity_Weyl_annihilations=detections,
                    embedding_sha256=sha256(embedding.tobytes()).hexdigest(),
                    cell=cell_metadata,code_words=words.tolist(),odd_degree_nonzero_entries=odd_counts,
                    structural_gate=bool(sum(odd_counts.values())),
                    identical_to_old_code=all(np.array_equal(a,b) for a,b in zip(planes,old)),
                    old_planes_sha256=old_metadata['planes_sha256'],
                    planes_sha256=sha256(b''.join(p.astype('<u8').tobytes() for p in planes)).hexdigest())
    return planes,metadata,original,changed


def seed():
    planes,metadata,_,_ = construct()
    return literal_seed(planes,metadata),planes,metadata
