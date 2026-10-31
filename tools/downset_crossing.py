#!/usr/bin/env python3
"""Exact elementary crossing matrices on binary word codes."""
from __future__ import annotations

from collections import defaultdict
from itertools import product

from search_auxiliary_coefficient_span import combine,reverse_rank


def one_plus_powers(degree):
    result = [1]
    for _ in range(degree):
        result.append(result[-1]^(result[-1] << 1))
    return result


def is_downset(code):
    code = set(code)
    return bool(code) and all(x^(1 << i) in code for x in code
                              for i in range(x.bit_length()) if x >> i & 1)


def coefficient_matrices(code,other=None):
    code = sorted(code)
    paired = other is not None
    other = sorted(other) if paired else code
    n,h = len(code),len(other)
    index,second = {x:i for i,x in enumerate(code)},{x:i for i,x in enumerate(other)}
    blocks = defaultdict(list)
    for x,y in product(code,other):
        blocks[x&y,x^y].append((x,y))
    degree = max(s.bit_count() for _,s in blocks)
    cuts = [[[0]*(n*h) for _ in range(degree+1)] for _ in range(2)]
    powers = one_plus_powers(degree)
    support = 0
    for (_,diff),pairs in blocks.items():
        distance = diff.bit_count()
        for (x,y),(v,u) in product(pairs,repeat=2):
            changed = (x^u).bit_count()
            polynomial = powers[distance-changed] << changed
            support += 1
            while polynomial:
                bit = polynomial & -polynomial
                k = bit.bit_length()-1
                polynomial ^= bit
                cuts[0][k][index[x]*h+second[y]] |= 1 << (second[u]*n+index[v])
                cuts[1][k][index[x]*h+second[u]] |= 1 << (second[y]*n+index[v])
    endpoint = [1 << (j*n+i) for i,j in product(range(n),range(h))]
    assert all(combine(cut,(1 << (degree+1))-1) == endpoint for cut in cuts)
    metadata = {'code_size':n,'degree':degree,'row_blocks':len(blocks),
                'nonzero_entry_polynomials':support,'endpoint_matrix_rows':2*n*h}
    if paired:
        metadata['other_code_size'] = h
    return cuts,metadata


def column_blocks(code,other=None):
    paired = other is not None
    other = other if paired else code
    assert is_downset(code) and is_downset(other)
    blocks = defaultdict(list)
    for x,u in product(sorted(code),sorted(other)):
        diff = x^u
        blocks[diff,x&diff].append(x&~diff)
    result = []
    for (diff,orientation),outside in sorted(blocks.items()):
        assert len(outside) == len(set(outside))
        outside = tuple(sorted(outside))
        assert is_downset(outside)
        if not paired:
            assert outside == tuple(sorted(blocks[diff,diff^orientation]))
        result.append((diff,orientation,outside))
    return result


def normalized_column_rank(blocks,mask):
    result = 0
    for diff,_,outside in blocks:
        shift = diff.bit_count()
        rows = [sum(1 << j for j,c in enumerate(outside)
                    if not b&c and mask >> (shift+b.bit_count()+c.bit_count()) & 1)
                for b in outside]
        result += reverse_rank(rows)
    return result


def verify_local_normal_form(outside):
    assert is_downset(outside)
    degree = 2*max(s.bit_count() for s in outside)
    powers = one_plus_powers(degree)
    identities = 0
    for b,c in product(outside,repeat=2):
        # Z H Z^T, with constant Boolean zeta Z; all subsets are present.
        normalized = 0
        for s,t in product(outside,repeat=2):
            if s&b == s and t&c == t:
                normalized ^= powers[(s^t).bit_count()]
        assert normalized == (0 if b&c else 1 << (b.bit_count()+c.bit_count()))
        # H=L(z) diag(z^(2|s|)) L(z)^T with polynomial unitriangular L.
        reconstructed = 0
        for s in outside:
            if s&b == s and s&c == s:
                reconstructed ^= powers[b.bit_count()+c.bit_count()-2*s.bit_count()] << (2*s.bit_count())
        assert reconstructed == powers[(b^c).bit_count()]
        identities += 2
    return identities


def literal_entry(x,y,u,v,length):
    result = 1
    for j in range(length):
        a,b,c,d = [(word >> j)&1 for word in (x,y,u,v)]
        identity,crossing = int(a == c and b == d),int(a == d and b == c)
        local = identity|((identity^crossing) << 1)
        result = (result if local & 1 else 0) ^ (result << 1 if local & 2 else 0)
    return result
