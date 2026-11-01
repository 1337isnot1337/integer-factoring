#!/usr/bin/env python3
"""F16 paired path encoders without a classical spin-invariance assumption."""
from __future__ import annotations

from functools import lru_cache
from itertools import product

from bounded_spin_crossing import binary_path
from check_fourier_multiplication_gates import BinaryField

FIELD = BinaryField(0b10011)
MULTIPLY = tuple(tuple(FIELD.multiply(a,b) for b in range(16)) for a in range(16))


def power(q,k):
    if not q:
        raise ValueError('Use a nonzero field element.')
    return FIELD.power(q,k % 15)


@lru_cache(None)
def splitting_blocks(q=8):
    denominator = q ^ power(q,-1)
    if not denominator:
        raise ValueError('The q-integer denominator vanishes.')

    def integer(k):
        return MULTIPLY[power(q,k) ^ power(q,-k)][power(denominator,-1)]

    result = {}
    for h in range(1,4):
        determinant = integer(h+1)
        if not determinant:
            raise ValueError('The retained splitting has zero determinant.')
        for m in range(1,h+1):
            a,b = 1,integer(m)
            c,d = power(q,-h+m-1),MULTIPLY[power(q,m)][integer(h-m+1)]
            assert MULTIPLY[a][d] ^ MULTIPLY[b][c] == determinant
            inverse = tuple(tuple(MULTIPLY[power(determinant,-1)][x] for x in row)
                            for row in ((d,b),(c,a)))
            result[h,m] = ((a,b),(c,d)),inverse
    return result


def step(height,direction,magnetic,bit,blocks,dual=False):
    sector = magnetic+bit
    following = sector-int(direction == -1)
    if not 0 <= following <= height+direction:
        return None
    scalar = 1
    if 0 < sector < height+1:
        primal,inverse = blocks[height,sector]
        branch = int(direction == -1)
        scalar = inverse[branch][bit] if dual else primal[bit][branch]
    return following,scalar


def path_entry(path,word,q=8,dual=False):
    height,magnetic,value = 0,0,1
    blocks = splitting_blocks(q)
    for i,direction in enumerate(path):
        transition = step(height,direction,magnetic,(word >> i)&1,blocks,dual)
        if transition is None:
            return 0
        magnetic,scalar = transition
        value = MULTIPLY[value][scalar]
        height += direction
    return value if magnetic == 0 else 0


def polynomial_entry(four_paths,q=8):
    """Coefficients in z, with primal X,Y and inverse-dual U,V legs."""
    length = len(four_paths[0])
    assert all(len(p) == length and sum(p) == 0 for p in four_paths)
    blocks = splitting_blocks(q)
    states,heights = {(0,0,0,0):[1]},[0]*4
    for i in range(length):
        directions = [p[i] for p in four_paths]
        following = {}
        for magnetic,polynomial in states.items():
            for bits in product(range(2),repeat=4):
                x,y,u,v = bits
                identity,crossing = int(x == u and y == v),int(x == v and y == u)
                local = (identity,identity ^ crossing)
                if not any(local):
                    continue
                target,scalar = [],1
                for leg,(h,s,r,b) in enumerate(zip(heights,directions,magnetic,bits)):
                    transition = step(h,s,r,b,blocks,dual=leg >= 2)
                    if transition is None:
                        scalar = 0
                        break
                    nr,c = transition
                    target.append(nr)
                    scalar = MULTIPLY[scalar][c]
                if not scalar:
                    continue
                out = following.setdefault(tuple(target),[0]*(i+2))
                for j,value in enumerate(polynomial):
                    value = MULTIPLY[scalar][value]
                    if local[0]:
                        out[j] ^= value
                    if local[1]:
                        out[j+1] ^= value
        states = following
        heights = [h+s for h,s in zip(heights,directions)]
    return states.get((0,0,0,0),[0]*(length+1))


def binary_entry(width,labels,q=8):
    return polynomial_entry(tuple(binary_path(width,label) for label in labels),q)
