#!/usr/bin/env python3
"""Coherent height-three spin paths and integral/F5 crossing constructors."""
from __future__ import annotations

from collections import defaultdict
from itertools import product

PRIME = 5
HEIGHT = 3


def paths(length,end_height=0):
    def walk(prefix,height):
        if len(prefix) == length:
            if height == end_height:
                yield prefix
            return
        for step in (1,-1):
            if 0 <= height+step <= HEIGHT:
                yield from walk(prefix+(step,),height+step)
    return tuple(walk((),0))


def suffix_counts(length,end_height=0):
    counts = [[int(h == end_height) for h in range(HEIGHT+1)]]
    for _ in range(length):
        counts.append([sum(counts[-1][j] for j in (h-1,h+1) if 0 <= j <= HEIGHT)
                       for h in range(HEIGHT+1)])
    return counts


def unrank_path(length,index,end_height=0):
    counts = suffix_counts(length,end_height)
    assert 0 <= index < counts[length][0]
    result,height = [],0
    for remaining in range(length-1,-1,-1):
        up = counts[remaining][height+1] if height < HEIGHT else 0
        if index < up:
            result.append(1)
            height += 1
        else:
            index -= up
            result.append(-1)
            height -= 1
            assert height >= 0
    assert height == end_height and index == 0
    return tuple(result)


def binary_path(width,label):
    assert width >= 1 and 0 <= label < 1 << width
    length = 4
    while suffix_counts(length)[length][0] < 1 << width:
        length += 2
    return unrank_path(length,label)


def magnetic_step(height,step,magnetic,bit):
    following = magnetic+bit-int(step == -1)
    if not 0 <= following <= height+step:
        return None
    scalar = 1 if step == 1 else height-magnetic if bit else -magnetic
    return following,scalar


def path_entry(path,word):
    height,magnetic,value = 0,0,1
    for i,step in enumerate(path):
        assert step in (-1,1) and 0 <= height+step <= HEIGHT
        transition = magnetic_step(height,step,magnetic,(word >> i)&1)
        if transition is None:
            return 0
        magnetic,scalar = transition
        value *= scalar
        height += step
    return value if magnetic == 0 else 0


def path_norm(path):
    height,result = 0,1
    for step in path:
        assert step in (-1,1) and 0 <= height+step <= HEIGHT
        if step == -1:
            result *= height*(height+1)
        height += step
    assert result % PRIME
    return result


def recursive_vector(path):
    # Independent integral highest-weight recursion, used for controls.
    vector,height = {0:1},0
    for n,step in enumerate(path):
        if step == -1:
            following = defaultdict(int)
            for word,value in vector.items():
                following[word|(1 << n)] += height*value
                for j in range(n):
                    if not word >> j & 1:
                        following[word|(1 << j)] -= value
            vector = {word:value for word,value in following.items() if value}
        height += step
    return vector


def _raw_polynomial_entry(four_paths,prime):
    length = len(four_paths[0])
    assert all(len(path) == length for path in four_paths)
    heights = [0]*4
    states = {(0,0,0,0):[1]}
    for i in range(length):
        steps = [path[i] for path in four_paths]
        following = {}
        for magnetic,polynomial in states.items():
            for bits in product(range(2),repeat=4):
                target,scalar = [],1
                for h,s,r,b in zip(heights,steps,magnetic,bits):
                    transition = magnetic_step(h,s,r,b)
                    if transition is None:
                        scalar = 0
                        break
                    nr,c = transition
                    target.append(nr)
                    scalar *= c
                    if prime is not None:
                        scalar %= prime
                if not scalar:
                    continue
                x,y,u,v = bits
                identity,crossing = int(x == u and y == v),int(x == v and y == u)
                if not identity and not crossing:
                    continue
                out = following.setdefault(tuple(target),[0]*(i+2))
                for j,c in enumerate(polynomial):
                    out[j] += scalar*identity*c
                    out[j+1] += scalar*(crossing-identity)*c
                    if prime is not None:
                        out[j] %= prime
                        out[j+1] %= prime
        states = following
        heights = [h+s for h,s in zip(heights,steps)]
        assert all(0 <= h <= HEIGHT for h in heights)
    result = states.get((0,0,0,0),[0]*(length+1))
    return result


def polynomial_entry(four_paths,normalize=True):
    result = _raw_polynomial_entry(four_paths,PRIME)
    if normalize:
        inverse = pow(path_norm(four_paths[2])*path_norm(four_paths[3]),-1,PRIME)
        result = [c*inverse % PRIME for c in result]
    return result


def integral_polynomial_entry(four_paths):
    result = _raw_polynomial_entry(four_paths,None)
    denominator = path_norm(four_paths[2])*path_norm(four_paths[3])
    scale = 12**len(four_paths[0])
    assert scale % denominator == 0
    return [c*(scale//denominator) for c in result]


def literal_entry(four_paths):
    vectors = [recursive_vector(path) for path in four_paths]
    length = len(four_paths[0])
    result = [0]*(length+1)
    for selected in product(*(tuple(v.items()) for v in vectors)):
        words = [word for word,_ in selected]
        scalar = 1
        for _,c in selected:
            scalar = scalar*c % PRIME
        polynomial = [scalar]
        for i in range(length):
            x,y,u,v = [(word >> i)&1 for word in words]
            identity,crossing = int(x == u and y == v),int(x == v and y == u)
            out = [0]*(len(polynomial)+1)
            for j,c in enumerate(polynomial):
                out[j] = (out[j]+identity*c) % PRIME
                out[j+1] = (out[j+1]+(crossing-identity)*c) % PRIME
            polynomial = out
        result = [(a+b) % PRIME for a,b in zip(result,polynomial)]
    inverse = pow(path_norm(four_paths[2])*path_norm(four_paths[3]),-1,PRIME)
    return [c*inverse % PRIME for c in result]
