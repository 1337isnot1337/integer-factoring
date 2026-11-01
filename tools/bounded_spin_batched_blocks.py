#!/usr/bin/env python3
"""Exact batched six-vertex coordinate blocks; requires NumPy and FLINT."""
from __future__ import annotations

from collections import defaultdict
from itertools import product
from math import comb

import numpy as np

from bounded_spin_crossing import paths,path_entry
from native_gf5_flint import multiply


def coefficient_tensor(n):
    ps = paths(n)
    words = [w for w in range(1 << n) if w.bit_count() == n//2]
    indices = {w:i for i,w in enumerate(words)}
    frame = np.array([[path_entry(p,w) % 5 for w in words] for p in ps],dtype=np.uint8)
    groups = defaultdict(list)
    for x,y in product(words,repeat=2):
        groups[x & y,x ^ y].append((x,y))
    r = n//2
    monomials = [(a,b,r-a-b) for a in range(r+1) for b in range(r-a+1)]
    slots = {m:i for i,m in enumerate(monomials)}
    size = len(ps)**2
    result = np.zeros((len(monomials),size,size),dtype=np.uint8)
    for s in range(0,n+1,2):
        blocks = [pairs for (_,difference),pairs in groups.items() if difference.bit_count() == s]
        local_size = comb(s,s//2)
        assert all(len(block) == local_size for block in blocks)
        pairs = [pair for block in blocks for pair in block]
        x = frame[:,[indices[x] for x,_ in pairs]]
        y = frame[:,[indices[y] for _,y in pairs]]
        a = (x[:,None,:]*y[None,:,:]).reshape(size,len(pairs)) % 5
        canonical = [w for w in range(1 << s) if w.bit_count() == s//2]
        mask = (1 << s)-1
        distances = np.array([[(x ^ (mask ^ y)).bit_count()//2 for y in canonical]
                              for x in canonical],dtype=np.uint8)
        for k in range(s//2+1):
            local = (distances == k).astype(np.uint8)
            left = multiply(a.reshape(size*len(blocks),local_size),local).reshape(size,len(pairs))
            result[slots[r-s//2,s//2-k,k]] = multiply(left,a.T)
    q = len(ps)
    return monomials,result.reshape(len(monomials),q,q,q,q).transpose(0,1,2,4,3)
