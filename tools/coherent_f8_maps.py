#!/usr/bin/env python3
"""Exact F8 coherent retractions and complete two-cut occupation profiles.

Vectors are sparse dictionaries {occupation_word: F8_coefficient}.
F8 uses x^3=x+1, encoded by integers0 through7. Pairings are bilinear.
"""
from __future__ import annotations
import ctypes
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import os, platform, shlex, shutil, subprocess, tempfile, time
from coherent_hot_compounds import height


def mul(a,b):
    z=0
    for _ in range(3):
        if b&1:z^=a
        b>>=1;a<<=1
        if a&8:a^=11
    return z


def axpy(x,a,y):
    z=dict(x)
    for k,v in y.items():
        z[k]=z.get(k,0)^mul(a,v)
        if not z[k]:del z[k]
    return z


def basis(vectors):
    pivots={}
    for v in vectors:
        v=dict(v)
        while v:
            k=max(v)
            if k in pivots:v=axpy(v,v[k],pivots[k])
            else:
                inv=(0,1,5,6,7,2,3,4)[v[k]];pivots[k]={i:mul(inv,a) for i,a in v.items()};break
    return pivots


def in_span(v,pivots):
    v=dict(v)
    while v:
        k=max(v)
        if k not in pivots:return False
        v=axpy(v,v[k],pivots[k])
    return True


def dot(x,y):
    result=0
    for k,v in x.items():result^=mul(v,y.get(k,0))
    return result


def validate_maps(embedding,retraction,modes):
    if not 1<=modes<=16 or not embedding or len(embedding)!=len(retraction):raise ValueError('Invalid map dimensions.')
    q=len(embedding);hs=[];charges=[]
    for p,l in zip(embedding,retraction):
        if not p or not l or any(not isinstance(w,int) or not 0<=w<1<<modes or (type(c)is not int or c not in range(1,8)) for v in (p,l) for w,c in v.items()):raise ValueError('Invalid sparse F8 vector.')
        words=set(p)|set(l);hh={height(w) for w in words};ch={w.bit_count() for w in words}
        if len(hh)!=1 or len(ch)!=1:raise ValueError('Maps must preserve height and occupation.')
        hs.append(hh.pop());charges.append(ch.pop())
    if any(dot(l,p)!=int(i==j) for i,l in enumerate(retraction) for j,p in enumerate(embedding)):raise ValueError('LP=I failed.')
    for vectors in (embedding,retraction):
        pivots=basis(vectors)
        if len(pivots)!=q:raise ValueError('Dependent map.')
        for v in vectors:
            for bit in range(modes):
                d={w^(1<<bit):c for w,c in v.items() if w>>bit&1}
                if not in_span(d,pivots):raise ValueError('Deletion closure failed.')
    return hs,charges



def tensor_vector(x,y,modes):return {a|(b<<modes):mul(c,d) for a,c in x.items() for b,d in y.items()}
