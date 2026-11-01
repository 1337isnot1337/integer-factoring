#!/usr/bin/env python3
"""Coherent six-vertex crossing with its complete monomial coefficient span."""
from __future__ import annotations

from collections import defaultdict
from math import comb

from bounded_spin_crossing import magnetic_step, path_norm

# Local tensor entries: a=f=d=1, b=e=t, c=w. The names refer to
# 0000,1111,1001; 0101,1010; and 0110, respectively.
VERTICES = (((0,0,0,0),0,0), ((1,1,1,1),0,0),
            ((1,0,0,1),0,0), ((0,1,0,1),1,0),
            ((1,0,1,0),1,0), ((0,1,1,0),0,1))


def monomial_entry(four_paths, prime=5, normalize=True):
    """Return coefficients in u^a v^b w^c, a+b+c=n/2, with u=1.

    Before the final weight-balance collapse, the finite-state contraction
    tracks t^(2b) w^c. No square root or division by a factorial is used.
    prime=None requests the unnormalized integer contraction.
    """
    n = len(four_paths[0])
    assert n % 2 == 0 and all(len(p) == n and sum(p) == 0 for p in four_paths)
    assert prime is not None or not normalize
    heights = [0]*4
    states = {((0,0,0,0),0,0):1}
    for i in range(n):
        steps = [p[i] for p in four_paths]
        following = defaultdict(int)
        for (magnetic, ell, c), value in states.items():
            for bits, dt, dw in VERTICES:
                target, scalar = [], value
                for h, s, r, bit in zip(heights, steps, magnetic, bits):
                    transition = magnetic_step(h,s,r,bit)
                    if transition is None:
                        scalar = 0
                        break
                    nr, coefficient = transition
                    target.append(nr)
                    scalar *= coefficient
                if scalar:
                    key = (tuple(target),ell+dt,c+dw)
                    following[key] += scalar
                    if prime is not None:
                        following[key] %= prime
        states = {key:value for key,value in following.items() if value}
        heights = [h+s for h,s in zip(heights,steps)]
    result = {}
    for (magnetic,ell,c),value in states.items():
        if magnetic != (0,0,0,0):
            continue
        assert ell % 2 == 0 and ell//2+c <= n//2
        result[n//2-ell//2-c,ell//2,c] = value
    if normalize:
        inverse = pow(path_norm(four_paths[2])*path_norm(four_paths[3]),-1,prime)
        result = {key:value*inverse % prime for key,value in result.items()}
    return {key:value for key,value in result.items() if value}


def unveil_coefficients(monomials, length, prime=5):
    """Substitute v=(1-z)^2, w=z^(n+1); disjoint windows expose all pieces."""
    result = [0]*(length//2*(length+1)+1)
    for (a,b,c),value in monomials.items():
        assert a+b+c == length//2
        for j in range(2*b+1):
            result[c*(length+1)+j] += value*comb(2*b,j)*(-1)**j
    return [value % prime for value in result] if prime is not None else result


def polynomial_entry(four_paths):
    return unveil_coefficients(monomial_entry(four_paths),len(four_paths[0]))


def integral_polynomial_entry(four_paths):
    n = len(four_paths[0])
    scale = 12**n
    denominator = path_norm(four_paths[2])*path_norm(four_paths[3])
    assert scale % denominator == 0
    raw = monomial_entry(four_paths,prime=None,normalize=False)
    return unveil_coefficients({m:c*(scale//denominator) for m,c in raw.items()},n,prime=None)
