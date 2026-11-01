#!/usr/bin/env python3
"""Thirty-state orthogonal path projector; sixty states for binary selection."""
from __future__ import annotations

from collections import defaultdict

from bounded_spin_crossing import binary_path, magnetic_step

STATES = tuple((h,r,s) for h in range(4) for r in range(h+1) for s in range(h+1))


def local_transitions(state, left_bit, right_bit):
    height,r,s = state
    for step in (1,-1):
        if not 0 <= height+step <= 3:
            continue
        left = magnetic_step(height,step,r,left_bit)
        right = magnetic_step(height,step,s,right_bit)
        if left is None or right is None:
            continue
        nr,a = left
        ns,b = right
        scalar = a*b % 5
        if step == -1:
            scalar = scalar*pow(height*(height+1),-1,5) % 5
        if scalar:
            yield step,(height+step,nr,ns),scalar


def projector_entry(length, left_word, right_word, last_path=None):
    """Sum v_P(left)v_P(right)/norm(P) in F5 over paths up to last_path.

    Without a cutoff all height-three return paths are included. With a
    cutoff, the path order is up first, as in the public unranking routine.
    """
    assert 0 <= left_word < 1 << length and 0 <= right_word < 1 << length
    if last_path is not None:
        assert len(last_path) == length and sum(last_path) == 0
    # comparison=0 means equal to the cutoff prefix; 1 means earlier.
    states = {((0,0,0),0):1}
    for i in range(length):
        following = defaultdict(int)
        for (state,comparison),value in states.items():
            for step,target,scalar in local_transitions(
                    state,(left_word >> i)&1,(right_word >> i)&1):
                nc = comparison
                if last_path is not None and comparison == 0 and step != last_path[i]:
                    if step == -1:  # Down follows up, so this path exceeds the cutoff.
                        continue
                    nc = 1
                key = target,nc
                following[key] = (following[key]+value*scalar) % 5
        states = {key:value for key,value in following.items() if value}
    return sum(states.get(((0,0,0),comparison),0) for comparison in (0,1)) % 5


def binary_projector_entry(width, left_word, right_word):
    last = binary_path(width,(1 << width)-1)
    return projector_entry(len(last),left_word,right_word,last)
