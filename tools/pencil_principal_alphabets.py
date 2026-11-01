#!/usr/bin/env python3
"""Matched principal restrictions of the strict five-letter F4 pencil."""
from itertools import product


def restrict(values, natural, alphabet_x, alphabet_y):
    ax, ay = list(alphabet_x), list(alphabet_y)
    for alphabet in (ax,ay):
        if not 2 <= len(alphabet) <= 5 or len(set(alphabet)) != len(alphabet) or any(type(v) is not int or not 0<=v<5 for v in alphabet):
            raise ValueError('Use two through five distinct original letters in each pair.')
    if len(values) != 625 or len(natural) != 33 or natural[-1] != 1:
        raise ValueError('Use the five-letter pencil and its unit-scale natural heights.')
    # Native tensor order is x+A*(y+B*(u+B*v)).
    restricted = [values[x+5*(y+5*(u+5*v))] for v,u,y,x in product(ax,ay,ay,ax)]
    alphabets = (ax,ay,ay,ax)
    groups = []
    for place in range(2):
        side = []
        for leg,alphabet in enumerate(alphabets):
            original = [0]+list(natural[place*16+4*leg:place*16+4*(leg+1)])
            chosen = [original[i]-original[alphabet[0]] for i in alphabet]
            side.append(chosen)
        groups.append(side)
    parameters = [v for side in groups for group in side for v in group[1:]]+[1]
    gx = [sum(groups[p][leg][i] for p in range(2) for leg in (0,3)) for i in range(len(ax))]
    gy = [sum(groups[p][leg][i] for p in range(2) for leg in (1,2)) for i in range(len(ay))]
    return {'q_x':len(ax),'q_u':len(ay),'values':restricted},parameters,gx,gy


def pressure_key(profiles, gx, gy):
    """Exact equivalence of inherited scalar curves, not of new free gauges."""
    minima = [min(p) for p in profiles]
    orders = sorted(tuple(v-m for v in p) for p,m in zip(profiles,minima))
    physical = sorted(tuple(sorted(v-min(g) for v in g)) for g in (gx,gy))
    gap = min(gx)+min(gy)-sum(minima)
    return tuple(orders),tuple(physical),gap
