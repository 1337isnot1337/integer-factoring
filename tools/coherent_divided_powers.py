#!/usr/bin/env python3
"""Graded dimensions of paired divided-power derivative spaces in char two."""
from collections import Counter
from itertools import combinations
from math import comb

from check_binary_coefficient_entropy import inclusion_formula


def disjoint_rank(n, a, b):
    """Rank of [A intersect B is empty], for |A|=a and |B|=b."""
    if not 0 <= a <= n or not 0 <= b <= n or a+b > n:
        raise ValueError('Use nonnegative subset sizes with a+b <= n.')
    return inclusion_formula(n, min(a, b), n-max(a, b))


def paired_top(modes, power):
    if modes < 2 or modes % 2 or not 0 <= power <= modes//2:
        raise ValueError('Use a positive even number of modes and a valid divided power.')
    pairs = [(1 << i) | (1 << (modes-1-i)) for i in range(modes//2)]
    return {sum(chosen): 1 for chosen in combinations(pairs, power)}


def graded_dimensions(modes, power):
    """Counts by (occupation,height), without enumerating physical states."""
    if modes < 2 or modes % 2 or not 0 <= power <= modes//2:
        raise ValueError('Use a positive even number of modes and a valid divided power.')
    pairs = modes//2
    # broken[b][h] counts choices of b pairs and one surviving mode per pair.
    broken = [Counter() for _ in range(power+1)]
    broken[0][0] = 1
    for i in range(pairs):
        for b in range(min(power, i+1), 0, -1):
            for h, multiplicity in broken[b-1].items():
                broken[b][h+2*i+1] += multiplicity
                broken[b][h+2*(modes-1-i)+1] += multiplicity
    result = Counter()
    for b in range(power+1):
        for deleted_pairs in range(power-b+1):
            surviving_pairs = power-b-deleted_pairs
            rank = disjoint_rank(pairs-b, deleted_pairs, surviving_pairs)
            for h, multiplicity in broken[b].items():
                result[(b+2*surviving_pairs, h+2*modes*surviving_pairs)] += rank*multiplicity
    return dict(sorted((key, value) for key, value in result.items() if value))
