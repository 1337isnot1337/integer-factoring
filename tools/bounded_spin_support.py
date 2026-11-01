#!/usr/bin/env python3
"""Physical word support and explicit factors for the extreme spin coefficient."""
from __future__ import annotations

from math import comb

from bounded_spin_crossing import path_entry, path_norm


def supported_word(length,word):
    height = 0
    for i in range(length):
        height += 1-2*((word >> i)&1)
        if not -3 <= height <= 3:
            return False
    return height == 0


def suffix_counts(length):
    counts = [[int(h == 0) for h in range(-3,4)]]
    for _ in range(length):
        counts.append([sum(counts[-1][j+3] for j in (h-1,h+1) if -3 <= j <= 3)
                       for h in range(-3,4)])
    return counts


def orbit_count(length):
    assert length > 0 and length % 2 == 0
    return suffix_counts(length-1)[-1][4]  # First bit zero gives height one.


def orbit_word(length,index):
    counts = suffix_counts(length-1)
    assert 0 <= index < counts[-1][4]
    word,height = 0,1
    for i in range(1,length):
        remaining = length-i-1
        zero = counts[remaining][height+4] if height < 3 else 0
        if index < zero:
            height += 1
        else:
            index -= zero
            word |= 1 << i
            height -= 1
            assert height >= -3
    assert height == 0 and index == 0
    return word


def unrank_subset(size,weight,index):
    assert 0 <= weight <= size and 0 <= index < comb(size,weight)
    result = []
    for position in range(size):
        chosen = comb(size-position-1,weight-1) if weight else 0
        if index < chosen:
            result.append(position)
            weight -= 1
        else:
            index -= chosen
    assert weight == 0 and index == 0
    return result


def corner_channel_count(length,b,c):
    assert b >= 0 and c >= 0 and b+c <= length//2
    return orbit_count(length)*comb(length,2*(b+c))*comb(2*(b+c),2*c)


def corner_words(length,b,c,channel):
    assert 0 <= channel < corner_channel_count(length,b,c)
    choices = comb(2*(b+c),2*c)
    cross_index,channel = channel % choices,channel//choices
    choices = comb(length,2*(b+c))
    positions_index,word_index = channel % choices,channel//choices
    positions = unrank_subset(length,2*(b+c),positions_index)
    cross = [positions[i] for i in unrank_subset(len(positions),2*c,cross_index)]
    j,l = sum(1 << i for i in positions),sum(1 << i for i in cross)
    word = orbit_word(length,word_index)
    return word,word^j,word^l,word^j^l


def vertex_factor(first_path,second_path,channel,dual_first=False,dual_second=False,scale=1):
    n = len(first_path)
    assert len(second_path) == n and channel >= 0
    if channel >= orbit_count(n):
        return 0
    word = orbit_word(n,channel)
    value = scale*path_entry(first_path,word)*path_entry(second_path,word) % 5
    if dual_first:
        value = value*pow(path_norm(first_path),-1,5) % 5
    if dual_second:
        value = value*pow(path_norm(second_path),-1,5) % 5
    return value
