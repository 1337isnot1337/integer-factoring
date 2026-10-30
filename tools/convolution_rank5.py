"""Exact rank of convolution by a 32-bit kernel on F2^5."""

from functools import lru_cache
from itertools import combinations


def anf_from_group(kernel, width=5):
    for j in range(width):
        step = 1 << j
        lower = ((1 << (1 << width))-1) // ((1 << step)+1)
        kernel ^= (kernel >> step) & lower
    return kernel


def group_from_anf(anf, width=5):
    # The upper-subset zeta transform is its own inverse in characteristic two.
    return anf_from_group(anf,width)


def alternating_rank_from_pattern(pattern, width=5):
    rows = [0]*width
    for bit,(i,j) in enumerate(combinations(range(width),2)):
        if pattern >> bit & 1:
            rows[i] |= 1 << j
            rows[j] |= 1 << i
    rank = 0
    for col in range(width):
        pivot = next((i for i in range(rank,width) if rows[i] >> col & 1),None)
        if pivot is None:
            continue
        rows[rank],rows[pivot] = rows[pivot],rows[rank]
        for i in range(rank+1,width):
            if rows[i] >> col & 1:
                rows[i] ^= rows[rank]
        rank += 1
    return rank


QUAD_MASKS = tuple((1 << i) | (1 << j) for i,j in combinations(range(5),2))
CUBIC_MASKS = tuple(31 ^ mask for mask in QUAD_MASKS)
ALT_RANKS = tuple(alternating_rank_from_pattern(i) for i in range(1024))
LINEAR_BITS = sum(1 << (1 << j) for j in range(5))
QUARTIC_BITS = sum(1 << (31 ^ (1 << j)) for j in range(5))


@lru_cache(maxsize=65536)
def multiply3(left,right):
    result = 0
    for a in range(8):
        if left >> a & 1:
            for b in range(8):
                if right >> b & 1 and not a & b:
                    result ^= 1 << (a | b)
    return result


def split_pair(anf, pair):
    bit_i = pair & -pair
    bit_j = pair ^ bit_i
    remaining = [1 << j for j in range(5) if (1 << j) not in (bit_i,bit_j)]
    pieces = [0]*4
    for small in range(8):
        original = sum(bit for j,bit in enumerate(remaining) if small >> j & 1)
        for flags in range(4):
            index = original | (bit_i if flags & 1 else 0) | (bit_j if flags & 2 else 0)
            pieces[flags] |= ((anf >> index) & 1) << small
    a,b,c,d = pieces
    assert d & 1
    return a ^ multiply3(multiply3(b,c),d)


def rank_anf5(anf):
    if anf & 1:
        return 32
    if anf & LINEAR_BITS:
        return 16
    quadratic = sum(((anf >> mask) & 1) << bit for bit,mask in enumerate(QUAD_MASKS))
    if quadratic:
        if ALT_RANKS[quadratic] == 4:
            return 12
        pair = QUAD_MASKS[(quadratic & -quadratic).bit_length()-1]
        residual = split_pair(anf,pair)
        assert residual in (0,128)
        return 8 + 2 * int(bool(residual))
    cubic = sum(((anf >> mask) & 1) << bit for bit,mask in enumerate(CUBIC_MASKS))
    if cubic:
        return 2 + ALT_RANKS[cubic]
    if anf & QUARTIC_BITS:
        return 2
    return int(bool(anf))


def convolution_rank5(kernel):
    return rank_anf5(anf_from_group(kernel))
