"""Constant-size local description of the bounded-weight crossing profiles.

Polynomial integers encode coefficients over F2. The derivation and the
permutation reduction are in research/hamming-tensor-power.md.
"""

from functools import lru_cache
from itertools import product

from modular_pfaffian import polynomial_product


# Vertex blocks are U,A,B,Y, two vertices each. The low bit of every part is 1.
LOW_FACTORS = {(0, 1): 6, (0, 2): 3, (0, 3): 2,
               (1, 2): 2, (1, 3): 3, (2, 3): 6}


def matchings(vertices):
    if not vertices:
        yield ()
        return
    first = vertices[0]
    for other in vertices[1:]:
        if first//2 == other//2:
            continue
        remaining = tuple(v for v in vertices[1:] if v != other)
        for tail in matchings(remaining):
            yield ((first, other),)+tail


def local_edge(i, j, bits):
    block, sheet = divmod(i, 2)
    other_block, other_sheet = divmod(j, 2)
    # Input bits have order x,y,u,v; the vertex blocks have order U,A,B,Y.
    values = (bits[2], bits[0], bits[3], bits[1])
    p, q = values[block] ^ sheet, values[other_block] ^ other_sheet
    pair = block, other_block
    if pair == (0, 1):
        return 1 << (p | q)
    if pair == (0, 2):
        return (1 if q == 0 else 3) if q <= p else 0
    if pair == (0, 3):
        return (1 << q) if p == q else 0
    if pair == (1, 2):
        return (1 << p) if p <= q else 0
    if pair == (1, 3):
        return 1 << (p & q)
    if pair == (2, 3):
        return polynomial_product(1 << q, 1 if p == 0 else 3) if not p & q else 0
    raise ValueError(pair)


def substitute(polynomial, weight):
    answer = 0
    while polynomial:
        bit = polynomial & -polynomial
        answer |= 1 << (weight*(bit.bit_length()-1))
        polynomial ^= bit
    return answer


@lru_cache(maxsize=1)
def local_terms():
    """Return (local 4-leg polynomial tensor, low-bit scalar) pairs."""
    terms = {}
    for matching in matchings(tuple(range(8))):
        scalar = 1
        for i, j in matching:
            scalar = polynomial_product(scalar, LOW_FACTORS[i//2, j//2])
        local = []
        for bits in product(range(2), repeat=4):
            entry = 1
            for i, j in matching:
                entry = polynomial_product(entry, local_edge(i, j, bits))
            local.append(entry)
        local = tuple(local)
        # The fixed leading bit has input bits 0000, for every label.
        if local[0]:
            terms[local] = terms.get(local, 0) ^ scalar
    return tuple((local, scalar) for local, scalar in terms.items() if scalar)


def tensor_profile(width, indices, weights=None):
    if weights is None:
        weights = (1,)*(width+2)
    assert width >= 1 and len(indices) == 4
    assert len(weights) == width+2 and weights[0] == 1 and min(weights) >= 1
    assert all(0 <= a < 1 << width for a in indices)
    answer = 0
    for local, scalar in local_terms():
        value = polynomial_product(scalar, substitute(local[0], weights[-1]))
        for bit in range(width):
            x, y, u, v = ((a >> bit) & 1 for a in indices)
            value = polynomial_product(value, substitute(
                local[(2*x+y)*4+2*u+v], weights[bit+1]))
        answer ^= value
    return answer


def permuted_label(width, label, permutation):
    """Permute all non-low bits, then choose the leading-zero representative.

    permutation[j] is the old bit placed in new non-low position j.
    """
    assert sorted(permutation) == list(range(width+1))
    changed = sum(((label >> old) & 1) << j for j, old in enumerate(permutation))
    if (changed >> width) & 1:
        changed ^= (1 << (width+1))-1
    return changed
