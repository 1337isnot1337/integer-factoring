#!/usr/bin/env python3
"""Exact extreme paired charges from integer height functions on a cube."""
from __future__ import annotations


def cube_edges(dimension):
    if not isinstance(dimension, int) or dimension < 1:
        raise ValueError('Use a positive integer cube dimension.')
    return tuple((state, state | (1 << bit)) for state in range(1 << dimension)
                 for bit in range(dimension) if not (state >> bit) & 1)


def height_functions(dimension):
    """Enumerate every anchored cube function with all edge differences +/-1."""
    cube_edges(dimension)
    order = sorted(range(1, 1 << dimension), key=lambda state: (state.bit_count(), state))
    values = [0]*(1 << dimension)

    def visit(position):
        if position == len(order):
            yield tuple(values)
            return
        state = order[position]
        predecessors = [values[state ^ (1 << bit)] for bit in range(dimension) if (state >> bit) & 1]
        choices = {predecessors[0]-1, predecessors[0]+1}
        for value in predecessors[1:]:
            choices.intersection_update((value-1, value+1))
        for value in sorted(choices):
            values[state] = value
            yield from visit(position+1)

    yield from visit(0)


def paired_charge_vertices(dimension):
    zero = (0,)*(1 << dimension)
    for height in height_functions(dimension):
        potential = tuple(2*value for value in height)
        yield potential, zero
        yield zero, potential


def parity_charge_vertex(dimension, signs):
    """Assign independent signs to odd vertices and zero to every even vertex.

    The returned potential already has edge-Lipschitz constant two. Its
    oscillation is at most four, so the strict typed lattice needs no
    endpoint pole larger than three, at any dimension.
    """
    cube_edges(dimension)
    odd = tuple(state for state in range(1 << dimension) if state.bit_count() % 2)
    if len(signs) != len(odd) or any(value not in (-1, 1) for value in signs):
        raise ValueError('Supply one sign, -1 or +1, for every odd cube vertex.')
    values = [0]*(1 << dimension)
    for state, sign in zip(odd, signs):
        values[state] = 2*sign
    return tuple(values)
