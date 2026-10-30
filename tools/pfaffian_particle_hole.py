#!/usr/bin/env python3
"""Exact particle-hole reduction and signed-charge Pfaffian cut ranks."""

from collections import defaultdict
from itertools import product

from check_modular_crossing import rank_two
from complement_shell_tensor import multiply


CHARGE = (0, -1, 1, 0)


def particle_hole_mask(parameter, state):
    """New bits: 1-y_low, 1-x_low, x_high, y_high."""
    return ((1 ^ (state & 1)) | ((1 ^ (parameter & 1)) << 1)
            | ((parameter & 2) << 1) | ((state & 2) << 2))


PARTICLE_HOLE = tuple(particle_hole_mask(x, y) for x, y in product(range(4), repeat=2))
INVERSE_PARTICLE_HOLE = tuple(PARTICLE_HOLE.index(mask) for mask in range(16))


def one_particle_matrix(p, q, r):
    """G=h tensor g, with det(g)=det(h)=p+r."""
    g = ((r, p), (1, 1))
    h = ((1, 1), (q, p ^ q ^ r))
    return tuple(tuple(multiply(h[a >> 1][b >> 1], g[a & 1][b & 1])
                       for b in range(4)) for a in range(4))


def local_exterior_table(p, q, r):
    """The universal local Pfaffian, as a polynomial formula even at p=r.

    Return the standard flattened (x,y,u,v) table. The one-particle and
    three-particle blocks are G; the middle block is Lambda^2(G)/(p+r).
    The displayed middle formula removes the apparent denominator.
    """
    delta, s = p ^ r, p ^ q ^ r
    p2, q2, r2 = (multiply(value, value) for value in (p, q, r))
    pr, qr, pq = multiply(p, r), multiply(q, r), multiply(p, q)
    six = (
        (1, 0, 1, 1, 0, 1),
        (0, r2, pr, pr, p2, 0),
        (q, r, r ^ q, p ^ q, p, s),
        (q, r, p ^ q, r ^ q, p, s),
        (0, 1, 1, 1, 1, 0),
        (q2, 0, qr ^ q2 ^ pq, qr ^ q2 ^ pq, 0, p2 ^ q2 ^ r2),
    )
    matrix = [[0] * 16 for _ in range(16)]
    matrix[0][0] = matrix[15][15] = delta
    g = one_particle_matrix(p, q, r)
    one = (1, 2, 4, 8)
    three = (7, 11, 13, 14)
    for a, b in product(range(4), repeat=2):
        matrix[one[a]][one[b]] = g[a][b]
        matrix[three[a]][three[b]] = g[a][b]
    two = (3, 5, 9, 6, 10, 12)
    for a, b in product(range(6), repeat=2):
        matrix[two[a]][two[b]] = six[a][b]
    return tuple(matrix[particle_hole_mask(x, y)][particle_hole_mask(v, u)]
                 for x, y, u, v in product(range(4), repeat=4))


def signed_charge_code(label):
    """Base-five digits charge+1; sums encode row blocks without carries."""
    result = 0
    for mask in reversed(label):
        result = 5 * result + CHARGE[mask] + 1
    return result


def compile_signed_blocks(parameters, states=None):
    parameters = tuple(parameters)
    states = parameters if states is None else tuple(states)
    if not parameters or not states:
        raise ValueError("both label collections must be nonempty")
    pairs = len(parameters[0])
    if any(len(label) != pairs for label in parameters + states):
        raise ValueError("all labels must use the same complementary mode pairs")
    if len(set(parameters)) != len(parameters) or len(set(states)) != len(states):
        raise ValueError("labels must be distinct")
    parameter_codes = [signed_charge_code(label) for label in parameters]
    state_codes = [signed_charge_code(label) for label in states]
    row_blocks, col_blocks = defaultdict(list), defaultdict(list)
    for x, a in enumerate(parameter_codes):
        for y, b in enumerate(state_codes):
            row_blocks[a + b].append((x, y))
            col_blocks[a - b].append((x, y))
    rows = tuple(row_blocks[key] for key in sorted(row_blocks))
    cols = tuple(col_blocks[key] for key in sorted(col_blocks))
    col_positions = [[None] * len(states) for _ in parameters]
    for block, labels in enumerate(cols):
        for position, (x, y) in enumerate(labels):
            col_positions[x][y] = block, position
    return {
        "parameters": parameters, "states": states, "row_blocks": rows,
        "column_blocks": cols, "column_positions": col_positions,
        "evaluated_entries": sum(len(block) ** 2 for block in rows),
        "dense_entries": len(parameters) ** 2 * len(states) ** 2,
    }


def coefficient_blocks(code, p=6, q=3, r=2):
    local = local_exterior_table(p, q, r)
    parameters, states = code["parameters"], code["states"]
    rows, cols = {}, {}
    cache = {}
    for block, labels in enumerate(code["row_blocks"]):
        for row, (x, y) in enumerate(labels):
            for col, (v, u) in enumerate(labels):
                factors = tuple(sorted(local[((a * 4 + b) * 4 + c) * 4 + d]
                                       for a, b, c, d in zip(parameters[x], states[y],
                                                             states[u], parameters[v])))
                if factors and factors[0] == 0:
                    continue
                value = cache.get(factors)
                if value is None:
                    value = 1
                    for factor in factors:
                        value = multiply(value, factor)
                    cache[factors] = value
                cblock, crow = code["column_positions"][x][u]
                cblock_check, ccol = code["column_positions"][v][y]
                assert cblock == cblock_check
                while value:
                    bit = value & -value
                    degree = bit.bit_length() - 1
                    if degree not in rows:
                        rows[degree] = [[0] * len(part) for part in code["row_blocks"]]
                        cols[degree] = [[0] * len(part) for part in code["column_blocks"]]
                    rows[degree][block][row] ^= 1 << col
                    cols[degree][cblock][crow] ^= 1 << ccol
                    value ^= bit
    degrees = sorted(rows)
    return degrees, [rows[k] for k in degrees], [cols[k] for k in degrees]


def rank_profile(parameters, states=None, p=6, q=3, r=2):
    code = compile_signed_blocks(parameters, states)
    degrees, rows, cols = coefficient_blocks(code, p, q, r)
    return {
        "parameter_labels": len(code["parameters"]), "state_labels": len(code["states"]),
        "row_block_count": len(code["row_blocks"]),
        "column_block_count": len(code["column_blocks"]),
        "largest_row_block": max(map(len, code["row_blocks"])),
        "largest_column_block": max(map(len, code["column_blocks"])),
        "evaluated_entries": code["evaluated_entries"], "dense_entries": code["dense_entries"],
        "profile": [[degree, sum(rank_two(block) for block in first),
                     sum(rank_two(block) for block in second)]
                    for degree, first, second in zip(degrees, rows, cols)],
    }
