#!/usr/bin/env python3
"""Finite polynomial representatives of regular F4 row-condition frames."""
from __future__ import annotations

from itertools import product

from f4_lattice_linear import MUL
from f4_repaired_identity_filters import INV,add,multiply
from search_binary_polynomial_flags import constant_row


def positions(targets,row):
    if tuple(targets)!=tuple(sorted(targets)):
        raise ValueError('Row targets must be in nondecreasing order.')
    return tuple((j,k) for j in range(row) for k in range(1,targets[row]-targets[j]))


def dimensions(targets):
    return tuple(len(positions(targets,row)) for row in range(len(targets)))


def row_catalogue(targets,row,field_order=4,base=None):
    if field_order not in (2,4):
        raise ValueError('Use F2 or F4 coefficients.')
    size=len(targets)
    if base is None:
        base=tuple(tuple(int(i==j) for j in range(size)) for i in range(size))
    cells=positions(targets,row)
    for values in product(range(field_order),repeat=len(cells)):
        vector=[list(pair) for pair in constant_row(base[row])]
        for (j,k),scalar in zip(cells,values):
            if not scalar:continue
            for column,entry in enumerate(base[j]):
                value=MUL[scalar][entry]
                vector[column][0]^=(value&1)<<k
                vector[column][1]^=(value>>1)<<k
        yield tuple(map(tuple,vector))


def truncated(value,precision):
    if any(power<0 for power in value):
        raise ValueError('A regular series may not have negative powers.')
    return {power:coefficient for power,coefficient in value.items()
            if power<precision and coefficient}


def unit_inverse(value,precision):
    constant=value.get(0,0)
    if not constant:
        raise ValueError('The formal diagonal pivot is not a unit.')
    result={0:INV[constant]}
    for power in range(1,precision):
        coefficient=0
        for k in range(1,power+1):
            coefficient^=MUL[value.get(k,0)][result.get(power-k,0)]
        coefficient=MUL[INV[constant]][coefficient]
        if coefficient:result[power]=coefficient
    return result


def normal_form(matrix,targets,precision=None):
    """Normalize W(0)=I modulo terms irrelevant to the target conditions.

    The input and output use dictionary polynomials. Only degrees below
    max(targets)-min(targets) can occur in the canonical representative.
    """
    size=len(targets)
    positions(targets,0)
    needed=max(1,max(targets)-min(targets))
    if precision is None:precision=needed
    if precision<needed:raise ValueError('Insufficient formal precision.')
    if len(matrix)!=size or any(len(row)!=size for row in matrix):
        raise ValueError('Expected a square frame with one row per target.')
    if any(matrix[i][j].get(0,0)!=int(i==j) for i in range(size) for j in range(size)):
        raise ValueError('Remove the chosen constant frame before normalization.')
    rows=[[truncated(entry,precision) for entry in row] for row in matrix]
    def mul(a,b):return truncated(multiply(a,b),precision)
    # Reverse Gaussian elimination uses only upper row operations and
    # diagonal units, all of which preserve sorted row conditions.
    for j in range(size-1,-1,-1):
        inverse=unit_inverse(rows[j][j],precision)
        rows[j]=[mul(inverse,entry) for entry in rows[j]]
        for i in range(j):
            coefficient=rows[i][j]
            if coefficient:
                rows[i]=[add(a,mul(coefficient,b)) for a,b in zip(rows[i],rows[j])]
    # Removing the high tail in column j can only change earlier columns.
    for i in range(1,size):
        for j in range(i-1,-1,-1):
            gap=targets[i]-targets[j]
            tail={power:value for power,value in rows[i][j].items() if power>=gap}
            if tail:rows[i]=[add(a,mul(tail,b)) for a,b in zip(rows[i],rows[j])]
    assert all(rows[i][i]=={0:1} for i in range(size))
    assert all(not rows[i][j] for i in range(size) for j in range(i+1,size))
    assert all(all(0<power<targets[i]-targets[j] for power in rows[i][j])
               for i in range(size) for j in range(i))
    return rows
