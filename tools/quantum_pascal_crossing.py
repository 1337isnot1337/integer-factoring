#!/usr/bin/env python3
"""Quantum-Pascal trace crossings over F2[t], with complete polynomial entries."""
from __future__ import annotations

from weighted_smith_oracle import polynomial_multiply


def quantum_binomials(n):
    """Rows [a choose b]_t of the Gaussian binomial triangle, reduced modulo two."""
    rows = [[1]]
    for level in range(1, n+1):
        previous = rows[-1]
        rows.append([(previous[i-1] if i else 0)
                     ^ ((previous[i] << i) if i < len(previous) else 0)
                     for i in range(level+1)])
    return rows


def construct(alphabet, vacuum=None):
    """Return the sparse tensor Tr(O_xy O_uv), its degree, and ladder entries.

    O_xy = X^x diag(vacuum) D_t^[y], with divided lowering coefficients
    [k choose y]_t. The default vacuum is [alphabet choose i]_t.
    Tensor keys outside x+u=y+v are zero. Endpoint normalization is unnecessary.
    """
    if alphabet < 2:
        raise ValueError('the alphabet must have at least two letters')
    binomials = quantum_binomials(alphabet)
    vacuum = binomials[alphabet][:-1] if vacuum is None else list(vacuum)
    if len(vacuum) != alphabet or vacuum[0] != 1 or any(v < 0 for v in vacuum):
        raise ValueError('supply alphabet binary polynomials with vacuum[0]=1')
    ladder = {(k, y): polynomial_multiply(binomials[k][y], vacuum[k-y])
              for k in range(alphabet) for y in range(k+1)}
    tensor = {}
    for x in range(alphabet):
        for y in range(alphabet):
            for u in range(alphabet):
                v = x+u-y
                if not 0 <= v < alphabet:
                    continue
                value = 0
                for k in range(v, alphabet):
                    j = k-v+u
                    if y <= j < alphabet:
                        value ^= polynomial_multiply(ladder[k, v], ladder[j, y])
                tensor[x, y, u, v] = value
    degree = max(value.bit_length()-1 for value in tensor.values())
    return tensor, degree, ladder


def cut_blocks(tensor, alphabet, cut):
    """Exact square grade blocks, retaining original physical row/column pairs."""
    if cut not in (0, 1):
        raise ValueError('choose cut zero or one')
    grades = range(1-alphabet, alphabet) if cut == 0 else range(2*alphabet-1)
    blocks = []
    for grade in grades:
        rows = [(a, b) for a in range(alphabet) for b in range(alphabet)
                if (a-b if cut == 0 else a+b) == grade]
        columns = [(a, b) for a in range(alphabet) for b in range(alphabet)
                   if (b-a if cut == 0 else a+b) == grade]
        matrix = [[tensor.get((a, b, c, e) if cut == 0 else (a, c, b, e), 0)
                   for c, e in columns] for a, b in rows]
        blocks.append((grade, rows, columns, matrix))
    return blocks


def evaluate(polynomial, value, field):
    answer = 0
    for i in range(polynomial.bit_length()-1, -1, -1):
        answer = field.mul(answer, value) ^ ((polynomial >> i) & 1)
    return answer
