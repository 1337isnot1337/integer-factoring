#!/usr/bin/env python3
"""Exact rank factorizations and whole-tensor coefficient descent controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

ROOT = Path(__file__).resolve().parent.parent


class Field:
    def __init__(self, modulus):
        self.modulus = modulus
        self.degree = modulus.bit_length() - 1
        self.order = 1 << self.degree

    def mul(self, a, b):
        value = 0
        while b:
            if b & 1:
                value ^= a
            b >>= 1
            a <<= 1
            if a & self.order:
                a ^= self.modulus
        return value

    def power(self, a, n):
        value = 1
        while n:
            if n & 1:
                value = self.mul(value, a)
            a = self.mul(a, a)
            n >>= 1
        return value

    def rank(self, matrix):
        a = [row[:] for row in matrix]
        rank = 0
        for j in range(len(a[0])):
            i = next((i for i in range(rank, len(a)) if a[i][j]), None)
            if i is None:
                continue
            a[rank], a[i] = a[i], a[rank]
            inverse = self.power(a[rank][j], self.order-2)
            a[rank] = [self.mul(v, inverse) for v in a[rank]]
            for i in range(rank+1, len(a)):
                coefficient = a[i][j]
                if coefficient:
                    a[i] = [v ^ self.mul(coefficient, u) for v, u in zip(a[i], a[rank])]
            rank += 1
        return rank

    def matmul(self, left, right):
        answer = [[0]*len(right[0]) for _ in left]
        for i, row in enumerate(left):
            for k, a in enumerate(row):
                for j, b in enumerate(right[k]):
                    answer[i][j] ^= self.mul(a, b)
        return answer


def factor_controls(field):
    rng = Random(64418 + field.order)
    base = Field(0b11)
    d = field.degree
    gram = [[field.mul(1 << i, 1 << j) & 1 for j in range(d)] for i in range(d)]
    assert base.rank(gram) == d
    sharp = [[field.mul(1 << i, 1 << j) for j in range(d)] for i in range(d)]
    assert field.rank(sharp) == 1
    cases = 0
    for rows, columns, inner in product((2,4,7), (3,5), (1,2,3)):
        for _ in range(3):
            left = [[rng.randrange(field.order) for _ in range(inner)] for _ in range(rows)]
            right = [[rng.randrange(field.order) for _ in range(columns)] for _ in range(inner)]
            matrix = field.matmul(left, right)
            projected = [[v & 1 for v in row] for row in matrix]
            coordinates = [[(v >> j) & 1 for v in row for j in range(d)] for row in left]
            expanded = [[field.mul(1 << j, right[k][c]) & 1 for c in range(columns)]
                        for k in range(inner) for j in range(d)]
            assert base.matmul(coordinates, expanded) == projected
            assert base.rank(projected) <= d*field.rank(matrix)
            cases += 1
    return dict(field=field.order, factorizations=cases, sharp_rank_increase=[1,d])


def tensor_controls():
    field, base = Field(0b111), Field(0b11)
    rng = Random(7319)
    local = []
    for v,u,y,x in product(range(2), repeat=4):
        a = rng.randrange(4)
        local.append([a, a ^ int(x == v and y == u)])
    records = []
    for length in (1,2,3):
        q = 1 << length
        coefficients = [[0]*q**4 for _ in range(length+1)]
        for v,u,y,x in product(range(q), repeat=4):
            polynomial = [1]
            for bit in range(length):
                index = ((x >> bit) & 1) + 2*((y >> bit) & 1) + 4*((u >> bit) & 1) + 8*((v >> bit) & 1)
                factors = local[index]
                result = [0]*(len(polynomial)+1)
                for i,a in enumerate(polynomial):
                    for j,b in enumerate(factors):
                        result[i+j] ^= field.mul(a,b)
                polynomial = result
            index = x+q*(y+q*(u+q*v))
            for degree,value in enumerate(polynomial):
                coefficients[degree][index] = value
            total = 0
            for value in polynomial:
                total ^= value & 1
            assert total == int(x == v and y == u)
        ranks = []
        for degree,values in enumerate(coefficients):
            cuts = []
            for cut in (0,1):
                matrix = [[0]*q**2 for _ in range(q**2)]
                for v,u,y,x in product(range(q), repeat=4):
                    row, column = (x+q*y,u+q*v) if cut == 0 else (x+q*u,y+q*v)
                    matrix[row][column] = values[x+q*(y+q*(u+q*v))]
                r_extension = field.rank(matrix)
                r_base = base.rank([[a & 1 for a in row] for row in matrix])
                assert r_base <= 2*r_extension
                cuts.append([r_extension,r_base])
            ranks.append(dict(degree=degree,rank_pairs=cuts))
        records.append(dict(length=length,alphabet=q,endpoint_entries=q**4,coefficients=ranks))
    assert (field.mul(2,2) & 1) != ((2 & 1)*(2 & 1))
    assert 1 ^ field.power(1,2) == 0
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    evidence = dict(status='exact controls for fixed-extension rank descent',
                    matrix_controls=[factor_controls(Field(m)) for m in (0b111,0b1011)],
                    tensor_controls=tensor_controls(),
                    nonmultiplicative_projection_checked=True,
                    trace_of_one_zero_control_checked=True,
                    source_sha256={str(Path(__file__).relative_to(ROOT)):sha256(Path(__file__).read_bytes()).hexdigest()})
    text = json.dumps(evidence,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text)


if __name__ == '__main__':
    main()
