#!/usr/bin/env python3
"""Finite controls for the partial-transpose rank and power-rank inequalities.

No search result is used as a proof of the arbitrary-field theorem.
"""
import argparse
import json
import random
from pathlib import Path


def rank_mod(matrix, p):
    a = [row[:] for row in matrix]
    r = 0
    for j in range(len(a[0]) if a else 0):
        pivot = next((i for i in range(r, len(a)) if a[i][j] % p), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(int(a[r][j]), -1, p)
        a[r] = [(v*inv) % p for v in a[r]]
        for i in range(r+1, len(a)):
            if a[i][j]:
                c = a[i][j]
                a[i] = [(x-c*y) % p for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def binary_rank(rows):
    pivots = {}
    for v in rows:
        while v:
            k = v.bit_length()-1
            if k not in pivots:
                pivots[k] = v
                break
            v ^= pivots[k]
    return len(pivots)


def partial_transpose(m, q):
    return [[m[q*x+v][q*u+y] for u in range(q) for v in range(q)]
            for x in range(q) for y in range(q)]


def realignment(m, q):
    return [[m[q*x+y][q*u+v] for y in range(q) for v in range(q)]
            for x in range(q) for u in range(q)]


def phi(m, q):
    return [[m[q*x+u][q*v+y] for u in range(q) for v in range(q)]
            for x in range(q) for y in range(q)]


def multiply(a, b, p):
    n = len(a)
    return [[sum(a[i][k]*b[k][j] for k in range(n)) % p
             for j in range(n)] for i in range(n)]


def check_identity(m, q):
    pt = partial_transpose(m, q)
    assert phi(m, q) == [[row[q*v+u] for u in range(q) for v in range(q)]
                         for row in pt]


def binary_q2_exhaustive():
    histogram = {}
    invertible = regular = 0
    for mask in range(1 << 16):
        m = [[(mask >> (4*i+j)) & 1 for j in range(4)] for i in range(4)]
        def br(a):
            return binary_rank([sum(v << j for j, v in enumerate(row)) for row in a])
        ordinary, osr, pt = br(m), br(realignment(m, 2)), br(partial_transpose(m, 2))
        assert pt <= ordinary*osr
        assert ordinary <= pt*osr
        check_identity(m, 2)
        histogram[(ordinary, osr, pt)] = histogram.get((ordinary, osr, pt), 0)+1
        if ordinary != 4:
            continue
        invertible += 1
        square = multiply(m, m, 2)
        r2 = br(phi(square, 2))
        assert br(realignment(square, 2)) <= pt*pt
        assert pt*pt*r2 >= 4
        s = [[v ^ (i == j) for j, v in enumerate(row)] for i, row in enumerate(m)]
        s2 = multiply(s, s, 2)
        s3 = multiply(s2, s, 2)
        s4 = multiply(s2, s2, 2)
        if not any(map(any, s4)) and any(map(any, s3)):
            regular += 1
    assert invertible == 20160 and regular == 2520
    return dict(all_binary_operators=65536, rank_triples=len(histogram),
                invertible_power_bounds_checked=invertible,
                regular_unipotent_controls=regular)


def finite_field_controls():
    rng = random.Random(20261030)
    output = []
    for p, q, trials in [(2,3,300), (2,4,200), (3,2,300), (3,3,300), (5,3,100)]:
        powers = 0
        for _ in range(trials):
            n = q*q
            m = [[rng.randrange(p) for _ in range(n)] for _ in range(n)]
            ordinary = rank_mod(m, p)
            osr = rank_mod(realignment(m, q), p)
            pt = rank_mod(partial_transpose(m, q), p)
            assert pt <= ordinary*osr and ordinary <= pt*osr
            check_identity(m, q)
            if ordinary == n:
                square = multiply(m, m, p)
                cube = multiply(square, m, p)
                r2, r3 = rank_mod(phi(square, q), p), rank_mod(phi(cube, q), p)
                assert rank_mod(realignment(square, q), p) <= pt*pt
                assert pt*pt*r2 >= n and pt*r2*r3 >= n
                powers += 1
        output.append(dict(prime=p, q=q, matrices=trials, invertible_power_controls=powers))
    return output


def sharp_square_control():
    q, n = 4, 16
    m = [[0]*n for _ in range(n)]
    for j in range(n):
        output = ((j << 1) & 15) | (j >> 3)
        m[output][j] = 1
    square = multiply(m, m, 2)
    transpose = [[int(i == q*(j % q)+j//q) for j in range(n)] for i in range(n)]
    assert square == transpose
    r, r2 = rank_mod(phi(m, q), 2), rank_mod(phi(square, q), 2)
    assert (r, r2) == (4, 1) and r*r*r2 == n
    return dict(q=q, first_choi_rank=r, second_choi_rank=r2,
                square_bound_is_equality=True,
                stronger_product_r_times_r2_bound_is_false=True)


def regular_rank_two_threshold():
    # B >= 2n+(n-1)(n/4-4) = n^2/4-9n/4+4.
    # For n>=25 its ratio to n^2 is at least 104/625.
    # 104>25^(10/7)>25^sqrt(2), certified with integers.
    assert 100 > 2*49
    assert 104**7 > 25**10
    rows = []
    for q in range(5, 33):
        n = q*q
        k = max(0, (n+3)//4-4)
        b = 2*n+(n-1)*k
        assert 4*b >= n*n-9*n+16
        assert 625*(n*n-9*n+16) >= 416*n*n
        assert b**7 > n**10
        rows.append(dict(q=q, zero_order_count=2, one_order_count_floor=k,
                         smith_cost_floor=b))
    return dict(all_q_at_least=5, base_integer_comparison=dict(left=104**7,right=25**10),
                finite_threshold_controls=rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = dict(status='arbitrary-field structural bound; finite controls only; no exponent improvement',
                  binary_q2=binary_q2_exhaustive(), finite_fields=finite_field_controls(),
                  sharp_square_control=sharp_square_control(),
                  regular_rank_two_threshold=regular_rank_two_threshold())
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
