#!/usr/bin/env python3
"""Exact controls for nilpotent-superoperator crossing families over F2.

No numerical rank or factorization package is required. All polynomial
entries use bit masks, with bit i the coefficient of z**i.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

def xorall(a):
 o=0
 for b in a:o^=b
 return o

def matmul(A,B,q):return [xorall(A[i*q+k]&B[k*q+j] for k in range(q)) for i in range(q) for j in range(q)]

def trace(A,q):return xorall(A[i*q+i] for i in range(q))

def inverse(A):
 n=len(A);r=[sum((x&1)<<j for j,x in enumerate(a))|(1<<(n+i)) for i,a in enumerate(A)]
 for j in range(n):
  p=next((i for i in range(j,n) if r[i]>>j&1),None)
  if p is None:raise ValueError('singular')
  r[j],r[p]=r[p],r[j]
  for i in range(n):
   if i!=j and r[i]>>j&1:r[i]^=r[j]
 return [[r[i]>>(n+j)&1 for j in range(n)] for i in range(n)]

def dual(basis,q):
 gram=[[trace(matmul(A,B,q),q) for B in basis] for A in basis]
 inv=inverse(gram)
 return [[xorall(inv[j][i]&basis[j][a] for j in range(q*q)) for a in range(q*q)] for i in range(q*q)]

def kron(A,B,q):return [[A[i*q+k]&B[j*q+l] for k in range(q) for l in range(q)] for i in range(q) for j in range(q)]

def construction(basis,q):
 n=q*q;L=n-1;D=dual(basis,q);K=[[0]*n for _ in range(n)]
 for i in range(n):
  for j in range(n-i):
   B=[D[L-j][a]^(D[L-j+1][a] if j else 0) for a in range(n)]
   R=kron(basis[i],B,q)
   for a in range(n):
    for b in range(n):K[a][b]^=R[a][b]<<(i+j)
 return K,D

def pmul(a,b):
 c=0
 while b:
  if b&1:c^=a
  a<<=1;b>>=1
 return c

def rank(rows):
 p={}
 for x in rows:
  while x:
   k=x.bit_length()-1
   if k in p:x^=p[k]
   else:p[k]=x;break
 return len(p)

def inverse_series(p,cap):
 assert p&1
 out=1;m=1
 while m<cap:
  m=min(2*m,cap)
  out=pmul(p,pmul(out,out))&((1<<m)-1)
 return out

def smith_series(K,cap):
 n=len(K);A=[row[:] for row in K];mask=(1<<cap)-1;ans=[]
 for k in range(n):
  best=None
  for i in range(k,n):
   for j in range(k,n):
    p=A[i][j]
    if p:
     v=(p&-p).bit_length()-1
     if best is None or v<best[0]:best=(v,i,j)
  if best is None:return ans+[cap]*(n-k)
  v,i,j=best;ans.append(v);A[k],A[i]=A[i],A[k]
  for r in A:r[k],r[j]=r[j],r[k]
  inv=inverse_series(A[k][k]>>v,cap)
  for i in range(k+1,n):
   if not A[i][k]:continue
   quotient=pmul(A[i][k]>>v,inv)&mask
   for j in range(k+1,n):A[i][j]^=pmul(quotient,A[k][j])&mask
   A[i][k]=0
 return ans

def row_infinity(K):
 n=len(K);L=n-1
 return [[sum(((x>>k)&1)<<(L-k) for k in range(L+1)) for x in row] for row in K]

def pdiv(a,b):
 assert b
 q=0;db=b.bit_length()
 while a and a.bit_length()>=db:
  d=a.bit_length()-db;q^=1<<d;a^=b<<d
 if a:raise ValueError(('nonexact division',a,b))
 return q

def det_poly(K):
 A=[row[:] for row in K];n=len(A);prev=1
 for k in range(n-1):
  p=next((i for i in range(k,n) if A[i][k]),None)
  if p is None:return 0
  A[k],A[p]=A[p],A[k];pivot=A[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):A[i][j]=pdiv(pmul(pivot,A[i][j])^pmul(A[i][k],A[k][j]),prev)
   A[i][k]=0
  prev=pivot
 return A[-1][-1]

def vector_mask(entries):
    return sum((v & 1) << i for i, v in enumerate(entries))


def apply_columns(columns, vector):
    out = 0
    while vector:
        bit = vector & -vector
        out ^= columns[bit.bit_length() - 1]
        vector ^= bit
    return out


def compose_columns(left, right):
    return [apply_columns(left, column) for column in right]


def phi(columns, q):
    """Phi(M)(x,y,u,v) = M(E_vy)_(x,u), in the row cut."""
    return [
        [
            (columns[v * q + y] >> (x * q + u)) & 1
            for u in range(q) for v in range(q)
        ]
        for x in range(q) for y in range(q)
    ]


def nilpotent_family(columns, q, index):
    n = q * q
    powers = [[1 << i for i in range(n)]]
    for _ in range(index):
        powers.append(compose_columns(columns, powers[-1]))
    assert not any(powers[-1])
    assert index == 1 or any(powers[-2])
    result = [[0] * n for _ in range(n)]
    for k in range(index):
        term = phi([a ^ b for a, b in zip(powers[k], powers[k + 1])], q)
        for i in range(n):
            for j in range(n):
                result[i][j] ^= term[i][j] << (index - 1 - k)
    return result


def column_cut(polynomial_rows, q):
    return [
        [polynomial_rows[x * q + y][u * q + v]
         for y in range(q) for v in range(q)]
        for x in range(q) for u in range(q)
    ]


def reverse_at_degree(matrix, degree):
    return [[sum(((v >> k) & 1) << (degree - k)
                 for k in range(degree + 1)) for v in row] for row in matrix]


def toeplitz_rank(matrix, precision):
    n = len(matrix)
    rows = []
    for a in range(precision):
        for i in range(n):
            rows.append(sum(
                ((matrix[i][j] >> (a - b)) & 1) << (b * n + j)
                for b in range(a + 1) for j in range(n)))
    return rank(rows)


def check_orders(matrix, orders, precision, toeplitz_precisions):
    actual = smith_series(matrix, precision)
    assert actual == orders, (actual, orders)
    for h in toeplitz_precisions:
        assert toeplitz_rank(matrix, h) == sum(max(0, h - v) for v in orders)
    return actual


def polynomial_remainder(a, b):
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def polynomial_gcd(a, b):
    while b:
        a, b = b, polynomial_remainder(a, b)
    return a


def irreducible(poly):
    degree = poly.bit_length() - 1
    value = 2
    # This slightly redundant all-degree test is exact and tiny here.
    for k in range(1, degree + 1):
        value = polynomial_remainder(pmul(value, value), poly)
        if k <= degree // 2 and polynomial_gcd(value ^ 2, poly) != 1:
            return False
    return value == 2


def make_q4_channel():
    q, n = 4, 16
    constants = {"P": 4681, "inverse": 37448, "A": 2122, "B": 33957}
    mats = {key: [(mask >> i) & 1 for i in range(n)]
            for key, mask in constants.items()}
    P, Q, A, B = (mats[key] for key in ("P", "inverse", "A", "B"))
    ident = [int(i // q == i % q) for i in range(n)]
    assert matmul(P, Q, q) == ident == matmul(Q, P, q)
    AB = matmul(A, B, q)
    assert AB == matmul(B, A, q)
    C = [a ^ b for a, b in zip(ident, AB)]

    def operator(X):
        XT = [X[j * q + i] for i in range(q) for j in range(q)]
        Y = [a ^ b for a, b in zip(matmul(XT, C, q),
                                   matmul(matmul(A, XT, q), B, q))]
        UX = matmul(matmul(P, Y, q), Q, q)
        return [a ^ b for a, b in zip(X, UX)]

    columns = [vector_mask(operator([int(i == j) for i in range(n)]))
               for j in range(n)]
    ident_mask = vector_mask(ident)
    assert apply_columns(columns, ident_mask) == 0
    assert all(trace([(column >> i) & 1 for i in range(n)], q) == 0
               for column in columns)
    chain = [1]
    for _ in range(n - 1):
        chain.append(apply_columns(columns, chain[-1]))
    assert chain[-1] == ident_mask
    assert rank(chain) == n
    assert apply_columns(columns, chain[-1]) == 0
    basis = [[(vector >> i) & 1 for i in range(n)] for vector in reversed(chain)]
    return constants, columns, basis


def verify():
    q2_basis = [[1, 0, 0, 1], [0, 1, 0, 0],
                [0, 0, 1, 0], [1, 0, 1, 0]]
    q2, _ = construction(q2_basis, 2)
    assert q2 == [[13, 6, 10, 12], [12, 5, 8, 10],
                  [12, 8, 5, 6], [0, 12, 12, 13]]
    assert det_poly(q2) == 1
    q2_row = check_orders(reverse_at_degree(q2, 3), [0, 0, 1, 11], 13,
                          [1, 2, 4, 11, 12])
    q2_col = check_orders(column_cut(q2, 2), [0, 4, 4, 4], 13, [1, 4, 5])

    constants, columns, basis = make_q4_channel()
    q4, dual_basis = construction(basis, 4)
    assert q4 == nilpotent_family(columns, 4, 16)
    assert dual_basis[-1] == [int(i // 4 == i % 4) for i in range(16)]
    for i in range(16):
        for j in range(16):
            assert (q4[i][j] & 1) == int(i == j)
            x, y = divmod(i, 4)
            u, v = divmod(j, 4)
            assert (q4[i][j].bit_count() & 1) == int(x == v and y == u)
    # Four independent cyclic-vector changes leave the entire family unchanged.
    chain = [vector_mask(v) for v in reversed(basis)]
    for shift in (1, 3, 8, 15):
        changed = [chain[0] ^ chain[shift]]
        for _ in range(15):
            changed.append(apply_columns(columns, changed[-1]))
        changed_basis = [[(v >> i) & 1 for i in range(16)]
                         for v in reversed(changed)]
        assert construction(changed_basis, 4)[0] == q4
    q4_row = check_orders(reverse_at_degree(q4, 15), [0, 0] + [1] * 12 + [9, 11],
                          241, [1, 2, 3, 9, 10, 11, 12])
    q4_col = check_orders(column_cut(q4, 4), [0] + [16] * 15, 241, [1, 2, 16, 17])
    assert det_poly(column_cut(q4, 4)) == 1 << 240
    row_determinant = det_poly(q4)
    factor_specs = [([0, 3, 4, 5, 7], 4),
                    ([0, 2, 4, 5, 6, 8, 9, 12, 14], 8),
                    ([0, 3, 5, 6, 8, 12, 17], 4)]
    factors = [(sum(1 << i for i in powers), mult) for powers, mult in factor_specs]
    product = 1
    for factor, multiplicity in factors:
        assert irreducible(factor)
        assert factor & 1 and (factor.bit_count() & 1)
        for _ in range(multiplicity):
            product = pmul(product, factor)
    assert product == row_determinant
    assert row_determinant.bit_length() - 1 == 208
    U = [column ^ (1 << i) for i, column in enumerate(columns)]
    r0 = rank([vector_mask(row) for row in phi(U, 4)])
    r2 = rank([vector_mask(row) for row in phi(compose_columns(U, U), 4)])
    assert (r0, r2) == (2, 14)
    assert q4_row.count(1) >= r2 - 2 * r0

    nonregular = []
    for q in (2, 3):
        n = q * q
        S = [(1 << (i * q + j)) ^ (1 << (j * q + i))
             for i in range(q) for j in range(q)]
        K = nilpotent_family(S, q, 2)
        twos = q * (q - 1) // 2
        expected = [0] * twos + [1] * q + [2] * twos
        check_orders(column_cut(K, q), expected, n + 1, [1, 2, 3])
        row_expected = [0] + [1] * (n - 1) if q % 2 else [0] + [1] * (n - 2) + [2]
        check_orders(reverse_at_degree(K, 1), row_expected, n + 1, [1, 2, 3])
        nonregular.append({"q": q, "index": 2, "row_infinity": row_expected, "jordan_blocks": [2] * twos + [1] * q,
                           "column_orders": expected})
    single_channel = []
    for q in (4, 6):
        n, half = q * q, q // 2
        P = [0] * n
        for i in range(half):
            P[i * q + half + i] = 1
            P[(half + i) * q + i] = 1
            if i + 1 < half:
                P[(half + i) * q + i + 1] = 1
        Q = sum(inverse([P[i * q:(i + 1) * q] for i in range(q)]), [])
        S = []
        for j in range(n):
            i0, j0 = divmod(j, q)
            XT = [int(i == j0 * q + i0) for i in range(n)]
            UX = matmul(matmul(P, XT, q), Q, q)
            S.append(vector_mask(UX) ^ (1 << j))
        power = [1 << i for i in range(n)]
        ranks = [n]
        for m in range(1, 2 * q + 1):
            power = compose_columns(S, power)
            ranks.append(rank(power))
            if not any(power):
                break
        else:
            raise AssertionError("single-channel control was not nilpotent")
        counts_at_least = [ranks[j - 1] - ranks[j] for j in range(1, m + 1)] + [0]
        blocks = [j for j in range(1, m + 1)
                  for _ in range(counts_at_least[j - 1] - counts_at_least[j])]
        col_expected = sorted(v for b in blocks for v in [m - b] + [m] * (b - 1))
        row_expected = [0] + [1] * (n - 2) + [2]
        K = nilpotent_family(S, q, m)
        check_orders(reverse_at_degree(K, m - 1), row_expected, n * (m - 1) + 1, [1, 2, 3])
        check_orders(column_cut(K, q), col_expected, n * (m - 1) + 1, [1, m, m + 1])
        single_channel.append({"q": q, "P_mask": vector_mask(P), "nilpotence_index": m,
                               "jordan_blocks": blocks, "row_infinity": row_expected,
                               "column_zero": col_expected})
    return {
        "single_transpose_controls": single_channel,
        "q2": {"basis": q2_basis, "row_infinity": q2_row, "column_zero": q2_col,
               "row_determinant_mask": 1},
        "q4": {**constants, "row_infinity": q4_row, "column_zero": q4_col,
               "row_determinant_degree": 208, "row_determinant_mask": row_determinant,
               "irreducible_factor_masks_and_multiplicities": factors,
               "rank_phi_U": r0, "rank_phi_U_squared": r2},
        "nonregular_controls": nonregular,
        "claims": {"new_factoring_exponent": False, "all_size_proof_in_note": True}
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify()
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload)
    print(json.dumps({"passed": True, "q2_row_infinity": result["q2"]["row_infinity"],
                      "q4_row_infinity": result["q4"]["row_infinity"],
                      "q4_column_zero": result["q4"]["column_zero"]}))


if __name__ == "__main__":
    main()
