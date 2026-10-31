#!/usr/bin/env python3
"""Exact binary controls for the all-field regular-pencil lemma."""
import argparse
import json
from pathlib import Path
from check_nilpotent_flag_crossings import inverse, kron, matmul, rank, trace, vector_mask


def binary_matrix(mask, q):
    return [(mask >> i) & 1 for i in range(q*q)]


def determinant_pencil_zero_2(A, C):
    return not (
        (A[0]*A[3] ^ A[1]*A[2])
        or (A[0]*C[3] ^ C[0]*A[3] ^ A[1]*C[2] ^ C[1]*A[2])
        or (C[0]*C[3] ^ C[1]*C[2])
    )


def matrix_rank(rows):
    return rank([vector_mask(row) for row in rows])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    mats = [binary_matrix(m, 2) for m in range(16)]
    pencils = tests = 0
    for A in mats:
        for C in mats:
            if not determinant_pencil_zero_2(A, C):
                continue
            pencils += 1
            for B in mats:
                for D in mats:
                    AB, CD = kron(A, B, 2), kron(C, D, 2)
                    assert matrix_rank([[a ^ b for a, b in zip(r, s)]
                                        for r, s in zip(AB, CD)]) < 4
                    tests += 1
    A = [1,0,0,0,0,0,0,0,1]
    C = [0,0,0,0,1,0,0,0,1]
    identity3 = [1,0,0,0,1,0,0,0,1]
    # Multiplication by t in F2[t]/(t^3+t+1), columns 2,4,3.
    D = [0,0,1,1,0,1,0,1,0]
    assert all(matrix_rank([M[3*i:3*i+3] for i in range(3)]) < 3
               for M in [A, C, [a ^ c for a, c in zip(A, C)]])
    AI, CD = kron(A, identity3, 3), kron(C, D, 3)
    extension_rank = matrix_rank([[a ^ b for a, b in zip(r, s)]
                                 for r, s in zip(AI, CD)])
    assert extension_rank == 9
    identity2 = [1,0,0,1]
    normalized = commuting = 0
    for P in mats:
        try:
            Q = inverse([P[:2], P[2:]])
        except ValueError:
            continue
        Q = [v for row in Q for v in row]
        for A in mats:
            for B in mats:
                AB, BA = matmul(A, B, 2), matmul(B, A, 2)
                Q0 = matmul([a ^ b for a, b in zip(identity2, AB)], Q, 2)
                P1, Q1 = matmul(P, A, 2), matmul(B, Q, 2)
                columns = [[(P[2*k+j]*Q0[2*i+l]) ^ (P1[2*k+j]*Q1[2*i+l])
                            for k in range(2) for l in range(2)]
                           for i in range(2) for j in range(2)]
                image_identity = [columns[0][i] ^ columns[3][i] for i in range(4)]
                assert image_identity == identity2
                preserves_trace = all(trace(columns[j], 2) == identity2[j] for j in range(4))
                assert preserves_trace == (AB == BA)
                normalized += 1
                commuting += preserves_trace
    result = dict(status='all-field lemma with exact finite controls; no score improvement',
                  singular_binary_pencils=pencils, singular_tensor_sums_checked=tests,
                  extension_example_tensor_rank=extension_rank,
                  normalized_binary_channels_checked=normalized,
                  trace_preserving_commuting_cases=commuting)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
