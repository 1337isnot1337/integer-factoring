#!/usr/bin/env python3
"""Independent small exact checks for the matrix-unit monodromy."""
from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import subprocess
import tempfile
import time


def direct(t, p, q, boundary):
    dp = [int(i == j) for i in range(4) for j in range(4)]
    for position in range(t):
        a, b = (p >> (2*position)) & 3, (q >> (2*position)) & 3
        following = dp.copy() if a == b else [0]*16
        for i in range(4):
            following[4*i+a] ^= dp[4*i+b] ^ (dp[4*i+b] << 1)
        dp = following
    answer = 0
    for a in range(4):
        if boundary >> a & 1:
            answer ^= dp[5*a]
    return answer


def walk(t, p, q, boundary):
    answer = 0
    for initial in range(4):
        if not (boundary >> initial & 1):
            continue
        state, a, b, valid = initial, 0, 0, True
        for position in range(t):
            u, v = (p >> (2*position)) & 3, (q >> (2*position)) & 3
            if u == v:
                a += u == state
            elif v != state:
                valid = False
                break
            else:
                state = u
                b += 1
        if valid and state == initial:
            polynomial = 1 << a
            for _ in range(b):
                polynomial ^= polynomial << 1
            answer ^= polynomial
    return answer


def rank(rows):
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length()-1
            if pivot not in basis:
                basis[pivot] = row
                break
            row ^= basis[pivot]
    return len(basis)


def zip_bits(t, x, y):
    return sum((((x >> i) & 1)+2*((y >> i) & 1)) << (2*i)
               for i in range(t))


def profile(t, boundary):
    n = 1 << t
    row_polynomials = [[walk(t, p, q, boundary) for q in range(n*n)]
                       for p in range(n*n)]
    column_polynomials = [[walk(t, zip_bits(t, x, y), zip_bits(t, v, u), boundary)
                           for v in range(n) for y in range(n)]
                          for x in range(n) for u in range(n)]
    profiles = []
    for polynomials in (row_polynomials, column_polynomials):
        profiles.append([rank([sum(((value >> degree) & 1) << j
                                    for j, value in enumerate(row))
                               for row in polynomials])
                         for degree in range(t+1)])
    return profiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--native', type=Path, default=None,
                        help='Use an existing native executable; otherwise compile the sibling C++ source in a temporary directory.')
    args = parser.parse_args()
    started = time.monotonic()
    native_workspace = tempfile.TemporaryDirectory(prefix='matrix-unit-monodromy-')
    native_directory = Path(native_workspace.name)
    if args.native is None:
        native_path = native_directory / 'matrix_unit_monodromy'
        subprocess.run(['c++', '-O2', '-std=c++17',
                        str(Path(__file__).with_name('matrix_unit_monodromy.cpp')),
                        '-o', str(native_path)], check=True)
    else:
        native_path = args.native.resolve()
    entry_checks = endpoint_checks = profile_checks = protected_checks = 0
    for boundary in (1, 14):
        for t in range(1, 4):
            for p, q in product(range(4**t), repeat=2):
                value = direct(t, p, q, boundary)
                assert value == walk(t, p, q, boundary)
                assert value.bit_count() % 2 == (p == q)
                entry_checks += 1
                endpoint_checks += 1
            rr, cc = profile(t, boundary)
            temporary = native_directory / f'profile-t{t}-b{boundary}.json'
            native = json.loads(subprocess.check_output(
                [str(native_path), str(t), str(temporary), str(boundary)], text=True))
            assert rr == native['row_ranks'] and cc == native['column_ranks']
            profile_checks += 2
    for t in range(1, 5):
        labels = list(product((2, 1, 3), repeat=t))
        encoded = [sum(a << (2*i) for i, a in enumerate(word)) for word in labels]
        for i, a in enumerate(labels):
            x = sum((v & 1) << j for j, v in enumerate(a))
            u = sum((v >> 1) << j for j, v in enumerate(a))
            for j, b in enumerate(labels):
                assert walk(t, encoded[i], encoded[j], 1) == (i == j)
                v = sum((w & 1) << k for k, w in enumerate(b))
                y = sum((w >> 1) << k for k, w in enumerate(b))
                value = walk(t, zip_bits(t, x, y), zip_bits(t, v, u), 1)
                if i == j:
                    assert value == 1
                elif i < j:
                    assert value == 0
                protected_checks += 2
    assert direct(2, 1, 4, 1) == 5  # (1+z)^2 at two odd-parity positions.
    output = {'status': 'passed', 'entry_algorithm_comparisons': entry_checks,
              'endpoint_checks': endpoint_checks, 'full_cut_rank_profile_comparisons': profile_checks,
              'protected_row_identity_and_column_triangular_checks': protected_checks,
              'outside_local_even_parity_span': True,
              'scope': 'Finite verification only; the full coefficient-span restriction is proved in the companion note.',
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'seconds': time.monotonic()-started}
    print(json.dumps(output, indent=2, sort_keys=True))
    if args.output:
        args.output.write_text(json.dumps(output, indent=2, sort_keys=True)+'\n')
    native_workspace.cleanup()


if __name__ == '__main__':
    main()
