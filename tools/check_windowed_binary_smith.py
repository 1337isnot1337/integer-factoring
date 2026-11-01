#!/usr/bin/env python3
"""Check finite-frame equivalence and full affine coefficient equations."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from random import Random

from full_alphabet_two_cut_lattices import FullAlphabetTwoCutLattice
from native_matrix_smith_f4 import MatrixOracle
from weyl_lattice_filter_synthesis import regular_smith_basis, matrix_multiply, polynomial_multiply
from windowed_binary_smith import windowed_regular_smith_basis
from windowed_full_alphabet_lattices import WindowedFullAlphabetTwoCutLattice


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng, matrices, compared = Random(0x01A195C2F), [], 0
    for n in (2, 3, 4):
        for trial in range(8):
            matrix = []
            for i in range(n):
                row = []
                for j in range(n):
                    q = rng.randrange(16) if rng.random() < .6 else 0
                    row.append(int(i == j) ^ q ^ (q << 1))
                matrix.append(row)
            full = regular_smith_basis(matrix)
            result = windowed_regular_smith_basis(matrix, 5)
            assert result[0] == full[0]
            mask = (1 << result[4]['output_precision'])-1
            for actual, original in zip(result[1:4], full[1:4]):
                assert actual == [[v & mask for v in row] for row in original]
                compared += n*n
            native = MatrixOracle([[(v, 0) for v in row] for row in matrix], n).profile()['orders']
            assert list(result[0]) == native
            matrices.append({'size': n, 'trial': trial, 'source': matrix,
                             'orders': result[0], 'window': result[4]})
    # Nonzero orders and dense polynomial units stress the precision lost
    # when the pivot is shifted down by its valuation.
    for n in (2, 3, 4):
        for trial in range(4):
            units = []
            for _ in range(2):
                unit = [[int(i == j) for j in range(n)] for i in range(n)]
                for _ in range(8):
                    i, j = rng.sample(range(n), 2)
                    p = rng.randrange(256)
                    unit[i] = [a ^ polynomial_multiply(p, b) for a, b in zip(unit[i], unit[j])]
                units.append(unit)
            expected = tuple(2+7*i for i in range(n))
            diagonal = [[(1 << expected[i]) if i == j else 0 for j in range(n)] for i in range(n)]
            matrix = matrix_multiply(matrix_multiply(units[0], diagonal), units[1])
            full = regular_smith_basis(matrix)
            result = windowed_regular_smith_basis(matrix, 5)
            assert result[0] == full[0] == expected
            mask = (1 << result[4]['output_precision'])-1
            for actual, original in zip(result[1:4], full[1:4]):
                assert actual == [[v & mask for v in row] for row in original]
                compared += n*n
            native = MatrixOracle([[(v, 0) for v in row] for row in matrix], n).profile()['orders']
            assert list(expected) == native
            matrices.append({'size': n, 'trial': trial, 'family': 'polynomial-unit diagonal',
                             'source': matrix, 'orders': result[0], 'window': result[4]})
    # All coefficient equations must agree, not merely the source profiles.
    affine = []
    for trial in range(6):
        q, degree = 2, 3
        values = []
        for index in range(q**4):
            x, y, u, v = ((index//q**i) % q for i in range(4))
            p = rng.randrange(1 << degree)
            values.append((int(x == v and y == u) ^ p ^ (p << 1), 0))
        heights = [[rng.randint(-2, 2) for _ in range(q)] for _ in range(8)]
        full = FullAlphabetTwoCutLattice(values, q, degree, heights, raise_cap=3)
        window = WindowedFullAlphabetTwoCutLattice(values, q, degree, heights, raise_cap=3)
        assert full.orders == window.orders and full.caps == window.caps
        assert full.equations == window.equations
        for _ in range(4):
            targets = [rng.randrange(cap+1) for cap in full.caps]
            assert full.solve(targets) == window.solve(targets)
        affine.append({'trial': trial, 'source': values, 'heights': heights,
                       'complete_equations': sum(len(row) for row in full.equations),
                       'window': window.frame_windows})
    # A large order forces a genuine precision retry; a singular input fails.
    retry = windowed_regular_smith_basis([[1, 0], [0, 1 << 100]], 4)
    assert retry[0] == (0, 100) and len(retry[4]['attempts']) > 1
    rejected = False
    try:
        windowed_regular_smith_basis([[1, 1], [1, 1]], 4)
    except ArithmeticError:
        rejected = True
    assert rejected
    root = Path(__file__).resolve().parents[1]
    names = ['tools/check_windowed_binary_smith.py', 'tools/windowed_binary_smith.py',
             'tools/full_alphabet_two_cut_lattices.py', 'tools/windowed_full_alphabet_lattices.py',
             'tools/weyl_lattice_filter_synthesis.py',
             'tools/native_matrix_smith_f4.py']
    report = {'status': 'pass', 'complete_small_matrix_controls': len(matrices),
              'frame_and_diagonal_coefficients_compared': compared,
              'complete_affine_equation_controls': affine, 'matrix_controls': matrices,
              'forced_retry': retry[4], 'singular_rejected': rejected,
              'source_sha256': {n: sha256((root/n).read_bytes()).hexdigest() for n in names},
              'scope': 'Exact low-coefficient equivalence to the original polynomial frames. '
              'No changed feasible coefficient fiber, stronger pressure or factoring bound.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({key: report[key] for key in
                      ('status', 'complete_small_matrix_controls', 'frame_and_diagonal_coefficients_compared')}))


if __name__ == '__main__':
    main()
