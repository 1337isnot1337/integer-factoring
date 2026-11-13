#!/usr/bin/env python3
"""Compare the two twelve-site class-sum T normalizations at ten sites.

The coarse (T,R) certificate uses the Kronecker class sum

    T_coarse = sum_i M_i kron M_i.

The coefficient table resolves T through entry functionals on quadruples
(a,b,c,d), which gives

    T_entry[(a,b),(c,d)] = sum_i M_i[a,d] M_i[b,c].

These differ by the path-pair swap on the column pair:
T_entry = T_coarse R.  At ten sites both matrices fit in memory.  Their
generalized nullities show that the swap changes the spectrum, so the
twelve-site coarse T-primary dimensions cannot be identified with the
coefficient-table T groups without a separate computation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
COARSE = ROOT / 'evidence/bounded-spin-coarse-twelve.json'
CUT = ROOT / 'evidence/bounded-spin-twelve-cut-trace-form.json'
sys.path.insert(0, str(ROOT / 'tools'))

from check_bounded_spin_transposition_subalgebra import (
    frame, replica_swap, transposition_matrices,
)
from native_gf5_flint import multiply, rank

PRIME = 5
LENGTH = 10
EXPECTED_COARSE = {0: 441, 1: 225, 2: 490, 3: 0, 4: 0}
EXPECTED_ENTRY = {0: 441, 1: 225, 2: 335, 3: 155, 4: 0}
TEN_PRIMARY_DIMENSIONS = [34, 1, 225, 300, 35, 406, 155]


def generalized_nullities(matrix: np.ndarray) -> dict[int, list[int]]:
    size = len(matrix)
    identity = np.eye(size, dtype=np.int64)
    result = {}
    for value in range(PRIME):
        shift = (matrix - value * identity) % PRIME
        power = identity.copy()
        values = []
        for _ in range(3):
            power = multiply(power, shift)
            values.append(int(size - rank(power)))
        result[value] = values
    return result


def check() -> dict:
    ps, words, vectors, inverse_norms = frame(LENGTH)
    matrices = transposition_matrices(LENGTH, words, vectors, inverse_norms)
    paths_count = len(ps)
    size = paths_count * paths_count
    coarse = np.zeros((size, size), dtype=np.int64)
    entry = np.zeros((size, size), dtype=np.int64)
    for matrix in matrices:
        m = matrix.astype(np.int64)
        coarse = (coarse + np.kron(m, m)) % PRIME
    swap = replica_swap(paths_count).astype(np.int64)
    entry = (coarse @ swap) % PRIME

    for a, b, c, d in ((0, 1, 2, 3), (3, 0, 1, 2), (2, 3, 0, 1)):
        column = np.zeros((size, 1), dtype=np.int64)
        column[a * paths_count + b, 0] = 1
        actual = int(multiply(entry, column)[c * paths_count + d, 0])
        expected = sum(int(matrix[a, d]) * int(matrix[b, c])
                       for matrix in matrices) % PRIME
        if actual != expected:
            raise ValueError('entry-convention comparison failed')

    coarse_nullities = generalized_nullities(coarse)
    entry_nullities = generalized_nullities(entry)
    coarse_single = {value: values[-1] for value, values in coarse_nullities.items()}
    entry_single = {value: values[-1] for value, values in entry_nullities.items()}
    if coarse_single != EXPECTED_COARSE:
        raise ValueError('unexpected coarse generalized nullities')
    if entry_single != EXPECTED_ENTRY:
        raise ValueError('unexpected entry-convention generalized nullities')
    if sum(EXPECTED_COARSE.values()) != size:
        raise ValueError('coarse nullities do not sum to the module dimension')
    if sum(EXPECTED_ENTRY.values()) != size:
        raise ValueError('entry nullities do not sum to the module dimension')
    if sum(TEN_PRIMARY_DIMENSIONS) != size:
        raise ValueError('primary dimensions do not sum to the module dimension')
    if entry_single[2] + entry_single[3] != coarse_single[2]:
        raise ValueError('the coarse T=2 block does not split as expected')

    coarse_twelve = json.loads(COARSE.read_text())
    joint = {(int(t), int(r)): int(dimension)
            for t, r, dimension in coarse_twelve['coarse_dimensions']}
    entry_twelve = {}
    for (value, sign), dimension in joint.items():
        entry_value = value if sign == 1 else (-value) % PRIME
        entry_twelve[entry_value] = entry_twelve.get(entry_value, 0) + dimension
    expected_twelve = {0: 1485, 1: 4500, 4: 1936}
    if entry_twelve != expected_twelve:
        raise ValueError('unexpected twelve-site entry-convention dimensions')
    cut = json.loads(CUT.read_text())
    row_group_sums = {}
    for index, (value, _) in enumerate([
            (0, 0), (1, 0), (1, 0), (4, 0), (4, 0), (1, 0), (0, 0), (1, 0)]):
        row_group_sums[value] = (row_group_sums.get(value, 0)
                                 + cut['idempotent_trace_residues']['row'][index]) % PRIME
    if {value: total % PRIME for value, total in row_group_sums.items()} != {
            value: dimension % PRIME for value, dimension in expected_twelve.items()}:
        raise ValueError('entry dimensions do not match the row trace group sums')

    return {
        'passed': True,
        'length': LENGTH,
        'prime': PRIME,
        'paths': paths_count,
        'module_dimension': size,
        'transposition_matrices': len(matrices),
        'coarse_T_generalized_nullities': {
            str(value): values for value, values in coarse_nullities.items()},
        'entry_T_generalized_nullities': {
            str(value): values for value, values in entry_nullities.items()},
        'coarse_T_primary_dimensions': {
            str(value): dimension for value, dimension in coarse_single.items()
            if dimension},
        'entry_T_primary_dimensions': {
            str(value): dimension for value, dimension in entry_single.items()
            if dimension},
        'ten_site_primary_dimensions': TEN_PRIMARY_DIMENSIONS,
        'coarse_T2_splits_as_entry_T2_plus_T3': [
            entry_single[2], entry_single[3]],
        'twelve_site_coarse_T_primary_dimensions': {
            str(value): sum(dimension for (t, _), dimension in joint.items() if t == value)
            for value in (0, 1, 4)},
        'twelve_site_entry_T_primary_dimensions': {
            str(value): dimension for value, dimension in expected_twelve.items()},
        'scope': (
            'Exact ten-site control and exact twelve-site layer comparison for the '
            'coarse Kronecker class sum '
            'T_coarse = sum M_i kron M_i with the coefficient-entry class sum '
            'T_entry = T_coarse R, where R swaps the column path pair. Both '
            'acts on the 1156-dimensional path-pair module. The spectra differ: '
            'T_coarse has generalized dimensions 441,225,490 at eigenvalues '
            '0,1,2, while T_entry has dimensions 441,225,335,155 at eigenvalues '
            '0,1,2,3. The published primary dimensions 34,1,225,300,35,406,155 '
            'split the coarse T=2 block as 335+155, matching T_entry. Therefore '
            'the twelve-site exact coarse T-primary dimensions 1485,3281,3155 '
            'are not the coefficient-table T group dimensions. Because '
            'T_entry = T_coarse R and R^2 = 1, the twelve-site entry '
            'dimensions are obtained by flipping the R sign: 1485, 4500, 1936 '
            'for T=0,1,4. Their residues 0,0,1 match the coefficient-table row '
            'trace group sums. Individual central-idempotent ranks and the '
            'two-cut pairing remain open; no cut rank or score change is '
            'claimed.'
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    text = json.dumps(result, indent=2)
    if args.output:
        args.output.write_text(text + '\n')
    print(text)


if __name__ == '__main__':
    main()
