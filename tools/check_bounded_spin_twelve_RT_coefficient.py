#!/usr/bin/env python3
"""Identify R,T inside the exact twelve-site coefficient table and its blocks."""
from __future__ import annotations

from hashlib import sha256
import itertools
import json
from math import comb
from pathlib import Path
import random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

from bounded_spin_crossing import paths
from bounded_spin_six_vertex import monomial_entry
from check_bounded_spin_transposition_subalgebra import (
    frame, transposition_matrices,
)

PRIME = 5
ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / 'evidence/bounded-spin-twelve-coefficient-table.json'
IDEMPOTENTS = ROOT / 'evidence/bounded-spin-twelve-coefficient-idempotents.json'
TEN_MODEL = ROOT / 'evidence/bounded-spin-ten-algebra-model.json'
OUTPUT = ROOT / 'evidence/bounded-spin-twelve-RT-coefficient.json'


def original_form(four_paths: tuple) -> list[int]:
    result = [0] * (len(four_paths[0]) + 1)
    for (_, b, c), value in monomial_entry(four_paths).items():
        for j in range(2 * b + 1):
            result[2 * c + j] = (result[2 * c + j]
                                 + value * comb(2 * b, j) * (-1) ** j) % PRIME
    return result


def rank_mod(rows: list[list[int]]) -> int:
    if not rows:
        return 0
    a = [[int(x) % PRIME for x in row] for row in rows]
    r = 0
    for column in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][column]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inverse = pow(int(a[r][column]), -1, PRIME)
        a[r] = [x * inverse % PRIME for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][column]:
                factor = a[i][column]
                a[i] = [(x - factor * y) % PRIME for x, y in zip(a[i], a[r])]
        r += 1
    return r


def solve_mod(rows: list[list[int]], rhs: list[int]) -> list[int]:
    a = [[int(x) % PRIME for x in row] + [int(y) % PRIME]
         for row, y in zip(rows, rhs)]
    row_count, column_count = len(a), len(a[0]) - 1
    r = 0
    pivots = []
    for column in range(column_count):
        pivot = next((i for i in range(r, row_count) if a[i][column]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inverse = pow(int(a[r][column]), -1, PRIME)
        a[r] = [x * inverse % PRIME for x in a[r]]
        for i in range(row_count):
            if i != r and a[i][column]:
                factor = a[i][column]
                a[i] = [(x - factor * y) % PRIME for x, y in zip(a[i], a[r])]
        pivots.append(column)
        r += 1
    for i in range(r, row_count):
        if not any(a[i][:column_count]) and a[i][-1]:
            raise ValueError('inconsistent system')
    solution = [0] * column_count
    for i, column in enumerate(pivots):
        solution[column] = a[i][-1] % PRIME
    return solution


def multiply(regular: np.ndarray, vector: np.ndarray) -> np.ndarray:
    return (regular @ vector) % PRIME


def product(regular: np.ndarray, left, right) -> np.ndarray:
    return (np.tensordot(np.array(left) % PRIME, regular, axes=(0, 0))
            @ np.array(right) % PRIME) % PRIME


def generated_dimension(regular: np.ndarray, unit, generators) -> int:
    space = [np.array(unit) % PRIME]
    frontier = [space[0]]
    while frontier:
        next_frontier = []
        for current in frontier:
            for generator in generators:
                candidate = product(regular, current, generator)
                if rank_mod(space + [candidate]) > len(space):
                    space.append(candidate)
                    next_frontier.append(candidate)
        frontier = next_frontier
    return len(space)


def tensor_entries(length: int) -> tuple[list, np.ndarray]:
    ps, words, vectors, inverse_norms = frame(length)
    matrices = transposition_matrices(length, words, vectors, inverse_norms)
    return ps, matrices


def class_sum_entry(matrices, quadruple) -> int:
    a, b, c, d = quadruple
    return sum(int(matrix[a, d]) * int(matrix[b, c])
               for matrix in matrices) % PRIME


def identity_entry(quadruple) -> int:
    a, b, c, d = quadruple
    return int(a == c and b == d)


def replica_entry(quadruple) -> int:
    a, b, c, d = quadruple
    return int(a == c and b == d)


def coordinates_for(entry_matrix, witness, matrices, evaluator) -> list[int]:
    rhs = [evaluator(w, matrices) for w in witness]
    return solve_mod(entry_matrix, rhs)


def ten_control() -> dict:
    model = json.loads(TEN_MODEL.read_text())
    basis = model['basis_indices']
    ps, matrices = tensor_entries(10)
    q = len(ps)
    candidates = []
    for a in range(q):
        for b in range(q):
            candidates.extend(((a, b, b, a), (a, b, a, b), (a, a, b, b)))
    witness = []
    for quadruple in candidates:
        trial = witness + [quadruple]
        entry_matrix = [[original_form(tuple(ps[i] for i in w))[k]
                         for k in basis] for w in trial]
        if rank_mod(entry_matrix) > len(witness):
            witness.append(quadruple)
        if len(witness) == len(basis):
            break
    entry_matrix = [[original_form(tuple(ps[i] for i in w))[k]
                     for k in basis] for w in witness]
    if rank_mod(entry_matrix) != len(basis):
        raise ValueError('ten-site witness does not separate the coefficient span')
    class_sum = coordinates_for(entry_matrix, witness, matrices,
                                lambda w, m: class_sum_entry(m, w))
    replica = coordinates_for(entry_matrix, witness, matrices,
                              lambda w, m: replica_entry(w))
    expected_class_sum = [0, 2, 0, 1, 0, 2, 2, 0]
    expected_replica = [1, 0, 0, 0, 0, 0, 0, 0]
    if class_sum != expected_class_sum or replica != expected_replica:
        raise ValueError('ten-site entry control disagrees with the published model')
    regular = np.array(model['regular'], dtype=np.int64)
    unit = model['unit']
    dimension = generated_dimension(
        regular, unit, [np.array(class_sum), np.array(replica)])
    if dimension != 6:
        raise ValueError('ten-site generated dimension is not six')
    return {
        'basis_indices': basis,
        'replica_swap_coordinates': replica,
        'class_sum_coordinates': class_sum,
        'generated_dimension': dimension,
        'matches_published_model': True,
    }


def holdout_checks(ps, entry_matrix, basis, coordinates, matrices) -> dict:
    q = len(ps)
    rng = random.Random(20261108)
    transpositions = [(i, j) for i in range(12) for j in range(i + 1, 12)]
    nonzero, zero = [], []
    for _ in range(256):
        p1, p2 = rng.randrange(q), rng.randrange(q)
        i, j = rng.choice(transpositions)
        matrix = matrices[transpositions.index((i, j))]
        choices_c = [k for k in range(q) if matrix[p1, k]]
        choices_d = [k for k in range(q) if matrix[p2, k]]
        if not choices_c or not choices_d:
            continue
        quadruple = (p1, p2, rng.choice(choices_c), rng.choice(choices_d))
        actual = class_sum_entry(matrices, quadruple)
        predicted = sum(
            int(coordinates[k])
            * original_form(tuple(ps[i] for i in quadruple))[basis[k]]
            for k in range(len(basis))) % PRIME
        if actual != predicted:
            raise ValueError('class-sum holdout does not lie in the coefficient span')
        if actual:
            nonzero.append(quadruple)
        else:
            zero.append(quadruple)
        if len(nonzero) >= 24 and len(zero) >= 8:
            break
    if len(nonzero) < 8:
        raise ValueError('not enough nonzero class-sum holdouts')
    return {
        'nonzero_entries_checked': len(nonzero),
        'zero_entries_checked': len(zero),
        'all_checked': True,
        'sample_nonzero': [list(w) for w in nonzero[:4]],
    }


def block_action(regular, transposition, block) -> dict:
    vectors = np.array(block['block_basis'], dtype=np.int64)
    frame = vectors.T
    images = [product(regular, transposition, frame[:, j])
              for j in range(frame.shape[1])]
    columns = [solve_mod(frame, image) for image in images]
    matrix = np.array(columns, dtype=np.int64).T % PRIME
    dimension = matrix.shape[0]
    if dimension == 1:
        trace = int(matrix[0, 0])
        determinant = int(matrix[0, 0])
        characteristic = [-trace % PRIME, 1]
    elif dimension == 2:
        trace = int(matrix[0, 0] + matrix[1, 1]) % PRIME
        determinant = int(matrix[0, 0] * matrix[1, 1]
                          - matrix[0, 1] * matrix[1, 0]) % PRIME
        characteristic = [determinant % PRIME, -trace % PRIME, 1]
    else:
        raise ValueError('unexpected block dimension')
    roots = []
    for value in range(PRIME):
        total = 0
        for coefficient in reversed(characteristic):
            total = (total * value + coefficient) % PRIME
        if total == 0:
            roots.append(value)
    return {
        'idempotent_index': block['index'],
        'kind': block['kind'],
        'regular_rank': block['regular_rank'],
        'action_matrix': [[int(x) for x in row] for row in matrix],
        'trace': trace,
        'determinant': determinant,
        'characteristic_coefficients_low_to_high': characteristic,
        'eigenvalues': roots,
    }


def main() -> None:
    table = json.loads(TABLE.read_text())
    blocks = json.loads(IDEMPOTENTS.read_text())
    if not table.get('passed') or not blocks.get('status') == 'pass':
        raise ValueError('input certificates are not marked as passed')
    if sha256(TABLE.read_bytes()).hexdigest() != blocks['table_sha256']:
        raise ValueError('idempotent certificate is not pinned to this table')
    basis = table['basis_indices']
    regular = np.array(table['regular'], dtype=np.int64)
    unit = np.array(table['unit'], dtype=np.int64)
    q = 89
    ps, matrices = tensor_entries(12)
    if len(ps) != q:
        raise ValueError('path count mismatch')
    index = {tuple(path): i for i, path in enumerate(ps)}
    witness = [tuple(index[tuple(path)] for path in quadruple)
               for quadruple in table['witness_quadruples']]
    entry_matrix = table['entry_matrix']
    if rank_mod(entry_matrix) != len(basis):
        raise ValueError('stored entry witness is not full rank')
    replica = coordinates_for(entry_matrix, witness, matrices,
                              lambda w, m: replica_entry(w))
    class_sum = coordinates_for(entry_matrix, witness, matrices,
                                lambda w, m: class_sum_entry(m, w))
    if replica != unit.tolist():
        raise ValueError('C0 is not the replica coefficient')
    holdouts = holdout_checks(ps, entry_matrix, basis, class_sum, matrices)

    powers = [unit]
    for _ in range(8):
        powers.append(product(regular, powers[-1], class_sum))
    for k in range(1, 9):
        if rank_mod(powers[:k + 1]) != min(k + 1, 5):
            raise ValueError('class-sum powers do not stop at dimension five')
    minimal_coefficients = solve_mod(np.array(powers[:5]).T, powers[5])
    if minimal_coefficients != [0, 4, 0, 2, 0]:
        raise ValueError('unexpected class-sum minimal relation')
    expected = (np.array(minimal_coefficients) @ np.array(powers[:5])) % PRIME
    if not np.array_equal(expected, powers[5]):
        raise ValueError('minimal relation does not reproduce T^5')

    extensions = {}
    coordinate_vectors = np.eye(len(basis), dtype=np.int64)
    for i, degree in enumerate(basis):
        extensions[str(degree)] = generated_dimension(
            regular, unit, [class_sum, coordinate_vectors[i]])
    completing_pairs = []
    for i, j in itertools.combinations(range(len(basis)), 2):
        if generated_dimension(regular, unit,
                               [class_sum, coordinate_vectors[i],
                                coordinate_vectors[j]]) == len(basis):
            completing_pairs.append([basis[i], basis[j]])
    if completing_pairs != [[2, 6], [3, 5], [3, 6]]:
        raise ValueError('unexpected completing pairs')

    separator = np.zeros(len(basis), dtype=np.int64)
    for degree in (3, 5):
        separator[basis.index(degree)] = 1
    if generated_dimension(
            regular, unit, [class_sum, separator]) != len(basis):
        raise ValueError('the 3+5 separator does not complete the algebra')
    if [3, 5] not in completing_pairs:
        raise ValueError('the 3+5 separator is not a completing pair')

    block_actions = [block_action(regular, class_sum, block)
                     for block in blocks['primitive_central_idempotents']]
    expected_actions = {
        0: (0, 0), 1: (1, 1), 2: (2, 1), 3: (4, 4),
        4: (3, 1), 5: (1, 1), 6: (0, 0), 7: (1, 1),
    }
    expected_eigenvalues = {
        0: [0], 1: [1], 2: [1], 3: [4],
        4: [4], 5: [1], 6: [0], 7: [1],
    }
    for action in block_actions:
        trace, determinant = expected_actions[action['idempotent_index']]
        if (action['trace'], action['determinant']) != (trace, determinant):
            raise ValueError('block action does not match the expected trace')
        if action['eigenvalues'] != expected_eigenvalues[action['idempotent_index']]:
            raise ValueError('block action has an unexpected eigenvalue')
        if action['kind'] == 'dual' and action['regular_rank'] != 2:
            raise ValueError('dual block does not have rank two')

    separator_actions = [block_action(regular, separator, block)
                         for block in blocks['primitive_central_idempotents']]
    expected_separator_matrices = {
        0: [[2]], 1: [[0]], 2: [[4, 4], [4, 0]], 3: [[3]],
        4: [[2, 4], [4, 3]], 5: [[3]], 6: [[4]], 7: [[1]],
    }
    for action in separator_actions:
        if action['action_matrix'] != expected_separator_matrices[
                action['idempotent_index']]:
            raise ValueError('separator has an unexpected block action')
    joint = {
        action['idempotent_index']: (
            expected_eigenvalues[action['idempotent_index']][0],
            action['eigenvalues'][0],
        )
        for action in separator_actions
    }
    expected_joint = {
        0: (0, 2), 1: (1, 0), 2: (1, 2), 3: (4, 3),
        4: (4, 0), 5: (1, 3), 6: (0, 4), 7: (1, 1),
    }
    if joint != expected_joint:
        raise ValueError('unexpected joint class-sum/separator characters')
    if len(set(joint.values())) != len(joint):
        raise ValueError('class sum and separator do not separate the primaries')
    class_sum_groups = {}
    for action in block_actions:
        class_sum_groups.setdefault(
            tuple(action['eigenvalues']), []).append(action['idempotent_index'])
    separator_group_actions = {
        index: separator_actions[index]['eigenvalues']
        for index in range(len(separator_actions))
    }
    for group in class_sum_groups.values():
        if len(group) > 1 and len({tuple(separator_group_actions[i])
                                   for i in group}) == 1:
            raise ValueError('separator does not refine a class-sum block')

    ten = ten_control()
    receipt = {
        'passed': True,
        'prime': PRIME,
        'length': 12,
        'basis_indices': basis,
        'table_sha256': sha256(TABLE.read_bytes()).hexdigest(),
        'replica_swap_coordinates': replica,
        'class_sum_coordinates': class_sum,
        'class_sum_powers': [power.tolist() for power in powers[:6]],
        'class_sum_subalgebra_dimension': 5,
        'class_sum_minimal_polynomial_coefficients_low_to_high':
            minimal_coefficients,
        'class_sum_minimal_polynomial': 'x*(x-1)^2*(x+1)^2',
        'extensions_by_one_coefficient': extensions,
        'completing_coefficient_pairs': completing_pairs,
        'central_idempotent_block_actions': block_actions,
        'separator_coefficients': separator.tolist(),
        'separator_completes_algebra': True,
        'separator_block_actions': separator_actions,
        'class_sum_separator_joint_characters': {
            str(index): [value[0], value[1]] for index, value in joint.items()},
        'class_sum_separator_separates_primaries': True,
        'class_sum_holdouts': holdouts,
        'ten_site_control': ten,
        'scope': (
            'Exact coefficient-algebra computations at twelve sites: the '
            'replica-swap unit C0 and transposition class sum T generate a '
            'five-dimensional commutative subalgebra, with T satisfying '
            'x(x-1)^2(x+1)^2. The eight primitive central idempotents see '
            'eigenvalues 0,1,4, with the two dual blocks carrying nontrivial '
            'nilpotent actions. The coefficient combination C3+C5 completes '
            'the full ten-dimensional algebra and, together with T, separates '
            'the eight central idempotents at the character level; this is a '
            'coefficient-algebra statement, not a split of the cut module '
            'coarse blocks. This is not the eight-dimensional R,T orbit on the '
            'cut module; no cut-module rank, two-cut pairing, rank bound, or '
            'factoring score change is claimed.'
        ),
        'source_sha256': {
            'bounded_spin_crossing.py': sha256(
                (ROOT / 'tools/bounded_spin_crossing.py').read_bytes()).hexdigest(),
            'bounded_spin_six_vertex.py': sha256(
                (ROOT / 'tools/bounded_spin_six_vertex.py').read_bytes()).hexdigest(),
            'check_bounded_spin_transposition_subalgebra.py': sha256(
                (ROOT / 'tools/check_bounded_spin_transposition_subalgebra.py')
                .read_bytes()).hexdigest(),
        },
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({
        'passed': True,
        'class_sum_coordinates': class_sum,
        'subalgebra_dimension': 5,
        'minimal_polynomial': 'x*(x-1)^2*(x+1)^2',
        'completing_pairs': completing_pairs,
        'separator_coefficients': separator.tolist(),
        'separator_joint_characters': {
            str(index): [value[0], value[1]] for index, value in joint.items()},
        'ten_site_dimension': ten['generated_dimension'],
        'holdouts': holdouts,
    }, indent=2))


if __name__ == '__main__':
    main()
