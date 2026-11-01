#!/usr/bin/env python3
"""Complete F4 rank-two all-subset Gram catalogue through six columns.

Sets contain distinct nonzero vectors. Repeated vectors already fail the
two-column condition. The accompanying dimension/counting proof excludes
larger sets; this finite catalogue is not a factoring-pressure claim.
"""
from hashlib import sha256
from itertools import combinations, permutations, product
import argparse
import json
from pathlib import Path

from coherent_hot_compounds_f4 import mul
from differential_pair_hilbert import hilbert
from orthogonal_differential_pairs_f4 import subset_gram_certificate


VECTORS = [(a, b) for a in range(4) for b in range(4) if a or b]


def outer(vector):
    a, b = vector
    return mul(a, a) | (mul(a, b) << 2) | (mul(b, b) << 4)


def gram_rank(encoded):
    if not encoded:
        return 0
    a, b, c = encoded & 3, (encoded >> 2) & 3, (encoded >> 4) & 3
    return 2 if mul(a, c) ^ mul(b, b) else 1


def valid(columns):
    # Each record is one previously scanned subset, with its exact outer sum,
    # one column spanning its first direction, and its column-space rank.
    records = [(0, None, 0)]
    for vector in columns:
        a, b = vector
        updates = []
        for encoded, first, column_rank in records:
            if not column_rank:
                new_rank = 1
            else:
                new_rank = 2 if column_rank == 2 or mul(first[0], b) ^ mul(first[1], a) else 1
            new_gram = encoded ^ outer(vector)
            if gram_rank(new_gram) != new_rank:
                return False
            updates.append((new_gram, vector if first is None else first, new_rank))
        records.extend(updates)
    return True


def direction(vector):
    a, b = vector
    return (1, mul(b, mul(a, a))) if a else (0, 1)


def canonical_rows(matrix):
    """Reduced row echelon form, retaining the specified column order."""
    A = [row[:] for row in matrix]
    cursor = 0
    for j in range(len(A[0])):
        pivot = next((i for i in range(cursor, len(A)) if A[i][j]), None)
        if pivot is None:
            continue
        A[pivot], A[cursor] = A[cursor], A[pivot]
        inverse = mul(A[cursor][j], A[cursor][j])
        A[cursor] = [mul(inverse, x) for x in A[cursor]]
        for i in range(len(A)):
            if i != cursor:
                coefficient = A[i][j]
                A[i] = [x ^ mul(coefficient, y) for x, y in zip(A[i], A[cursor])]
        cursor += 1
        if cursor == len(A):
            break
    return tuple(tuple(row) for row in A[:cursor])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, default=Path('evidence/tensor-f4-rank-two-catalogue.json'))
    args = ap.parse_args()
    counts, retained = {}, {}
    attempted = 0
    for size in range(1, 7):
        good = []
        for indices in combinations(range(15), size):
            attempted += 1
            columns = [VECTORS[i] for i in indices]
            if valid(columns):
                good.append(columns)
        counts[size] = len(good)
        retained[size] = good
    assert counts == {1: 15, 2: 105, 3: 390, 4: 450, 5: 18, 6: 0}

    largest = {tuple(columns) for columns in retained[5]}
    C = [[1, 1, 2, 2, 0], [0, 1, 3, 1, 3]]
    seed = list(zip(*C))
    assert valid(seed)
    independent_checks = 0
    for columns in retained[5]:
        matrix = [[v[i] for v in columns] for i in (0, 1)]
        assert subset_gram_certificate(matrix)['all_subsets_nondegenerate']
        assert len({direction(v) for v in columns}) == 5
        independent_checks += 32
    orbit = set()
    group_elements = 0
    for a, b, c, d in product(range(4), repeat=4):
        if not (mul(a, d) ^ mul(b, c)):
            continue
        group_elements += 1
        image = tuple(sorted((mul(a, x) ^ mul(b, y), mul(c, x) ^ mul(d, y)) for x, y in seed))
        orbit.add(image)
    assert group_elements == 180 and orbit == largest

    ordered_kernels = set()
    for columns in retained[5]:
        for permutation in permutations(columns):
            ordered_kernels.add(canonical_rows([[v[i] for v in permutation] for i in (0, 1)]))
    assert len(ordered_kernels) == 12
    frobenius_orbits = set()
    for rows in ordered_kernels:
        conjugate = canonical_rows([[mul(x, x) for x in row] for row in rows])
        assert conjugate in ordered_kernels and conjugate != rows
        frobenius_orbits.add(min(rows, conjugate))
    assert len(frobenius_orbits) == 6

    # Independently enumerate all symmetric rank-one matrices and the
    # rank-one outer-product parametrization used in the six-column proof.
    rank_one_matrices = {value for value in range(64) if gram_rank(value) == 1}
    assert len(rank_one_matrices) == 15
    assert rank_one_matrices == {outer(v) for v in VECTORS}
    pairs = [(i, 9 - i) for i in range(5)]
    graded, subset_ranks = hilbert(pairs, C, 10)
    assert sum(graded.values()) == 320
    assert all(rank == min(2, mask.bit_count()) for mask, rank in enumerate(subset_ranks))

    report = {'scope': __doc__, 'passed': True, 'nonzero_vectors': VECTORS,
              'complete_sets_examined': attempted, 'valid_sets_by_size': counts,
              'largest_sets': retained[5], 'independent_subset_gram_checks': independent_checks,
              'gl2_group_elements': group_elements, 'largest_set_orbit_size': len(orbit),
              'ordered_kernels_modulo_row_operations': len(ordered_kernels),
              'frobenius_profile_classes': len(frobenius_orbits),
              'frobenius_profile_representatives': sorted(frobenius_orbits),
              'symmetric_nonzero_rank_one_matrices': len(rank_one_matrices),
              'explicit_coefficient_matrix': C, 'pairs': pairs,
              'full_ten_mode_kernel_dimension': sum(graded.values()),
              'weighted_hilbert_cells': [{'height': h, 'occupation': k, 'dimension': value}
                                         for (h, k), value in sorted(graded.items())],
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest()
                                for n in (Path(__file__).name, 'coherent_hot_compounds_f4.py',
                                          'orthogonal_differential_pairs_f4.py', 'differential_pair_hilbert.py')}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in
                      ('nonzero_vectors', 'largest_sets', 'weighted_hilbert_cells', 'source_sha256')}))


if __name__ == '__main__':
    main()
