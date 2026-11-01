#!/usr/bin/env python3
"""Exact independent-minor checks for deterministic orthogonal scalings."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path
from random import Random
import time

from tensor_orthogonal_scaling import TowerField, gram, orthogonal_scaling, rank


def determinant(A, F):
    n = len(A)
    value = 0
    for permutation in permutations(range(n)):
        term = 1
        for i, j in enumerate(permutation):
            term = F.mul(term, A[i][j])
        value ^= term
    return value


def add_poly(P, Q):
    R = dict(P)
    for degree, coefficient in Q.items():
        R[degree] = R.get(degree, 0) ^ coefficient
        if not R[degree]:
            del R[degree]
    return R


def multiply_poly(P, Q, F):
    R = {}
    for i, a in P.items():
        for j, b in Q.items():
            R[i + j] = R.get(i + j, 0) ^ F.mul(a, b)
            if not R[i + j]:
                del R[i + j]
    return R


def polynomial_determinant(A, F):
    result = {}
    for permutation in permutations(range(len(A))):
        term = {0: 1}
        for i, j in enumerate(permutation):
            term = multiply_poly(term, A[i][j], F)
        result = add_poly(result, term)
    return result


def basis_polynomial(C, columns, F, assigned_weights, degree_bound):
    selected = []
    for row in C:
        restricted = [row[i] for i in columns]
        if rank(selected + [restricted], F) > len(selected):
            selected.append(restricted)
    s = len(selected)
    Q, weights = {}, []
    for local in combinations(range(len(columns)), s):
        coefficient = determinant([[row[i] for i in local] for row in selected], F)
        if coefficient:
            degree = sum(assigned_weights[columns[i]] for i in local)
            Q[degree] = Q.get(degree, 0) ^ coefficient
            if not Q[degree]:
                del Q[degree]
            weights.append(degree)
    assert Q and weights.count(min(weights)) == 1
    assert min(Q) == min(weights)
    assert max(Q) - min(Q) <= degree_bound
    assert max(Q) <= degree_bound
    # Direct formal determinant, independent of the maximal-minor sum.
    G = []
    for row in selected:
        out = []
        for other in selected:
            entry = {}
            for j, column in enumerate(columns):
                coefficient = F.mul(row[j], other[j])
                if coefficient:
                    exponent = 2 * assigned_weights[column]
                    entry[exponent] = entry.get(exponent, 0) ^ coefficient
                    if not entry[exponent]:
                        del entry[exponent]
            out.append(entry)
        G.append(out)
    assert polynomial_determinant(G, F) == {2 * i: F.mul(a, a) for i, a in Q.items()}
    return s, Q, min(weights)


def slow_mul(a, b, level):
    """Four-product quotient multiplication, separate from the three-product code."""
    if level == 0:
        return a & b
    half = 1 << (level - 1)
    mask = (1 << half) - 1
    a0, a1, b0, b1 = a & mask, a >> half, b & mask, b >> half
    p00 = slow_mul(a0, b0, level - 1)
    p01 = slow_mul(a0, b1, level - 1)
    p10 = slow_mul(a1, b0, level - 1)
    p11 = slow_mul(a1, b1, level - 1)
    return (p00 ^ slow_mul(p11, 1 << (half - 1), level - 1)) | ((p01 ^ p10 ^ p11) << half)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, default=Path('evidence/tensor-orthogonal-scaling-controls.json'))
    args = ap.parse_args()
    start = time.monotonic()
    rng = Random(20261101)
    multiplication_checks = 0
    tower_checks = []
    for level in range(5):
        F = TowerField(level)
        pairs = list(product(range(F.order), repeat=2)) if level <= 2 else [
            (rng.randrange(F.order), rng.randrange(F.order)) for _ in range(160)]
        for a, b in pairs:
            assert F.mul(a, b) == slow_mul(a, b, level)
            multiplication_checks += 1
        witness = 1 << (F.bits - 1)
        assert F.absolute_trace(witness) == 1
        samples = list(range(F.order)) if level <= 2 else [rng.randrange(F.order) for _ in range(32)]
        for a in samples:
            assert F.pow(a, F.order) == a
            if a:
                assert F.mul(a, F.inv(a)) == 1
        tower_checks.append({'level': level, 'absolute_degree': F.bits,
                             'trace_one_element': witness, 'field_samples': len(samples)})

    cases = []
    # All 2-by-4 binary coefficient matrices, including loops and dependent rows.
    for bits in range(256):
        cases.append((0, [[(bits >> (4 * i + j)) & 1 for j in range(4)] for i in range(2)]))
    # The five projective points over F4, with all 4-subsets and sampled orders.
    points = [(1, 0), (1, 1), (1, 2), (1, 3), (0, 1)]
    for omitted in range(5):
        columns = [point for i, point in enumerate(points) if i != omitted]
        for order in list(permutations(range(4)))[::4]:
            cases.append((1, [[columns[j][i] for j in order] for i in range(2)]))
    # Rank-three, repeated-column, zero-column and redundant-row F4 controls.
    for m in (3, 4, 5, 6):
        for _ in range(5):
            cases.append((1, [[rng.randrange(4) for _ in range(m)] for _ in range(3)]))

    subset_checks = 0
    maximum_extension_bits = 0
    examples = []
    for case_index, (base_level, C) in enumerate(cases):
        result = orthogonal_scaling(C, base_level)
        base = TowerField(base_level)
        E = TowerField(result['extension_level'])
        maximum_extension_bits = max(maximum_extension_bits, E.bits)
        d, theta = result['relative_degree'], result['theta']
        if d > 1:
            assert E.pow(theta, 1 << (base.bits * d)) == theta
            assert E.pow(theta, 1 << (base.bits * (d // 2))) != theta
        for mask in range(1 << len(C[0])):
            columns = [i for i in range(len(C[0])) if mask >> i & 1]
            expected_rank, Q, minimum = basis_polynomial(C, columns, base,
                                                        result['weights'], result['degree_bound'])
            value = 0
            for exponent, coefficient in Q.items():
                value ^= E.mul(coefficient, E.pow(theta, exponent - minimum))
            assert value != 0
            scaled = [[row[i] for i in columns] for row in result['scaled_matrix']]
            assert rank(scaled, E) == expected_rank
            assert rank(gram(scaled, E), E) == expected_rank
            subset_checks += 1
        if case_index in (0, 255, 256, len(cases) - 1):
            examples.append({'coefficient_matrix': C, 'construction': result})

    rank_one_controls = []
    for m in range(1, 10):
        result = orthogonal_scaling([[1] * m])
        pivots = {}
        for element in result['scales']:
            v = element
            while v:
                pivot = v.bit_length() - 1
                if pivot in pivots:
                    v ^= pivots[pivot]
                else:
                    pivots[pivot] = v
                    break
            assert v
        assert len(pivots) == m
        rank_one_controls.append({'columns': m, 'absolute_field_degree': result['extension_bits'],
                                  'independent_binary_vectors': len(pivots)})

    corank_one_controls = []
    corank_one_subset_checks = 0
    full_rank_controls = []
    for r in range(1, 8):
        identity = [[int(i == j) for j in range(r)] for i in range(r)]
        unchanged = orthogonal_scaling(identity, 1)
        assert unchanged['relative_degree'] == 1 and unchanged['scales'] == [1] * r
        full_rank_controls.append({'rank': r, 'relative_degree': 1})
        for _ in range(3):
            C = [row + [rng.randrange(4)] for row in identity]
            result = orthogonal_scaling(C, 1)
            assert result['corank'] == 1 and result['degree_bound'] == 1
            assert result['relative_degree'] == 2
            E = TowerField(result['extension_level'])
            for mask in range(1 << (r + 1)):
                selected = [[x for j, x in enumerate(row) if mask >> j & 1]
                            for row in result['scaled_matrix']]
                assert rank(selected, E) == rank(gram(selected, E), E)
                corank_one_subset_checks += 1
            corank_one_controls.append({'matrix': C, 'construction': result})

    report = {'scope': __doc__, 'passed': True, 'complete_binary_2_by_4_matrices': 256,
              'coefficient_matrices': len(cases), 'all_subset_gram_checks': subset_checks,
              'independent_formal_cauchy_binet_checks': subset_checks,
              'unique_minimum_basis_checks': subset_checks,
              'independent_tower_multiplication_checks': multiplication_checks,
              'tower_checks': tower_checks, 'rank_one_controls': rank_one_controls,
              'corank_one_controls': corank_one_controls,
              'corank_one_subset_gram_checks': corank_one_subset_checks,
              'full_column_rank_controls': full_rank_controls,
              'maximum_absolute_extension_degree': maximum_extension_bits, 'examples': examples,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest()
                                for n in (Path(__file__).name, 'tensor_orthogonal_scaling.py')},
              'seconds': time.monotonic() - start}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('tower_checks', 'rank_one_controls', 'corank_one_controls', 'full_column_rank_controls', 'examples', 'source_sha256')}))


if __name__ == '__main__':
    main()
