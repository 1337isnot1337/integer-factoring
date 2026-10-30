#!/usr/bin/env python3
"""Exact complete-piece recombinations and structured binary label screens.

The profile evaluator is reused unchanged from PR #22, commit 735b3dd.
Finite rank improvements do not establish an asymptotic factoring improvement.
"""

import argparse
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path

from check_modular_crossing import regular_partitions, rank_two
from check_partition_channels import rank_mod
from modular_pfaffian import final_piece_profile, polynomial_product


CROSSING_SYMMETRIES = tuple(
    p for p in permutations(range(4))
    if {frozenset((p[0], p[3])), frozenset((p[1], p[2]))}
    == {frozenset((0, 3)), frozenset((1, 2))})


def pair_labels(width, spacing=8, offset=0):
    """Two equal-sum choices per block; weight O(width^2) for fixed spacing."""
    assert width >= 1 and spacing >= 8 and offset >= 0 and spacing % 8 == offset % 8 == 0
    return tuple(tuple(sorted(
        part for j in range(width)
        for part in ((spacing*j + offset + 1, spacing*j + offset + 7)
                     if (x >> j) & 1 else
                     (spacing*j + offset + 3, spacing*j + offset + 5))))
        for x in range(1 << width))


def encode(rows, dimension):
    return sum(row << (i * dimension**2) for i, row in enumerate(rows))


def decode(vector, dimension):
    mask = (1 << dimension**2) - 1
    return [(vector >> (i * dimension**2)) & mask
            for i in range(dimension**2)]


def permute(rows, dimension, permutation):
    answer = [0] * dimension**2
    for x, y, u, v in product(range(dimension), repeat=4):
        coordinates = (x, y, u, v)
        a, b, c, d = (coordinates[k] for k in permutation)
        answer[x*dimension+y] |= (
            (rows[a*dimension+b] >> (c*dimension+d)) & 1
        ) << (u*dimension+v)
    return answer


def ranks(rows, dimension):
    return (rank_two(rows), rank_two(permute(rows, dimension, (0, 2, 1, 3))))


def crossing(dimension):
    return [1 << (y*dimension+x)
            for x in range(dimension) for y in range(dimension)]


def profile(labels):
    dimension = len(labels)
    weight = sum(labels[0])
    assert len(set(labels)) == dimension
    assert all(sum(label) == weight and len(set(label)) == len(label)
               and all(part > 0 and part % 2 for part in label) for label in labels)
    matrices = {}
    for x, y, u, v in product(range(dimension), repeat=4):
        polynomial = final_piece_profile(labels[x], labels[v], labels[u], labels[y])
        assert polynomial.bit_count() % 2 == int(x == v and y == u)
        while polynomial:
            bit = polynomial & -polynomial
            shift = bit.bit_length() - 1 - 2*weight
            matrices.setdefault(shift, [0] * dimension**2)[x*dimension+y] |= (
                1 << (u*dimension+v))
            polynomial ^= bit
    total = [0] * dimension**2
    for rows in matrices.values():
        total = [a ^ b for a, b in zip(total, rows)]
    assert total == crossing(dimension)
    polynomial_product.cache_clear()
    return matrices


def insert(basis, vector, mask):
    while vector:
        pivot = vector.bit_length() - 1
        if pivot not in basis:
            basis[pivot] = (vector, mask)
            return True
        other, other_mask = basis[pivot]
        vector ^= other
        mask ^= other_mask
    return False


def represent(basis, target):
    mask = 0
    while target:
        pivot = target.bit_length() - 1
        if pivot not in basis:
            return None
        vector, other_mask = basis[pivot]
        target ^= vector
        mask ^= other_mask
    return mask


def candidates(matrices, dimension, pairs):
    answer = {}
    for shift, rows in sorted(matrices.items()):
        for permutation in CROSSING_SYMMETRIES:
            vector = encode(permute(rows, dimension, permutation), dimension)
            if vector:
                answer.setdefault(vector, ((shift, permutation),))
    if pairs:
        original = tuple(answer.items())
        for i, (left, left_description) in enumerate(original):
            for right, right_description in original[:i]:
                vector = left ^ right
                if vector:
                    answer.setdefault(vector, left_description + right_description)
    return answer


def select(candidate_descriptions, dimension):
    """Minimum feasible rank-product cutoff within exactly this candidate list."""
    costs = {}
    for vector in candidate_descriptions:
        row_rank, column_rank = ranks(decode(vector, dimension), dimension)
        costs[vector] = (row_rank * column_rank, row_rank, column_rank)
    ordered = sorted(costs, key=lambda vector: (costs[vector][0], vector))
    basis, selected = {}, []
    target = encode(crossing(dimension), dimension)
    for vector in ordered:
        if insert(basis, vector, 1 << len(selected)):
            selected.append(vector)
        witness = represent(basis, target)
        if witness is not None:
            cutoff = costs[vector][0]
            lower_basis = {}
            for other in ordered:
                if costs[other][0] >= cutoff:
                    break
                insert(lower_basis, other, 0)
            assert represent(lower_basis, target) is None
            chosen = [v for i, v in enumerate(selected) if (witness >> i) & 1]
            total = 0
            for v in chosen:
                total ^= v
            assert total == target
            assert chosen and (total ^ chosen[0]) != target
            record = [{
                'row_rank': costs[v][1], 'column_rank': costs[v][2],
                'terms': [{'shift': s, 'permutation': list(p)}
                          for s, p in candidate_descriptions[v]],
            } for v in chosen]
            return {
                'candidate_count': len(candidate_descriptions),
                'minimum_cutoff_in_this_list': cutoff,
                'strictly_lower_cutoff_span_rank': len(lower_basis),
                'accepted_span_rank': len(basis),
                'witness': record, 'witness_piece_count': len(chosen),
                'omission_negative_control': True,
            }
    raise AssertionError('candidate list does not span the crossing')


def replay(record, matrices, dimension):
    total = [0] * dimension**2
    for piece in record['witness']:
        value = [0] * dimension**2
        for term in piece['terms']:
            rows = permute(matrices[term['shift']], dimension, term['permutation'])
            value = [a ^ b for a, b in zip(value, rows)]
        assert ranks(value, dimension) == (piece['row_rank'], piece['column_rank'])
        # Independent list-of-entries elimination checks the packed-bit ranks.
        for rows, expected in ((value, piece['row_rank']),
                               (permute(value, dimension, (0, 2, 1, 3)),
                                piece['column_rank'])):
            dense = [[(row >> j) & 1 for j in range(dimension**2)] for row in rows]
            assert rank_mod(dense, 2) == expected
        total = [a ^ b for a, b in zip(total, value)]
    assert total == crossing(dimension)


def regular_case(weight, limit=None):
    labels = regular_partitions(2, weight)
    if limit:
        labels = labels[:limit]
    matrices = profile(labels)
    dimension = len(labels)
    single = select(candidates(matrices, dimension, False), dimension)
    pair = select(candidates(matrices, dimension, True), dimension)
    replay(single, matrices, dimension)
    replay(pair, matrices, dimension)
    return {'weight': weight, 'dimension': dimension, 'labels': labels,
            'single_orbit': single, 'pair_orbit': pair}


def structured_case(width, spacing, offset):
    labels = pair_labels(width, spacing, offset)
    matrices = profile(labels)
    dimension = len(labels)
    records = [{'shift': s, 'row_rank': r, 'column_rank': c}
               for s, rows in sorted(matrices.items())
               for r, c in [ranks(rows, dimension)]]
    return {'width': width, 'spacing': spacing, 'offset': offset,
            'weight': sum(labels[0]), 'dimension': dimension,
            'max_rank_product': max(p['row_rank']*p['column_rank'] for p in records),
            'max_width_sum': max((p['row_rank']-1).bit_length()
                                 + (p['column_rank']-1).bit_length() for p in records),
            'ranks': records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    # Independent source recursion checks the complete profile at weight eight.
    from check_modular_crossing import normalized_piece
    labels = regular_partitions(2, 8)
    matrices = profile(labels)
    source_entries = 0
    for shift in range(-16, 17):
        source = normalized_piece(2, 8, shift, labels)
        rows = matrices.get(shift, [0]*4)
        for x, y, u, v in product(range(2), repeat=4):
            assert source[x*2+y][v*2+u] == ((rows[x*2+y] >> (u*2+v)) & 1)
            source_entries += 1
    regular = [regular_case(h) for h in (8, 12, 16)]
    if args.full:
        regular += [regular_case(20), regular_case(24, 8)]
    settings = ((8, 0), (16, 0), (32, 0), (8, 8), (8, 32))
    structured = [structured_case(width, spacing, offset)
                  for width in range(1, 4 if args.full else 3)
                  for spacing, offset in settings]
    result = {
        'status': 'exact finite checks; no asymptotic improvement established',
        'full': args.full, 'source_profile_entries': source_entries,
        'profile_source_commit': '735b3ddb37d2e63dc469cf833687355e42097dd4',
        'sha256': {name: sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
                   for name in ('modular_pfaffian.py', Path(__file__).name)},
        'regular_cases': regular, 'structured_cases': structured,
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({
        'status': result['status'], 'full': args.full,
        'source_profile_entries': source_entries,
        'regular_cutoffs': [(case['weight'], case['dimension'],
                             case['single_orbit']['minimum_cutoff_in_this_list'],
                             case['pair_orbit']['minimum_cutoff_in_this_list']) for case in regular],
        'structured_cases': len(structured),
    }, indent=2))


if __name__ == '__main__':
    main()
