#!/usr/bin/env python3
"""Independent complete polynomial minors for aligned Gaussian complements."""
from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import random

from check_block_smith_crossings import binary_rank
from replay_source_provenance import resolve_source_hashes
from saturated_dual_kernel_directions import literal_directions
from search_f4_recentered_crossings import digest


def multiply(first, second):
    result = 0
    while second:
        bit = second & -second
        second ^= bit
        result ^= first << (bit.bit_length()-1)
    return result


def strict_matrix():
    """Closed strict kernel D=[[uG,w^16 I],[I,uH]], u=1+w."""
    result = [[0]*16 for _ in range(16)]
    for i in range(8):
        for j in range(8):
            g, h = 1 << (2*(i & j)), 1 << (2*(7-(i | j)))
            result[i][j] = g ^ (g << 1)
            result[8+i][8+j] = h ^ (h << 1)
            if i == j:
                result[i][8+j] = 1 << 16
                result[8+i][j] = 1
    return result


def minor(matrix, rows, columns):
    rr = [i for i in range(16) if (rows >> i) & 1]
    cc = [i for i in range(16) if (columns >> i) & 1]
    assert len(rr) == len(cc)

    @lru_cache(None)
    def determinant(mask):
        if not mask:
            return 1
        row = len(rr)-mask.bit_count()
        result = 0
        remaining = mask
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            coefficient = matrix[rr[row]][cc[bit.bit_length()-1]]
            if coefficient:
                result ^= multiply(coefficient, determinant(mask ^ bit))
        return result

    return determinant((1 << len(rr))-1)


def expand_square(polynomial):
    result = 0
    while polynomial:
        bit = polynomial & -polynomial
        polynomial ^= bit
        result |= 1 << (2*(bit.bit_length()-1))
    return result


def run(path):
    record = json.loads(path.read_text())
    provenance, unresolved = resolve_source_hashes(record['source_sha256'])
    assert not unresolved, unresolved
    code = record['code_masks']
    assert code == [15, 77, 43, 23, 51]
    assert digest(record['reference_values']) == record['reference_digest']
    matrix = strict_matrix()
    baseline_checks = 0
    for index, (polynomial, high) in enumerate(record['reference_values']):
        x, y, u, v = ((index//5**leg) % 5 for leg in range(4))
        actual = minor(matrix, code[v] | (code[y] << 8), code[u] | (code[x] << 8))
        assert high == 0 and actual == expand_square(polynomial) << 16
        baseline_checks += 1
    expected_complements = {sum(1 << i for i in positions) for positions in combinations(range(8), 4)}-set(code)
    assert expected_complements == {row['complement_mask'] for row in record['raw_slabs']} | set(record['zero_complement_masks'])
    rng = random.Random(2026103201)
    sampled_minors = endpoint_entries = normalization_identities = direct_entries = 0
    for item in record['raw_slabs']:
        raw, normalized = item['raw_slab_in_z'], item['normalized_slab']
        assert len(raw) == len(normalized) == 125
        assert max(value.bit_length()-1 for value in raw) == item['raw_degree']
        assert max(value.bit_length()-1 for value in normalized) == item['normalized_degree']
        for old, new in zip(raw, normalized):
            assert old.bit_count() % 2 == new.bit_count() % 2 == 0
            assert multiply(item['entry_gcd'], new) == old ^ (old << 1)
            endpoint_entries += 1
            normalization_identities += 1
        for index in rng.sample(range(125), 2):
            y, u, v = ((index//5**leg) % 5 for leg in range(3))
            actual = minor(matrix, code[v] | (code[y] << 8), code[u] | (item['complement_mask'] << 8))
            assert actual == expand_square(raw[index])
            sampled_minors += 1
        direct = all(not (value & 255) for value in raw) and item['raw_degree'] <= 26
        assert direct == item['direct_degree18_direction']
        if direct:
            for value in raw:
                assert (value >> 8).bit_length() <= 19 and (value >> 8).bit_count() % 2 == 0
                direct_entries += 1
    selected = [item for item in record['raw_slabs'] if item['normalized_degree'] <= record['degree']]
    assert [item['complement_mask'] for item in selected] == [item['complement_mask'] for item in record['primitives']]
    assert [item['normalized_slab'] for item in selected] == [item['normalized_slab'] for item in record['primitives']]
    directions = literal_directions(record)
    encoded = []
    for direction in directions:
        assert len(direction) == 625 and all(high == 0 and low.bit_length() <= 19 and low.bit_count() % 2 == 0
                                             for low, high in direction)
        encoded.append(sum(low << (19*entry) for entry, (low, high) in enumerate(direction)))
    rank = binary_rank(encoded)
    return {'record': str(path), 'record_sha256': sha256(path.read_bytes()).hexdigest(), 'status': 'passed',
            'source_provenance': provenance,
            'independent_full_baseline_minors': baseline_checks,
            'independent_complementary_minors': sampled_minors,
            'literal_endpoint_zero_slab_entries': endpoint_entries,
            'exact_polynomial_normalization_identities': normalization_identities,
            'direct_degree18_slab_entries': direct_entries,
            'candidate_directions': len(directions), 'independent_direction_rank_over_F2_and_F4': rank}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run(args.record)
    result['checker_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
