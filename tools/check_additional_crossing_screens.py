#!/usr/bin/env python3
"""Independent small audits of archived group-control and bipartite screens."""

import argparse
import hashlib
from itertools import product
import json
from math import isqrt, sqrt
from pathlib import Path

from check_modular_crossing import rank_two
from complement_shell_tensor import multiply
from explore_complement_edge_pencils import all_length_labels, minimum_decomposition
from explore_bipartite_pfaffian_codes import coefficient_basis, invariant_partition_space
from pfaffian_particle_hole import local_exterior_table


ROOT = Path(__file__).resolve().parents[1]


def compose(first, second):
    return tuple(first[second[i]] for i in range(4))


def square_action_group():
    """Build D8 from actual permutations of a square, not its native cocycle."""
    identity = tuple(range(4))
    first = tuple((-i) % 4 for i in range(4))
    second = tuple((1 - i) % 4 for i in range(4))
    central = tuple((i + 2) % 4 for i in range(4))
    actions = []
    for value in range(8):
        action = identity
        for bit, generator in ((1, second), (0, first), (2, central)):
            if value >> bit & 1:
                action = compose(action, generator)
        actions.append(action)
    assert len(set(actions)) == 8
    positions = {action: i for i, action in enumerate(actions)}
    table = [[positions[compose(x, y)] for y in actions] for x in actions]
    inverse = [next(y for y in range(8) if table[x][y] == 0) for x in range(8)]
    return table, inverse, [positions[first], positions[second], positions[central]]


def d8_coefficient_matrices(weights, controls, indices):
    table, inverse, generators = square_action_group()
    kernels = []
    for control in controls:
        values = {0: 1}
        for bit, (generator, weight) in enumerate(zip(generators, weights)):
            if not (control >> bit & 1):
                continue
            following = {}
            for value, polynomial in values.items():
                following[value] = following.get(value, 0) ^ (polynomial << weight)
                other = table[value][generator]
                following[other] = following.get(other, 0) ^ polynomial ^ (polynomial << weight)
            values = {value: polynomial for value, polynomial in following.items() if polynomial}
        kernels.append(values)
    rows = [[0] * 64 for _ in indices]
    columns = [[0] * 64 for _ in indices]
    # Enumerate all four legs. This does not reuse the native orbit loops.
    entries = 0
    for x, y, u, v in product(range(8), repeat=4):
        source = table[x][inverse[y]]
        if source != table[v][inverse[u]]:
            continue
        displacement = table[inverse[y]][u]
        value = kernels[source].get(displacement, 0)
        for i, degree in enumerate(indices):
            if value >> degree & 1:
                rows[i][8 * x + y] ^= 1 << (8 * v + u)
                columns[i][8 * x + u] ^= 1 << (8 * y + v)
        entries += 1
    assert entries == 512
    return rows, columns


def combination_ranks(basis):
    values = [0] * len(basis[0])
    ranks = [0] * (1 << len(basis))
    previous = 0
    for step in range(1, len(ranks)):
        gray = step ^ (step >> 1)
        changed = (gray ^ previous).bit_length() - 1
        values = [a ^ b for a, b in zip(values, basis[changed])]
        ranks[gray] = rank_two(values)
        previous = gray
    return ranks


def verify_d8():
    results = []
    for key in ('111', '112', '124'):
        path = ROOT / f'evidence/arbitrary-D8-controls-{key}.jsonl'
        events = [json.loads(line) for line in path.read_text().splitlines()]
        start, final = events[0], events[-1]
        assert final['event'] == 'complete'
        assert final['permutations_tested'] == 40320
        assert sum(count for lower, count in final['integer_lower_histogram']) == 40320
        assert min(lower for lower, count in final['integer_lower_histogram'] if count) == final['global_integer_lower']
        row, column = d8_coefficient_matrices(start['weights'], final['best_control_word_permutation'], start['basis_indices'])
        rr, cc = combination_ranks(row), combination_ranks(column)
        assert rr == start['row_ranks'] and cc == final['best_column_ranks']
        target = start['target']
        assert rr[target] == cc[target] == 64
        lower, _ = minimum_decomposition([isqrt(r * c) for r, c in zip(rr, cc)], target=target)
        cost, _ = minimum_decomposition([sqrt(r * c) for r, c in zip(rr, cc)], target=target)
        assert lower == final['global_integer_lower'] == 50
        assert abs(cost - final['best_cost']) < 1e-10
        assert lower ** 70 > 8 ** 99
        witness = 0
        for mask in final['best_pieces']:
            witness ^= mask
        assert witness == target
        results.append({'weights': start['weights'], 'dimension': start['dimension'],
                        'independent_dense_combinations': len(rr), 'integer_lower': lower,
                        'best_cost': cost, 'native_log_sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    return results


def verify_bipartite_dense():
    records = []
    for pairs, length, index in ((2, 2, 0), (2, 2, 1), (3, 3, 18)):
        labels = all_length_labels(pairs, length)
        count, unrank = invariant_partition_space(labels)
        first, second = unrank(index)
        report = coefficient_basis(labels, first, second)
        d, stride = len(first), length + 1
        degrees = [a + stride * b for a, b, c in report['basis_monomials_p_q_r']]
        local = local_exterior_table(2, 1 << stride, 1)
        rows = [[0] * (d * d) for _ in degrees]
        columns = [[0] * (d * d) for _ in degrees]
        for x, y, u, v in product(range(d), repeat=4):
            polynomial = 1
            for a, b, c, e in zip(labels[first[x]], labels[second[y]], labels[second[u]], labels[first[v]]):
                polynomial = multiply(polynomial, local[((a * 4 + b) * 4 + c) * 4 + e])
            assert polynomial & 1 == int(x == v and y == u)
            assert not (polynomial >> length) & 1
            assert not (polynomial >> (length * stride)) & 1
            for i, degree in enumerate(degrees):
                if polynomial >> degree & 1:
                    rows[i][x * d + y] ^= 1 << (v * d + u)
                    columns[i][x * d + u] ^= 1 << (v * d + y)
        full = [combination_ranks(rows), combination_ranks(columns)]
        for cut, key in enumerate(('basis_row_blocks', 'basis_column_blocks')):
            blocks = [[list(map(lambda v: int(v, 16), block)) for block in tensor] for tensor in report[key]]
            blocked = [0] * (1 << len(blocks))
            for block in range(len(report['block_sizes'])):
                values = combination_ranks([tensor[block] for tensor in blocks])
                blocked = [a + b for a, b in zip(blocked, values)]
            assert full[cut] == blocked
        lower, _ = minimum_decomposition([isqrt(r * c) for r, c in zip(*full)], target=1)
        assert lower ** 70 > d ** 99
        records.append({'pairs': pairs, 'length': length, 'partition_index': index,
                        'partitions_available': count, 'labels': d,
                        'independent_dense_entries': d ** 4, 'all_binary_combinations': len(full[0]),
                        'integer_lower': lower})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    sources = ['explore_arbitrary_group_controls.cpp', 'explore_bipartite_pfaffian_codes.py',
               'explore_mixed_radix_group_crossings.py']
    report = {'status': 'pass',
              'scope': 'independent winning-profile and small dense audits; native exhaustive ranges recorded separately',
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'source_sha256': {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in sources},
              'D8_square_action_checks': verify_d8(),
              'bipartite_exterior_dense_checks': verify_bipartite_dense()}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps(report, sort_keys=True))


if __name__ == '__main__':
    main()
