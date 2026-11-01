#!/usr/bin/env python3
"""Enumerate all seven-column F4 Gram systems containing a fixed basis.

Every rank-three system can be normalized to contain that basis by a row
change. The search is exhaustive for unordered nonzero column sets in
this normalization, not a list of distinct GL3 orbits or physical kernels.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path
import time

from check_tensor_orthogonal_capacity import Field, subset_grams


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/tensor-f4-rank-three-catalogue.json'))
    args = parser.parse_args()
    started = time.monotonic()
    field = Field(2, 0b111)
    mul = field.table
    vectors = [(a, b, c) for a, b, c in product(range(4), repeat=3)]
    basis = [vectors.index(v) for v in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
    positions = [(i, j) for i in range(3) for j in range(i, 3)]
    outers = [sum(mul[v[i]][v[j]] << (2 * k) for k, (i, j) in enumerate(positions))
              for v in vectors]
    gram_ranks = []
    for encoded in range(1 << 12):
        matrix = [[0] * 3 for _ in range(3)]
        for k, (i, j) in enumerate(positions):
            matrix[i][j] = matrix[j][i] = encoded >> (2 * k) & 3
        gram_ranks.append(field.rank(matrix))
    pair_independent = [[any(mul[vectors[a][i]][vectors[b][j]]
                             ^ mul[vectors[a][j]][vectors[b][i]]
                             for i in range(3) for j in range(i))
                         for b in range(64)] for a in range(64)]

    def determinant(a, b, c):
        x, y, z = vectors[a], vectors[b], vectors[c]
        return (mul[x[0]][mul[y[1]][z[2]] ^ mul[y[2]][z[1]]]
                ^ mul[x[1]][mul[y[0]][z[2]] ^ mul[y[2]][z[0]]]
                ^ mul[x[2]][mul[y[0]][z[1]] ^ mul[y[1]][z[0]]])

    def extend(records, column):
        updates = []
        for gram, rank, first, second in records:
            new_gram = gram ^ outers[column]
            if rank == 0:
                new_rank, a, b = 1, column, 0
            elif rank == 1 and pair_independent[first][column]:
                new_rank, a, b = 2, first, column
            elif rank == 2 and determinant(first, second, column):
                new_rank, a, b = 3, first, second
            else:
                new_rank, a, b = rank, first, second
            if gram_ranks[new_gram] != new_rank:
                return None
            updates.append((new_gram, new_rank, a, b))
        return records + updates

    records = [(0, 0, 0, 0)]
    for column in basis:
        records = extend(records, column)
        assert records is not None
    candidates = [v for v in range(1, 64) if v not in basis
                  and extend(records, v) is not None]
    assert len(candidates) == 47
    counts = Counter({0: 1})
    attempts = Counter()
    component_counts, hilbert_counts = Counter(), Counter()
    representatives = {}
    largest_sets = []
    excluded_completions = 0

    def classify(columns):
        parent = list(range(3))

        def find(i):
            while parent[i] != i:
                i = parent[i]
            return i

        for column in columns:
            support = [i for i, x in enumerate(vectors[column]) if x]
            for i in support[1:]:
                parent[find(i)] = find(support[0])
        component_ranks = Counter(find(i) for i in range(3))
        component_widths = Counter(find(next(i for i, x in enumerate(vectors[v]) if x))
                                   for v in columns)
        return tuple(sorted((rank, component_widths[key])
                            for key, rank in component_ranks.items()))

    def visit(selected, state, start):
        nonlocal excluded_completions
        depth = len(selected)
        if depth == 4:
            columns = basis + selected
            kind = classify(columns)
            dimension = sum(1 << (7 - rank) for _, rank, _, _ in state)
            component_counts[kind] += 1
            hilbert_counts[dimension] += 1
            representatives.setdefault(str(kind), [[vectors[v][i] for v in columns] for i in range(3)])
            largest_sets.append(selected[:])
            return
        for index in range(start, len(candidates) - (3 - depth)):
            attempts[depth + 1] += 1
            column = candidates[index]
            new_state = extend(state, column)
            if new_state is not None:
                counts[depth + 1] += 1
                visit(selected + [column], new_state, index + 1)
            else:
                excluded_completions += comb(len(candidates) - index - 1, 3 - depth)

    visit([], records, 0)
    assert len(largest_sets) == 12
    assert excluded_completions + len(largest_sets) == comb(47, 4)
    assert component_counts == {((1, 2), (2, 5)): 12}
    independent_subset_checks = 0
    for extra in largest_sets:
        matrix = [[vectors[v][i] for v in basis + extra] for i in range(3)]
        independent_subset_checks += len(subset_grams(field, matrix))
    report = {'passed': True, 'scope': __doc__, 'basis_columns': [vectors[v] for v in basis],
              'candidate_extra_columns': [vectors[v] for v in candidates],
              'prefix_counts': dict(sorted(counts.items())),
              'extension_attempts': dict(sorted(attempts.items())),
              'all_four_extra_choices_before_single_column_filter': comb(60, 4),
              'all_four_extra_choices_after_single_column_filter': comb(47, 4),
              'completions_excluded_by_bad_prefix': excluded_completions,
              'independent_subset_gram_checks': independent_subset_checks,
              'component_types': [{'components_rank_width': kind, 'sets': count}
                                  for kind, count in sorted(component_counts.items())],
              'full_fourteen_mode_kernel_dimensions': dict(sorted(hilbert_counts.items())),
              'representatives_by_component_type': representatives,
              'largest_extra_column_indices': largest_sets,
              'vector_indexing': vectors, 'elapsed_seconds': time.monotonic() - started,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, 'check_tensor_orthogonal_capacity.py')}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: report[key] for key in ('prefix_counts', 'extension_attempts',
                      'component_types', 'full_fourteen_mode_kernel_dimensions', 'elapsed_seconds')}))


if __name__ == '__main__':
    main()
