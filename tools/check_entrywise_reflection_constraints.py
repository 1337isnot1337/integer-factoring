#!/usr/bin/env python3
"""Exhaustive small-field controls and exact full-tensor jet constraints."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random

from entrywise_reflection_constraints import field_dot, local_affine, solve_target
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import add, digest, scale
from screen_affine_finite_places import bitplane_rank, rows
from search_entrywise_reflection_jets import basis_catalogue, bell_jets


def packed(values):
    return (sum((value & 1) << index for index, value in enumerate(values)),
            sum((value >> 1) << index for index, value in enumerate(values)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(312463)
    local_records, compared = [], 0
    for trial in range(36):
        width = 2 + trial % 4
        equations = [(packed([rng.randrange(4) for _ in range(width)]), rng.randrange(4))
                     for _ in range(1 + trial % 7)]
        actual = {packed(values) for values in product(range(4), repeat=width)
                  if all(field_dot(vector, packed(values)) == rhs for vector, rhs in equations)}
        compared += 4 ** width
        result = local_affine(equations, width)
        if result is None:
            assert not actual
            generated = set()
        else:
            particular, kernel, rank = result
            generated = set()
            for coefficients in product(range(4), repeat=len(kernel)):
                value = particular
                for scalar, direction in zip(coefficients, kernel):
                    value = add(value, scale(direction, scalar))
                generated.add(value)
            assert len(generated) == 4 ** (width - rank)
            assert generated == actual
        local_records.append({'trial': trial, 'width': width, 'equations': len(equations),
                              'solution_count': len(actual)})
    base = tensor(4)
    matrices = basis_catalogue(base)
    targets = [
        (24, 100, 0, 'exact', [], []),
        (16, 100, 0, 'exact', [], []),
        (24, 8, 16, 'exact', [(2, 0, (0, 1), (2, 3), 2), (3, 1, (0, 1), (2, 3), 2)], []),
        (16, 0, 0, 'exact', [(2, 0, (0, 1), (2, 3), 8), (3, 1, (0, 1), (2, 3), 8)], []),
        (16, 0, 0, 'exact', [], [bell_jets(matrices['source_Bell'], 0, 2, 6),
                                 bell_jets(matrices['identity'], 1, 3, 6)]),
        (24, 8, 16, 'zero-pattern', [], [bell_jets(matrices['source_Bell'], 0, 2, 4),
                                        bell_jets(matrices['source_Bell'], 1, 3, 4)]),
    ]
    records = []
    for cap, zero, infinity, mode, conditions, prescribed in targets:
        result = solve_target(base, cap, conditions, zero, infinity, 16, mode, prescribed)
        seed = result.pop('seed', None)
        record = {'cap': cap, 'zero_cutoff': zero, 'infinity_cutoff': infinity,
                  'mode': mode, 'rectangle_conditions': conditions,
                  'prescribed_orders': [[item['point'], len(item['values'])] for item in prescribed],
                  **result}
        if seed is not None:
            finite_ranks = [[bitplane_rank(rows(seed['values'], cut, point)) for cut in (0, 1)]
                            for point in (2, 3)]
            record.update(actual_degree=seed['degree'], seed_sha256=digest(seed['values']),
                          endpoint_entries=256, derivative_entries=256, paired_trace_entries=32,
                          finite_ranks=finite_ranks)
        records.append(record)
    assert records[0]['feasible'] and records[0]['seed_sha256'] == digest(base['values'])
    assert not records[1]['feasible']
    assert records[2]['feasible'] and records[3]['feasible'] and records[4]['feasible']
    assert records[2]['finite_ranks'] == records[3]['finite_ranks'] == [[4, 4], [4, 4]]
    assert records[4]['finite_ranks'][0][0] == records[4]['finite_ranks'][1][1] == 1
    names = (Path(__file__).name, 'entrywise_reflection_constraints.py',
             'search_entrywise_reflection_jets.py', 'solve_companion_affine_degree.py',
             'reflection_dual_rail_gaussian.py', 'screen_affine_finite_places.py')
    report = {'status': 'passed', 'scope': 'Complete exact local solution spaces and full endpoint, derivative, paired-trace and Hasse constraints; no pressure certificate.',
              'exhaustive_local_systems': len(local_records), 'brute_force_assignments': compared,
              'local_systems': local_records, 'full_tensor_cases': records,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ('source_sha256', 'local_systems', 'full_tensor_cases')}))


if __name__ == '__main__':
    main()
