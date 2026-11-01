#!/usr/bin/env python3
"""Check nonorthogonal literal cells and exhaustive projective enumeration."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

import tensor_hot_compounds_f16 as field
from check_tensor_literal_hot_pfaffians_f16 import check_map
from tensor_literal_rank_one_cell_search import all_cases, cell_code, orbit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/tensor-literal-rank-one-cells-controls.json'))
    args = parser.parse_args()
    started = time.monotonic()
    enumerations = []
    for size in (2, 4, 16):
        projective = set()
        for p in product(range(size), repeat=2):
            if p == (0, 0):
                continue
            first = next(x for x in p if x)
            inv = next(x for x in range(1, size) if field.mul(x, first) == 1)
            projective.add(tuple(field.mul(x, inv) for x in p))
        expected = {(p, left) for p in projective for left in product(range(size), repeat=2)
                    if field.mul(p[0], left[0]) ^ field.mul(p[1], left[1]) == 1}
        generated = set()
        for representative in all_cases(size):
            for case in orbit(representative):
                code = cell_code(10, 9, case)
                generated.add((tuple(code['physical_cell_embedding']), tuple(code['physical_cell_retraction'])))
        assert generated == expected
        enumerations.append({'field_size': size, 'projective_lines': len(projective),
                             'retraction_pairs': len(expected),
                             'frobenius_representatives': len(all_cases(size))})
    cases = []
    for height, case in [(8, (0, 1, 2)), (9, (0, 1, 2)), (9, (0, 2, 4))]:
        code = cell_code(10, height, case)
        assert not code['orthogonal']
        assert all(not value['deletion_closed'] for value in code['deletion_diagnostics'].values())
        code['retraction_kind'] = 'nonorthogonal'
        cases.append(check_map(code, f'height-{height}-case-{case}'))
    names = [Path(__file__).name, 'tensor_literal_rank_one_cell_search.py',
             'check_tensor_literal_hot_pfaffians_f16.py', 'reflection_dual_rail_gaussian.py',
             'tensor_literal_hot_pfaffians_f16.py', 'tensor_literal_hot_pfaffians_f16.cpp',
             'tensor_hot_boundary_retractions.py', 'tensor_hot_compounds_f16_fast.cpp',
             'tensor_hot_compounds_f16.py', 'tensor_hot_compounds_f16.cpp',
             'tensor_orthogonal_scaling.py', 'check_tensor_hot_compounds_f16.py',
             'coherent_hot_compounds_f4.cpp', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp']
    report = {'passed': True, 'enumeration_checks': enumerations, 'cases': cases,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'seconds': time.monotonic()-started}
    for key in ['literal_polynomial_entries', 'independent_toeplitz_rank_checks',
                'capped_order_comparisons', 'endpoint_checks', 'infinity_reciprocity_checks']:
        report[key] = sum(case[key] for case in cases)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('cases', 'source_sha256')}), flush=True)


if __name__ == '__main__':
    main()
