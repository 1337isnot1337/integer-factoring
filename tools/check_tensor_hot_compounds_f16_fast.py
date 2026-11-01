#!/usr/bin/env python3
"""Independent polynomial and Smith-profile controls for the F16 accelerator.

The binary determinant path is checked against literal permutation sums.
Both complete and capped profiles are compared to the existing constructor;
the small complete profiles also undergo independent binary Toeplitz ranks.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import random
import time

import check_tensor_hot_compounds_f16 as controls
import tensor_hot_compounds_f16 as base
import tensor_hot_compounds_f16_fast as fast


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/tensor-hot-compounds-f16-fast-controls.json'))
    args = parser.parse_args()
    started = time.monotonic()
    rng = random.Random(20261101)
    minor_checks = []
    for modes in (4, 8, 16):
        for charge in range(7):
            if charge > 2 * modes:
                continue
            for _ in range(12):
                rows = sorted(rng.sample(range(2 * modes), charge))
                columns = sorted(rng.sample(range(2 * modes), charge))
                rmask = sum(1 << i for i in rows)
                cmask = sum(1 << i for i in columns)
                matrix = [[(3 << ((i % modes) ^ (j % modes))) if i // modes == j // modes
                           else int(i % modes == j % modes) for j in columns] for i in rows]
                expected = controls.binary_determinant(matrix)
                actual = fast.binary_minor(modes, rmask, cmask)
                assert actual == expected, (modes, rows, columns)
                bound = sum(i % modes for i in rows + columns) + charge
                minor_checks.append({'modes': modes, 'row_mask': rmask, 'column_mask': cmask,
                                     'degree_bound': bound, 'polynomial': actual,
                                     'path': 'subset' if bound < 64 else 'existing fallback'})

    original = base.compound_block
    comparisons = {'calls': 0, 'orders': 0, 'literal_entries': 0}

    def checked(*arguments, **keywords):
        reference = original(*arguments, **keywords)
        accelerated = fast.compound_block(*arguments, **keywords)
        assert {k:v for k,v in reference.items() if k != 'precision_attempts'} == {k:v for k,v in accelerated.items() if k != 'precision_attempts'}
        comparisons['calls'] += 1
        comparisons['orders'] += len(reference['orders_t'])
        if 'matrix' in reference:
            comparisons['literal_entries'] += len(reference['orders_t']) ** 2
        return accelerated

    base.compound_block = checked
    cases = []
    try:
        for coefficient in (2, 4, 5, 8):
            embedding = [{0: 1}] + [{1 << i: 1} for i in range(4)] + [{9: 1, 6: coefficient}]
            scale = base.INV[base.dot(embedding[-1], embedding[-1])]
            left = embedding[:-1] + [{word: base.mul(scale, value)
                                     for word, value in embedding[-1].items()}]
            cases.append(controls.one_case(embedding, left, 4, f'coefficient_{coefficient}'))
        code = base.differential_kernel([{9: 1, 6: 4}], 4, 8)
        left = base.orthogonal_retraction(code)
        assert left is not None
        cases.append(controls.one_case(code['embedding'], left, 4, 'differential_kernel'))
    finally:
        base.compound_block = original

    names = [Path(__file__).name, 'tensor_hot_compounds_f16_fast.py',
             'tensor_hot_compounds_f16_fast.cpp', 'check_tensor_hot_compounds_f16.py',
             'tensor_hot_compounds_f16.py', 'tensor_hot_compounds_f16.cpp',
             'tensor_orthogonal_scaling.py', 'coherent_hot_compounds_f4.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp']
    report = {'passed': True, 'scope': __doc__, 'minor_checks': minor_checks,
              'permutation_minor_checks': len(minor_checks),
              'subset_path_checks': sum(x['path'] == 'subset' for x in minor_checks),
              'fallback_path_checks': sum(x['path'] == 'existing fallback' for x in minor_checks),
              'native_constructor_comparisons': comparisons,
              'independent_toeplitz_rank_checks': sum(x['toeplitz_rank_checks'] for x in cases),
              'small_profile_cases': cases,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'seconds': time.monotonic() - started}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ('minor_checks', 'small_profile_cases', 'source_sha256')}))


if __name__ == '__main__':
    main()
