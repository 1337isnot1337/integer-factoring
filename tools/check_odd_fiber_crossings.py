#!/usr/bin/env python3
"""Independent controls for coherent odd-fiber physical crossing codes."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_physical_split_crossings import check_case, literal_retraction
from physical_split_crossings import retract
from reflection_dual_rail_gaussian import tensor
from search_physical_code_gauges import odd_fiber_maps


def check():
    cases = []
    small = tensor(4)
    maps = odd_fiber_maps([(0, 1, 2), (3,)], 4)
    cases.append(check_case(small, maps, maps, [0, 8]))
    large = tensor(8)
    for groups, energy in [
        ([(0, 1, 3), (4,), (8,)], [0, 8, 24]),
        ([(0, 1, 3), (4, 5, 7), (8, 9, 11)], [0, 8, 24]),
        ([(0, 1, 3), (2,), (4,), (8,), (12,)], [0, 0, 4, 8, 12]),
    ]:
        maps = odd_fiber_maps(groups, 16)
        cases.append(check_case(large, maps, maps, energy))
    triple = (7, 11, 15)
    groups = [triple] + [(symbol,) for symbol in range(16) if symbol not in triple]
    maps = odd_fiber_maps(groups, 16)
    projected = retract(large, maps, maps)
    assert projected['values'] == literal_retraction(large, maps, maps, maps, maps)
    rejected = 0
    for groups in [[], [(0, 1)], [(0, 1, 2), (2,)], [(16,)]]:
        try:
            odd_fiber_maps(groups, 16)
        except ValueError:
            rejected += 1
    assert rejected == 4
    sources = [Path(__file__), Path(__file__).with_name('search_physical_code_gauges.py'),
               Path(__file__).with_name('physical_split_crossings.py'),
               Path(__file__).with_name('reflection_dual_rail_gaussian.py'),
               Path(__file__).with_name('check_physical_split_crossings.py'),
               Path(__file__).with_name('check_block_smith_crossings.py')]
    return {'status': 'pass', 'cases': cases,
            'endpoint_entries': sum(case['endpoint_entries'] for case in cases) + len(projected['values']),
            'independent_profiles': 2 * len(cases),
            'toeplitz_checks': sum(profile['toeplitz_checks'] for case in cases
                                   for profile in case['profiles']),
            'larger_literal_contraction': {'dimensions': projected['dimensions'],
                                           'triple': triple, 'entries': len(projected['values'])},
            'invalid_fibers_rejected': rejected,
            'source_sha256': {path.name: sha256(path.read_bytes()).hexdigest() for path in sources},
            'scope': 'Exact odd-fiber split maps, literal full contractions, complete endpoints, and independent binary Smith/Toeplitz controls; numerical optimization is separate.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = check()
    rendered = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
