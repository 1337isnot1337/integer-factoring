#!/usr/bin/env python3
"""Exact cross-controls for iteration on arbitrary F4 reference tensors."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random
import time

from entrywise_reflection_constraints import field_dot, hasse_row
from f4_two_cut_smith_fibers import F4FullAlphabetLattice, CompressedF4Fiber
from finite_precision_smith_frames import build_lattice
from full_alphabet_two_cut_lattices import CompressedFullAlphabetFiber
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('evidence/finite-fiber-basis-input.json'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started, rng = time.monotonic(), random.Random(3917)
    h, zero = [0, 6, 2, 8], [0]*4
    heights = [zero, zero, h, h, h, h, zero, zero]
    source = tensor(4)
    generic = F4FullAlphabetLattice(source['values'], 4, 24, heights, 4)
    binary = build_lattice(source['values'], 4, 24, heights, 4)
    assert generic.orders == binary.orders
    equations = 0
    for first, second in zip(generic.equations, binary.equations):
        assert first == tuple((power, (row, 0), rhs) for power, row, rhs in second)
        equations += len(first)
    floors = tuple(max(0, value-4) for value in generic.orders)
    first = CompressedF4Fiber(generic, floors)
    second = CompressedFullAlphabetFiber(binary, floors)
    assert first.variables == second.variables
    controls = []
    for targets in (floors, generic.orders):
        pivots = first.solve(targets)
        assert pivots is not None
        values = first.sample(pivots, rng)
        profiles = generic.check_witness(values, targets)
        controls.append({'reference': 'binary reflection', 'tensor_sha256': digest(values),
                         'targets': targets, 'profiles': profiles})
    original = args.input.read_bytes()
    for case in json.loads(original)['cases']:
        seed = case['seed']
        assert digest(seed['values']) == case['tensor_sha256']
        generic = F4FullAlphabetLattice(seed['values'], 4, seed['degree'], heights, 4)
        floors = tuple(max(0, value-4) for value in generic.orders)
        fiber = CompressedF4Fiber(generic, floors)
        for targets in (generic.orders, floors):
            pivots = fiber.solve(targets)
            assert pivots is not None
            witness = fiber.sample(pivots, rng)
            requests = [(rng.randrange(256), rng.choice((0, 2, 3)), rng.randrange(3)) for _ in range(24)]
            equations_at_points = []
            for entry, point, order in requests:
                value = field_dot(hasse_row(point, order, seed['degree']), witness[entry])
                equations_at_points.append(fiber.jet_equation(entry, point, order, value))
            refined = fiber.solve(targets, equations_at_points)
            assert refined is not None
            values = fiber.sample(refined, rng)
            for entry, point, order in requests:
                row = hasse_row(point, order, seed['degree'])
                assert field_dot(row, values[entry]) == field_dot(row, witness[entry])
            profiles = generic.check_witness(values, targets)
            controls.append({'reference': case['label'], 'tensor_sha256': digest(values),
                             'fiber_variables': fiber.variables, 'finite_jet_equations': len(requests),
                             'targets': targets, 'profiles': profiles})
    names = (Path(__file__).name, 'f4_two_cut_smith_fibers.py', 'f4_regular_smith_frames.py',
             'finite_jet_smith_fibers.py', 'finite_precision_smith_frames.py',
             'full_alphabet_two_cut_lattices.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp')
    report = {'scope': __doc__, 'passed': True, 'input_sha256': sha256(original).hexdigest(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'binary_equations_agree': equations, 'binary_common_fiber_dimension': first.variables,
              'full_tensor_frame_replays': len(controls), 'endpoint_entries': 256*len(controls),
              'finite_jet_equalities': sum(row.get('finite_jet_equations', 0) for row in controls),
              'complete_native_profiles': 2*len(controls), 'controls': controls,
              'elapsed_seconds': time.monotonic()-started}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
