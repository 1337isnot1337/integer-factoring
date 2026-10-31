#!/usr/bin/env python3
"""Independent finite-field and literal full-tensor checks for colored row cells."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random

from colored_row_fibers import ColoredRowFiber, insert_equation, sample_assignment
from typed_charge_lattice_synthesis import TypedChargeLattice, natural_heights
from search_full_alphabet_two_cut_lattices import root_for_profiles


def literal_multiply(a, b):
    polynomial = 0
    for i in range(2):
        for j in range(2):
            if ((a >> i) & 1) and ((b >> j) & 1):
                polynomial ^= 1 << (i+j)
    if polynomial & 4:
        polynomial ^= 7
    return polynomial


def literal_satisfies(equations, assignment):
    for coefficients, rhs in equations:
        total = 0
        for a, b in zip(coefficients, assignment):
            total ^= literal_multiply(a, b)
        if total != rhs:
            return False
    return True


def linear_controls():
    rng = random.Random(318081)
    cases = solutions = samples = 0
    for variables in range(1, 5):
        for _ in range(40):
            equations = [(tuple(rng.randrange(4) for _ in range(variables)), rng.randrange(4))
                         for _ in range(rng.randrange(1, 2*variables+2))]
            brute = [assignment for assignment in product(range(4), repeat=variables)
                     if literal_satisfies(equations, assignment)]
            pivots = {}; consistent = True
            for coefficients, rhs in equations:
                low = sum((c & 1) << i for i, c in enumerate(coefficients))
                high = sum((c >> 1) << i for i, c in enumerate(coefficients))
                if not insert_equation(pivots, low, high, rhs):
                    consistent = False; break
            assert consistent == bool(brute)
            if consistent:
                assert len(brute) == 4**(variables-len(pivots))
                for _ in range(16):
                    low, high = sample_assignment(pivots, variables, rng)
                    assignment = tuple(((low >> i) & 1) | (((high >> i) & 1) << 1)
                                       for i in range(variables))
                    assert literal_satisfies(equations, assignment)
                    samples += 1
            cases += 1; solutions += len(brute)
    return {'complete_enumerated_affine_systems': cases,
            'literal_solutions_counted': solutions, 'sampled_assignments_checked': samples}


def tensor_controls():
    rng = random.Random(318082)
    heights, _ = natural_heights(4)
    lattice = TypedChargeLattice(4, heights, charge_weights=[0]*4, raise_cap=8)
    floors = tuple(max(0, value-6) for value in lattice.orders)
    prepared = ColoredRowFiber(lattice, floors, ()).prepared
    cases = []
    for case in range(8):
        operations = []
        for side in range(2):
            for _ in range(case % 4+1):
                i, j = rng.sample(range(lattice.size), 2)
                power = max(0, lattice.orders[side*lattice.size+i]-lattice.orders[side*lattice.size+j])
                operations.append((side, i, j, power+rng.randrange(3), 2+(case % 2)))
        fiber = ColoredRowFiber(lattice, floors, operations, prepared=prepared)
        targets = list(lattice.orders)
        for row in rng.sample(range(32), 16):
            targets[row] = floors[row]
        pivots = fiber.solve(targets)
        assert pivots is not None
        for row in rng.sample(range(32), 32):
            while targets[row] < lattice.caps[row]:
                proposal = fiber.increment(pivots, row, targets[row])
                if proposal is None:
                    break
                pivots = proposal; targets[row] += 1
        witness = fiber.sample(pivots, rng)
        profiles = fiber.check_witness(witness, targets)
        assert lattice.check_witness(witness, floors) == profiles
        assert any(high for low, high in witness)
        cases.append({'operations': operations, 'targets': targets, 'profiles': profiles,
            'solution_dimension': fiber.variables-len(pivots),
            'witness_sha256': sha256(json.dumps(witness, separators=(',', ':')).encode()).hexdigest()})
    return {'common_fiber_variables': len(prepared[1]), 'full_untruncated_tensor_replays': len(cases), 'cases': cases}


def replay_record(path):
    data = json.loads(path.read_text())
    for name, digest in data['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == digest, name
    args = data['arguments']
    values = None if data['reference'] is None else tuple((pair[0], 0) for pair in data['reference']['values'])
    lattice = TypedChargeLattice(4, data['heights_infinity_then_zero'],
        raise_cap=args['upper_window'], charge_weights=data['charge_weights'],
        negative_boundary=args['negative_boundary'], reference_values=values)
    prepared = ColoredRowFiber(lattice, data['row_target_floors'], ()).prepared
    assert lattice.variables == data['variables'] and len(prepared[1]) == data['fiber_variables']
    assert lattice.profiles(lattice.reference_values) == tuple(map(tuple, data['reference_profiles']))
    for witness in data['retained_witnesses']:
        operations = data['frames'][witness['frame']]['operations']
        fiber = ColoredRowFiber(lattice, data['row_target_floors'], operations, prepared=prepared)
        profiles = fiber.check_witness(witness['values'], witness['targets'])
        assert profiles == tuple(map(tuple, witness['profiles']))
        cost = root_for_profiles(profiles, *data['physical_paired_energies'])
        assert abs(cost['rho']-witness['rho']) < 1e-9
        assert abs(cost['critical_pressure']-witness['critical_pressure']) < 1e-8
    return {'path': str(path), 'completed_requested_trials': data['completed_requested_trials'],
        'completed_trials': data['completed_trials'], 'complete_frames': data['completed_frames'],
        'replayed_witnesses': len(data['retained_witnesses']), 'best_recorded_rho': data['best_recorded_rho']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--search-records', type=Path, nargs='*', default=[])
    args = parser.parse_args()
    result = {'status': 'Exact F4 affine systems and literal, untruncated row-frame products verified.',
        'linear': linear_controls(), 'tensor': tensor_controls(),
        'search_replays': [replay_record(path) for path in args.search_records]}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'linear': result['linear'],
        'full_untruncated_tensor_replays': result['tensor']['full_untruncated_tensor_replays'],
        'search_replays': result['search_replays']}))


if __name__ == '__main__':
    main()
