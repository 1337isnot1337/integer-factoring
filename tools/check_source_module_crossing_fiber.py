#!/usr/bin/env python3
"""Independent complete-equation and literal-image tests for module shift reuse."""

from hashlib import sha256
import json
from pathlib import Path
import random

from f4_recentered_lattices import RecenteredCrossingLattice, dot, multiply_matrix
from restricted_f4_crossing_fiber import RestrictedCrossingFiber
from search_f4_recentered_crossings import five_letter_source, increment
from source_module_crossing_fiber import SourceModuleCrossingFiber


def small_source(q, degree):
    base = []
    for index in range(q ** 4):
        x, y, u, v = ((index // q ** leg) % q for leg in range(4))
        r, e = int(x == y and u == v), int(x == v and y == u)
        base.append((r | ((r ^ e) << 1), 0))
    primitives = []
    for i in range(3):
        row = [0] * q ** 3
        row[i] = 3
        row[i + q] ^= 3 << i
        primitives.append({'normalized_slab': row, 'normalized_degree': i + 1})
    descriptors = [(i, first, shift) for i, primitive in enumerate(primitives)
                   for first in range(q) for shift in range(degree - primitive['normalized_degree'] + 1)]
    source = {'dimensions': [q] * 4, 'degree': degree,
              'primitives': primitives, 'direction_descriptors': descriptors}
    directions = []
    for i, first, shift in descriptors:
        tensor = [(0, 0)] * q ** 4
        for j, value in enumerate(primitives[i]['normalized_slab']):
            tensor[first + q * j] = value << shift, 0
        directions.append(tensor)
    return base, source, directions


def literal_residuals(fiber, low, high):
    values = fiber.seed(low, high)
    checked = 0
    for side in range(2):
        lattice = fiber.full_lattice
        matrix = lattice.local_matrix(values, side)
        left, right = lattice.frames[side]
        precision = lattice.precisions[side]
        actual = multiply_matrix(multiply_matrix(left, matrix, precision), right, precision)
        for a in range(fiber.size):
            equations = fiber.equations[side * fiber.size + a]
            complete = {}
            for power in range(fiber.caps[side * fiber.size + a]):
                for b in range(fiber.size):
                    value = (actual[a][b] >> (2 * power)) & 3
                    complete[power, b] = value
            positions = fiber.equation_positions[side * fiber.size + a]
            assert len(positions) == len(equations)
            for position, (_, row_low, row_high, rhs) in zip(positions, equations):
                assert complete.pop(position) == dot(row_low, row_high, low, high) ^ rhs
            assert all(value == 0 for value in complete.values())
            checked += fiber.caps[side * fiber.size + a] * fiber.size
    assert fiber.contains(values)
    return checked


def verify():
    rng = random.Random(90761)
    comparisons = target_controls = literal_entries = 0
    for q, degree in ((2, 4), (3, 5)):
        values, source, directions = small_source(q, degree)
        zero = (0,) * q
        h = tuple(2 * i - 1 for i in range(q))
        heights = (zero, h, h, zero, h, zero, zero, h)
        old = RestrictedCrossingFiber(RecenteredCrossingLattice(values, q, degree, heights, 3), directions)
        new = SourceModuleCrossingFiber(values, source, heights, 3)
        assert old.orders == new.orders and old.caps == new.caps
        assert old.equations == new.equations
        comparisons += sum(len(row) for row in old.equations)
        for _ in range(8):
            low, high = rng.getrandbits(new.variables), rng.getrandbits(new.variables)
            assert old.seed(low, high) == new.seed(low, high)
            literal_entries += literal_residuals(new, low, high)
            targets = [rng.randrange(cap + 1) for cap in new.caps]
            assert old.solve(targets) == new.solve(targets)
            target_controls += 1

    path = Path('evidence/aligned-gaussian-module-directions-q5.json')
    source = json.loads(path.read_text())
    baseline = five_letter_source()[0]
    half, zero = (0, 6, 3, 0, 4), (0,) * 5
    heights = (zero, zero, half, half, half, half, zero, zero)
    fiber = SourceModuleCrossingFiber(baseline, source, heights, 5)
    reference = json.loads(Path('evidence/aligned-gaussian-module-lattice-two-zero-d18.json').read_text())
    assert fiber.orders == tuple(reference['reference_targets'])
    assert fiber.caps == tuple(reference['ceiling_targets'])
    pivots = fiber.solve(fiber.orders)
    assert pivots is not None and len(pivots) == fiber.variables == 2120
    raises = [row for row in range(50) if increment(fiber, pivots, row, fiber.orders[row]) is not None]
    assert raises == reference['complete_single_row_raises'] == []
    for _ in range(3):
        low, high = rng.getrandbits(fiber.variables), rng.getrandbits(fiber.variables)
        literal_entries += literal_residuals(fiber, low, high)
    retained = 0
    for witness in reference['retained_witnesses']:
        assert fiber.check_witness(witness['values'], witness['targets']) == tuple(map(tuple, witness['profiles']))
        retained += 1
    return {'status': 'exact primitive-shift crossing fiber certified',
            'source_record_sha256': sha256(path.read_bytes()).hexdigest(),
            'small_complete_projected_equations_compared': comparisons,
            'small_target_systems_compared': target_controls,
            'literal_full_F4_frame_coefficients_checked': literal_entries,
            'aligned_module_variables': fiber.variables,
            'primitive_frame_products': fiber.primitive_frame_images,
            'baseline_preserving_dimension': 0,
            'single_row_raises_feasible': raises,
            'retained_original_fiber_witnesses_replayed': retained,
            'factoring_bound_improved': False}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
