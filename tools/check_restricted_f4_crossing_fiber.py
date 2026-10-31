#!/usr/bin/env python3
"""Independent image-rank and literal substitution controls for prescribed fibers."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random

from check_block_smith_crossings import binary_rank
from check_f4_recentered_lattices import random_crossing
from f4_recentered_lattices import RecenteredCrossingLattice, RecenteredFiber, scale_row
from restricted_f4_crossing_fiber import RestrictedCrossingFiber


def run():
    rng = random.Random(2026103198)
    counts = {'independent_binary_image_ranks': 0, 'literal_span_memberships': 0,
              'feasible_witnesses': 0, 'complete_parameter_equation_comparisons': 0,
              'negative_controls': 0}
    cases = []
    for q, degree in ((2, 2), (2, 3), (3, 2)):
        base = random_crossing(q, degree, rng)
        heights = [tuple(rng.randrange(-1, 2) for _ in range(q)) for _ in range(8)]
        lattice = RecenteredCrossingLattice(base, q, degree, heights, raise_cap=3)
        directions = []
        for _ in range(7):
            values = random_crossing(q, degree, rng)
            directions.append(tuple((a ^ c, b ^ d) for (a, b), (c, d) in zip(values, base)))
        directions.append(tuple((0, 0) for _ in base))
        directions.append(tuple((b, a ^ b) for a, b in directions[0]))
        restricted = RestrictedCrossingFiber(lattice, directions)
        binary_rows = []
        for direction in directions:
            a = sum(low << (entry*(degree+1)) for entry, (low, high) in enumerate(direction))
            b = sum(high << (entry*(degree+1)) for entry, (low, high) in enumerate(direction))
            width = q**4*(degree+1)
            binary_rows.extend((a | (b << width), b | ((a ^ b) << width)))
        assert binary_rank(binary_rows) == 2*restricted.variables
        counts['independent_binary_image_ranks'] += 1
        assert restricted.variables <= 7
        floors = tuple(max(0, order-2) for order in lattice.orders)
        compressed = RecenteredFiber(restricted, floors)
        for trial in range(20):
            a, b = rng.getrandbits(restricted.variables), rng.getrandbits(restricted.variables)
            values = restricted.seed(a, b)
            literal = [list(pair) for pair in base]
            for i, source in enumerate(restricted.selected_direction_indices):
                coefficient = ((a >> i) & 1) | (((b >> i) & 1) << 1)
                for j, pair in enumerate(directions[source]):
                    x, y = scale_row(*pair, coefficient)
                    literal[j][0] ^= x
                    literal[j][1] ^= y
            assert tuple(map(tuple, literal)) == values
            assert restricted.contains(values)
            counts['literal_span_memberships'] += 1
            targets = lattice.orders if trial == 0 else tuple(
                rng.randrange(floor, min(cap, order+1)+1)
                for floor, cap, order in zip(floors, lattice.caps, lattice.orders))
            direct, small = restricted.solve(targets), compressed.solve(targets)
            assert (direct is None) == (small is None)
            if direct is not None:
                assert restricted.variables-len(direct) == compressed.variables-len(small)
                restricted.check_witness(compressed.sample(small, rng), targets)
                counts['feasible_witnesses'] += 1
        outside = list(base)
        outside[0] = (outside[0][0] ^ 1, outside[0][1])
        assert not restricted.contains(outside)
        counts['negative_controls'] += 1
        if q == 2:
            whole = []
            for entry in range(q**4):
                for power in range(degree):
                    direction = [(0, 0)]*(q**4)
                    direction[entry] = ((1 << power) ^ (1 << (power+1)), 0)
                    whole.append(tuple(direction))
            full = RestrictedCrossingFiber(lattice, whole)
            assert full.variables == lattice.variables
            assert full.equations == lattice.equations
            counts['complete_parameter_equation_comparisons'] += 1
        cases.append({'alphabet': q, 'degree': degree, 'input_directions': len(directions),
                      'image_rank': restricted.variables, 'common_floor_rank': compressed.variables})
    return {'scope': 'Exact finite image ranks, complete projected equations, and literal prescribed-family substitutions.',
            'checks': counts, 'cases': cases,
            'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                              for name in ('restricted_f4_crossing_fiber.py', 'f4_recentered_lattices.py', Path(__file__).name)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result['checks']))


if __name__ == '__main__':
    main()
