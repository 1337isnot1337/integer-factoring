#!/usr/bin/env python3
"""Compare the cached elimination oracle with exhaustive determinantal ideals."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import random

from f16_weighted_local_oracle import F16WeightedLocalOracle
from fused_weighted_place_profiles import FusedPlaceOracle


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    generator = random.Random(311715)
    records = []
    for m, theta in ((2, 2), (3, 4), (4, 4), (5, 4)):
        old = FusedPlaceOracle(m, theta)
        q = old.alphabet
        nonlinear = next((p for p in old.points if p is not None and len(p) > 2), None)
        points = [(0, 1), None] + ([nonlinear] if nonlinear is not None else [])
        engines = [[F16WeightedLocalOracle(block['matrix'], old.degree) for block in blocks]
                   for blocks in old.blocks]
        comparisons = 0
        for point in points:
            for cut in (0, 1):
                row_legs, column_legs = ((0, 1), (2, 3)) if cut == 0 else ((0, 2), (1, 3))
                for trial in range(3):
                    heights = [Fraction(generator.randrange(-4, 5), 2) for _ in range(4 * q)]
                    scale = Fraction(3, 2)
                    orders = []
                    for block, engine in zip(old.blocks[cut], engines[cut]):
                        labels = block['labels']
                        row_weights = [sum(heights[leg * q + labels[i][0][leg]] for leg in row_legs)
                                       for i in range(engine.size)]
                        column_weights = [sum(heights[leg * q + labels[0][j][leg]] for leg in column_legs)
                                          for j in range(engine.size)]
                        orders.extend(engine.profile(row_weights, column_weights, point, scale)['orders'])
                    assert tuple(sorted(orders)) == old.profile(cut, point, heights, scale)
                    comparisons += 1
        records.append({'m': m, 'theta': theta, 'alphabet': q,
                        'places': [list(point) if point is not None else None for point in points],
                        'full_weighted_profile_comparisons': comparisons})
    directory = Path(__file__).parent
    sources = (Path(__file__).name, 'f16_weighted_local_oracle.py', 'fused_weighted_place_profiles.py',
               'check_fused_closed_profiles.py', 'check_fused_multipoint_profiles.py',
               'fused_six_vertex_crossings.py', 'f16_identity_operator_profiles.py', 'f4_lattice_linear.py')
    result = {'scope': 'Exact cached fraction-free local Smith profiles versus exhaustive nonzero minors, with negative rational physical weights and nonlinear places.',
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in sources},
              'seed': 311715, 'cases': records,
              'total_full_weighted_profile_comparisons': sum(row['full_weighted_profile_comparisons'] for row in records)}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'comparisons': result['total_full_weighted_profile_comparisons']}))


if __name__ == '__main__':
    main()
