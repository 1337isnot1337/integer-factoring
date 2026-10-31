#!/usr/bin/env python3
"""Independent finite-field controls for packed F4 lattice equations."""
import argparse
from itertools import product
import json
from pathlib import Path
import random

from f4_lattice_linear import (
    MUL, dot, eliminate, encode, endpoint_remainders, f16_multiply,
    f16_scale_f4, full_support_lift, verify_assignment,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = random.Random(0xF4162026)
    for a, b, c in product(range(16), repeat=3):
        assert f16_multiply(a, b ^ c) == f16_multiply(a, b) ^ f16_multiply(a, c)
        assert f16_multiply(f16_multiply(a, b), c) == f16_multiply(a, f16_multiply(b, c))
    assert all(any(f16_multiply(a, b) == 1 for b in range(1, 16)) for a in range(1, 16))
    assert all(f16_multiply(a, b) == MUL[a][b] for a, b in product(range(4), repeat=2))
    systems = assignments = f4_feasible = f16_feasible = 0
    for size in (1, 4, 8):
        for variables in range(1, 6):
            for sample in range(24):
                raw = [tuple(rng.randrange(4) for _ in range(variables))
                       for _ in range(rng.randrange(variables + 2))]
                equations = [encode(row) for row in raw]
                masks = [rng.randrange(1, 1 << variables) for _ in range(size)]
                pivots = eliminate(equations)
                remainders = endpoint_remainders(masks, pivots)
                possible = False
                support = [False] * size
                nullity_count = 0
                for values in product(range(4), repeat=variables):
                    assignments += 1
                    scalar = encode(values)
                    independent = []
                    for row in raw:
                        total = 0
                        for a, b in zip(row, values):
                            total ^= MUL[a][b]
                        independent.append(total)
                    assert [dot(row, scalar) for row in equations] == independent
                    if any(independent):
                        continue
                    nullity_count += 1
                    endpoints = [dot((mask, 0), scalar) for mask in masks]
                    assert endpoints == [dot(row, scalar) for row in remainders]
                    possible |= all(endpoints)
                    support = [old or bool(value) for old, value in zip(support, endpoints)]
                assert nullity_count == 4 ** (variables - len(pivots))
                for field_size in (4, 16):
                    result = full_support_lift(remainders, pivots, field_size)
                    assert bool(result) == (possible if field_size == 4 else all(support))
                    if result:
                        endpoints, assignment, _ = result
                        assert verify_assignment(equations, assignment)
                        actual = []
                        for mask in masks:
                            value = 0
                            for column, coefficient in assignment.items():
                                if (mask >> column) & 1:
                                    value ^= coefficient
                            actual.append(value)
                        assert actual == list(endpoints) and all(actual)
                        assert all(0 < value < field_size for value in assignment.values())
                        if field_size == 4:
                            f4_feasible += 1
                        else:
                            f16_feasible += 1
                systems += 1
    # Seven projective F2 forms obstruct F4 full support but permit F16.
    masks = [1, 2, 3, 4, 5, 6, 7, 1]
    remainders = endpoint_remainders(masks, {})
    assert full_support_lift(remainders, {}, 4) is None
    assert full_support_lift(remainders, {}, 16)
    # At q coordinates, a particular greedy support chain can get stuck
    # even though a full-support vector exists.  Use the complete fallback.
    boundary = tuple(encode(row) for row in ((1, 0, 0, 1), (0, 1, 0, 2),
                                            (0, 0, 1, 3), (0, 0, 0, 1)))
    assert full_support_lift(boundary, {}, 4)
    result = {'status': 'All exact controls passed.', 'random_systems': systems,
              'complete_F4_assignments': assignments, 'F4_feasible_systems': f4_feasible,
              'F16_feasible_systems': f16_feasible,
              'F16_associativity_and_distributivity_triples': 16 ** 3,
              'four_coordinate_greedy_boundary_regression': True,
              'seven_form_F4_obstruction_and_F16_lift': True}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
