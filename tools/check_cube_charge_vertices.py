#!/usr/bin/env python3
"""Exact cube-height catalogue and independent tight-graph vertex controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from cube_charge_vertices import cube_edges, height_functions, paired_charge_vertices, parity_charge_vertex


def tight_connected(values, edges):
    reached = {0}
    while True:
        old = len(reached)
        for u, v in edges:
            if abs(values[u]-values[v]) == 1 and (u in reached or v in reached):
                reached.update((u, v))
        if len(reached) == old:
            return len(reached) == len(values)


def exhaustive_integer_vertices(dimension):
    """Independent full integer box, including unsaturated/parity-mixed values."""
    edges = cube_edges(dimension)
    ranges = [range(-state.bit_count(), state.bit_count()+1) for state in range(1, 1 << dimension)]
    vertices = set(); feasible = checked = 0
    for tail in product(*ranges):
        values = (0,)+tail; checked += 1
        if any(abs(values[u]-values[v]) > 1 for u, v in edges):
            continue
        feasible += 1
        if tight_connected(values, edges):
            vertices.add(values)
    return vertices, checked, feasible


def additive(values, dimension):
    weights = [values[1 << bit] for bit in range(dimension)]
    return all(value == sum(weights[bit] for bit in range(dimension) if (state >> bit) & 1)
               for state, value in enumerate(values))


def check():
    rows = []; total_edges = total_parity = integer_boxes = 0
    expected = {1: 2, 2: 6, 3: 38, 4: 990}
    for dimension in range(1, 5):
        edges = cube_edges(dimension)
        heights = tuple(height_functions(dimension))
        assert len(heights) == len(set(heights)) == expected[dimension]
        for values in heights:
            assert values[0] == 0
            assert all(abs(values[u]-values[v]) == 1 for u, v in edges)
            assert tight_connected(values, edges)
            assert all((value-state.bit_count()) % 2 == 0 for state, value in enumerate(values))
            total_edges += len(edges)
        brute = None
        if dimension <= 3:
            vertices, checked, feasible = exhaustive_integer_vertices(dimension)
            assert vertices == set(heights)
            integer_boxes += checked
            brute = {'complete_integer_box_points': checked, 'feasible_points': feasible, 'vertices': len(vertices)}
        catalogue = tuple(paired_charge_vertices(dimension))
        assert len(catalogue) == len(set(catalogue)) == 2*len(heights)
        additive_count = sum(additive(values, dimension) for values in heights)
        assert additive_count == 1 << dimension
        parity = []
        for signs in product((-1, 1), repeat=1 << (dimension-1)):
            values = parity_charge_vertex(dimension, signs)
            assert tuple(value//2 for value in values) in set(heights)
            assert max(values)-min(values) <= 4
            pole = max(0, max(values)-min(values)-1)
            assert pole <= 3
            parity.append(values); total_parity += 1
        assert len(parity) == len(set(parity)) == 1 << (1 << (dimension-1))
        rows.append({'cube_dimension': dimension, 'alphabet_size': 1 << dimension,
            'anchored_unit_height_vertices': len(heights), 'paired_charge_vertices': len(catalogue),
            'additive_paired_vertices': 2*additive_count,
            'nonlinear_paired_vertices': 2*(len(heights)-additive_count),
            'independent_integer_box': brute, 'parity_vertices_on_one_species': len(parity),
            'catalogue_sha256': sha256(json.dumps(heights, separators=(',', ':')).encode()).hexdigest(),
            'height_catalogue': heights})
    # Recover the previously certified complete four-symbol list exactly.
    expected_square = {(0, 2, 2, 4), (0, 2, -2, 0), (0, 2, 2, 0)}
    expected_square |= {tuple(-value for value in row) for row in expected_square}
    assert {tuple(2*value for value in row) for row in height_functions(2)} == expected_square
    return {'status': 'Complete cube charge vertices and bounded-pole parity subfamilies; no pressure improvement claim.',
        'all_edge_saturation_checks': total_edges, 'independent_integer_box_points': integer_boxes,
        'literal_parity_family_vertices': total_parity, 'cases': rows,
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in (Path(__file__).name, 'cube_charge_vertices.py')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = check()
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ('status', 'all_edge_saturation_checks',
        'independent_integer_box_points', 'literal_parity_family_vertices')}))
    print(json.dumps([{key: value for key, value in row.items() if key not in ('height_catalogue', 'catalogue_sha256')}
                      for row in result['cases']]))


if __name__ == '__main__':
    main()
