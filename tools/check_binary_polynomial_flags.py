#!/usr/bin/env python3
"""Literal polynomial controls for nonconstant two-cut row flags."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random

from check_binary_two_cut_lattices import pmul, polynomial, literal_profiles, endpoint
from f4_lattice_linear import dot
from search_binary_f4_flags import flags
from search_binary_polynomial_flags import (PolynomialFlagSearch, polynomial_row,
    constant_row, regular_frame, lower_jet_frames, intersect)
from search_binary_two_cut_lattices import flatten, weighted_monomial


def padd(first, second):
    result = dict(first)
    for key, value in second.items():
        result[key] = result.get(key, 0) ^ value
        if not result[key]:
            del result[key]
    return result


def matrix_product(first, second):
    out = [[{} for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for j in range(4):
            for k in range(4):
                out[i][j] = padd(out[i][j], pmul(first[i][k], second[k][j]))
    return out


def binary_polynomial(value):
    return {j: 1 for j in range(value.bit_length()) if (value >> j) & 1}


def literal_transformation(engine, seed, side, frame):
    lattice = engine.lattice
    source = [[{} for _ in range(4)] for _ in range(4)]
    for entry, pair in enumerate(seed['values']):
        i, j = flatten(entry, side)
        for power, value in polynomial(pair).items():
            exponent = weighted_monomial(entry, power, side, lattice.degree,
                                         lattice.heights).bit_length()-1
            source[i][j][exponent] = value
    left, right = lattice.bases[side]
    left = [[binary_polynomial(v) for v in row] for row in left]
    right = [[binary_polynomial(v) for v in row] for row in right]
    frame = [[polynomial(pair) for pair in row] for row in frame]
    return matrix_product(frame, matrix_product(left, matrix_product(source, right)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = random.Random(311636)
    engine = PolynomialFlagSearch(12, 4)
    frame_set = flags()
    coefficient_rows = endpoints = full_profiles = constant_controls = 0
    for trial in range(60):
        base = rng.choice(frame_set)
        degree = 1+trial%4
        frame = tuple(polynomial_row(row, tuple((k, tuple(rng.randrange(4) for _ in range(4)))
                          for k in range(1, degree+1))) for row in base)
        assert regular_frame(frame)
        seed = engine.sample({}, rng)
        augmented = engine.assignment(seed)
        endpoint(seed['values']); endpoints += 16
        for side in range(2):
            literal = literal_transformation(engine, seed, side, frame)
            for i, vector in enumerate(frame):
                expected = Counter((power, value) for entry in literal[i] for power, value in entry.items())
                actual = Counter((power, value) for power, equation in engine.polynomial_equations(side, vector)
                                 if (value := dot(equation, augmented)))
                assert actual == expected
                coefficient_rows += 1
        if trial < 12:
            assert engine.lattice.profiles(seed) == literal_profiles(seed['values'])
            full_profiles += 2
    for base in frame_set[::281]:
        frame = tuple(constant_row(row) for row in base)
        for side in range(2):
            for row, vector in zip(base, frame):
                assert engine.row(side, row) == engine.polynomial_equations(side, vector)
                constant_controls += 1
    identity = tuple(constant_row(tuple(int(i == j) for j in range(4))) for i in range(4))
    targets = engine.lattice.orders
    common = engine.common_system(targets)
    first = engine.polynomial_system(0, identity, targets[:4], common)
    second = engine.polynomial_system(1, identity, targets[4:], common)
    pivots = intersect(first, second, engine.lattice.variables)
    assert pivots is not None and engine.lattice.variables-len(pivots) == 1
    free = next(i for i in range(engine.lattice.variables) if i not in pivots)
    all_line_values = []
    # Force each of the four free F4 values instead of relying on a sampler.
    from f4_lattice_linear import encode, dot as fdot
    for color in range(4):
        assignment = [0]*engine.lattice.variables; assignment[free] = color
        planes = encode(assignment)
        for pivot, (a, b, rhs) in sorted(pivots.items()):
            value = fdot((a, b), planes) ^ rhs
            planes = planes[0] ^ ((value & 1) << pivot), planes[1] ^ ((value >> 1) << pivot)
        seed = engine.lattice.seed(*planes)
        profiles = engine.polynomial_verify(seed, (identity, identity), targets)
        assert profiles == literal_profiles(seed['values'])
        endpoint(seed['values']);endpoints += 16;full_profiles += 2
        all_line_values.append(seed['values'])
    # A jet at the max-target minus min-target threshold cannot change a
    # divisibility constraint on a source already satisfying the entry floor.
    tail_jet_controls = 0
    seed = engine.lattice.seed(0)
    for side in range(2):
        local = targets[4*side:4*side+4]
        depth = max(local)-min(local)
        for _ in range(10):
            frame = tuple(polynomial_row(tuple(int(i == j) for j in range(4)),
                         ((depth, tuple(rng.randrange(4) for _ in range(4))),)) for i in range(4))
            baseline = literal_transformation(engine, seed, side, identity)
            changed = literal_transformation(engine, seed, side, frame)
            for i in range(4):
                for j in range(4):
                    assert {k:v for k,v in baseline[i][j].items() if k<local[i]} == {
                           k:v for k,v in changed[i][j].items() if k<local[i]}
            tail_jet_controls += 1
    singular = tuple(constant_row((0, 0, 0, 0)) for _ in range(4))
    try:
        engine.polynomial_system(0, singular, targets[:4])
    except ValueError:
        pass
    else:
        raise AssertionError('Singular constant frame was accepted.')
    manifest_path = Path('evidence/binary-polynomial-flags-first-jet.json')
    manifest = json.loads(manifest_path.read_text())
    assert manifest['completed'] and len(manifest['variants']) == 59
    assert manifest['frames_per_cut'] == 4096
    assert len(set(lower_jet_frames())) == 4096
    feasible = [v for v in manifest['variants'] if v['feasible_pairs']]
    assert len(feasible) == 1 and feasible[0]['target'] == list(targets)
    assert feasible[0]['feasible_pairs'] == 1
    for witness in manifest['witnesses']:
        frames = tuple(tuple(tuple(tuple(v) for v in row) for row in frame) for frame in witness['frames'])
        seed = engine.lattice.seed(0); seed['values'] = witness['values']
        assert engine.polynomial_verify(seed, frames, witness['target']) == literal_profiles(seed['values'])
        assert tuple(map(tuple, witness['values'])) in [tuple(map(tuple, values)) for values in all_line_values]
        full_profiles += 2
    result = {'status': 'Literal polynomial, endpoint, complete minor and finite-manifest controls passed.',
              'literal_transformed_rows': coefficient_rows,
              'constant_frame_equation_controls': constant_controls,
              'literal_endpoint_coefficients': endpoints,
              'all_minor_cut_profiles': full_profiles,
              'full_affine_line_points': 4, 'irrelevant_tail_jet_controls': tail_jet_controls,
              'finite_search': {'targets': 59, 'nonconstant_frames_per_cut': 4096,
                  'feasible_targets': 1, 'feasible_frame_pairs': 1,
                  'scope': 'Complete I+zN catalogue with strictly lower triangular N; other jets and constant flags remain open.'},
              'source_sha256': {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                  for name in ('search_binary_polynomial_flags.py', Path(__file__).name)}}
    if args.output: args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
