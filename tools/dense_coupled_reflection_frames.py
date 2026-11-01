#!/usr/bin/env python3
"""Exact rank-one square-zero Laurent frames on full physical alphabets.

The contractions act directly in the original physical coordinates, using
the checked Bell factors only to conjugate the two x/v vectors. Every
four-leg interaction is retained until the complete pole test.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

from physical_split_crossings import F4_PRODUCT
from reflection_dual_rail_gaussian import tensor, physical_bit_energies
from reflection_laurent_physical_frames import (
    add, scale, degree, valuation, bell_normalize, coupled_admissible_poles,
)


def dot(left, right):
    value = 0
    for a, b in zip(left, right):
        value ^= F4_PRODUCT[a][b]
    return value


def matvec(matrix, vector):
    return [dot(row, vector) for row in matrix]


def transpose(matrix):
    return list(map(list, zip(*matrix)))


def projective_vectors(size):
    """Exactly one vector per F4 projective point: first nonzero entry is one."""
    for first in range(size):
        for tail in product(range(4), repeat=size - first - 1):
            yield [0] * first + [1] + list(tail)


def rank_one_generators(size):
    """Every nonzero rank-one square-zero matrix u v^T exactly once."""
    for left in projective_vectors(size):
        pivot = next(i for i, value in enumerate(left) if value)
        free = [i for i in range(size) if i != pivot]
        for coordinates in product(range(4), repeat=size - 1):
            right = [0] * size
            for i, value in zip(free, coordinates):
                right[i] = value
            right[pivot] = dot(left, right)
            if any(right):
                assert dot(left, right) == 0
                yield left, right


def apply_rank_one(values, size, leg, left, right):
    """Contract a leg with right, then expand it with left: O(q^4) work."""
    stride = size ** leg
    contracted = {}
    for index, value in values.items():
        symbol = index // stride % size
        coefficient = right[symbol]
        if coefficient:
            key = index - symbol * stride
            updated = add(contracted.get(key, (0, 0)), scale(value, coefficient))
            if updated != (0, 0):
                contracted[key] = updated
            elif key in contracted:
                del contracted[key]
    output = {}
    for symbol, coefficient in enumerate(left):
        if coefficient:
            for key, value in contracted.items():
                output[key + symbol * stride] = scale(value, coefficient)
    return output


def original_vectors(left, right, factors):
    ux = matvec(factors['first_factor'], left)
    vx = matvec(transpose(factors['inverse_factor']), right)
    return ((ux, vx), (right, left), (left, right), (vx, ux))


def rank_one_terms(seed, left, right, factors):
    """All f^0,...,f^4 coefficients for the restored coupled frame."""
    size = seed['q_x']
    if seed['q_u'] != size or len(left) != size or len(right) != size:
        raise ValueError('Use vectors on the complete common physical alphabet.')
    if not any(left) or not any(right) or dot(left, right):
        raise ValueError('Require nonzero u,v with v^T u=0.')
    vectors = original_vectors(left, right, factors)
    terms = [{i: value for i, value in enumerate(seed['values']) if value != (0, 0)}]
    terms += [{} for _ in range(4)]
    for leg, (u, v) in enumerate(vectors):
        for power in range(leg, -1, -1):
            added = apply_rank_one(terms[power], size, leg, u, v)
            for index, value in added.items():
                updated = add(terms[power + 1].get(index, (0, 0)), value)
                if updated != (0, 0):
                    terms[power + 1][index] = updated
                elif index in terms[power + 1]:
                    del terms[power + 1][index]
    return terms


def transformed_from_terms(seed, terms, pole):
    if pole < 1:
        raise ValueError('Use a positive zero-pole order.')
    shift, mask = 4 * pole, (1 << (4 * pole)) - 1
    values = [(0, 0)] * len(seed['values'])
    for power, entries in enumerate(terms):
        for index, (low, high) in entries.items():
            a, b = values[index]
            offset = shift - power * pole
            values[index] = a ^ (low << offset), b ^ (high << offset)
    if any((a | b) & mask for a, b in values):
        raise ArithmeticError('The complete four-leg Laurent tensor has a pole.')
    values = [(a >> shift, b >> shift) for a, b in values]
    size = seed['q_x']
    for index, (a, b) in enumerate(values):
        x, y, u, v = ((index // size ** leg) % size for leg in range(4))
        assert (a.bit_count() % 2, b.bit_count() % 2) == (int(x == v and y == u), 0)
    actual_degree = max(map(degree, values))
    assert actual_degree == seed['degree']
    return dict(seed, values=values, degree=actual_degree,
                monomial_content=min(valuation(value) for value in values if value != (0, 0)))


def filtration_record(left, right, factors, energy):
    vectors = original_vectors(left, right, factors)
    intervals = []
    for u, v in vectors:
        ui = [energy[i] for i, value in enumerate(u) if value]
        vi = [energy[i] for i, value in enumerate(v) if value]
        intervals.append([min(ui) - max(vi), max(ui) - min(vi)])
    return {'height_difference_ranges': intervals,
            'zero_active_minimum': min(intervals[0][0], intervals[1][0]),
            'infinity_active_minimum': min(intervals[2][0], intervals[3][0]),
            'total_support': sum(sum(bool(value) for value in u + v) for u, v in vectors)}


def digest(values):
    return sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest()


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(data, separators=(',', ':')) + '\n')
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval', type=int, default=4)
    parser.add_argument('--max-pole', type=int, default=8)
    parser.add_argument('--limit', type=int, default=5355)
    parser.add_argument('--seconds', type=float, default=600)
    parser.add_argument('--original-coordinate', action='append', default=[],
                        help='restrict to G=Q^-1 E_i,j Q, written i:j')
    parser.add_argument('--output', type=Path, default=Path('/tmp/dense-coupled-reflection-frames.json'))
    args = parser.parse_args()
    seed = tensor(args.interval)
    _, factors = bell_normalize(seed)
    size = seed['q_x']
    weights = physical_bit_energies(args.interval)
    energy = [sum(weight for bit, weight in enumerate(weights) if symbol >> bit & 1)
              for symbol in range(size)]
    if args.original_coordinate:
        generators = []
        for value in args.original_coordinate:
            i, j = map(int, value.split(':'))
            if not 0 <= i < size or not 0 <= j < size or i == j:
                raise ValueError('Use two distinct in-range letters i:j.')
            left = [row[i] for row in factors['inverse_factor']]
            right = list(factors['first_factor'][j])
            generators.append((left, right))
    else:
        generators = rank_one_generators(size)
    report = {'scope': 'exact square-zero rank-one full-letter Laurent frames; no pressure claim',
              'interval': args.interval, 'alphabet': size, 'degree': seed['degree'],
              'max_pole_checked': args.max_pole, 'energy': energy, 'bell_factors': factors,
              'requested_seconds': args.seconds, 'generator_limit': args.limit,
              'source_seed_sha256': digest(seed['values']), 'records': []}
    begun, total, candidates = time.monotonic(), 0, 0
    for left, right in generators:
        if total >= args.limit or time.monotonic() - begun >= args.seconds:
            break
        terms = rank_one_terms(seed, left, right, factors)
        allowed = coupled_admissible_poles(terms, args.max_pole)
        record = {'index': total, 'u': left, 'v': right, 'pole_orders': allowed,
                  'filtration': filtration_record(left, right, factors, energy)}
        if allowed:
            cases = []
            for pole in allowed:
                result = transformed_from_terms(seed, terms, pole)
                cases.append({'pole': pole, 'degree': result['degree'],
                              'monomial_content': result['monomial_content'],
                              'seed_sha256': digest(result['values'])})
            record['candidates'] = cases
            candidates += len(cases)
        report['records'].append(record)
        total += 1
        if total % 50 == 0:
            report.update(generators_checked=total, polynomial_frames=candidates,
                          elapsed_seconds=time.monotonic() - begun)
            save(args.output, report)
            print(json.dumps({'generators': total, 'polynomial_frames': candidates,
                              'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    report.update(generators_checked=total, polynomial_frames=candidates,
                  elapsed_seconds=time.monotonic() - begun,
                  full_q4_rank_one_catalogue=(size == 4 and not args.original_coordinate and total == 5355))
    save(args.output, report)
    print(json.dumps({key: report[key] for key in ('generators_checked', 'polynomial_frames',
                                                  'elapsed_seconds', 'full_q4_rank_one_catalogue')}))


if __name__ == '__main__':
    main()
