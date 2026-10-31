"""Polynomial vectors in the rational span of exact source kernel slabs.

The bounded-space calculation is complete: pivot coordinates of a
polynomial vector are themselves polynomial, and every remaining
coordinate is constrained by exact polynomial division. No search
over rational denominators is truncated.
"""

from pathlib import Path
import argparse
import json

from check_rectangular_defect_support import Rat, ZERO, ONE, pmul, pdivmod, pgcd


def rational_row_basis(slabs):
    indices = [i for i, row in enumerate(slabs) if any(row)]
    work = [[Rat(value) for value in slabs[i]] for i in indices]
    transformation = [[ONE if i == j else ZERO for j in range(len(indices))]
                      for i in range(len(indices))]
    pivots, position = [], 0
    for column in range(len(slabs[0])):
        row = next((i for i in range(position, len(work)) if work[i][column]), None)
        if row is None:
            continue
        work[position], work[row] = work[row], work[position]
        transformation[position], transformation[row] = transformation[row], transformation[position]
        pivot = work[position][column]
        work[position] = [value / pivot for value in work[position]]
        transformation[position] = [value / pivot for value in transformation[position]]
        for i in range(len(work)):
            if i == position or not work[i][column]:
                continue
            multiplier = work[i][column]
            work[i] = [a + multiplier * b for a, b in zip(work[i], work[position])]
            transformation[i] = [a + multiplier * b for a, b in
                                 zip(transformation[i], transformation[position])]
        pivots.append(column)
        position += 1
        if position == len(work):
            break
    return {'source_indices': indices, 'rank': position, 'pivot_columns': pivots,
            'rref': [[value.pair() for value in row] for row in work[:position]],
            'transformation': [[value.pair() for value in row]
                               for row in transformation[:position]]}


def column_denominators(record):
    rows, rank = record['rref'], record['rank']
    result = []
    for column in range(len(rows[0])):
        denominator = 1
        for i in range(rank):
            divisor = rows[i][column][1]
            quotient, remainder = pdivmod(divisor, pgcd(denominator, divisor))
            assert remainder == 0
            denominator = pmul(denominator, quotient)
        numerators = []
        for i in range(rank):
            quotient, remainder = pdivmod(denominator, rows[i][column][1])
            assert remainder == 0
            numerators.append(pmul(rows[i][column][0], quotient))
        result.append((denominator, numerators))
    return result


def insert_binary(pivots, row):
    while row:
        pivot = row.bit_length() - 1
        if pivot not in pivots:
            pivots[pivot] = row
            return True
        row ^= pivots[pivot]
    return False


def bounded_polynomial_space(record, maximum_degree):
    if maximum_degree < 0:
        raise ValueError('Use a nonnegative polynomial degree bound.')
    rank = record['rank']
    width = maximum_degree + 1
    variables = rank * width
    columns = column_denominators(record)
    pivots = {}
    for denominator, numerators in columns:
        equations = {}
        remainder_width = denominator.bit_length() - 1
        for i, numerator in enumerate(numerators):
            for power in range(width):
                quotient, remainder = pdivmod(numerator << power, denominator)
                bad = remainder | ((quotient >> width) << remainder_width)
                variable = i * width + power
                while bad:
                    bit = bad & -bad
                    bad ^= bit
                    coefficient = bit.bit_length() - 1
                    equations[coefficient] = equations.get(coefficient, 0) ^ (1 << variable)
        for equation in equations.values():
            insert_binary(pivots, equation)
    free = [i for i in range(variables) if i not in pivots]
    assignments, slabs = [], []
    for variable in free:
        assignment = 1 << variable
        for pivot, row in sorted(pivots.items()):
            assignment |= ((row & assignment).bit_count() & 1) << pivot
        assignments.append(assignment)
        slab = []
        for denominator, numerators in columns:
            quotient = remainder = 0
            active = assignment
            while active:
                bit = active & -active
                active ^= bit
                i, power = divmod(bit.bit_length() - 1, width)
                q, r = pdivmod(numerators[i] << power, denominator)
                quotient ^= q
                remainder ^= r
            assert not remainder and quotient.bit_length() <= width
            slab.append(quotient)
        slabs.append(slab)
    return {'degree_bound': maximum_degree, 'dimension': len(free),
            'equation_pivots': [[pivot, row] for pivot, row in sorted(pivots.items())],
            'free_coordinates': free, 'assignments': assignments, 'slabs': slabs}


def row_degree(row):
    return max(value.bit_length() - 1 for value in row)


def encode(row, width):
    return sum(value << (width * j) for j, value in enumerate(row))


def minimal_generators(record, maximum_degree=18):
    generators, dimensions = [], []
    for degree in range(maximum_degree + 1):
        space = bounded_polynomial_space(record, degree)
        dimensions.append(space['dimension'])
        span = {}
        for generator in generators:
            for power in range(degree - row_degree(generator) + 1):
                assert insert_binary(span, encode([value << power for value in generator], degree + 1))
        for row in space['slabs']:
            if insert_binary(span, encode(row, degree + 1)):
                generators.append(row)
        assert len(span) == space['dimension']
        if len(generators) == record['rank']:
            break
    leading = [sum(((value >> row_degree(row)) & 1) << j for j, value in enumerate(row))
               for row in generators]
    leading_pivots = {}
    for row in leading:
        assert insert_binary(leading_pivots, row), 'the displayed rows are not degree reduced'
    return {'generators': generators, 'degrees': list(map(row_degree, generators)),
            'bounded_dimensions': dimensions, 'leading_rank': len(leading_pivots)}


def direction_record(source, degree=18):
    rational = rational_row_basis(source['kernel_slabs'])
    minimal = minimal_generators(rational, max(degree, 16))
    if len(minimal['generators']) != rational['rank']:
        raise ArithmeticError('The degree window did not expose a full polynomial module basis.')
    maximum = max(minimal['degrees'])
    saturation_space = bounded_polynomial_space(rational, max(0, maximum - 1))
    predicted = sum(max(0, maximum - value) for value in minimal['degrees'])
    assert saturation_space['dimension'] == predicted
    primitives = []
    for index, (row, source_degree) in enumerate(zip(minimal['generators'], minimal['degrees'])):
        normalized = [value ^ (value << 1) for value in row]
        primitives.append({'source_module_basis': index, 'source_degree': source_degree,
                           'normalized_degree': source_degree + 1,
                           'normalized_slab': normalized})
    descriptors = [(i, logical, shift)
                   for i, item in enumerate(primitives)
                   for logical in range(source['dimensions'][0])
                   for shift in range(max(0, degree - item['normalized_degree'] + 1))]
    return {'constructor': 'source_slab_rational_module_directions',
            'dimensions': source['dimensions'], 'degree': degree,
            'rational_row_certificate': rational, 'minimal_basis': minimal,
            'saturation_certificate': {
                'maximum_generator_degree': maximum,
                'checked_degree_bound': max(0, maximum - 1),
                'bounded_polynomial_dimension': saturation_space['dimension'],
                'shifted_generator_dimension': predicted,
                'equation_rank': rational['rank'] * max(1, maximum) - saturation_space['dimension'],
                'complete_polynomial_module': True},
            'primitives': primitives, 'direction_descriptors': descriptors,
            'direction_count_before_rank_reduction': len(descriptors),
            'scope': 'Exact rational-span saturation; each displayed polynomial direction vanishes at one.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--degree', type=int, default=18)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = direction_record(json.loads(args.source.read_text()), args.degree)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'rank': result['rational_row_certificate']['rank'],
                      'minimal_degrees': result['minimal_basis']['degrees'],
                      'directions': result['direction_count_before_rank_reduction']}))


if __name__ == '__main__':
    main()
