#!/usr/bin/env python3
"""Apply constant physical core maps to cached Laurent coordinate matrices.

This module performs exact algebra on the supplied coefficient windows. It
does not certify the endpoint of arbitrary supplied maps: callers must use
a valid flag family or an exact C F W target identity before promising rank.
"""
from __future__ import annotations

from collections import Counter
import time

from check_rectangular_f4_repair import binary_mul, binary_power, binary_divide, scale, times_binary, multiply
from check_rectangular_identity_transports import checked_matrix
from check_rectangular_defect_support import pgcd


CORE = (1, 0, 2, 3)


def core_height(index):
    return (0, 0, 1, 1)[index % 4]+(0, 0, 1)[index // 4]


def physical(place, M, letter, state):
    mask = CORE[state]
    low, top = mask & 1, (mask >> 1)&1
    if place == 'infinity':
        if letter == 'x':
            return M*(low+top)
        if letter == 'y':
            return -M*(low+top)
        return (1-M)*low+(M-1)*top
    if letter == 'x':
        return (1-2*M)*low-top
    if letter == 'y':
        return low+(2*M-1)*top
    return 0


def normalize(value, shift):
    union = value[0] | value[1]
    if not union:
        return (0, 0), 0
    order = (union & -union).bit_length()-1
    return (value[0] >> order, value[1] >> order), shift+order


def add_entry(left, left_shift, right, right_shift):
    if not (right[0] | right[1]):
        return left, left_shift
    if not (left[0] | left[1]):
        return right, right_shift
    shift = min(left_shift, right_shift)
    value = tuple((left[p] << (left_shift-shift)) ^ (right[p] << (right_shift-shift)) for p in (0, 1))
    return normalize(value, shift)


def coordinate_lookup(states):
    seen, lookup, positions = Counter(), {}, []
    for index, state in enumerate(states):
        position = seen[state]
        seen[state] += 1
        lookup[state, position] = index
        positions.append(position)
    if len(set(seen.values())) != 1:
        raise ValueError('Core states must have a common vacuum multiplicity.')
    return lookup, positions


def replacements(states, letters, lookup, positions, target, source):
    mask = target_bits = source_bits = 0
    for letter, state in target.items():
        if letter in letters:
            offset = 2*letters.index(letter)
            mask |= 3 << offset
            target_bits |= CORE[state] << offset
            source_bits |= CORE[source[letter]] << offset
    return [(i, lookup[((state & ~mask)|source_bits, positions[i])])
            for i, state in enumerate(states) if state & mask == target_bits]


def trim(matrix, shifts, cutoff):
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            if not (value[0] | value[1]):
                continue
            length = cutoff-shifts[i][j]
            if length <= 0:
                matrix[i][j], shifts[i][j] = (0, 0), 0
            elif max(value).bit_length() > length:
                mask = (1 << length)-1
                matrix[i][j], shifts[i][j] = normalize((value[0]&mask, value[1]&mask), shifts[i][j])


def transform_coordinate_matrix(compiled, data, V=None, W=None):
    """Return the cached matrix after output V and input W core dressing.

    V is indexed by target/source (v,y), W by source/target (x,u).
    Vacuum multiplicity coordinates are preserved. For inexact input, the
    returned known_below is a uniform absolute cutoff, including zero cells.
    """
    started = time.monotonic()
    source_matrix, source_shifts, original_metadata = compiled
    metadata = dict(original_metadata)
    if metadata['place'] != data['place'] or metadata['M'] != data['M']:
        raise ValueError('The compiled matrix and support data do not match.')
    n, place, M = len(source_matrix), data['place'], data['M']
    matrix = [[tuple(value) for value in row] for row in source_matrix]
    shifts = [row[:] for row in source_shifts]
    exact = bool(metadata['exact_polynomial_matrix'])
    nonzero_shifts = [shifts[i][j] for i, row in enumerate(matrix) for j, value in enumerate(row) if value[0] | value[1]]
    input_window = int(metadata.get('profile_precision', 512))
    cutoff = int(metadata['known_below'])
    if not exact and nonzero_shifts:
        cutoff = min(cutoff, min(nonzero_shifts)+input_window)
    initial_cutoff = cutoff
    if not exact:
        trim(matrix, shifts, cutoff)
    row_states, column_states = metadata['row_core_states'], metadata['column_core_states']
    row_letters, column_letters = metadata['row_core_letters'], metadata['column_core_letters']
    row_lookup, row_positions = coordinate_lookup(row_states)
    column_lookup, column_positions = coordinate_lookup(column_states)
    all_letters = set(row_letters) | set(column_letters)
    denominator = int(metadata['unit_denominator_hex'], 16)
    steps = []
    for side, supplied, letters in (('output', V, ('v', 'y')), ('input', W, ('x', 'u'))):
        if supplied is None:
            continue
        supplied = checked_matrix(supplied)
        terms = []
        for r, line in enumerate(supplied):
            for s, value in enumerate(line):
                coefficient = value ^ int(r == s)
                if not coefficient:
                    continue
                target_index, source_index = (r, s) if side == 'output' else (s, r)
                target = {letters[0]: target_index % 4, letters[1]: target_index // 4}
                source = {letters[0]: source_index % 4, letters[1]: source_index // 4}
                dh = (core_height(source_index)-core_height(target_index) if side == 'output'
                      else core_height(target_index)-core_height(source_index))
                dz = sum(physical(place, M, letter, target[letter])-physical(place, M, letter, source[letter]) for letter in letters)
                if place == 'infinity':
                    dz -= dh
                terms.append((coefficient, dh, dz, target, source))
        if not terms:
            continue
        if not set(letters) <= all_letters:
            raise ValueError(f'The supplied coordinate bases do not contain all {side} core letters.')
        padding = max(0, -min(term[1] for term in terms))
        old_matrix, old_shifts = matrix, shifts
        unit = binary_power(3, padding)
        matrix = [[times_binary(value, unit) if padding else value for value in row] for row in old_matrix]
        shifts = [row[:] for row in old_shifts]
        for coefficient, dh, dz, target, source in terms:
            factor = binary_power(3, padding+dh)
            rows = replacements(row_states, row_letters, row_lookup, row_positions, target, source)
            columns = replacements(column_states, column_letters, column_lookup, column_positions, target, source)
            for i, source_i in rows:
                for j, source_j in columns:
                    value = old_matrix[source_i][source_j]
                    if not (value[0] | value[1]):
                        continue
                    value = scale(times_binary(value, factor), coefficient)
                    matrix[i][j], shifts[i][j] = add_entry(matrix[i][j], shifts[i][j], value, old_shifts[source_i][source_j]+dz)
        shift_floor = min(0, *(term[2] for term in terms))
        cutoff += shift_floor
        denominator = binary_mul(denominator, unit)
        if not exact:
            trim(matrix, shifts, cutoff)
        steps.append({'side': side, 'nonidentity_entries': len(terms), 'unit_denominator_padding': padding,
                      'minimum_laurent_shift': shift_floor, 'known_below_after_stage': cutoff})
    maximum_degree = max((max(value).bit_length()-1 for row in matrix for value in row), default=-1)
    precision = max(128, input_window, 1 << max(0, maximum_degree).bit_length())
    if precision > 8192:
        raise ArithmeticError('The transformed exact coefficient window exceeds the supported8192-term profile limit.')
    minimum = min((shifts[i][j] for i, row in enumerate(matrix) for j, value in enumerate(row) if value[0] | value[1]), default=0)
    metadata.update({'unit_denominator_hex': hex(denominator), 'unit_denominator_degree': denominator.bit_length()-1,
                     'profile_precision': precision, 'minimum_valuation': minimum,
                     'maximum_normalized_degree': maximum_degree,
                     'exact_polynomial_matrix': exact, 'known_below': cutoff,
                     'coordinate_transform_seconds': time.monotonic()-started,
                     'coordinate_transform_precision': {'input_exact': exact, 'input_uniform_known_below': initial_cutoff,
                                                        'steps': steps, 'output_uniform_known_below': cutoff},
                     'endpoint_validation_required': True})
    return matrix, shifts, metadata


def linear_combination_coordinate_matrices(terms):
    """Combine (compiled, [(F4 coefficient, local-z exponent), ...]) terms.

    Scalar exponents refer to z at the selected place, so a w-power has
    opposite sign at infinity. All matrices must use the identical physical
    core coordinate bases. Unit denominators and unknown tails are retained.
    """
    started = time.monotonic()
    prepared = []
    for compiled, scalar_terms in terms:
        coefficients = {}
        for coefficient, exponent in scalar_terms:
            if coefficient not in (0, 1, 2, 3) or not isinstance(exponent, int):
                raise ValueError('A scalar term needs an F4 coefficient and an integer local exponent.')
            coefficients[exponent] = coefficients.get(exponent, 0) ^ coefficient
        coefficients = {exponent: coefficient for exponent, coefficient in coefficients.items() if coefficient}
        if not coefficients:
            continue
        base = min(coefficients)
        polynomial = tuple(sum(((coefficient >> p)&1) << (exponent-base)
                               for exponent, coefficient in coefficients.items()) for p in (0, 1))
        prepared.append((compiled, polynomial, base))
    if not prepared:
        raise ValueError('At least one nonzero scalar term is required.')
    first = prepared[0][0]
    metadata = dict(first[2])
    n = len(first[0])
    keys = ('M', 'place', 'full_core', 'row_core_states', 'column_core_states',
            'row_core_letters', 'column_core_letters', 'row_charges', 'column_charges')
    denominator, exact, cutoffs, precisions = 1, True, [], []
    for (matrix, shifts, meta), polynomial, base in prepared:
        if len(matrix) != n or any(meta.get(key) != metadata.get(key) for key in keys):
            raise ValueError('Scalar combinations require identical physical core coordinate bases.')
        current = int(meta['unit_denominator_hex'], 16)
        if not current & 1:
            raise ValueError('Only valuation-unit denominators may be cleared.')
        denominator = binary_mul(denominator, binary_divide(current, pgcd(denominator, current)))
        precision = int(meta.get('profile_precision', 512))
        precisions.append(precision)
        if not meta['exact_polynomial_matrix']:
            exact = False
            nonzero = [shifts[i][j] for i, row in enumerate(matrix) for j, value in enumerate(row) if value[0] | value[1]]
            known = int(meta['known_below'])
            if nonzero:
                known = min(known, min(nonzero)+precision)
            cutoffs.append(known+base)
    matrix, shifts = [[(0, 0)]*n for _ in range(n)], [[0]*n for _ in range(n)]
    for (source, source_shifts, meta), polynomial, base in prepared:
        unit = binary_divide(denominator, int(meta['unit_denominator_hex'], 16))
        for i, row in enumerate(source):
            for j, value in enumerate(row):
                if not (value[0] | value[1]):
                    continue
                value = multiply(times_binary(value, unit), polynomial)
                matrix[i][j], shifts[i][j] = add_entry(matrix[i][j], shifts[i][j], value, source_shifts[i][j]+base)
    cutoff = min(cutoffs) if cutoffs else min(int(compiled[2]['known_below'])+base for compiled, _, base in prepared)
    if not exact:
        trim(matrix, shifts, cutoff)
    degree = max((max(value).bit_length()-1 for row in matrix for value in row), default=-1)
    precision = max(128, *precisions, 1 << max(0, degree).bit_length())
    if precision > 8192:
        raise ArithmeticError('The scalar combination exceeds the supported8192-term profile limit.')
    metadata.update({'unit_denominator_hex': hex(denominator), 'unit_denominator_degree': denominator.bit_length()-1,
                     'profile_precision': precision, 'maximum_normalized_degree': degree,
                     'minimum_valuation': min((shifts[i][j] for i, row in enumerate(matrix) for j, value in enumerate(row) if value[0] | value[1]), default=0),
                     'exact_polynomial_matrix': exact, 'known_below': cutoff,
                     'scalar_combination_seconds': time.monotonic()-started,
                     'scalar_combination_terms': len(prepared), 'endpoint_validation_required': True})
    return matrix, shifts, metadata
