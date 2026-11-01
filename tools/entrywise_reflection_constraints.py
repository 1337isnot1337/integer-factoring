#!/usr/bin/env python3
"""Direct exact polynomial-entry constraints followed by paired trace coupling.

The first stage solves each entry's finite Hasse and coefficient equations.
Only entries with x=v or y=u occur in the second stage, whose equations are
the two complete polynomial partial traces.  No sampled evaluation basis
or physical-orbit saturation assumption is used to certify an output.
"""
from __future__ import annotations

from math import comb

from physical_split_crossings import F4_PRODUCT
from reflection_laurent_physical_frames import add, degree, scale
from solve_companion_affine_degree import Span, coefficient, scaled, witness_entries


def field_dot(left, right):
    a, b = left
    c, d = right
    low = ((a & c).bit_count() ^ (b & d).bit_count()) & 1
    high = ((a & d).bit_count() ^ (b & c).bit_count() ^ (b & d).bit_count()) & 1
    return low | (high << 1)


def local_affine(equations, width):
    """Return one coefficient vector and a complete local kernel basis."""
    rows = {}
    for vector, rhs in equations:
        a, b = vector
        while a or b:
            pivot = (a | b).bit_length() - 1
            value = coefficient((a, b), pivot)
            if pivot not in rows:
                inverse = (0, 1, 3, 2)[value]
                rows[pivot] = scaled((a, b), inverse), F4_PRODUCT[inverse][rhs]
                break
            stored, answer = rows[pivot]
            first, second = scaled(stored, value)
            a, b = a ^ first, b ^ second
            rhs ^= F4_PRODUCT[value][answer]
        else:
            if rhs:
                return None

    def solution(free=None):
        result = (0, 0) if free is None else (1 << free, 0)
        for pivot, (vector, rhs) in sorted(rows.items()):
            value = field_dot(vector, result) ^ (rhs if free is None else 0)
            result = (result[0] ^ ((value & 1) << pivot),
                      result[1] ^ ((value >> 1) << pivot))
        return result

    particular = solution()
    kernel = [solution(free) for free in range(width) if free not in rows]
    assert all(field_dot(vector, particular) == rhs for vector, rhs in equations)
    assert all(field_dot(vector, direction) == 0 for direction in kernel for vector, _ in equations)
    return particular, kernel, len(rows)


def hasse_row(point, order, cap):
    powers = [1]
    for _ in range(cap):
        powers.append(F4_PRODUCT[powers[-1]][point])
    a = b = 0
    for exponent in range(order, cap + 1):
        if exponent & order == order:
            value = powers[exponent - order]
            a |= (value & 1) << exponent
            b |= (value >> 1) << exponent
    return a, b


def trace_column(entry, polynomial, width):
    x, y, u, v = ((entry // 4 ** leg) % 4 for leg in range(4))
    result = (0, 0)
    if x == v:
        result = tuple(value << (width * (4 * y + u)) for value in polynomial)
    if y == u:
        second = tuple(value << (width * (16 + 4 * x + v)) for value in polynomial)
        result = add(result, second)
    return result


def trace_target(power, width):
    a = 0
    for pair in range(2):
        for letter in range(4):
            a |= (1 | (1 << power)) << (width * (16 * pair + 5 * letter))
    return a, 0


def entries_for_rectangle(condition):
    point, cut, first, second, order = condition
    selected_rows, selected_columns = {5 * i for i in first}, {5 * i for i in second}
    result = []
    for entry in range(256):
        x, y, u, v = ((entry // 4 ** leg) % 4 for leg in range(4))
        row, column = (4 * x + y, 4 * u + v) if cut == 0 else (4 * x + u, 4 * y + v)
        if row not in selected_rows and column not in selected_columns:
            result.append(entry)
    return result


def solve_target(base, cap, conditions, zero_cutoff=0, infinity_cutoff=0,
                 trace_power=16, protection_mode='exact', prescribed_jets=()):
    if base['q_x'] != 4 or base['q_u'] != 4 or len(base['values']) != 256:
        raise ValueError('This direct coordinate solver currently uses four letters.')
    if cap < 2 or not 0 < trace_power <= cap:
        return {'feasible': False, 'reason': 'The requested trace exceeds the degree cap.'}
    if protection_mode not in ('exact', 'zero-pattern'):
        raise ValueError('Choose exact or zero-pattern coefficient protection.')
    width, source_degree = cap + 1, base['degree']
    h = (0, 6, 2, 8)
    local_conditions = [[] for _ in range(256)]
    finite_equations = 0
    for condition in conditions:
        point, cut, first, second, order = condition
        if point not in (0, 2, 3) or cut not in (0, 1) or order < 1:
            raise ValueError('Use a non-endpoint F4 place and a positive Hasse order.')
        for entry in entries_for_rectangle(condition):
            local_conditions[entry].extend((hasse_row(point, jet, cap), 0) for jet in range(order))
            finite_equations += order
    for prescribed in prescribed_jets:
        point = prescribed['point']
        if point not in (0, 2, 3):
            raise ValueError('Use a non-endpoint F4 place for prescribed jets.')
        for jet, values in enumerate(prescribed['values']):
            if len(values) != 256 or any(value not in range(4) for value in values):
                raise ValueError('Every prescribed coefficient tensor has 256 F4 entries.')
            row = hasse_row(point, jet, cap)
            for entry, value in enumerate(values):
                local_conditions[entry].append((row, value))
                finite_equations += 1
    particular, columns, local_nullity, protected = [], [], 0, 0
    endpoint_row = ((1 << width) - 1, 0)
    derivative_row = (sum(1 << exponent for exponent in range(1, width, 2)), 0)
    for entry, source in enumerate(base['values']):
        x, y, u, v = ((entry // 4 ** leg) % 4 for leg in range(4))
        equations = [(endpoint_row, int(x == v and y == u)), (derivative_row, 0)] + local_conditions[entry]
        for exponent in range(max(source_degree, cap) + 1):
            low = exponent + h[x] + h[y] < zero_cutoff
            high = source_degree - exponent + h[u] + h[v] < infinity_cutoff
            if not (low or high):
                continue
            rhs = coefficient(source, exponent)
            if protection_mode == 'zero-pattern' and rhs:
                continue
            vector = (1 << exponent, 0) if exponent <= cap else (0, 0)
            equations.append((vector, rhs))
            protected += 1
        solved = local_affine(equations, width)
        if solved is None:
            return {'feasible': False, 'reason': 'An entrywise polynomial system is inconsistent.',
                    'first_inconsistent_entry': entry, 'protected_equations': protected,
                    'finite_equations': finite_equations}
        value, kernel, _ = solved
        particular.append(value)
        local_nullity += len(kernel)
        if x == v or y == u:
            columns.extend((entry, vector) for vector in kernel)
    target = trace_target(trace_power, width)
    for entry, polynomial in enumerate(particular):
        target = add(target, trace_column(entry, polynomial, width))
    system = Span(witnesses=True)
    for index, (entry, polynomial) in enumerate(columns):
        system.insert(trace_column(entry, polynomial, width), index)
    possible, witness, residual = system.solve(target)
    result = {'feasible': possible, 'protected_equations': protected,
              'finite_equations': finite_equations, 'local_free_parameters': local_nullity,
              'trace_coupled_parameters': len(columns), 'trace_constraint_rank': len(system.rows),
              'affine_solution_dimension': local_nullity - len(system.rows) if possible else None}
    if not possible:
        result.update(reason='The remaining paired trace equations are inconsistent.',
                      trace_residual_bitplanes=[hex(value) for value in residual])
        return result
    values = list(particular)
    for index, scalar in witness_entries(witness):
        entry, polynomial = columns[index]
        values[entry] = add(values[entry], scale(polynomial, scalar))
    seed = dict(base, values=values, degree=max(map(degree, values)),
                requested_degree_bound=cap, construction='entrywise Hasse constraints and paired traces')
    seed.pop('scalar_w_power', None)
    verify_target(base, seed, cap, conditions, zero_cutoff, infinity_cutoff,
                  trace_power, protection_mode, prescribed_jets)
    result['seed'] = seed
    return result


def verify_target(base, seed, cap, conditions, zero_cutoff, infinity_cutoff,
                  trace_power, protection_mode, prescribed_jets=()):
    """Replay output identities using direct coefficient and binomial sums."""
    assert max(map(degree, seed['values'])) <= cap
    h = (0, 6, 2, 8)
    total_trace = (0, 0)
    for entry, value in enumerate(seed['values']):
        x, y, u, v = ((entry // 4 ** leg) % 4 for leg in range(4))
        assert tuple(part.bit_count() % 2 for part in value) == (int(x == v and y == u), 0)
        assert all(sum((part >> r) & 1 for r in range(1, cap + 1, 2)) % 2 == 0 for part in value)
        total_trace = add(total_trace, trace_column(entry, value, cap + 1))
        for exponent in range(max(base['degree'], cap) + 1):
            protected = (exponent + h[x] + h[y] < zero_cutoff
                         or base['degree'] - exponent + h[u] + h[v] < infinity_cutoff)
            wanted = coefficient(base['values'][entry], exponent)
            if protected and (protection_mode == 'exact' or not wanted):
                assert coefficient(value, exponent) == wanted
    assert total_trace == trace_target(trace_power, cap + 1)
    for condition in conditions:
        point, cut, first, second, order = condition
        powers = [1]
        for _ in range(cap):
            powers.append(F4_PRODUCT[powers[-1]][point])
        for entry in entries_for_rectangle(condition):
            for jet in range(order):
                total = 0
                for exponent in range(jet, cap + 1):
                    if comb(exponent, jet) % 2:
                        total ^= F4_PRODUCT[coefficient(seed['values'][entry], exponent)][powers[exponent - jet]]
                assert total == 0
    for prescribed in prescribed_jets:
        point = prescribed['point']
        powers = [1]
        for _ in range(cap):
            powers.append(F4_PRODUCT[powers[-1]][point])
        for jet, wanted in enumerate(prescribed['values']):
            for entry in range(256):
                total = 0
                for exponent in range(jet, cap + 1):
                    if comb(exponent, jet) % 2:
                        total ^= F4_PRODUCT[coefficient(seed['values'][entry], exponent)][powers[exponent - jet]]
                assert total == wanted[entry]
