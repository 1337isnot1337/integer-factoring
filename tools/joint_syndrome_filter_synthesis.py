#!/usr/bin/env python3
"""Joint exact polynomial synthesis on a binary syndrome alphabet."""

from functools import lru_cache
from heapq import heappop, heappush
from math import log
import random

from boolean_multiplier_profiles import code_group_signatures_for_types, signature_patterns
from polynomial_state_filter_synthesis import residue_patterns
from weyl_lattice_filter_synthesis import endpoint_image


@lru_cache(maxsize=128)
def coefficient_tables(rank, column_types, spacing, degrees):
    if spacing < 1 or len(degrees) != 1 << rank or any(value < 0 for value in degrees):
        raise ValueError('Use a positive spacing and one nonnegative degree cap per syndrome.')
    signatures = code_group_signatures_for_types(rank, column_types)
    patterns = tuple(signature_patterns(signatures, rank, syndrome)
                     for syndrome in range(1 << rank))
    variables = tuple((syndrome, degree) for syndrome, cap in enumerate(degrees)
                      for degree in range(cap + 1))
    endpoint_masks = [0] * (1 << rank)
    maximum = max(degrees) + spacing * rank
    equations = [[0] * (maximum + 1) for _ in signatures]
    for column, (syndrome, degree) in enumerate(variables):
        bit = 1 << column
        endpoint_masks[syndrome] |= bit
        for state, pattern in enumerate(patterns[syndrome]):
            for layer in range(rank + 1):
                if pattern & (1 << layer):
                    equations[state][degree + spacing * (rank - layer)] ^= bit
    return variables, tuple(endpoint_masks), tuple(map(tuple, equations)), patterns


def reduce_row(row, pivots):
    for pivot in sorted(pivots, reverse=True):
        if row & (1 << pivot):
            row ^= pivots[pivot]
    return row


def add_equation(row, pivots, endpoint_remainders):
    reduced = reduce_row(row, pivots)
    if not reduced:
        return endpoint_remainders
    pivot = reduced.bit_length() - 1
    remainders = tuple(value ^ reduced if value & (1 << pivot) else value
                       for value in endpoint_remainders)
    if not all(remainders):
        return None
    pivots[pivot] = reduced
    return remainders


def transform(polynomials, rank, column_types, spacing):
    signatures = code_group_signatures_for_types(rank, column_types)
    patterns = tuple(signature_patterns(signatures, rank, syndrome)
                     for syndrome in range(1 << rank))
    output = []
    for state in range(len(signatures)):
        value = {}
        for syndrome, polynomial in enumerate(polynomials):
            pattern = patterns[syndrome][state]
            for layer in range(rank + 1):
                if pattern & (1 << layer):
                    shift = spacing * (rank - layer)
                    for degree, coefficient in polynomial.items():
                        at = degree + shift
                        value[at] = value.get(at, 0) ^ coefficient
        output.append({degree: coefficient for degree, coefficient in value.items() if coefficient})
    return tuple(output)


def construct_witness(rank, column_types, spacing, degrees, targets, pivots):
    variables, masks, _, _ = coefficient_tables(rank, column_types, spacing, tuple(degrees))
    image = endpoint_image(variables, pivots, masks)
    support = 0
    for endpoint, _ in image:
        support |= endpoint
    if support != (1 << (1 << rank)) - 1:
        return None
    field_size = 1 << (rank + 1)
    endpoints = [0] * (1 << rank)
    coefficients = {}
    for vector, lift in image:
        union = tuple(index for index in range(len(endpoints)) if endpoints[index] or vector & (1 << index))
        if all(endpoints[index] for index in union):
            continue
        scalar = next(value for value in range(field_size)
                      if all(endpoints[index] ^ (value if vector & (1 << index) else 0)
                             for index in union))
        for index in range(len(endpoints)):
            endpoints[index] ^= scalar if vector & (1 << index) else 0
        for index, variable in enumerate(variables):
            if lift & (1 << index):
                coefficients[variable] = coefficients.get(variable, 0) ^ scalar
                if not coefficients[variable]:
                    del coefficients[variable]
        if all(endpoints):
            break
    assert all(endpoints)
    f = tuple({degree: coefficient for (syndrome, degree), coefficient in coefficients.items()
               if syndrome == index} for index in range(1 << rank))
    delta = transform(f, rank, column_types, spacing)
    assert all(delta)
    raw_degrees, raw_orders = tuple(map(max, f)), tuple(map(min, delta))
    assert all(actual <= bound for actual, bound in zip(raw_degrees, degrees))
    assert all(actual >= bound for actual, bound in zip(raw_orders, targets))
    common = min(map(min, f))
    f = tuple({degree - common: coefficient for degree, coefficient in value.items()} for value in f)
    delta = tuple({degree - common: coefficient for degree, coefficient in value.items()} for value in delta)
    return {
        'rank': rank, 'column_types': column_types, 'spacing': spacing,
        'coefficient_field_order': field_size, 'degree_bounds': tuple(degrees),
        'order_bounds': tuple(targets), 'removed_common_monomial_order': common,
        'actual_degrees': tuple(degree - common for degree in raw_degrees),
        'actual_orders': tuple(value - common for value in raw_orders),
        'endpoint_values': tuple(endpoints), 'f': f, 'delta': delta,
        'constraint_rank': len(pivots), 'coefficient_variables': len(variables),
        'endpoint_image_dimension': len(image),
    }


def synthesize(rank, column_types, spacing, degrees, targets, construct=True):
    column_types, degrees, targets = tuple(column_types), tuple(degrees), tuple(targets)
    variables, masks, equations, _ = coefficient_tables(rank, column_types, spacing, degrees)
    if len(targets) != len(equations) or any(value < 0 for value in targets):
        raise ValueError('Use one nonnegative target per parity-group state.')
    pivots, remainders = {}, masks
    for target, rows in zip(targets, equations):
        if target > len(rows):
            return None
        for row in rows[:target]:
            update = add_equation(row, pivots, remainders)
            if update is None:
                return None
            remainders = update
    if not construct:
        return True
    return construct_witness(rank, column_types, spacing, degrees, targets, pivots)


def greedy_synthesize(rank, column_types, spacing, degrees, probabilities, tau,
                      seed=0, priority_noise=0.0, initial_targets=None):
    """A weighted target walk; feasibility and the final witness are exact."""
    column_types, degrees = tuple(column_types), tuple(degrees)
    variables, masks, equations, _ = coefficient_tables(rank, column_types, spacing, degrees)
    if len(probabilities) != len(equations) or tau <= 0:
        raise ValueError('Supply one physical probability per parity-group state and positive tau.')
    generator = random.Random(seed)
    priorities = [log(value) if value > 0 else -float('inf') for value in probabilities]
    priorities = [value + priority_noise * generator.gauss(0, 1) for value in priorities]
    pivots, remainders = {}, masks
    targets = [0] * len(equations) if initial_targets is None else list(initial_targets)
    if len(targets) != len(equations) or any(value < 0 for value in targets):
        raise ValueError('Initial targets must be nonnegative and cover every output.')
    for target, rows in zip(targets, equations):
        if target >= len(rows):
            return None
        for row in rows[:target]:
            update = add_equation(row, pivots, remainders)
            if update is None:
                return None
            remainders = update
    queue = [(-priority + 2 * tau * targets[state] / spacing, state)
             for state, priority in enumerate(priorities)]
    from heapq import heapify
    heapify(queue)
    attempts = accepted = 0
    while queue:
        _, state = heappop(queue)
        target = targets[state]
        assert target < len(equations[state])
        attempts += 1
        update = add_equation(equations[state][target], pivots, remainders)
        if update is None:
            continue
        remainders = update
        targets[state] += 1
        accepted += 1
        heappush(queue, (-priorities[state] + 2 * tau * targets[state] / spacing, state))
    witness = construct_witness(rank, column_types, spacing, degrees, targets, pivots)
    assert witness
    witness.update(target_attempts=attempts, accepted_targets=accepted,
                   priority_seed=seed, priority_noise=priority_noise)
    return witness


def basis_release_targets(rank, column_types, spacing, degrees, targets):
    """Explicit feasible weakenings that release a coefficient-space direction.

    A dual vector evaluates to one on one echelon-basis row and zero on all
    other basis rows. For each output, keep only the prefix of constraints
    annihilating that vector. The remaining row space is a proper subspace.
    """
    _, masks, equations, _ = coefficient_tables(rank, tuple(column_types), spacing, tuple(degrees))
    pivots, remainders = {}, masks
    for target, rows in zip(targets, equations):
        for row in rows[:target]:
            update = add_equation(row, pivots, remainders)
            if update is None:
                raise ValueError('The source targets are infeasible.')
            remainders = update
    result = []
    for selected in sorted(pivots):
        vector = 0
        for pivot in sorted(pivots):
            required = int(pivot == selected)
            if ((pivots[pivot] & vector).bit_count() & 1) ^ required:
                vector ^= 1 << pivot
        assert all(((row & vector).bit_count() & 1) == int(pivot == selected)
                   for pivot, row in pivots.items())
        weakened = tuple(next((degree for degree, row in enumerate(rows[:target])
                               if (row & vector).bit_count() & 1), target)
                         for target, rows in zip(targets, equations))
        assert any(new < old for new, old in zip(weakened, targets))
        result.append({'released_pivot': selected, 'dual_vector': vector,
                       'order_bounds': weakened})
    return tuple(result)


@lru_cache(maxsize=4096)
def local_endpoint_image(rank, column_types, degrees, targets):
    """One residue class; degree -1 means that a coordinate is absent."""
    signatures = code_group_signatures_for_types(rank, column_types)
    patterns = tuple(signature_patterns(signatures, rank, syndrome)
                     for syndrome in range(1 << rank))
    variables = tuple((syndrome, degree) for syndrome, cap in enumerate(degrees)
                      for degree in range(max(0, cap + 1)))
    masks = [0] * (1 << rank)
    equations = {}
    for column, (syndrome, degree) in enumerate(variables):
        bit = 1 << column
        masks[syndrome] |= bit
        for state, pattern in enumerate(patterns[syndrome]):
            for layer in range(rank + 1):
                if pattern & (1 << layer):
                    at = degree + rank - layer
                    if at < targets[state]:
                        key = state, at
                        equations[key] = equations.get(key, 0) ^ bit
    pivots = {}
    for equation in equations.values():
        while equation:
            pivot = equation.bit_length() - 1
            if pivot in pivots:
                equation ^= pivots[pivot]
            else:
                pivots[pivot] = equation
                break
    return variables, endpoint_image(variables, pivots, masks), len(pivots)


def synthesize_compressed(rank, column_types, spacing, degrees, targets, construct=True):
    """Exact arbitrary-spacing feasibility with finitely many residue patterns."""
    column_types, degrees, targets = tuple(column_types), tuple(degrees), tuple(targets)
    if (spacing < 1 or len(degrees) != 1 << rank
            or len(targets) != 1 << len(column_types)
            or any(value < 0 for value in degrees + targets)):
        raise ValueError('Use nonnegative degree/order bounds and positive spacing.')
    patterns = residue_patterns(spacing, degrees, targets)
    basis = {}
    total_rank = maximum_local_variables = 0
    for index, (residue, local_degrees, local_targets) in enumerate(patterns):
        variables, image, constraint_rank = local_endpoint_image(rank, column_types, local_degrees, local_targets)
        next_residue = patterns[index + 1][0] if index + 1 < len(patterns) else spacing
        total_rank += (next_residue - residue) * constraint_rank
        maximum_local_variables = max(maximum_local_variables, len(variables))
        for endpoint, coefficients in image:
            lift = {(syndrome, residue + spacing * degree)
                    for column, (syndrome, degree) in enumerate(variables)
                    if coefficients & (1 << column)} if construct else set()
            while endpoint:
                pivot = endpoint.bit_length() - 1
                if pivot in basis:
                    old_endpoint, old_lift = basis[pivot]
                    endpoint ^= old_endpoint
                    if construct:
                        lift.symmetric_difference_update(old_lift)
                else:
                    basis[pivot] = endpoint, lift
                    break
    support = 0
    for endpoint, _ in basis.values():
        support |= endpoint
    if support != (1 << (1 << rank)) - 1:
        return None
    if not construct:
        return True
    field_size = 1 << (rank + 1)
    endpoints = [0] * (1 << rank)
    coefficients = {}
    for vector, lift in basis.values():
        union = tuple(index for index in range(len(endpoints)) if endpoints[index] or vector & (1 << index))
        if all(endpoints[index] for index in union):
            continue
        scalar = next(value for value in range(field_size)
                      if all(endpoints[index] ^ (value if vector & (1 << index) else 0)
                             for index in union))
        for index in range(len(endpoints)):
            endpoints[index] ^= scalar if vector & (1 << index) else 0
        for variable in lift:
            coefficients[variable] = coefficients.get(variable, 0) ^ scalar
            if not coefficients[variable]:
                del coefficients[variable]
        if all(endpoints):
            break
    assert all(endpoints)
    f = tuple({degree: coefficient for (syndrome, degree), coefficient in coefficients.items()
               if syndrome == index} for index in range(1 << rank))
    delta = transform(f, rank, column_types, spacing)
    assert all(delta)
    raw_degrees, raw_orders = tuple(map(max, f)), tuple(map(min, delta))
    assert all(actual <= bound for actual, bound in zip(raw_degrees, degrees))
    assert all(actual >= bound for actual, bound in zip(raw_orders, targets))
    common = min(map(min, f))
    return {
        'rank': rank, 'column_types': column_types, 'spacing': spacing,
        'coefficient_field_order': field_size, 'degree_bounds': degrees, 'order_bounds': targets,
        'removed_common_monomial_order': common,
        'actual_degrees': tuple(degree - common for degree in raw_degrees),
        'actual_orders': tuple(value - common for value in raw_orders),
        'endpoint_values': tuple(endpoints),
        'f': tuple({degree - common: coefficient for degree, coefficient in value.items()} for value in f),
        'delta': tuple({degree - common: coefficient for degree, coefficient in value.items()} for value in delta),
        'constraint_rank': total_rank, 'coefficient_variables': sum(degree + 1 for degree in degrees),
        'maximum_local_coefficient_variables': maximum_local_variables,
        'endpoint_image_dimension': len(basis), 'residue_pattern_count': len(patterns),
    }
