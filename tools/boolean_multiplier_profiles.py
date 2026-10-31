#!/usr/bin/env python3
"""Exact finite automata for one-support Boolean polynomial multipliers.

The homogeneous-layer pattern encodes g_0(k),...,g_d(k) as binary digits.
All polynomial arithmetic here is over F_2. Continuous cutoff profiles are
realized by a single nonbinary F_4 coefficient, as explained in the note.
"""

def order(polynomial):
    if not polynomial:
        raise ValueError('The zero polynomial has no finite order.')
    return (polynomial & -polynomial).bit_length() - 1


def multiply(left, right):
    result = 0
    while right:
        low = right & -right
        result ^= left << (low.bit_length() - 1)
        right ^= low
    return result


def rational_target(pattern):
    """Return a,A for S(X)=X^a/A(X), with A(0)=1."""
    if pattern <= 0:
        raise ValueError('The zero layer pattern has no rational target.')
    degree = pattern.bit_length() - 1
    denominator = sum(((pattern >> bit) & 1) << (degree - bit)
                      for bit in range(degree + 1))
    assert denominator & 1
    return degree, denominator


def series_prefix(pattern, length):
    degree, denominator = rational_target(pattern)
    remainder, result = 1 << degree, 0
    for position in range(length):
        coefficient = remainder & 1
        result |= coefficient << position
        remainder = (remainder ^ (denominator if coefficient else 0)) >> 1
    return result


def layer_pattern_excess(pattern, multiplier):
    """Return ord(beta_k)-|k|, independent of the actual occupation |k|."""
    if not pattern:
        return 0
    degree, denominator = rational_target(pattern)
    polynomial = (1 << degree) ^ multiply(multiplier, denominator)
    return order(polynomial) - degree


def multiplier_profile(multiplier, maximum_boolean_degree):
    if not multiplier or multiplier.bit_count() % 2:
        raise ValueError('Use a nonzero multiplier with R(1)=0.')
    return {
        'multiplier': multiplier,
        'degree': multiplier.bit_length() - 1,
        'excesses': tuple(layer_pattern_excess(pattern, multiplier)
                          for pattern in range(1 << (maximum_boolean_degree + 1))),
    }


def finite_degree_representatives(maximum_boolean_degree, maximum_degree):
    """Cover all endpoint-trivial binary R of degree <=maximum_degree.

    A branch with no remaining rational target has fixed orders, so only one
    feasible endpoint-trivial completion of that branch is needed.
    """
    targets = tuple(series_prefix(pattern, maximum_degree + 1)
                    for pattern in range(1, 1 << (maximum_boolean_degree + 1)))
    representatives = set()
    for degree in range(1, maximum_degree + 1):
        mask = (1 << (degree + 1)) - 1
        initial = tuple(sorted({target & mask for target in targets}))

        def visit(prefix, depth, active):
            if not active:
                candidate = prefix | (1 << degree)
                if candidate.bit_count() % 2:
                    if depth >= degree:
                        return
                    candidate ^= 1 << depth
                representatives.add(candidate)
                return
            if depth > degree:
                if prefix.bit_count() % 2 == 0:
                    representatives.add(prefix)
                return
            for coefficient in ((1,) if depth == degree else (0, 1)):
                visit(prefix | (coefficient << depth), depth + 1,
                      tuple(target for target in active
                            if ((target >> depth) & 1) == coefficient))

        visit(0, 0, initial)
    return tuple(sorted(representatives))


def minimal_divergence_multiplier(pattern, divergence):
    """Minimum-degree R(1)=0 whose first divergence from S is at divergence."""
    candidate = series_prefix(pattern, divergence + 1) ^ (1 << divergence)
    if candidate.bit_count() % 2:
        candidate ^= 1 << (divergence + 1)
    return candidate


def divergence_tail(pattern, start):
    """Finite startup and periodic degree offsets for all divergences >=start.

    Each row records D=h+offset. A zero candidate denotes the available
    baseline. The finite-state proof tracks the rational-series remainder,
    the prefix parity, and the distance to its last nonzero coefficient.
    """
    degree, denominator = rational_target(pattern)
    remainder, parity, gap = 1 << degree, 0, None

    def step(remainder, parity, gap):
        coefficient = remainder & 1
        next_remainder = (remainder ^ (denominator if coefficient else 0)) >> 1
        next_gap = 1 if coefficient else (None if gap is None else gap + 1)
        if not next_remainder:
            next_gap = 0  # Polynomial targets are monomials, with odd endpoint.
        return next_remainder, parity ^ coefficient, next_gap

    for _ in range(start):
        remainder, parity, gap = step(remainder, parity, gap)
    seen, rows = {}, []
    position = start
    while (remainder, parity, gap) not in seen:
        seen[remainder, parity, gap] = position
        coefficient = remainder & 1
        if parity == coefficient:
            offset = 1
        elif coefficient == 0:
            offset = 0
        else:
            offset = None if gap is None else -gap
        rows.append({'divergence': position, 'degree_offset': offset})
        remainder, parity, gap = step(remainder, parity, gap)
        position += 1
        if position > start + (1 << (degree + 2)) * (degree + 2) + degree + 4:
            raise ArithmeticError('Unexpected rational-series state count.')
    cycle_start = seen[remainder, parity, gap]
    period = position - cycle_start
    return {
        'startup': tuple(row for row in rows if row['divergence'] < cycle_start),
        'cycle': tuple(row for row in rows if row['divergence'] >= cycle_start),
        'cycle_start': cycle_start, 'period': period,
    }


def complete_multiplier_automaton(maximum_boolean_degree):
    """Finite exact lower-envelope description, with no multiplier degree cap.

    fixed_profiles cover branches that leave every target. continuous_profiles
    cover nonbinary/nonintegral first divergences. Singleton binary tails are
    finite startup rows and arithmetic progressions with periodic degree gaps.
    Profiles with D>=h in a singleton tail are dominated by its continuous ray.
    """
    count = 1 << (maximum_boolean_degree + 1)
    length = 2 * maximum_boolean_degree + 2
    targets = (0,) + tuple(series_prefix(pattern, length) for pattern in range(1, count))
    degrees = tuple(0 if not pattern else pattern.bit_length() - 1 for pattern in range(count))
    fixed = {}
    continuous = []
    singleton = []

    def coefficients(prefix, depth, active):
        mask = (1 << depth) - 1
        active_set = set(active)
        outside = []
        for pattern in range(count):
            if pattern == 0 or pattern in active_set:
                outside.append(0)
            else:
                mismatch = (targets[pattern] ^ prefix) & mask
                assert mismatch
                outside.append(degrees[pattern] - order(mismatch))
        return tuple(outside)

    def visit(prefix, depth, active):
        if depth:
            outside = coefficients(prefix, depth, active)
            continuous.append({
                'prefix': prefix, 'prefix_length': depth,
                'lower_degree': depth - 1,
                'upper_degree': None if len(active) == 1 else depth,
                'active_patterns': tuple(active), 'outside_exponents': outside,
            })
            if len(active) == 1:
                pattern = active[0]
                tail = divergence_tail(pattern, depth)
                singleton.append({
                    'pattern': pattern, 'prefix': prefix, 'prefix_length': depth,
                    'outside_exponents': outside,
                    'period': tail['period'], 'cycle_start': tail['cycle_start'],
                    'startup': tuple(row for row in tail['startup']
                                     if row['degree_offset'] is not None
                                     and row['degree_offset'] < 0),
                    'cycle': tuple(row for row in tail['cycle']
                                   if row['degree_offset'] is not None
                                   and row['degree_offset'] < 0),
                    'full_tail': tail,
                })
                return
        if depth >= length:
            raise ArithmeticError('Distinct rational targets failed to separate.')
        for coefficient in (0, 1):
            new_prefix = prefix | (coefficient << depth)
            child = tuple(pattern for pattern in active
                          if ((targets[pattern] >> depth) & 1) == coefficient)
            if child:
                visit(new_prefix, depth + 1, child)
            else:
                candidate = new_prefix
                if candidate.bit_count() % 2:
                    candidate ^= 1 << (depth + 1)
                if candidate:
                    profile = multiplier_profile(candidate, maximum_boolean_degree)
                    key = profile['degree'], profile['excesses']
                    fixed.setdefault(key, profile)

    visit(0, 0, tuple(range(1, count)))
    return {
        'maximum_boolean_degree': maximum_boolean_degree,
        'pattern_count': count,
        'fixed_profiles': tuple(fixed.values()),
        'continuous_profiles': tuple(continuous),
        'singleton_tails': tuple(singleton),
    }


def homogeneous_layer_patterns(truth, modes):
    """Independent Boolean Mobius transform followed by degree-wise evaluation."""
    if len(truth) != 1 << modes:
        raise ValueError('The truth table has the wrong size.')
    coefficients = list(truth)
    for bit in range(modes):
        for state in range(1 << modes):
            if state & (1 << bit):
                coefficients[state] ^= coefficients[state ^ (1 << bit)]
    result = [coefficient << state.bit_count()
              for state, coefficient in enumerate(coefficients)]
    for bit in range(modes):
        for state in range(1 << modes):
            if state & (1 << bit):
                result[state] ^= result[state ^ (1 << bit)]
    return tuple(result)


def group_translate(element, group_element):
    translated = 0
    while element:
        low = element & -element
        translated ^= 1 << ((low.bit_length() - 1) ^ group_element)
        element ^= low
    return translated


def code_group_signatures(rank):
    """All parity-type signatures in F_2[(Z/2)^rank], for rank <=4."""
    if not 1 <= rank <= 4:
        raise ValueError('This explicit group-state representation supports rank 1..4.')
    return code_group_signatures_for_types(rank, tuple(range(1, 1 << rank)))


def code_group_signatures_for_types(rank, column_types):
    """Parity-group signatures for specified distinct nonzero column types."""
    if len(set(column_types)) != len(column_types) or any(
            not 0 < value < (1 << rank) for value in column_types):
        raise ValueError('Use distinct nonzero columns in the specified binary group.')
    if len(column_types) > 16:
        raise ValueError('Use at most sixteen explicit parity types.')
    types = len(column_types)
    result = [(1,) + (0,) * rank]
    for state in range(1, 1 << types):
        low = state & -state
        group_element = column_types[low.bit_length() - 1]
        old = result[state ^ low]
        result.append((old[0],) + tuple(
            old[degree] ^ old[degree - 1] ^ group_translate(old[degree - 1], group_element)
            for degree in range(1, rank + 1)))
    return tuple(result)


def signature_patterns(signatures, rank, syndrome):
    return tuple(sum(((signature[degree] >> syndrome) & 1) << degree
                     for degree in range(rank + 1)) for signature in signatures)
