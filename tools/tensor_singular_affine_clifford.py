#!/usr/bin/env python3
"""Exact nonzero Smith profiles for singular one-mode/Gaussian circuits.

Regularize local matrices by epsilon I, encode the lexicographic valuation
with a proved Kronecker bound, and eliminate a 2(n+1)-row Clifford matrix.
The algorithm uses exact polynomial divisions and has no series rank promise.
Polynomial degrees, including expanded binary-encoded exponents, are charged.
"""
from tensor_affine_clifford_profiles import (
    ZERO, ONE, add, multiply, valuation, matrix_product, local_numerator,
)


def degree(value):
    return (value[0] | value[1]).bit_length() - 1


def scalar_product(a, b):
    a0, a1, b0, b1 = a & 1, a >> 1, b & 1, b >> 1
    return ((a0*b0) ^ (a1*b1)) | (((a0*b1) ^ (a1*b0) ^ (a1*b1)) << 1)


def scaled(value, coefficient):
    a, b = value
    if coefficient == 0:
        return ZERO
    if coefficient == 1:
        return value
    if coefficient == 2:
        return b, a ^ b
    if coefficient == 3:
        return a ^ b, a
    raise ValueError('An F4 coefficient must be 0, 1, 2 or 3.')


def exact_divide(numerator, denominator):
    """Polynomial long division, rejecting a nonzero remainder."""
    if denominator == ZERO:
        raise ZeroDivisionError('Polynomial division by zero.')
    if denominator == ONE or numerator == ZERO:
        return numerator
    result = ZERO
    remaining = numerator
    d = degree(denominator)
    leading = ((denominator[0] >> d) & 1) | (((denominator[1] >> d) & 1) << 1)
    inverse = scalar_product(leading, leading)
    while degree(remaining) >= d:
        r = degree(remaining)
        c = ((remaining[0] >> r) & 1) | (((remaining[1] >> r) & 1) << 1)
        c = scalar_product(c, inverse)
        power = r - d
        term = scaled(ONE, c)
        result = add(result, (term[0] << power, term[1] << power))
        subtract = scaled(denominator, c)
        remaining = add(remaining, (subtract[0] << power, subtract[1] << power))
    if remaining != ZERO:
        raise ArithmeticError('The Bareiss polynomial division was not exact.')
    return result


def polynomial_smith_orders(matrix):
    """All finite t-adic Smith orders, with generic rank found exactly.

    Complete minimum-valuation pivoting in fraction-free Bareiss elimination
    is local Smith elimination. Each pivot is a determinantal minor and its
    valuation minus the preceding pivot valuation is the next Smith order.
    """
    if not matrix:
        return []
    width = len(matrix[0])
    if any(len(row) != width for row in matrix):
        raise ValueError('A rectangular polynomial matrix is required.')
    current = [list(row) for row in matrix]
    previous, previous_order, orders = ONE, 0, []
    while current and current[0]:
        candidates = ((valuation(value), i, j) for i, row in enumerate(current)
                      for j, value in enumerate(row) if value != ZERO)
        found = min(candidates, default=None)
        if found is None:
            break
        _, i, j = found
        current[0], current[i] = current[i], current[0]
        for row in current:
            row[0], row[j] = row[j], row[0]
        pivot = current[0][0]
        pivot_order = valuation(pivot)
        orders.append(pivot_order - previous_order)
        following = []
        for row in current[1:]:
            following.append([
                exact_divide(add(multiply(pivot, row[k]),
                                 multiply(row[0], current[0][k])), previous)
                for k in range(1, len(current[0]))])
        previous, previous_order, current = pivot, pivot_order, following
    if orders != sorted(orders):
        raise ArithmeticError('Minimum-valuation elimination lost Smith ordering.')
    return orders


def _polynomial(value):
    if (not isinstance(value, (tuple, list)) or len(value) != 2
            or any(type(x) is not int or x < 0 for x in value)):
        raise ValueError('Supply nonnegative integer F4 bit planes.')
    return tuple(value)


def _factors(n, factors):
    if type(n) is not int or n < 1:
        raise ValueError('The number of physical modes must be a positive integer.')
    output = []
    numerator_degree = denominator_degree = 0
    for factor in factors:
        kind = factor['kind']
        if kind == 'local':
            i = factor['mode']
            if type(i) is not int or not 0 <= i < n:
                raise ValueError('A local mode index is out of range.')
            matrix = factor['matrix']
            if len(matrix) != 2 or any(len(row) != 2 for row in matrix):
                raise ValueError('A local matrix must be two by two.')
            matrix = [[_polynomial(x) for x in row] for row in matrix]
            d = max(0, max(degree(x) for row in matrix for x in row))
            numerator_degree += 2*d
            denominator_degree += 2*d
            output.append({'kind': kind, 'mode': i, 'matrix': matrix})
        elif kind in ('create_pair', 'annihilate_pair', 'mode_shear'):
            i, j = factor['modes']
            if any(type(x) is not int or not 0 <= x < n for x in (i, j)) or i == j:
                raise ValueError('A Gaussian shear requires two distinct valid modes.')
            coefficient = _polynomial(factor['coefficient'])
            numerator_degree += max(0, degree(coefficient))
            output.append({'kind': kind, 'modes': [i, j], 'coefficient': coefficient})
        elif kind == 'diagonal':
            weights = factor['weights']
            if len(weights) != n or any(type(x) is not int for x in weights):
                raise ValueError('Supply one integer exterior weight per physical mode.')
            bound = max(map(abs, weights), default=0)
            numerator_degree += 2*bound
            denominator_degree += bound
            output.append({'kind': kind, 'weights': list(weights)})
        elif kind == 'scalar':
            power = factor['power']
            if type(power) is not int:
                raise ValueError('The scalar Laurent power must be an integer.')
            output.append({'kind': kind, 'power': power})
        else:
            raise ValueError('Unknown circuit factor kind.')
    return output, numerator_degree, denominator_degree


def _identity(n, diagonal=ONE):
    return [[diagonal if i == j else ZERO for j in range(n)] for i in range(n)]


def singular_profile(n, factors, *, expand=False):
    """Return a compact exact profile. Factors act successively on the left.

    Supported gates are arbitrary polynomial one-mode matrices, Gaussian pair
    and mode shears, exterior Laurent diagonals and Laurent scalar powers.
    Setting expand=True explicitly requests the possibly exponential multiset.
    """
    factors, numerator_degree, denominator_degree = _factors(n, factors)
    modes, size = n + 1, 2*(n + 1)
    bound = size*numerator_degree + denominator_degree
    weight = 2*bound + 1
    epsilon = (1 << weight, 0)
    numerator = _identity(size)
    denominator_order = 0
    mean_twice = [0, 0]
    for factor in factors:
        kind = factor['kind']
        if kind == 'scalar':
            mean_twice[1] += 2*factor['power']
            continue
        if kind == 'local':
            original = factor['matrix']
            determinant = add(multiply(original[0][0], original[1][1]),
                              multiply(original[0][1], original[1][0]))
            trace = add(original[0][0], original[1][1])
            for power, coefficient in enumerate((determinant, trace, ONE)):
                if coefficient != ZERO:
                    mean_twice[0] += power
                    mean_twice[1] += valuation(coefficient)
                    break
            local = [list(row) for row in original]
            local[0][0] = add(local[0][0], epsilon)
            local[1][1] = add(local[1][1], epsilon)
            block, denominator = local_numerator(local)
            gate = _identity(size, denominator)
            positions = [factor['mode'], modes + factor['mode'], n, modes+n]
            for i, a in enumerate(positions):
                for j, b in enumerate(positions):
                    gate[a][b] = block[i][j]
            denominator_order += valuation(denominator)
        elif kind == 'diagonal':
            weights = factor['weights'] + [0]
            clearing = max(map(abs, weights))
            gate = [[ZERO]*size for _ in range(size)]
            for i, value in enumerate(weights):
                gate[i][i] = 1 << (clearing+value), 0
                gate[modes+i][modes+i] = 1 << (clearing-value), 0
            denominator_order += clearing
            mean_twice[1] += sum(weights)
        else:
            gate = _identity(size)
            i, j = factor['modes']
            coefficient = factor['coefficient']
            if kind == 'create_pair':
                gate[j][modes+i] = gate[i][modes+j] = coefficient
            elif kind == 'annihilate_pair':
                gate[modes+j][i] = gate[modes+i][j] = coefficient
            else:
                gate[i][j] = gate[modes+j][modes+i] = coefficient
        numerator = matrix_product(gate, numerator)
    numerator_orders = polynomial_smith_orders(numerator)
    if len(numerator_orders) != size:
        raise ArithmeticError('The regularized Clifford matrix lost its proved full rank.')
    encoded = [value-denominator_order for value in numerator_orders]
    signed = []
    for value in encoded:
        first = (value + weight//2)//weight
        second = value - weight*first
        if abs(second) > bound:
            raise ArithmeticError('Kronecker decoding exceeded its determinantal degree bound.')
        signed.append((first, second))
    if signed != sorted(signed) or signed != [(-a, -b) for a, b in signed[::-1]]:
        raise ArithmeticError('The complete quadratic Clifford profile is not paired.')
    gaps = signed[modes:]
    if not gaps or gaps[0] != (0, 0):
        raise ArithmeticError('The exact auxiliary zero gap is missing.')
    gaps = gaps[1:]
    difference = [mean_twice[k] - sum(gap[k] for gap in gaps) for k in (0, 1)]
    if any(value % 2 for value in difference):
        raise ArithmeticError('The determinant normalization gave a nonintegral minimum.')
    minimum = tuple(value//2 for value in difference)
    if minimum[0] < 0:
        raise ArithmeticError('An epsilon-polynomial Fock matrix acquired an epsilon pole.')
    surviving = [gap[1] for gap in gaps if gap[0] == 0]
    if any(value < 0 for value in surviving):
        raise ArithmeticError('A surviving lexicographically nonnegative gap became negative.')
    rank = 0 if minimum[0] else 1 << len(surviving)
    result = {
        'physical_modes': n, 'clifford_matrix_order': size,
        'rank': rank, 'minimum_nonzero_order': minimum[1] if rank else None,
        'nonzero_binomial_gaps': surviving if rank else [],
        'regularized_minimum_pair': list(minimum),
        'regularized_binomial_gaps': [list(gap) for gap in gaps],
        'clifford_signed_pairs': [list(pair) for pair in signed],
        'determinant_mean_twice_pair': mean_twice,
        'kronecker_weight': weight, 'signed_second_coordinate_bound': bound,
        'bivariate_numerator_t_degree_bound': numerator_degree,
        'bivariate_denominator_t_degree_bound': denominator_degree,
        'encoded_numerator_degree': max(degree(x) for row in numerator for x in row),
        'encoded_common_denominator_order': denominator_order,
        'exact_polynomial_rank': len(numerator_orders),
        'auxiliary_smith_multiplicity': 2,
    }
    if expand:
        orders = [minimum[1]] if rank else []
        for gap in surviving if rank else []:
            orders += [value+gap for value in orders]
        result['nonzero_orders'] = sorted(orders)
    return result
