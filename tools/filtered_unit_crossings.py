#!/usr/bin/env python3
"""Filtered unit curves in a characteristic-two Boolean group algebra.

GF(4)=F2[w]/(w^2+w+1) is used for exact finite controls. Scalars 0 and 1
give the binary subfamily. All constructor and rank operations are exact.
"""


def multiply_field(a, b):
    a0,a1=a&1,a>>1
    b0,b1=b&1,b>>1
    return (a0*b0 ^ a1*b1) | ((a0*b1 ^ a1*b0 ^ a1*b1) << 1)


def parameter_positions(width):
    if not isinstance(width, int) or width < 1:
        raise ValueError('The binary width must be a positive integer')
    return [(mask, degree) for mask in range(1, 1 << width)
            for degree in range(1, mask.bit_count())]


def coefficient_polynomials(width, state, parameters):
    """Return f_J(u) as a (u degree) by (nilpotent monomial) table."""
    positions = parameter_positions(width)
    dimension = 1 << width
    if not isinstance(state, int) or not 0 <= state < dimension:
        raise ValueError('State is outside the label space')
    if len(parameters) != len(positions) or any(type(x) is not int or not 0 <= x < 4 for x in parameters):
        raise ValueError('Supply one GF(4) scalar for every free coefficient')
    answer = [[0] * dimension for _ in range(width + 1)]
    answer[0][0] = 1
    for mask in range(1, dimension):
        if mask & state == mask:
            answer[mask.bit_count()][mask] = 1
    for (mask, degree), value in zip(positions, parameters):
        answer[degree][mask] ^= value
        answer[mask.bit_count()][mask] ^= value
    return answer


def group_kernel(width, state, parameters):
    """Return coefficients of z^width g_state P_state((1+z)/z)."""
    dimension = 1 << width
    coefficients = coefficient_polynomials(width, state, parameters)
    answer = [[0] * dimension for _ in range(width + 1)]
    for degree, row in enumerate(coefficients):
        for mask, scalar in enumerate(row):
            if not scalar:
                continue
            subset = mask
            while True:
                for power in range(degree + 1):
                    if power & degree == power:
                        answer[width - degree + power][subset ^ state] ^= scalar
                if subset == 0:
                    break
                subset = (subset - 1) & mask
    return answer


def tensor_kernels(width, parameters):
    dimension = 1 << width
    if len(parameters) != dimension:
        raise ValueError('Supply a parameter row for each logical state')
    return [group_kernel(width, state, row) for state, row in enumerate(parameters)]


def rank(matrix):
    a = [row[:] for row in matrix]
    if not a:
        return 0
    rows, columns = len(a), len(a[0])
    result = 0
    for column in range(columns):
        pivot = next((i for i in range(result, rows) if a[i][column]), None)
        if pivot is None:
            continue
        a[result], a[pivot] = a[pivot], a[result]
        inverse = 1 if a[result][column] == 1 else 5 - a[result][column]
        a[result] = [multiply_field(value, inverse) for value in a[result]]
        for i in range(result + 1, rows):
            scalar = a[i][column]
            if scalar:
                a[i] = [value ^ multiply_field(scalar, other) for value, other in zip(a[i], a[result])]
        result += 1
        if result == rows:
            break
    return result


def convolution_rank(values):
    dimension = len(values)
    return rank([[values[i ^ j] for j in range(dimension)] for i in range(dimension)])


def coefficient_ranks(kernels):
    dimension = len(kernels)
    degree = len(kernels[0]) - 1
    rows, columns = [], []
    for power in range(degree + 1):
        rows.append(sum(convolution_rank(kernels[state][power]) for state in range(dimension)))
        columns.append(sum(convolution_rank([kernels[state][power][word] for state in range(dimension)]) for word in range(dimension)))
    return rows, columns


def dense_coefficient(kernels, power, cut):
    """The literal tensor at (x,y)|(u,v), or (x,u)|(y,v)."""
    dimension = len(kernels)
    left, right = ((0,1),(2,3)) if cut == 0 else ((0,2),(1,3))
    matrix = []
    for a in range(dimension):
        for b in range(dimension):
            row = []
            for c in range(dimension):
                for d in range(dimension):
                    letters = [0] * 4
                    letters[left[0]], letters[left[1]] = a,b
                    letters[right[0]], letters[right[1]] = c,d
                    x,y,u,v = letters
                    row.append(kernels[x ^ y][power][x ^ u] if x ^ y == u ^ v else 0)
            matrix.append(row)
    return matrix
