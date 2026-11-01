"""Constant linear physical restrictions of binary elementary crossings.

Vectors and covectors over F2 are represented by integer bit masks.
An embedding P and retraction L with LP=I preserve both the identity
at zero and the literal swap endpoint at one.
"""

from functools import lru_cache
from itertools import product

from check_block_smith_crossings import reverse_polynomial, smith_profile
from check_cancellation_module_automata import pmul


def kronecker_vector(left, right, width):
    result = 0
    while left:
        index = (left & -left).bit_length()-1
        result |= right << (width*index)
        left &= left-1
    return result


@lru_cache(maxsize=None)
def elementary_coefficient_columns(length):
    """The coefficient matrices of the literal tensor, with columns packed."""
    if not isinstance(length, int) or length < 1:
        raise ValueError('The number of elementary sites must be positive.')
    alphabet = 1 << length
    columns = [[0]*(alphabet*alphabet) for _ in range(length+1)]
    for x, y, u, v in product(range(alphabet), repeat=4):
        polynomial = 1
        for bit in range(length):
            xb, yb, ub, vb = ((word >> bit) & 1 for word in (x, y, u, v))
            straight = int(xb == ub and yb == vb)
            crossing = int(xb == vb and yb == ub)
            polynomial = pmul(polynomial, straight | ((straight ^ crossing) << 1))
        for degree in range(length+1):
            if polynomial >> degree & 1:
                columns[degree][alphabet*u+v] |= 1 << (alphabet*x+y)
    return tuple(map(tuple, columns))


def linear_restriction_matrices(length, embedding, retraction):
    """Construct both polynomial cut matrices for a supplied constant retraction."""
    columns = elementary_coefficient_columns(length)
    alphabet, dimension = 1 << length, len(embedding)
    if not 1 <= dimension <= alphabet or len(retraction) != dimension:
        raise ValueError('The embedding and retraction need equal valid dimensions.')
    if any(not isinstance(value, int) or not 0 <= value < (1 << alphabet)
           for value in tuple(embedding)+tuple(retraction)):
        raise ValueError('Every vector or covector must fit the ambient physical space.')
    if any((retraction[i] & embedding[j]).bit_count() % 2 != int(i == j)
           for i, j in product(range(dimension), repeat=2)):
        raise ValueError('The supplied maps do not satisfy LP=I.')
    output_rows = [kronecker_vector(left, right, alphabet)
                   for left, right in product(retraction, repeat=2)]
    input_columns = [kronecker_vector(left, right, alphabet)
                     for left, right in product(embedding, repeat=2)]
    size = dimension*dimension
    row = [[0]*size for _ in range(size)]
    for j, vector in enumerate(input_columns):
        for degree, coefficients in enumerate(columns):
            image, remaining = 0, vector
            while remaining:
                index = (remaining & -remaining).bit_length()-1
                image ^= coefficients[index]
                remaining &= remaining-1
            for i, covector in enumerate(output_rows):
                if (image & covector).bit_count() % 2:
                    row[i][j] |= 1 << degree
    column = [[0]*size for _ in range(size)]
    for x, y, u, v in product(range(dimension), repeat=4):
        polynomial = row[dimension*x+y][dimension*u+v]
        assert polynomial & 1 == int(x == u and y == v)
        assert polynomial.bit_count() % 2 == int(x == v and y == u)
        column[dimension*x+u][dimension*y+v] = polynomial
    return {'row': row, 'column': column, 'degree_bound': length,
            'ambient_alphabet': alphabet, 'physical_alphabet': dimension}


def linear_restriction_profiles(length, embedding, retraction):
    """Return the two good-orientation profiles, with an exact degree bound."""
    matrices = linear_restriction_matrices(length, embedding, retraction)
    reversed_row = [[reverse_polynomial(value, length) for value in line]
                    for line in matrices['row']]
    row_orders, row_precision = smith_profile(reversed_row, length)
    column_orders, column_precision = smith_profile(matrices['column'], length)
    return {
        'row_infinity_orders': row_orders, 'column_zero_orders': column_orders,
        'row_degree_scores': tuple(length-order for order in row_orders),
        'row_precision': row_precision, 'column_precision': column_precision,
        'degree_bound': length, 'physical_alphabet': len(embedding),
    }
