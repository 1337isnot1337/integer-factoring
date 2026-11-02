"""Project supplied generators using certified high two-power row roots.

The ambient finite abelian group, canonical identifiers, and its size bound
are promises. Root identities and the required precision are checked here;
this module does not acquire roots or compute a group order.
"""

from check_positive_relation_kernel import determinant
from odd_quotient_relations import word_value


def adjugate(rows):
    size = len(rows)
    return [[(-1)**(i+j)*determinant([
        [entry for column, entry in enumerate(row) if column != i]
        for index, row in enumerate(rows) if index != j])
        for j in range(size)] for i in range(size)]


def two_power(value, depth, multiply):
    """Raise to 2**depth using depth squarings, without materializing 2**depth."""
    for _ in range(depth):
        value = multiply(value, value)
    return value


def project_from_high_roots(generators, rows, roots, root_depth,
                            ambient_order_bit_bound, identity, multiply, inverse):
    """Return generators of the full two-primary subgroup of <generators>.

    All generators AND roots must lie in one finite abelian ambient group
    of order at most 2**ambient_order_bit_bound. This is a promise, not a
    black-box size test. For d=det(rows)=2**s*u, require root_depth>=b+s
    and check roots[i]**(2**root_depth)==word_value(generators, rows[i]).
    The returned elements have order dividing 2**s. Cost charges the
    NUMERICAL root_depth, not merely its binary encoding length.
    """
    size = len(generators)
    if (not size or len(rows) != size or len(roots) != size
            or any(len(row) != size for row in rows)
            or any(not isinstance(x, int) for row in rows for x in row)):
        raise ValueError('A nonempty square integer matrix and one root per row are required')
    if (not isinstance(root_depth, int) or root_depth < 0
            or not isinstance(ambient_order_bit_bound, int) or ambient_order_bit_bound < 0):
        raise ValueError('Nonnegative integer depth and ambient size bounds are required')
    det = determinant(rows)
    if not det:
        raise ValueError('The row matrix must be nonsingular')
    odd, valuation = det, 0
    while odd % 2 == 0:
        odd //= 2
        valuation += 1
    if root_depth < ambient_order_bit_bound+valuation:
        raise ValueError('Root depth must cover the ambient bound plus v2(det R)')
    for row, root in zip(rows, roots):
        if two_power(root, root_depth, multiply) != word_value(
                generators, row, identity, multiply, inverse):
            raise ValueError('A supplied root does not satisfy its row identity')
    projected = []
    for generator, row in zip(generators, adjugate(rows)):
        root_word = word_value(roots, row, identity, multiply, inverse)
        correction = two_power(root_word, root_depth-valuation, multiply)
        original = word_value([generator], [odd], identity, multiply, inverse)
        projected.append(multiply(original, inverse(correction)))
    return projected
