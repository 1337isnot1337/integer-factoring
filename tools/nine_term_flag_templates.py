"""Sparse nine-term templates from the first rank-nine two-copy image vector."""

from itertools import product


COMPONENTS = ((0, 0, 1, 1, 0), (1, 0, 0, 0, 1), (0, 1, 1, 1, 0))
COORDINATES = (0, 2, 4)
HEIGHTS = (0, 5, 3, 1, 4)
PAIRS = tuple(product(range(3), repeat=2))
ANTITONE = tuple(PAIRS.index((2 - i, 2 - j)) for i, j in PAIRS)


def tensor_vector(first, second):
    return tuple(a * b for a in first for b in second)


def component_data():
    result = []
    for i, j in PAIRS:
        left = tensor_vector(COMPONENTS[i], COMPONENTS[j])
        coordinate = 5 * COORDINATES[i] + COORDINATES[j]
        right = tuple(int(k == coordinate) for k in range(25))
        minimum = min(HEIGHTS[k // 5] + HEIGHTS[k % 5] for k, value in enumerate(left) if value)
        height = HEIGHTS[COORDINATES[i]] + HEIGHTS[COORDINATES[j]]
        result.append({'pair': (i, j), 'left': left, 'coordinate_vector': right,
                       'minimum_left_height': minimum, 'coordinate_height': height,
                       'balanced_height_sum': minimum + height})
    return result


def one_plus_power(power):
    """The exact binary polynomial (1+z)^power."""
    result = 1
    for _ in range(power):
        result ^= result << 1
    return result


def updates(pairing=None, zero_orders=21, coefficients=None, shape='binomial'):
    """Return nine updates accepted by LowRankTwoCopyFrames.profiles.

    Each infinity leading order is -10, relative order six. The requested
    zero_orders are literal minimum entry orders of the individual terms.
    For the antitone pairing, every integer16..21 is feasible for every term.
    shape is 'binomial' for z^ell+z^d, or 'endpoint_power' for
    z^ell(1+z)^(d-ell). All vectors use the retained source's label order.
    """
    pairing = ANTITONE if pairing is None else tuple(pairing)
    if sorted(pairing) != list(range(9)):
        raise ValueError('Supply a permutation of the nine component pairs.')
    if isinstance(zero_orders, int):
        zero_orders = (zero_orders,) * 9
    if len(zero_orders) != 9:
        raise ValueError('Supply one zero leading order per matched term.')
    coefficients = (1,) * 9 if coefficients is None else tuple(coefficients)
    if len(coefficients) != 9 or any(value not in (1, 2, 3) for value in coefficients):
        raise ValueError('Supply nine nonzero F4 coefficients.')
    if shape not in ('binomial', 'endpoint_power'):
        raise ValueError('Use binomial or endpoint_power endpoint-zero polynomials.')
    data = component_data()
    result = []
    for i, j in enumerate(pairing):
        a, b = data[i], data[j]
        degree = 10 + a['coordinate_height'] + b['coordinate_height']
        lower = zero_orders[i] - a['minimum_left_height'] - b['minimum_left_height']
        if not isinstance(lower, int) or not 0 <= lower < degree <= 36:
            raise ValueError('The requested zero order is incompatible with the infinity degree and endpoint.')
        polynomial = ((1 << lower) ^ (1 << degree)) if shape == 'binomial' else one_plus_power(degree - lower) << lower
        assert polynomial.bit_count() % 2 == 0
        result.append({'vectors': (a['left'], b['left'], a['coordinate_vector'], b['coordinate_vector']),
                       'polynomial': polynomial, 'coefficient': coefficients[i],
                       'component_pair': (i, j), 'lowest_polynomial_degree': lower,
                       'highest_polynomial_degree': degree,
                       'minimum_entry_order_infinity': -10,
                       'minimum_entry_order_zero': zero_orders[i]})
    return result
