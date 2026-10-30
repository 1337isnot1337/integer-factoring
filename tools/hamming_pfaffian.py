"""A linear-degree crossing polynomial from the folded characteristic-two matrix.

This changes the polynomial entries, not the grading of the source operators.
At z=1 the matrix agrees with the source matrix on the equal-sum pair labels.
See research/hamming-pfaffian.md for the construction and full-span cost proof.
"""

from itertools import product

from modular_pfaffian import one_plus_z_power, polynomial_pfaffian, polynomial_product


def labels(width):
    assert width >= 1
    h = 1 << (width+2)
    return tuple((2*x+1, h-2*x-1) for x in range(1 << width))


def folded_matrix(x, v, u, y, weights=None):
    def degree(n):
        return n.bit_count() if weights is None else sum(
            weight for j, weight in enumerate(weights) if (n >> j) & 1)
    def oneplus(n):
        if weights is None:
            return one_plus_z_power(n.bit_count())
        answer = 1
        for j, weight in enumerate(weights):
            if (n >> j) & 1:
                answer = polynomial_product(answer, 1 | (1 << weight))
        return answer
    matrix = [[0]*8 for _ in range(8)]
    def put(i, j, value):
        matrix[i][j] = matrix[j][i] = value
    for i, t in enumerate(u):
        for j, n in enumerate(x):
            put(i, 2+j, 3 << degree(n | t))
        for j, m in enumerate(v):
            if ((m-1) & ~(t-1)) == 0:
                put(i, 4+j, oneplus(m))
        for j, other in enumerate(y):
            if t == other:
                put(i, 6+j, 1 << degree(t))
    for i, n in enumerate(x):
        for j, m in enumerate(v):
            if ((n-1) & ~(m-1)) == 0:
                put(2+i, 4+j, 1 << degree(n))
        for j, t in enumerate(y):
            put(2+i, 6+j, 3 << (degree(n & t)-1))
    for i, m in enumerate(v):
        for j, t in enumerate(y):
            if (m & t) == 1:
                put(4+i, 6+j, oneplus(m) << degree(t))
    return matrix


def profile(x, v, u, y, weights=None):
    return polynomial_pfaffian(folded_matrix(x, v, u, y, weights))


def matrices(width, weights=None):
    if weights is not None:
        assert len(weights) == width+2 and weights[0] == 1 and all(w >= 1 for w in weights)
    bound = 4*(width+3) if weights is None else 4*(2+sum(weights[1:]))
    domain = labels(width)
    dimension = len(domain)
    result = {}
    for x, y, u, v in product(range(dimension), repeat=4):
        polynomial = profile(domain[x], domain[v], domain[u], domain[y], weights)
        assert polynomial.bit_length() <= bound+1
        assert polynomial.bit_count() % 2 == int(x == v and y == u)
        while polynomial:
            bit = polynomial & -polynomial
            degree = bit.bit_length()-1
            result.setdefault(degree, [0]*dimension**2)[x*dimension+y] |= (
                1 << (u*dimension+v))
            polynomial ^= bit
    return result


def full_span_candidates(width, weights=None):
    """Uniform 2^O(width) full-span enumeration; very large constants."""
    from check_modular_recombination import candidates, insert
    original = list(candidates(matrices(width, weights), 1 << width, False).items())
    basis = {}
    for i, (vector, description) in enumerate(original):
        insert(basis, vector, 1 << i)
    generators = list(basis.values())
    vector, description_mask = 0, 0
    for index in range(1, 1 << len(generators)):
        bit = (index & -index).bit_length()-1
        other, other_mask = generators[bit]
        vector ^= other
        description_mask ^= other_mask
        yield vector, tuple(term for i, (original_vector, terms) in enumerate(original)
                            if (description_mask >> i) & 1 for term in terms)
