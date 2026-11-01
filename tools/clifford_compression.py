"""Exact low-order column jets of constant elementary-tensor restrictions.

All scalar arithmetic is over F2. Matrices are tuples of packed columns;
Gram coefficient matrices are lists of packed rows. The routines keep
the actual coefficients, including cancellations between monomials.
"""

from functools import lru_cache


def bit_rank(vectors):
    pivots = {}
    for value in vectors:
        while value:
            pivot = value.bit_length()-1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def apply_columns(matrix, vector):
    answer = 0
    while vector:
        index = (vector & -vector).bit_length()-1
        answer ^= matrix[index]
        vector &= vector-1
    return answer


def multiply_columns(left, right):
    return tuple(apply_columns(left, column) for column in right)


def add_columns(left, right):
    return tuple(a ^ b for a, b in zip(left, right))


@lru_cache(maxsize=None)
def clifford_generators(length):
    """The single-site E01,E10 operators, in that order at each site."""
    if not isinstance(length, int) or length < 1:
        raise ValueError('The site count must be a positive integer.')
    size = 1 << length
    generators = []
    for site in range(length):
        stride = 1 << site
        generators.append(tuple(1 << (index ^ stride) if index & stride else 0
                                for index in range(size)))
        generators.append(tuple(0 if index & stride else 1 << (index ^ stride)
                                for index in range(size)))
    return tuple(generators)


@lru_cache(maxsize=None)
def clifford_monomials(length):
    """The ordered squarefree PBW basis, indexed by a 2*length bit mask."""
    generators = clifford_generators(length)
    identity = tuple(1 << index for index in range(1 << length))
    monomials = [identity]
    for mask in range(1, 1 << (2*length)):
        index = (mask & -mask).bit_length()-1
        monomials.append(multiply_columns(generators[index], monomials[mask ^ (1 << index)]))
    return tuple(monomials)


def validate_retraction(length, embedding, retraction):
    size, dimension = 1 << length, len(embedding)
    if not 1 <= dimension <= size or len(retraction) != dimension:
        raise ValueError('The embedding and retraction dimensions must agree.')
    if any(not isinstance(value, int) or not 0 <= value < (1 << size)
           for value in tuple(embedding)+tuple(retraction)):
        raise ValueError('The maps must fit the declared ambient vector space.')
    for i, covector in enumerate(retraction):
        for j, vector in enumerate(embedding):
            if (covector & vector).bit_count() % 2 != int(i == j):
                raise ValueError('The maps must satisfy LP=I.')


def compressed_operator(matrix, embedding, retraction):
    dimension = len(embedding)
    answer = 0
    for column, vector in enumerate(embedding):
        image = apply_columns(matrix, vector)
        for row, covector in enumerate(retraction):
            if (image & covector).bit_count() % 2:
                answer |= 1 << (dimension*row+column)
    return answer


def compressed_clifford_basis(length, embedding, retraction):
    validate_retraction(length, embedding, retraction)
    return tuple(compressed_operator(matrix, embedding, retraction)
                 for matrix in clifford_monomials(length))


def paired_mask(mask, length):
    answer = 0
    for site in range(length):
        part = (mask >> (2*site)) & 3
        answer |= ({0: 0, 1: 2, 2: 1, 3: 3}[part] << (2*site))
    return answer


def column_jet_coefficients(length, embedding, retraction, precision):
    """Return K(z) mod z^precision using its exact Clifford Gram factorization."""
    if not isinstance(precision, int) or precision < 1:
        raise ValueError('The jet precision must be positive.')
    compressed = compressed_clifford_basis(length, embedding, retraction)
    dimension = len(embedding)**2
    coefficients = [[0]*dimension for _ in range(precision)]
    expansions = {}

    def expand(mask, cap):
        key = mask, cap
        if key in expansions:
            return expansions[key]
        inactive = [site for site in range(length) if not (mask >> (2*site)) & 3]
        output = [0]*(cap+1)
        for choices in range(1 << len(inactive)):
            degree = choices.bit_count()
            if degree > cap:
                continue
            extended = mask
            for index, site in enumerate(inactive):
                if choices >> index & 1:
                    extended |= 3 << (2*site)
            output[degree] ^= compressed[extended]
        expansions[key] = output
        return output

    for mask in range(1 << (2*length)):
        weight = mask.bit_count()
        if weight >= precision:
            continue
        cap = precision-weight-1
        left, right = expand(mask, cap), expand(paired_mask(mask, length), cap)
        for a, vector in enumerate(left):
            for b, covector in enumerate(right[:cap-a+1]):
                degree = weight+a+b
                remaining = vector
                while remaining:
                    row = (remaining & -remaining).bit_length()-1
                    coefficients[degree][row] ^= covector
                    remaining &= remaining-1
    return coefficients


def _polynomial_product(left, right, mask):
    answer = 0
    while right:
        if right & 1:
            answer ^= left
        left <<= 1
        right >>= 1
    return answer & mask


def _polynomial_power(base, exponent, precision):
    mask, answer = (1 << precision)-1, 1
    while exponent:
        if exponent & 1:
            answer = _polynomial_product(answer, base, mask)
        base = _polynomial_product(base, base, mask)
        exponent >>= 1
    return answer


def row_jet_coefficients(length, embedding, retraction, precision):
    """Return the degree-length row reversal mod z^precision.

The regular frame uses z/(1+z); all rational units are expanded
exactly to the declared precision. Separate left and right physical
maps are retained, including nonorthogonal retractions.
"""
    if not isinstance(precision, int) or precision < 1:
        raise ValueError('The jet precision must be positive.')
    validate_retraction(length, embedding, retraction)
    monomials = clifford_monomials(length)
    left_basis, right_basis = [], []
    for operator in monomials:
        tilted = tuple(reversed(operator))  # Right multiplication by J^tensor length.
        left_basis.append(compressed_operator(tilted, retraction, retraction))
        right_basis.append(compressed_operator(tilted, embedding, embedding))
    dimension, modulus = len(embedding)**2, (1 << precision)-1
    coefficients = [[0]*dimension for _ in range(precision)]
    ratio = modulus ^ 1  # z/(1+z) modulo z^precision.
    ratio_powers = [_polynomial_power(ratio, exponent, precision)
                    for exponent in range(length+1)]

    def expand(mask, cap, basis):
        inactive = [site for site in range(length) if not (mask >> (2*site)) & 3]
        output = [0]*(cap+1)
        for choices in range(1 << len(inactive)):
            extended = mask
            for index, site in enumerate(inactive):
                if choices >> index & 1:
                    extended |= 3 << (2*site)
            polynomial = ratio_powers[choices.bit_count()]
            for degree in range(cap+1):
                if polynomial >> degree & 1:
                    output[degree] ^= basis[extended]
        return output

    for mask in range(1 << (2*length)):
        weight = mask.bit_count()
        if weight >= precision:
            continue
        cap = precision-weight-1
        left, right = expand(mask, cap, left_basis), expand(mask, cap, right_basis)
        unit_exponent = length-weight
        unit_base = 3 if unit_exponent >= 0 else modulus
        unit = _polynomial_power(unit_base, abs(unit_exponent), precision)
        for a, vector in enumerate(left):
            for b, covector in enumerate(right[:cap-a+1]):
                for c in range(cap-a-b+1):
                    if not (unit >> c & 1):
                        continue
                    degree = weight+a+b+c
                    remaining = vector
                    while remaining:
                        row = (remaining & -remaining).bit_length()-1
                        coefficients[degree][row] ^= covector
                        remaining &= remaining-1
    return coefficients


def jet_order_multiplicities(coefficients):
    """Resolve all Smith multiplicities strictly below the supplied precision."""
    size = len(coefficients[0])
    if any(len(matrix) != size for matrix in coefficients):
        raise ValueError('All coefficient matrices must have the same row size.')
    ranks = [0]
    for precision in range(1, len(coefficients)+1):
        rows = []
        for block in range(precision):
            for row in range(size):
                packed = 0
                for degree in range(block+1):
                    packed ^= coefficients[degree][row] << (size*(block-degree))
                rows.append(packed)
        ranks.append(bit_rank(rows))
    counts = [ranks[1]]
    counts.extend(ranks[k+1]-2*ranks[k]+ranks[k-1] for k in range(1, len(coefficients)))
    return tuple(counts), tuple(ranks)
