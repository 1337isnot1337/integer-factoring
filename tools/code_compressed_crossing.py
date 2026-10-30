"""Binary-code dual compression of the elementary crossing pencil over F2."""

from itertools import combinations


def codewords(generators):
    values = [0]
    for row in generators:
        values += [value ^ row for value in values]
    assert len(set(values)) == len(values)
    return tuple(values)


def oneplus_power(exponent):
    answer, bit = 1, 0
    while exponent:
        if exponent & 1:
            answer ^= answer << (1 << bit)
        exponent >>= 1
        bit += 1
    return answer


def rref_subspaces(ambient_width, dimension):
    """Every binary linear subspace exactly once, using least-bit pivots."""
    for pivots in combinations(range(ambient_width), dimension):
        free = tuple((i, j) for i, pivot in enumerate(pivots)
                     for j in range(pivot+1, ambient_width) if j not in pivots)
        for bits in range(1 << len(free)):
            rows = [1 << pivot for pivot in pivots]
            for bit, (i, j) in enumerate(free):
                if bits >> bit & 1:
                    rows[i] |= 1 << j
            yield tuple(rows), pivots


def representatives(ambient_width, pivots):
    remaining = [j for j in range(ambient_width) if j not in pivots]
    return tuple(sum(((label >> i) & 1) << j for i, j in enumerate(remaining))
                 for label in range(1 << len(remaining)))


def entry(ambient_width, code, code_dimension, x, y, u, v):
    """Coefficient bit polynomial after dividing the common z^code_dimension."""
    if x ^ y ^ u ^ v:
        return 0
    difference = u ^ v
    complement = ((1 << ambient_width)-1) ^ difference
    target = (x ^ v) & difference
    result = 0
    for shift in code:
        if shift & difference == target:
            exponent = ((x ^ u ^ shift) & complement).bit_count()
            result ^= oneplus_power(exponent)
    result <<= difference.bit_count()
    assert result & ((1 << code_dimension)-1) == 0
    return result >> code_dimension


def literal_pencil_entry(ambient_width, x, y, u, v):
    """Independent product of the original one-bit R+z(E+R) entries."""
    result = 1
    for j in range(ambient_width):
        a, b, c, d = ((word >> j) & 1 for word in (x, y, u, v))
        r = int(a == b and c == d)
        e = int(a == d and b == c)
        factor = r | ((e ^ r) << 1)
        if not factor:
            return 0
        if factor == 2:
            result <<= 1
        elif factor == 3:
            result ^= result << 1
    return result


def matrices(ambient_width, generators, pivots):
    code = codewords(generators)
    reps = representatives(ambient_width, pivots)
    dimension = len(reps)
    width = ambient_width - len(generators)
    rows = [[0] * dimension**2 for _ in range(width+1)]
    cols = [[0] * dimension**2 for _ in range(width+1)]
    for xi, x in enumerate(reps):
        for yi, y in enumerate(reps):
            for ui, u in enumerate(reps):
                # Only the additive-conservation column can be nonzero.
                vi = xi ^ yi ^ ui
                v = reps[vi]
                polynomial = entry(ambient_width, code, len(generators), x,y,u,v)
                assert polynomial.bit_length() <= width+1
                assert polynomial.bit_count() % 2 == int(x == v and y == u)
                while polynomial:
                    bit = polynomial & -polynomial
                    degree = bit.bit_length()-1
                    rows[degree][xi*dimension+yi] |= 1 << (ui*dimension+vi)
                    cols[degree][xi*dimension+ui] |= 1 << (yi*dimension+vi)
                    polynomial ^= bit
    return rows, cols
