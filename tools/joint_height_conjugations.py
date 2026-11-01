"""Constant joint physical conjugations of the four-letter reflection seed.

Matrices act on the pair of physical row legs, preserve total height,
and commute with swap. The complete binary family has 1536 members.
"""

from functools import lru_cache
from itertools import product

from check_block_smith_crossings import binary_rank, smith_profile
from check_triangular_auxiliary_row_fiber import HEIGHT, LABELS, reflection_seed
from matroid_smith_channels import binary_inverse


def xor_sum(values):
    result = 0
    for value in values:
        result ^= value
    return result


def support(mask):
    return tuple(i for i in range(16) if mask >> i & 1)


@lru_cache(maxsize=1)
def source_matrix():
    matrix = [[0] * 16 for _ in range(16)]
    for (x, y, u, v), value in zip(LABELS, reflection_seed()):
        matrix[4 * x + y][4 * u + v] = value
    return tuple(map(tuple, matrix))


def middle_centralizers():
    """All invertible matrices commuting with two binary swap blocks."""
    result = []
    for mask in range(256):
        rows = [0] * 4
        for i, j in product(range(2), repeat=2):
            first = (mask >> (2 * (2 * i + j))) & 1
            second = (mask >> (2 * (2 * i + j) + 1)) & 1
            rows[2 * i] |= (first << (2 * j)) | (second << (2 * j + 1))
            rows[2 * i + 1] |= (second << (2 * j)) | (first << (2 * j + 1))
        if binary_rank(rows) == 4:
            result.append(tuple(rows))
    assert len(result) == 96
    return tuple(result)


def joint_height_centralizers():
    """Yield all 1536 invertible binary height-preserving swap centralizers."""
    middle = (3, 12, 6, 9)
    pairs = ((2, 8), (1, 4), (11, 14), (7, 13))
    for index, block in enumerate(middle_centralizers()):
        for flags in range(16):
            rows = [1 << i for i in range(16)]
            for i, ambient in enumerate(middle):
                rows[ambient] = xor_sum(1 << middle[j] for j in range(4)
                                       if block[i] >> j & 1)
            for j, (a, b) in enumerate(pairs):
                if flags >> j & 1:
                    rows[a], rows[b] = rows[b], rows[a]
            rows = tuple(rows)
            for i, row in enumerate(rows):
                assert all(HEIGHT[i // 4] + HEIGHT[i % 4]
                           == HEIGHT[k // 4] + HEIGHT[k % 4] for k in support(row))
                assert xor_sum(1 << (4 * (k % 4) + k // 4) for k in support(row)) \
                    == rows[4 * (i % 4) + i // 4]
            yield index, flags, rows, binary_inverse(rows)


def conjugate(matrix, left, inverse):
    temporary = [[xor_sum(matrix[k][j] for k in support(mask)) for j in range(16)]
                 for mask in left]
    columns = [tuple(k for k in range(16) if inverse[k] >> j & 1) for j in range(16)]
    return tuple(tuple(xor_sum(temporary[i][k] for k in columns[j])
                       for j in range(16)) for i in range(16))


def restored_endpoints(matrix):
    """Return the two regular endpoint matrices with restored charges."""
    zero, infinity = [[0] * 16 for _ in range(16)], [[0] * 16 for _ in range(16)]
    for x, y, u, v in product(range(4), repeat=4):
        value = matrix[4 * x + y][4 * u + v]
        assert value.bit_count() % 2 == int(x == v and y == u)
        shifted = value << (HEIGHT[x] + HEIGHT[y])
        assert shifted & 15 == 0
        zero[4 * x + y][4 * u + v] = shifted >> 4
        height = HEIGHT[u] + HEIGHT[v] + 4
        transformed = 0
        for degree in range(value.bit_length()):
            if value >> degree & 1:
                assert height >= degree
                transformed ^= 1 << (height - degree)
        infinity[4 * x + u][4 * y + v] = transformed
    return tuple(map(tuple, zero)), tuple(map(tuple, infinity))


def exact_profile(matrix):
    degree = max(value.bit_length() - 1 for row in matrix for value in row)
    return smith_profile(matrix, degree)[0]


def leading_infinity_operator(matrix):
    """Undo the opposite-cut reshuffle on its constant coefficient."""
    return tuple(sum(((matrix[i][j] >> (HEIGHT[j // 4] + HEIGHT[j % 4] + 4)) & 1) << j
                     for j in range(16)) for i in range(16))


def binary_product(left, right):
    return tuple(xor_sum(right[k] for k in support(row)) for row in left)
