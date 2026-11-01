"""Exact homogeneous physical retractions of the two-site reflection seed.

The endpoint gauges are inherited from the four-letter reflection tensor.
Embedding and retraction maps mix only equal-height physical words. All
scalar arithmetic is over F2; polynomials are packed coefficient masks.
"""

from functools import lru_cache
from itertools import product

from check_block_smith_crossings import smith_profile
from check_cancellation_module_automata import pmul
from check_triangular_auxiliary_row_fiber import HEIGHT, LABELS, reflection_seed


WORD_HEIGHTS = tuple(HEIGHT[word & 3] + HEIGHT[word >> 2] for word in range(16))
MIDDLE_FIBER = (3, 6, 9, 12)


@lru_cache(maxsize=1)
def two_site_endpoint_matrix():
    """The two ambient gauged endpoint matrices coincide after restoration.

    At zero the common scalar z^8 is removed. The resulting total letter
    energy is four times the two-site word height.
    """
    values = dict(zip(LABELS, reflection_seed()))
    local = [[0] * 16 for _ in range(16)]
    for x, y, u, v in product(range(4), repeat=4):
        value = values[x, y, u, v]
        zero = value << (HEIGHT[x] + HEIGHT[y])
        assert zero & 15 == 0
        zero >>= 4
        # Infinity uses the opposite cut. The physical labels y,u are
        # exchanged before comparing its regular coefficient matrix.
        opposite = values[x, u, y, v]
        height = HEIGHT[y] + HEIGHT[v] + 4
        infinity = 0
        for degree in range(opposite.bit_length()):
            if opposite >> degree & 1:
                assert height >= degree
                infinity ^= 1 << (height - degree)
        assert zero == infinity
        local[4 * x + y][4 * u + v] = zero
    matrix = []
    for x, y in product(range(16), repeat=2):
        row = []
        for u, v in product(range(16), repeat=2):
            value = pmul(local[4 * (x & 3) + (y & 3)][4 * (u & 3) + (v & 3)],
                         local[4 * (x >> 2) + (y >> 2)][4 * (u >> 2) + (v >> 2)])
            assert not value or (WORD_HEIGHTS[x] + WORD_HEIGHTS[y]
                                 + WORD_HEIGHTS[u] + WORD_HEIGHTS[v]) % 2 == 0
            assert value.bit_count() % 2 == int(x == v and y == u)
            row.append(value)
        matrix.append(tuple(row))
    return tuple(matrix)


def middle_fiber_retraction(vector, functional):
    """Keep every other height fiber and a rank-three part of the middle one."""
    if (not isinstance(vector, int) or not isinstance(functional, int)
            or not 1 <= vector < 16 or not 1 <= functional < 16
            or (vector & functional).bit_count() % 2 != 1):
        raise ValueError("Use nonzero four-bit masks with functional(vector)=1.")
    pivot = (functional & -functional).bit_length() - 1
    coordinates = tuple(index for index in range(16) if index != MIDDLE_FIBER[pivot])
    embedding, retraction = [], []
    for index in coordinates:
        if index not in MIDDLE_FIBER:
            embedding.append(1 << index)
            retraction.append(1 << index)
            continue
        local = MIDDLE_FIBER.index(index)
        embedding.append((1 << index) ^ (((functional >> local) & 1) << MIDDLE_FIBER[pivot]))
        covector = 1 << index
        if vector >> local & 1:
            for j, ambient in enumerate(MIDDLE_FIBER):
                if functional >> j & 1:
                    covector ^= 1 << ambient
        retraction.append(covector)
    return tuple(embedding), tuple(retraction), tuple(WORD_HEIGHTS[index] for index in coordinates)


def _support(mask):
    return tuple(index for index in range(16) if mask >> index & 1)


def restricted_endpoint_profiles(embedding, retraction, heights, *, keep_blocks=False):
    """Compute both exact profiles for a homogeneous constant retraction."""
    dimension = len(embedding)
    if (not 1 <= dimension <= 16 or len(retraction) != dimension
            or len(heights) != dimension):
        raise ValueError("Embedding, retraction and heights must have equal valid dimensions.")
    if any(not isinstance(mask, int) or not 0 <= mask < (1 << 16)
           for mask in tuple(embedding) + tuple(retraction)):
        raise ValueError("Each vector and covector must fit the sixteen-letter space.")
    if any((retraction[i] & embedding[j]).bit_count() % 2 != int(i == j)
           for i, j in product(range(dimension), repeat=2)):
        raise ValueError("The maps must satisfy LP=I.")
    inputs, outputs = tuple(map(_support, embedding)), tuple(map(_support, retraction))
    for height, first, second in zip(heights, inputs, outputs):
        if any(WORD_HEIGHTS[index] != height for index in first + second):
            raise ValueError("The maps must preserve each physical letter height.")
    groups = {}
    for a, b in product(range(dimension), repeat=2):
        groups.setdefault((heights[a] + heights[b]) % 2, []).append(dimension * a + b)

    def tensor_support(first, second):
        return tuple(16 * a + b for a in first for b in second)

    left_row = [tensor_support(a, b) for a, b in product(outputs, repeat=2)]
    right_row = [tensor_support(a, b) for a, b in product(inputs, repeat=2)]
    column = [tensor_support(a, b) for a in outputs for b in inputs]
    source = two_site_endpoint_matrix()
    profiles, saved_blocks = [], []
    for left, right in ((left_row, right_row), (column, column)):
        intermediate = []
        for row in source:
            changed = []
            for support in right:
                value = 0
                for index in support:
                    value ^= row[index]
                changed.append(value)
            intermediate.append(changed)
        orders, blocks = [], []
        for group in groups.values():
            block = []
            for i in group:
                row = []
                for j in group:
                    value = 0
                    for index in left[i]:
                        value ^= intermediate[index][j]
                    assert value.bit_count() % 2 == int(
                        i // dimension == j % dimension and i % dimension == j // dimension)
                    row.append(value)
                block.append(row)
            degree = max(value.bit_length() - 1 for row in block for value in row)
            resolved, precision = smith_profile(block, degree)
            orders.extend(resolved)
            if keep_blocks:
                blocks.append({"indices": group, "matrix": block,
                               "orders": resolved, "precision": precision})
        assert len(orders) == dimension**2
        profiles.append(tuple(sorted(orders)))
        saved_blocks.append(blocks)
    result = {"zero_orders": profiles[0], "infinity_orders": profiles[1],
              "letter_heights": tuple(heights), "type_energies": tuple(4 * height for height in heights)}
    if keep_blocks:
        result["blocks"] = saved_blocks
    return result
