"""A finite relabeling seed and its square-zero four-piece block model."""

from functools import lru_cache
from itertools import permutations, product

from check_hamming_joint_lower import source_basis
from check_modular_recombination import decode, encode, insert, represent


def relabel(vector, dimension, permutation):
    rows = decode(vector, dimension)
    result = [0] * dimension**2
    for x, y, u, v in product(range(dimension), repeat=4):
        value = (
            rows[permutation[x] * dimension + permutation[y]]
            >> (permutation[u] * dimension + permutation[v])
        ) & 1
        result[x * dimension + y] |= value << (u * dimension + v)
    return encode(result, dimension)


@lru_cache(maxsize=1)
def seed_data():
    original = source_basis()
    images = {}
    for permutation in permutations(range(4)):
        for j, vector in enumerate(original):
            images.setdefault(relabel(vector, 4, permutation), (permutation, j))
    vectors = list(images)
    basis = {}
    for j, vector in enumerate(vectors):
        insert(basis, vector, 1 << j)
    assert len(vectors) == 196 and len(basis) == 102
    xor_rows = [0] * 16
    for x, y, u, v in product(range(4), repeat=4):
        xor_rows[4 * x + y] |= int(x ^ y ^ u ^ v == 0) << (4 * u + v)
    xor_tensor = encode(xor_rows, 4)
    xor_mask = represent(basis, xor_tensor)
    assert xor_mask is not None
    pieces = [
        vectors[0],
        vectors[1] ^ vectors[0],
        vectors[19],
        vectors[49] ^ vectors[9],
        vectors[69] ^ vectors[9],
        xor_tensor,
    ]
    masks = [1, 3, 1 << 19, (1 << 49) ^ (1 << 9), (1 << 69) ^ (1 << 9), xor_mask]
    return original, tuple(images.items()), pieces, masks


def action(matrix, vector):
    return sum(((row & vector).bit_count() % 2) << j for j, row in enumerate(matrix))


def compose(left, right):
    result = []
    for row in left:
        value = 0
        while row:
            bit = row & -row
            value ^= right[bit.bit_length() - 1]
            row ^= bit
        result.append(value)
    return result


def normalized(vector):
    return [
        sum(
            ((row >> (4 * y + x)) & 1) << (4 * x + y)
            for x, y in product(range(4), repeat=2)
        )
        for row in decode(vector, 4)
    ]


def block_basis():
    """Build a constant F2 basis for R,C,A,X, where A is all-distinct."""
    _, _, pieces, _ = seed_data()
    identity = [1 << i for i in range(16)]
    r, c, xor_tensor = map(normalized, (pieces[0], pieces[1], pieces[5]))
    n = [a ^ b for a, b in zip(c, identity)]
    l = [a ^ b for a, b in zip(n, r)]
    image_basis, domain = {}, []
    for j in range(16):
        column = action(l, 1 << j)
        if insert(image_basis, column, 1 << len(domain)):
            domain.append(1 << j)
    assert len(domain) == 7
    small = []
    for matrix in (r, xor_tensor):
        columns = [represent(image_basis, action(matrix, v)) for v in domain]
        assert all(v is not None for v in columns)
        small.append(
            [sum(((columns[j] >> i) & 1) << j for j in range(7)) for i in range(7)]
        )
    p, q = small
    assert compose(p, p) == p and compose(p, q) == compose(q, p) == p
    assert compose(q, q) == p
    special = next(action(p, 1 << j) for j in range(7) if action(p, 1 << j))
    chains, independent = [], {}
    for j in range(7):
        v = (1 << j) ^ action(p, 1 << j)
        u = action(q, v)
        if insert(independent, u, 0):
            chains.append((u, v))
    assert len(chains) == 3

    def lift(coordinate):
        return sum(v for j, v in enumerate(domain) if coordinate >> j & 1)

    s = lift(special)
    basis = [action(l, s), s]
    for u, v in chains:
        u, v = lift(u), lift(v)
        basis += [action(l, u), action(l, v), u, v]
    spanning = {}
    for j, vector in enumerate(basis):
        assert insert(spanning, vector, 1 << j)
    for vector in range(1, 1 << 16):
        if action(n, vector) == action(r, vector) == action(xor_tensor, vector) == 0:
            if insert(spanning, vector, 1 << len(basis)):
                basis.append(vector)
                if len(basis) == 16:
                    break
    assert len(basis) == 16
    return basis
