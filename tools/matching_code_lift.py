#!/usr/bin/env python3
"""A multigraded characteristic-two crossing on uniform set partitions.

The original one-variable folded Wick formula is due to the parallel PR22
work, revision c12da897fcc6081685e8c663217ce9083db87b4a. This module gives its
division-free bit-coordinate lift and exact permutation-orbit rank tables.
See research/matching-code-lift.md. No asymptotic rank bound is asserted.
"""

from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations

from check_modular_crossing import column_bits, rank_two
from fermionic_piece_polynomial import polynomial_product


def uniform_partitions(vertices, block_size):
    """Unordered set partitions, encoded as sorted tuples of block bit masks."""
    if not vertices:
        yield ()
        return
    if len(vertices) % block_size:
        raise ValueError("block size must divide the number of vertices")
    first = vertices[0]
    for tail in combinations(vertices[1:], block_size - 1):
        block = (first,) + tail
        mask = sum(1 << i for i in block)
        rest = tuple(i for i in vertices if i not in block)
        for partition in uniform_partitions(rest, block_size):
            yield tuple(sorted((mask,) + partition))


def labels(vertices, block_size=2):
    return tuple(uniform_partitions(tuple(range(vertices)), block_size))


def modes(label):
    return tuple(1 + 2 * block for block in label)


def decode(vertices, exponent):
    """Decode L first, then radix-five high coordinates.

    The L radix 2*vertices+1 is valid for every block_size >= 2.
    """
    low, exponent = exponent % (2 * vertices + 1), exponent // (2 * vertices + 1)
    high = []
    for _ in range(vertices):
        high.append(exponent % 5)
        exponent //= 5
    if exponent:
        raise ValueError("exponent outside the proved degree box")
    return low, tuple(high)


def edge_polynomial(vertices, low, ones, binomials=0):
    exponent = low + (2 * vertices + 1) * sum(
        5 ** i for i in range(vertices) if (ones >> i) & 1)
    result = 1 << exponent
    for i in range(vertices):
        if (binomials >> i) & 1:
            result ^= result << ((2 * vertices + 1) * 5 ** i)
    return result


def lifted_polynomial(vertices, x, v, u, y):
    """Pfaffian of the six lifted blocks, ordered U,A,B,Y.

    All four arguments are uniform set partitions of the same vertex set.
    Monomial coefficients are packed into an integer using decode's radices.
    The proved degree box prevents carries between distinct coordinates.
    """
    groups = ([('U', e) for e in u] + [('A', e) for e in x]
              + [('B', e) for e in v] + [('Y', e) for e in y])
    count = len(groups)
    edges = {}
    for i, (first, e) in enumerate(groups):
        for j in range(i + 1, count):
            second, f = groups[j]
            pair = first, second
            value = 0
            if pair == ('U', 'A'):
                value = edge_polynomial(vertices, 1, e | f)
                value ^= value << 1
            elif pair == ('U', 'B') and e == f:
                value = edge_polynomial(vertices, 0, 0, e)
                value ^= value << 1
            elif pair in (('U', 'Y'), ('A', 'B')) and e == f:
                value = edge_polynomial(vertices, 1, e)
            elif pair == ('A', 'Y'):
                value = edge_polynomial(vertices, 0, e & f)
                value ^= value << 1
            elif pair == ('B', 'Y') and not (e & f):
                value = edge_polynomial(vertices, 1, f, e)
                value ^= value << 1
            if value:
                edges[i, j] = value

    @lru_cache(None)
    def pfaffian(remaining):
        if not remaining:
            return 1
        first_bit = remaining & -remaining
        first = first_bit.bit_length() - 1
        rest = remaining ^ first_bit
        value = 0
        candidates = rest
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            other = bit.bit_length() - 1
            edge = edges.get((first, other), 0)
            if edge:
                value ^= polynomial_product(edge, pfaffian(rest ^ bit))
        return value

    return pfaffian((1 << count) - 1)


def specialize(vertices, polynomial):
    """L=z, Z_i=z^(2^(i+1)); collisions are added in F_2."""
    answer = 0
    while polynomial:
        bit = polynomial & -polynomial
        polynomial ^= bit
        low, high = decode(vertices, bit.bit_length() - 1)
        exponent = low + sum(digit << (i + 1) for i, digit in enumerate(high))
        answer ^= 1 << exponent
    return answer


def permutation_action(labs, permutation):
    index = {label: i for i, label in enumerate(labs)}
    action = []
    for label in labs:
        image = tuple(sorted(sum(1 << permutation[i] for i in range(len(permutation))
                                 if (block >> i) & 1) for block in label))
        action.append(index[image])
    return action


def orbit_tables(vertices, block_size=2, verify_entry=None):
    """Exact matrices for one canonical multidegree per permutation orbit.

    Full coefficient matrices are unnecessary: permuting the bit positions
    permutes both flattenings' rows and columns. This is a finite optimizer,
    not the uniform exp(O(vertices)) entry-generation algorithm.
    """
    labs = labels(vertices, block_size)
    dimension = len(labs)
    perms = tuple(permutations(range(vertices)))
    actions = [permutation_action(labs, permutation) for permutation in perms]
    generators = []
    for i in range(vertices - 1):
        permutation = list(range(vertices))
        permutation[i], permutation[i + 1] = permutation[i + 1], permutation[i]
        generators.append(actions[perms.index(tuple(permutation))])

    patterns = tuple(combinations_with_replacement(range(5), vertices))
    charge_actions = {}
    powers = [5 ** i for i in range(vertices)]
    for ci, pattern in enumerate(patterns):
        for pi, permutation in enumerate(perms):
            original = sum(pattern[permutation[i]] * powers[i] for i in range(vertices))
            if original not in charge_actions:
                charge_actions[original] = [ci, []]
            charge_actions[original][1].append(pi)

    seen = bytearray(dimension ** 4)
    representatives = []
    for code in range(dimension ** 4):
        if seen[code]:
            continue
        number = code
        quartet = []
        for _ in range(4):
            quartet.append(number % dimension)
            number //= dimension
        quartet = tuple(quartet)  # x,y,v,u
        representatives.append(quartet)
        pending = [quartet]
        seen[code] = 1
        while pending:
            current = pending.pop()
            for action in generators:
                target = tuple(action[q] for q in current)
                target_code = sum(target[i] * dimension ** i for i in range(4))
                if not seen[target_code]:
                    seen[target_code] = 1
                    pending.append(target)

    matrices = {}
    maximum_monomials = 0
    for xi, yi, vi, ui in representatives:
        x, y, v, u = (labs[i] for i in (xi, yi, vi, ui))
        polynomial = lifted_polynomial(vertices, x, v, u, y)
        assert polynomial.bit_count() % 2 == int(x == v and y == u)
        if verify_entry is not None:
            verify_entry(vertices, x, v, u, y, polynomial)
        maximum_monomials = max(maximum_monomials, polynomial.bit_count())
        while polynomial:
            bit = polynomial & -polynomial
            polynomial ^= bit
            exponent = bit.bit_length() - 1
            low, charge = exponent % (2 * vertices + 1), exponent // (2 * vertices + 1)
            ci, compatible = charge_actions[charge]
            matrix = matrices.setdefault((low, ci), [0] * dimension ** 2)
            for pi in compatible:
                action = actions[pi]
                row = action[xi] * dimension + action[yi]
                column = action[vi] * dimension + action[ui]
                # Equivalent entries are copied, not added as new summands.
                matrix[row] |= 1 << column
        polynomial_product.cache_clear()

    profile = [(low, patterns[ci], rank_two(matrix), rank_two(column_bits(matrix, dimension)))
               for (low, ci), matrix in sorted(matrices.items())]
    maximum = max(row * column for _, _, row, column in profile)
    return {
        'vertices': vertices, 'block_size': block_size,
        'dimension': dimension, 'weight': sum(modes(labs[0])),
        'quartet_orbits': len(representatives),
        'degree_orbits': len(profile), 'maximum_monomials': maximum_monomials,
        'maximum_rank_product': maximum,
        'worst': [entry for entry in profile if entry[2] * entry[3] == maximum],
        'profile': profile,
    }, {(low, patterns[ci]): matrix for (low, ci), matrix in matrices.items()}
