#!/usr/bin/env python3
"""Exact rank tools for the full coefficient tensor span of the CNOT path.

Coefficient indices use base-four digits a,b,c,d, least significant site first.
The row-rank formula is valid over every characteristic-two field.  The fast
implementation below uses F2 packed rows; small extension-field controls use
the separate, deliberately direct field implementation in verify().
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
import random
from pathlib import Path

# Original local row matrices, represented by packed binary rows.
C0 = [1, 2, 4, 8]
C1 = [0, 0, 12, 12]
C2 = [0, 12, 10, 6]
C3 = [0, 10, 0, 10]
ROW = [[x ^ y for x, y in zip(C0, C2)], C2, C1, C3]
# Independent constant column-basis change from the module certificate.
COLUMN = [[11, 12, 8, 8], [11, 8, 8, 8], [0, 12, 12, 0], [2, 2, 0, 0]]
COLUMN_DIAGONAL = [[0, 5, 3, 7], [0, 2, 4, 0], [1, 0, 0, 0], [0, 0, 0, 8]]
# In the u,b,c,d coefficient basis, where u=a+b.
COLUMN_ARROW = [[1, 0, 0, 0], [8, 0, 8, 3], [0, 2, 0, 0], [0, 0, 4, 0]]


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            p = row.bit_length() - 1
            if p in pivots:
                row ^= pivots[p]
            else:
                pivots[p] = row
                break
    return len(pivots)


def kron(A, B):
    m = len(B)
    result = []
    for a in A:
        for b in B:
            result.append(sum(b << (m * j) for j in range(len(A)) if a >> j & 1))
    return result


def digits(word, t):
    return [(word >> (2 * i)) & 3 for i in range(t)]


def tensor_matrix(coefficients, t, local=ROW):
    """Literal matrix, useful for small exact checks and column ranks."""
    n = 4 ** t
    result = [0] * n
    for word in coefficients:
        term = [1]
        for d in reversed(digits(word, t)):
            term = kron(term, local[d])
        result = [a ^ b for a, b in zip(result, term)]
    return result


def row_blocks(coefficients, t):
    """Return the exact nonzero direct-sum blocks, indexed by active sites.

    A block on S has multiplicity 2**(t-len(S)).  Within it, b,c,d act as
    I_2,E_01,E_10.  An a at an active site, or b/c/d at an inactive site,
    annihilates that block, so every coefficient word contributes to one S.
    """
    blocks = {}
    for word in coefficients:
        ds = digits(word, t)
        active = [i for i, d in enumerate(ds) if d]
        mask = sum(1 << i for i in active)
        size = 1 << len(active)
        rows = blocks.setdefault(mask, [0] * size)
        row_fixed = col_fixed = 0
        free = []
        for j, i in enumerate(active):
            if ds[i] == 1:
                free.append(j)
            elif ds[i] == 2:
                col_fixed |= 1 << j
            else:
                row_fixed |= 1 << j
        for assignments in range(1 << len(free)):
            common = sum(((assignments >> k) & 1) << j for k, j in enumerate(free))
            rows[row_fixed | common] ^= 1 << (col_fixed | common)
    return blocks


def row_rank(coefficients, t):
    return sum((1 << (t - mask.bit_count())) * rank(rows)
               for mask, rows in row_blocks(coefficients, t).items())


def column_character_bound(coefficients, t):
    """Rank of the simultaneous diagonal quotient, counted with multiplicity.

    The three characters are p=a+b, c, d, with multiplicities 2,1,1.
    This is a lower bound, not generally the full column rank.
    """
    evaluations = set()
    for word in coefficients:
        character = 0
        for i, d in enumerate(digits(word, t)):
            character |= (0 if d < 2 else d - 1) << (2 * i)
        if character in evaluations:
            evaluations.remove(character)
        else:
            evaluations.add(character)
    return sum(1 << sum(d == 0 for d in digits(character, t))
               for character in evaluations)


def ranks(coefficients, t):
    return row_rank(coefficients, t), column_rank(coefficients, t)


def column_rank(coefficients, t):
    # On W^tensor(t), the independent two-sided diagonalization is exact.
    if all(all(digits(word, t)) for word in coefficients):
        return sum(1 << digits(word, t).count(1) for word in coefficients)
    u_coefficients = swap_a_u_basis(coefficients, t)
    b_words = [word for word in u_coefficients if 1 in digits(word, t)]
    if not b_words:
        return len(u_coefficients)
    if len(u_coefficients) == 1:
        return 1 << digits(b_words[0], t).count(1)
    return rank(tensor_matrix(coefficients, t, COLUMN_DIAGONAL))


def swap_a_u_basis(coefficients, t):
    """The involutive coefficient substitution a=u+b, equivalently u=a+b."""
    result = set()
    for word in coefficients:
        replacements = coefficient_product([{0, 1} if d == 0 else {d} for d in digits(word, t)])
        result.symmetric_difference_update(replacements)
    return result


def low_column_rank(u_coefficients, t):
    """Exact 0/1/2 rank classification in the u,b,c,d basis; None means >2."""
    b_words = [word for word in u_coefficients if 1 in digits(word, t)]
    if not b_words:
        return len(u_coefficients) if len(u_coefficients) <= 2 else None
    if len(u_coefficients) == 1 and digits(b_words[0], t).count(1) == 1:
        return 2
    return None


def two_product_row_rank(left, right, t):
    """Row rank of nonzero multiples of two distinct u/c/d product words."""
    A, B = digits(left, t), digits(right, t)
    assert left != right and 1 not in A and 1 not in B
    ka, kb = A.count(0), B.count(0)
    if any({a, b} == {2, 3} for a, b in zip(A, B)):
        return 4 ** ka + 4 ** kb
    common = sum(a == b == 0 for a, b in zip(A, B))
    only_a, only_b = ka - common, kb - common
    if not only_a or not only_b:
        return 4 ** (common + max(only_a, only_b))
    return 4 ** common * (4 ** only_a + 4 ** only_b - 2)


def column_rank_two_pieces(t):
    """All nonzero F2 column-rank-two pieces, returned in the a,b,c,d basis."""
    words = [sum(d << (2 * i) for i, d in enumerate(ds))
             for ds in itertools.product([0, 2, 3], repeat=t)]
    for i, left in enumerate(words):
        for right in words[:i]:
            yield swap_a_u_basis({left, right}, t), two_product_row_rank(left, right, t)
    for site in range(t):
        for ds in itertools.product([0, 2, 3], repeat=t - 1):
            full = list(ds)
            full.insert(site, 1)
            word = sum(d << (2 * i) for i, d in enumerate(full))
            yield swap_a_u_basis({word}, t), 2 * 4 ** full.count(0)


def coefficient_product(local_choices):
    """Tensor product of local coefficient supports in the a,b,c,d basis."""
    words = {0}
    for i, choices in enumerate(local_choices):
        words = {word | (d << (2 * i)) for word in words for d in choices}
    return words


def cut(M):
    return [sum(((M[x * 2 + y] >> (u * 2 + v)) & 1) << (y * 2 + v)
                for v in range(2) for y in range(2))
            for x in range(2) for u in range(2)]


def binary_inverse(rows):
    n = len(rows)
    a = [row | (1 << (n + i)) for i, row in enumerate(rows)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i] >> j & 1)
        a[j], a[pivot] = a[pivot], a[j]
        for i in range(n):
            if i != j and a[i] >> j & 1:
                a[i] ^= a[j]
    return [row >> n for row in a]


def multiply(A, B):
    return [xor_all(B[j] for j in range(len(B)) if a >> j & 1) for a in A]


def xor_all(values):
    result = 0
    for value in values:
        result ^= value
    return result


def conjugate(matrix, columns):
    basis = [sum(((column >> i) & 1) << j for j, column in enumerate(columns))
             for i in range(len(columns))]
    return multiply(multiply(binary_inverse(basis), matrix), basis)


def field_multiply(a, b, modulus=7):
    result = 0
    degree = modulus.bit_length() - 1
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a >> degree:
            a ^= modulus
    return result


def field_rank(matrix, modulus=7):
    a = [row[:] for row in matrix]
    q = 1 << (modulus.bit_length() - 1)
    r = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = next(x for x in range(1, q) if field_multiply(a[r][j], x, modulus) == 1)
        a[r] = [field_multiply(x, inv, modulus) for x in a[r]]
        for i in range(r + 1, len(a)):
            if a[i][j]:
                c = a[i][j]
                a[i] = [x ^ field_multiply(c, y, modulus) for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def field_control(coefficients, t):
    """Independent dense F4 check of the block identity, not a fast engine."""
    n = 4 ** t
    literal = [[0] * n for _ in range(n)]
    block_matrices = {}
    for word, value in coefficients.items():
        tensor = tensor_matrix({word}, t)
        for i, row in enumerate(tensor):
            for j in range(n):
                if row >> j & 1:
                    literal[i][j] ^= value
        for mask, rows in row_blocks({word}, t).items():
            size = len(rows)
            block = block_matrices.setdefault(mask, [[0] * size for _ in range(size)])
            for i, row in enumerate(rows):
                for j in range(size):
                    if row >> j & 1:
                        block[i][j] ^= value
    direct = field_rank(literal)
    predicted = sum((1 << (t - mask.bit_count())) * field_rank(matrix)
                    for mask, matrix in block_matrices.items())
    assert direct == predicted
    return direct


def verify():
    expected_row = [[1, 2, 0, 0], [0, 0, 4, 8], [0, 0, 8, 0], [0, 0, 0, 4]]
    for original, normal, column in zip(ROW, expected_row, COLUMN):
        assert conjugate(original, [1, 14, 12, 10]) == normal
        assert conjugate(cut(original), [6, 9, 1, 2]) == column
    # Independent two-sided diagonalization of the all-active subspace.
    L, R = kron([3, 1], [3, 1]), kron([3, 2], [3, 2])
    diagonal_controls = [[0, 7, 7, 7], [0, 2, 4, 0], [1, 0, 0, 0], [0, 0, 0, 8]]
    for original, expected in zip([C0, C2, C1, C3], diagonal_controls):
        assert multiply(multiply(L, cut(original)), R) == expected
    arrow_left, arrow_right = [2, 1, 10, 6], [2, 8, 11, 4]
    for original, expected in zip(diagonal_controls, COLUMN_ARROW):
        assert multiply(multiply(arrow_left, original), arrow_right) == expected
    rng = random.Random(0xC007)
    counts = []
    for t, count in [(1, 16), (2, 512), (3, 128)]:
        for k in range(count):
            mask = k if t == 1 else rng.getrandbits(4 ** t)
            coeff = {word for word in range(4 ** t) if mask >> word & 1}
            direct = rank(tensor_matrix(coeff, t))
            assert row_rank(coeff, t) == direct
            column_direct = rank(tensor_matrix(coeff, t, COLUMN))
            assert column_rank(coeff, t) == column_direct
            assert column_character_bound(coeff, t) <= column_direct
        counts.append({'tensor_exponent': t, 'binary_checks': count})
    field_controls = []
    for t, count in [(1, 256), (2, 64), (3, 8)]:
        for k in range(count):
            if t == 1:
                coeff = {word: (k >> (2 * word)) & 3 for word in range(4)}
            else:
                coeff = {word: rng.randrange(4) for word in range(4 ** t)}
            field_control({word: value for word, value in coeff.items() if value}, t)
        field_controls.append({'tensor_exponent': t, 'F4_checks': count})
    local_rank_one_counts = [0, 0]
    for values in itertools.product(range(4), repeat=4):
        a, b, c, d = values
        matrices = []
        for local in [ROW, COLUMN]:
            matrices.append([[xor_all(value for value, base in zip(values, local)
                                      if base[i] >> j & 1)
                              for j in range(4)] for i in range(4)])
        rr, cr = map(field_rank, matrices)
        row_condition = a == 0 and any(values) and field_multiply(b, b) == field_multiply(c, d)
        column_condition = ((a == b and a != 0 and c == d == 0)
                            or (a == b == d == 0 and c != 0)
                            or (a == b == c == 0 and d != 0))
        assert (rr == 1) == row_condition
        assert (cr == 1) == column_condition
        local_rank_one_counts[0] += rr == 1
        local_rank_one_counts[1] += cr == 1
    assert local_rank_one_counts == [15, 9]
    # Complete binary square, using literal matrices and a Gray-code traversal.
    t = 2
    basis = [tensor_matrix({word}, t, COLUMN_ARROW) for word in range(16)]
    current, previous = [0] * 16, 0
    square_counts = {1: 0, 2: 0}
    for i in range(1, 1 << 16):
        mask = i ^ (i >> 1)
        changed = (mask ^ previous).bit_length() - 1
        current = [x ^ y for x, y in zip(current, basis[changed])]
        previous = mask
        direct = rank(current)
        coeff = {word for word in range(16) if mask >> word & 1}
        predicted = low_column_rank(coeff, t)
        assert predicted == (direct if direct <= 2 else None)
        if direct <= 2:
            square_counts[direct] += 1
    assert square_counts == {1: 9, 2: 42}
    rank_two_profiles = []
    for t in range(1, 5):
        histogram = {}
        for coeff, predicted_row in column_rank_two_pieces(t):
            assert row_rank(coeff, t) == predicted_row
            if t <= 2:
                assert rank(tensor_matrix(coeff, t, COLUMN_DIAGONAL)) == 2
            histogram[predicted_row] = histogram.get(predicted_row, 0) + 1
        assert sum(histogram.values()) == math.comb(3 ** t, 2) + t * 3 ** (t - 1)
        rank_two_profiles.append({'tensor_exponent': t, 'row_rank_histogram': histogram,
                                  'column_rank': 2})
    rank_one_counts = []
    for t in range(1, 5):
        for choices in itertools.product(({2}, {3}, {1, 2, 3}), repeat=t):
            coeff = coefficient_product(choices)
            rr, cr = ranks(coeff, t)
            v_count = sum(len(choice) == 3 for choice in choices)
            assert rr == 1 and cr == 4 ** v_count
        for choices in itertools.product(({0, 1}, {2}, {3}), repeat=t):
            coeff = coefficient_product(choices)
            rr, cr = ranks(coeff, t)
            u_count = sum(len(choice) == 2 for choice in choices)
            assert cr == 1 and rr == 4 ** u_count
        rank_one_counts.append({'tensor_exponent': t, 'rank_one_pieces_in_each_cut': 3 ** t})
    conic_profiles = []
    for t in range(1, 5):
        profile = []
        for degree in range(2 * t + 1):
            coeff = {sum(d << (2 * i) for i, d in enumerate(word))
                     for word in itertools.product([1, 2, 3], repeat=t)
                     if sum({1: 1, 2: 0, 3: 2}[d] for d in word) == degree}
            rr, cr = ranks(coeff, t)
            predicted = [min(degree + 1, 2 * t - degree + 1, t + 1), math.comb(2 * t, degree)]
            assert [rr, cr] == predicted
            assert cr == column_character_bound(coeff, t)
            profile.append(predicted)
        conic_profiles.append({'tensor_exponent': t, 'coefficient_rank_pairs': profile})
    # Exact arithmetic used by the full eight-label all-field theorem.
    assert 128**2 * 32 > 19**2 * 33**2
    assert 63**2 > 11**2 * 32
    assert 27 > 3**2
    assert 19**70 > 2**297
    assert 99**2 > 2 * 70**2
    q = 8
    flip = [[int(i == (j % q) * q + j // q) for j in range(q*q)] for i in range(q*q)]
    shifted_flip_ranks = []
    for scalar in range(4):
        shifted = [[entry ^ (scalar if i == j else 0) for j, entry in enumerate(row)]
                   for i, row in enumerate(flip)]
        r = field_rank(shifted)
        assert r == (28 if scalar == 1 else 64)
        shifted_flip_ranks.append(r)
    return {'row_basis': 'a=C0+C2, b=C2, c=C1, d=C3',
            'local_row_matrices': ROW, 'local_column_normal_matrices': COLUMN,
            'binary_direct_checks': counts, 'extension_field_checks': field_controls,
            'rank_one_product_controls': rank_one_counts,
            'local_F4_rank_one_counts': {'row': 15, 'column': 9},
            'all_active_two_sided_column_forms_u_b_c_d': diagonal_controls,
            'column_arrow_left': arrow_left, 'column_arrow_right': arrow_right,
            'column_arrow_forms_u_b_c_d': COLUMN_ARROW,
            'complete_binary_square_low_column_rank_counts': square_counts,
            'column_rank_two_product_profiles': rank_two_profiles,
            'rank_one_conic_coefficient_profiles': conic_profiles,
            'q8_shifted_flip_ranks_F4': shifted_flip_ranks,
            'q8_full_span_cost_theorem': 'S>19>8^(99/70)>8^sqrt(2) over every characteristic-two field',
            'scope': 'exact rank engine and all-field structural theorems; no exponent improvement'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify()
    rendered = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end='')
