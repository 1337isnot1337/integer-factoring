#!/usr/bin/env python3
"""Exact local symbols and dense rank checks for fixed-target Feistel quotients."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random
from time import perf_counter

from check_frame_filtration import binary_rank, omega_rank


def literal_pencil(width, perm):
    """Expand original four-leg entries, independently of the rank quotient.

    This small routine is reproduced from this session's independently
    checked nonlinear-shear search helper; it uses no block-rank shortcut.
    """
    dimension = 1 << width
    rows = [[0]*dimension**2 for _ in range(width+1)]
    cols = [[0]*dimension**2 for _ in range(width+1)]
    for x, y, u, v in product(range(dimension), repeat=4):
        polynomial = 1
        for bit in range(width):
            a, b, c, d = ((perm[label] >> bit) & 1 for label in (x, y, u, v))
            ee = int(a == d and b == c)
            cc = int(a == c and b == d)
            local = cc | ((ee ^ cc) << 1)
            polynomial = ((polynomial if local & 1 else 0)
                          ^ (polynomial << 1 if local & 2 else 0))
        while polynomial:
            low = polynomial & -polynomial
            degree = low.bit_length()-1
            rows[degree][x*dimension+y] |= 1 << (u*dimension+v)
            cols[degree][x*dimension+u] |= 1 << (y*dimension+v)
            polynomial ^= low
    return list(zip(rows, cols))


def tensor_xor(left, right):
    return tuple([a ^ b for a, b in zip(x, y)] for x, y in zip(left, right))


def subsets(mask):
    result = []
    while True:
        result.append(mask)
        if not mask:
            return result
        mask = (mask-1) & result[0]


def multiply(columns, vector):
    result = 0
    while vector:
        bit = vector & -vector
        result ^= columns[bit.bit_length()-1]
        vector ^= bit
    return result


def compose(left, right):
    return [multiply(left, column) for column in right]


def augmentation_permutation(pi, targets):
    """Conjugate a label permutation by the target monomial basis matrix."""
    dimension = len(pi)
    mask = (1 << targets)-1
    basis = [sum(1 << ((index & ~mask) | subset) for subset in subsets(index & mask))
             for index in range(dimension)]
    assert compose(basis, basis) == [1 << i for i in range(dimension)]
    permutation = [1 << pi[i] for i in range(dimension)]
    return compose(basis, compose(permutation, basis))


def first_nontrivial_symbol(columns, targets):
    mask = (1 << targets)-1
    for source, column in enumerate(columns):
        error = column ^ (1 << source)
        while error:
            bit = error & -error
            destination = bit.bit_length()-1
            if (destination & mask).bit_count() <= (source & mask).bit_count():
                return [source, destination]
            error ^= bit
    return None


def tensor_sum(pencil, mask):
    dimension = len(pencil[0][0])
    result = ([0]*dimension, [0]*dimension)
    for degree, tensor in enumerate(pencil):
        if mask >> degree & 1:
            result = tensor_xor(result, tensor)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    start = perf_counter()
    # Independent local calculation in basis (1,a,b,ab).
    basis = [sum(1 << j for j in subsets(i)) for i in range(4)]
    swap = [1, 4, 2, 8]
    identity = [1, 2, 4, 8]
    r_standard = [9, 0, 0, 9]
    r = compose(basis, compose(r_standard, basis))
    assert r == [14, 14, 14, 0]
    n = [a ^ b for a, b in zip(swap, identity)]
    symbol = [sum(1 << j for j in range(4)
                  if (column >> j & 1) and j.bit_count() == i.bit_count())
              for i, column in enumerate(r)]
    assert symbol == n == [0, 6, 6, 0]
    u = [a ^ b for a, b in zip(swap, r)]
    assert compose(u, u) == identity
    assert compose(u, r) == r
    assert compose(r, r) == [0]*4
    assert binary_rank(r) == binary_rank(n) == 1

    rng = Random(730172)
    structural_models = 0
    dense_models = 0
    rank_comparisons = 0
    for width in range(2, 7):
        dimension = 1 << width
        for targets in range(1, width):
            for trial in range(8):
                values = [rng.randrange(1 << targets) for _ in range(1 << (width-targets))]
                pi = [x ^ values[x >> targets] for x in range(dimension)]
                assert sorted(pi) == list(range(dimension))
                assert first_nontrivial_symbol(augmentation_permutation(pi, targets), targets) is None
                structural_models += 1
                if width > 4 or trial >= 2:
                    continue
                original = literal_pencil(width, list(range(dimension)))
                sheared = literal_pencil(width, pi)
                assert [binary_rank(row) for row, _ in original] == [omega_rank(width, i) for i in range(width+1)]
                for sample in range(64):
                    a, b = [rng.randrange(1 << (width+1)) for _ in range(2)]
                    actual = tensor_xor(tensor_sum(original, a), tensor_sum(sheared, b))
                    quotient = tensor_sum(original, a ^ b)
                    for matrix, small in zip(actual, quotient):
                        assert binary_rank(matrix) >= binary_rank(small)
                        rank_comparisons += 1
                dense_models += 1

    # A mixture of three unrelated, arbitrary truth-table functions uses the
    # same coefficient quotient; no quadratic-degree assumption enters.
    mixed_models = 0
    for width, targets in [(3, 1), (4, 1), (4, 2), (4, 3)]:
        dimension = 1 << width
        original = literal_pencil(width, list(range(dimension)))
        pencils = []
        for family in range(3):
            values = [rng.randrange(1 << targets) for _ in range(1 << (width-targets))]
            pi = [x ^ values[x >> targets] for x in range(dimension)]
            pencils.append(literal_pencil(width, pi))
        for sample in range(128):
            actual = ([0]*dimension**2, [0]*dimension**2)
            coefficient_mask = 0
            for pencil in pencils:
                mask = rng.randrange(1 << (width+1))
                coefficient_mask ^= mask
                actual = tensor_xor(actual, tensor_sum(pencil, mask))
            quotient = tensor_sum(original, coefficient_mask)
            for matrix, small in zip(actual, quotient):
                assert binary_rank(matrix) >= binary_rank(small)
                rank_comparisons += 1
        mixed_models += 1

    # The different-target two-round control fails the common-symbol premise.
    negative_pi = []
    for label in range(16):
        first = label ^ ((label >> 1 & 1) & (label >> 2 & 1))
        second = first ^ (((first & 1) & (first >> 3 & 1)) << 1)
        negative_pi.append(second)
    witness = first_nontrivial_symbol(augmentation_permutation(negative_pi, 2), 2)
    assert witness is not None
    negative_inverse = [0]*16
    for i, j in enumerate(negative_pi):
        negative_inverse[j] = i
    negative_original = literal_pencil(4, list(range(16)))
    negative_relabeled = literal_pencil(4, negative_inverse)
    negative_tensor = tensor_xor(tensor_xor(negative_original[0], negative_original[1]),
                                 negative_relabeled[1])
    negative_ranks = list(map(binary_rank, negative_tensor))
    assert negative_ranks == [254, 18]
    assert list(map(binary_rank, negative_original[0])) == [256, 1]
    source_profile = [[omega_rank(7, i), omega_rank(7, 7-i)] for i in range(8)]
    assert source_profile[3] == [1820, 576]
    assert min(r for r, _ in source_profile[:4]) == 1820
    assert min(c for _, c in source_profile[3:]) == 576
    assert 1820*576 > 956**2
    # 577/408 > sqrt(2), and 128^(577/408) < 956, by integers only.
    assert 577**2 > 2*408**2
    assert 128**577 < 956**408
    result = {'status':'PASS', 'scope':'Fixed-target elementary-pencil Feistel coefficient spans only.',
              'local_column_augmentation_columns':r, 'local_column_symbol_columns':symbol,
              'structural_arbitrary_truth_table_models':structural_models,
              'dense_two_pencil_models':dense_models, 'dense_three_pencil_models':mixed_models,
              'individual_cut_rank_comparisons':rank_comparisons,
              'width7_elementary_rank_profiles':source_profile,
              'width7_forced_piece_rank_product':1820*576,
              'different_target_negative_symbol_witness':witness,
              'different_target_negative_rank_quotient':{
                  'combination':'B_0+B_1+B_1^pi', 'actual_ranks':negative_ranks,
                  'proposed_elementary_quotient_ranks':[256, 1]},
              'sha256':{name:sha256((Path(__file__).resolve().parents[1]/name).read_bytes()).hexdigest()
                        for name in ['tools/check_feistel_filtration.py',
                                     'research/fixed-target-feistel-barrier.md',
                                     'tools/check_frame_filtration.py']},
              'seconds':perf_counter()-start}
    text = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text)


if __name__ == '__main__':
    main()
