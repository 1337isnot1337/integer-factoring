#!/usr/bin/env python3
"""Replay exact code-section identities and the small complete rank screens.

The optional native checkpoint audit verifies coverage and checksums; replay
check_code_kernel_width5.cpp to independently regenerate all 2^26 states.
"""

import argparse
import hashlib
import json
import random
from collections import Counter
from pathlib import Path
from time import perf_counter

from code_compressed_crossing import literal_pencil_entry
from code_kernel_crossing import (base_kernels, convolution_rank,
                                 matrices_from_kernels, rank_profile,
                                 toggle, transpose, xor_shift)
from convolution_rank5 import (ALT_RANKS, QUAD_MASKS, convolution_rank5,
                               group_from_anf)
from general_section_crossing import (constant_representatives,
                                     general_kernels, slope_family,
                                     validate_slopes)


def matrix_rank(rows):
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return len(basis)


def polynomial_product(a, b):
    answer = 0
    while b:
        bit = b & -b
        answer ^= a << (bit.bit_length()-1)
        b ^= bit
    return answer


def determinant_polynomial(rows):
    """Leibniz expansion, only used for matrices of order at most three."""
    from itertools import permutations
    answer = 0
    for permutation in permutations(range(len(rows))):
        term = 1
        for i, j in enumerate(permutation):
            term = polynomial_product(term, rows[i][j])
        answer ^= term
    return answer


def literal_section(width, physical, quartet):
    """Independent sum of original local pencils, without normal forms."""
    ambient = len(physical)
    representatives = [sum(((a & label).bit_count() & 1) << j
                           for j, (_, a) in enumerate(physical))
                       for label in range(1 << width)]
    code = []
    for word in range(1 << ambient):
        projection = 0
        for j, (p, _) in enumerate(physical):
            if word >> j & 1:
                projection ^= p
        if projection == 0:
            code.append(word)
    assert len(code) == 1 << (ambient-width)
    x, y, u, v = (representatives[a] for a in quartet)
    result = 0
    for a in code:
        for b in code:
            result ^= literal_pencil_entry(ambient, x ^ a, y ^ b, u, v)
    divisor = ambient-width
    assert result & ((1 << divisor)-1) == 0
    return result >> divisor


def describe_physical(width, physical):
    slopes = [0] * (1 << width)
    constants = 0
    for p, a in physical:
        if p:
            slopes[p] ^= a
            constants ^= 1 << (p-1)
    validate_slopes(slopes, width)
    return slopes, constants


def kernel_entry(kernels, quartet):
    x, y, u, v = quartet
    if x ^ y ^ u ^ v:
        return 0
    return sum(((row[x ^ y] >> (x ^ u)) & 1) << degree
               for degree, row in enumerate(kernels))


def check_augmentations(kernels, slopes, constants, width):
    dimension = 1 << width
    for s in range(dimension):
        polynomial = sum((rows[s].bit_count() & 1) << degree
                         for degree, rows in enumerate(kernels))
        assert polynomial == 1 << width
    matrix = [[0] * width for _ in range(width)]
    for p in range(1, dimension):
        for i in range(width):
            if p >> i & 1:
                matrix[i][0:width] = [value ^ ((slopes[p] >> j) & 1)
                                     for j, value in enumerate(matrix[i])]
    determinant = determinant_polynomial([
        [((2 if i == j else 0) ^ (3 if matrix[i][j] else 0))
         for j in range(width)] for i in range(width)])
    for w in range(dimension):
        polynomial = sum((sum((row[s] >> w) & 1 for s in range(dimension)) & 1)
                         << degree for degree, row in enumerate(kernels))
        assert polynomial == determinant
    for s in range(dimension):
        selected = [p for p in range(1, dimension)
                    if ((constants >> (p-1)) & 1)
                    ^ ((slopes[p] & s).bit_count() & 1)]
        c = [0] * width
        for p in selected:
            for i in range(width):
                if p >> i & 1:
                    c[i] ^= p
        expected = ((1 << dimension)-1) if matrix_rank(c) == width else 0
        assert kernels[0][s] == expected


def check_identities(rng):
    counts = Counter()
    for width in range(1, 4):
        dimension = 1 << width
        # Physical basis columns plus extra columns. Arbitrary section rows
        # on extra coordinates, with basis-coordinate rows solving Pi J=I.
        for case in range(12):
            extras = [(rng.randrange(dimension), rng.randrange(dimension))
                      for _ in range(case % 3)]
            physical = [(1 << i, 1 << i) for i in range(width)]
            for p, a in extras:
                for i in range(width):
                    if p >> i & 1:
                        physical[i] = (1 << i, physical[i][1] ^ a)
            physical += extras
            slopes, constants = describe_physical(width, physical)
            kernels = general_kernels(width, slopes, constants)
            check_augmentations(kernels, slopes, constants, width)
            counts['augmentation_and_determinant_profiles'] += 1
            for _ in range(32):
                quartet = tuple(rng.randrange(dimension) for _ in range(4))
                assert literal_section(width, physical, quartet) == kernel_entry(kernels, quartet)
                counts['literal_section_quartets'] += 1
            dense_rows, dense_cols = matrices_from_kernels(kernels, width)
            assert rank_profile(kernels, width) == [
                (matrix_rank(r), matrix_rank(c)) for r, c in zip(dense_rows, dense_cols)]
            counts['dense_rank_profiles'] += 1
            for a in range(dimension):
                shifted = constants
                for p in range(1, dimension):
                    if (slopes[p] & a).bit_count() & 1:
                        shifted ^= 1 << (p-1)
                other = general_kernels(width, slopes, shifted)
                for degree in range(width+1):
                    for s in range(dimension):
                        assert other[degree][s] == xor_shift(kernels[degree][s ^ a], a, width)
                assert rank_profile(other, width) == rank_profile(kernels, width)
                counts['translation_gauges'] += 1
        # Coordinate-section toggles agree with the general-section formula.
        slopes = slope_family(width, 0)
        base_constants = sum(1 << ((1 << j)-1) for j in range(width))
        kernels = base_kernels(width)
        assert kernels == general_kernels(width, slopes, base_constants)
        constants = base_constants
        for _ in range(32):
            p = rng.randrange(1, dimension)
            kernels = toggle(kernels, p, width)
            constants ^= 1 << (p-1)
            assert kernels == general_kernels(width, slopes, constants)
            counts['coordinate_toggle_profiles'] += 1
        # The augmentation formula holds without Pi J=I too.
        for _ in range(32):
            slopes = [0] + [rng.randrange(dimension) for _ in range(dimension-1)]
            constants = rng.randrange(1 << (dimension-1))
            kernels = general_kernels(width, slopes, constants)
            check_augmentations(kernels, slopes, constants, width)
            counts['arbitrary_slope_profiles'] += 1
    return dict(counts)


def check_classifier(rng, full):
    counts = Counter()
    for _ in range(10000):
        kernel = rng.getrandbits(32)
        assert convolution_rank5(kernel) == convolution_rank(5, kernel)
        counts['random_width_five'] += 1
    if full:
        higher = [m for m in range(32) if m.bit_count() >= 3]
        # One representative of each possible quadratic rank, with two
        # additional rank-two controls; every higher-degree perturbation.
        for pattern in (0, 1, 3, 11, 129):
            base = sum(1 << mask for bit, mask in enumerate(QUAD_MASKS)
                       if pattern >> bit & 1)
            for value in range(1 << len(higher)):
                anf = base | sum(((value >> j) & 1) << m
                                 for j, m in enumerate(higher))
                kernel = group_from_anf(anf)
                assert convolution_rank5(kernel) == convolution_rank(5, kernel)
                counts[f'complete_quadratic_pattern_{pattern}_rank_{ALT_RANKS[pattern]}'] += 1
            convolution_rank.cache_clear()
        for kernel in range(1 << 16):
            assert convolution_rank5(kernel) == 2 * convolution_rank(4, kernel)
            counts['embedded_width_four'] += 1
    convolution_rank.cache_clear()
    return dict(counts)


def screen_coordinate(width):
    parities = [p for p in range(1, 1 << width) if p.bit_count() > 1]
    kernels = base_kernels(width)
    baseline = rank_profile(kernels, width)
    best = max(r*c for r, c in baseline)
    histogram = Counter()
    digest = hashlib.sha256()
    failures = 0
    for index in range(1 << len(parities)):
        if index:
            bit = (index & -index).bit_length()-1
            kernels = toggle(kernels, parities[bit], width)
        profile = rank_profile(kernels, width)
        digest.update(json.dumps(profile, separators=(',', ':')).encode())
        q = max(r*c for r, c in profile)
        histogram[q] += 1
        best = min(best, q)
        failures += int(any(r < br or c < bc for (r, c), (br, bc) in zip(profile, baseline)))
    return dict(width=width, cases=1 << len(parities), optimum=best,
                baseline=baseline, coefficient_dominance_failures=failures,
                maximum_product_histogram=dict(sorted(histogram.items())),
                profile_digest=digest.hexdigest())


def screen_general(width):
    dimension = 1 << width
    free_count = dimension-1-width
    histogram = Counter()
    digest = hashlib.sha256()
    for free in range(1 << (width*free_count)):
        slopes = slope_family(width, free)
        validate_slopes(slopes, width)
        for constants in constant_representatives(slopes, width):
            kernels = general_kernels(width, slopes, constants)
            profile = rank_profile(kernels, width)
            digest.update(json.dumps(profile, separators=(',', ':')).encode())
            histogram[max(r*c for r, c in profile)] += 1
    assert sum(histogram.values()) == 1 << ((width+1)*free_count)
    return dict(width=width, cases=sum(histogram.values()), optimum=min(histogram),
                maximum_product_histogram=dict(sorted(histogram.items())),
                profile_digest=digest.hexdigest())


def audit_native(prefix):
    chunks = [json.loads(Path(f'{prefix}-chunk-{i}.json').read_text()) for i in range(64)]
    witnesses = [0] * 6
    checksum = 0
    for i, chunk in enumerate(chunks):
        assert chunk['chunk'] == i
        assert chunk['begin'] == i << 20
        assert chunk['end'] == (i+1) << 20
        assert chunk['covered'] == 1 << 20
        assert chunk['threshold'] == 15312
        assert chunk['violations'] == 0
        assert sum(chunk['witness_degree_counts']) == 1 << 20
        witnesses = [a+b for a, b in zip(witnesses, chunk['witness_degree_counts'])]
        checksum = (checksum+chunk['checksum_mod_2_64']) % (1 << 64)
    return dict(width=5, scope='coordinate sections; all ambient lengths',
                chunks=64, classes=1 << 26, optimum=15312, violations=0,
                witness_degree_counts=witnesses, checksum_mod_2_64=checksum,
                elapsed_seconds=max(c['elapsed_seconds'] for c in chunks),
                checkpoint_note='Coverage audit; native source regenerates every state.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full-classifier', action='store_true')
    parser.add_argument('--native-prefix', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    start = perf_counter()
    rng = random.Random(261030)
    result = {'identities': check_identities(rng),
              'classifier': check_classifier(rng, args.full_classifier),
              'coordinate_sections': [screen_coordinate(width) for width in range(1, 5)],
              'general_sections': [screen_general(width) for width in range(1, 4)]}
    if args.native_prefix:
        result['native_width_five'] = audit_native(args.native_prefix)
    result['source_sha256'] = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                             for name in ('code_compressed_crossing.py', 'code_kernel_crossing.py',
                                          'general_section_crossing.py', 'convolution_rank5.py',
                                          'check_code_sections.py', 'check_code_kernel_width5.cpp')}
    result['seconds'] = perf_counter()-start
    text = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text)


if __name__ == '__main__':
    main()
