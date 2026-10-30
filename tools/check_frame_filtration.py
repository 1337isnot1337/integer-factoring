#!/usr/bin/env python3
"""Exact two-cut rank quotient for mixed affine-selector crossing frames."""

import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import comb, log2
from pathlib import Path
import random
from time import perf_counter

from code_kernel_crossing import convolution_rank, rank_profile, transpose, xor_shift
from general_section_crossing import general_kernels, slope_family, validate_slopes


def binary_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length()-1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return len(pivots)


def group_to_square_zero(kernel, width):
    for j in range(width):
        step = 1 << j
        lower = ((1 << (1 << width))-1)//((1 << step)+1)
        kernel ^= (kernel >> step) & lower
    return kernel


def polynomial_coordinates(kernel, width):
    """Rows B encode coefficients of Boolean S_B times square-zero N_A."""
    rows = [group_to_square_zero(row, width) for row in kernel]
    for j in range(width):
        for s in range(1 << width):
            if s >> j & 1:
                rows[s] ^= rows[s ^ (1 << j)]
    return rows


def quotient_flags(kernel, width):
    rows = polynomial_coordinates(kernel, width)
    flags = 0
    for degree in range(width+1):
        representative = (1 << degree)-1
        flags |= ((rows[representative] >> representative) & 1) << degree
    for b, row in enumerate(rows):
        while row:
            bit = row & -row
            a = bit.bit_length()-1
            assert a.bit_count() >= b.bit_count()
            if a.bit_count() == b.bit_count():
                assert a == b
            row ^= bit
        assert ((rows[b] >> b) & 1) == ((flags >> b.bit_count()) & 1)
    return flags


def canonical_kernels(width):
    """Independent coefficients of product_j (1+z S_j N_j)."""
    dimension = 1 << width
    result = [[0]*dimension for _ in range(width+1)]
    for s in range(dimension):
        for w in range(dimension):
            if w & ~s:
                continue
            remaining = s.bit_count()-w.bit_count()
            for degree in range(w.bit_count(), s.bit_count()+1):
                needed = degree-w.bit_count()
                if not needed & ~remaining:
                    result[degree][s] |= 1 << w
    return result


def combine(kernels, flags):
    result = [0]*len(kernels[0])
    for i, kernel in enumerate(kernels):
        if flags >> i & 1:
            result = [a ^ b for a, b in zip(result, kernel)]
    return result


def reverse_flags(flags, width):
    return sum(((flags >> i) & 1) << (width-i) for i in range(width+1))


@lru_cache(None)
def exterior_rank(width, flags):
    monomials = [a | (a << width) for a in range(1 << width)
                 if flags >> a.bit_count() & 1]
    rows = []
    for source in range(1 << (2*width)):
        row = 0
        for monomial in monomials:
            if not source & monomial:
                row ^= 1 << (source | monomial)
        rows.append(row)
    return binary_rank(rows)


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def inclusion_formula(n, low, high):
    gap = high-low
    bound = min(low, n-high)
    return sum(choose(n, a)-choose(n, a-1) for a in range(bound+1)
               if not gap & (bound-a))


def inclusion_direct(n, low, high):
    inputs = [sum(1 << j for j in subset) for subset in combinations(range(n), low)]
    outputs = [sum(1 << j for j in subset) for subset in combinations(range(n), high)]
    return binary_rank(sum((source & target == source) << j for j, target in enumerate(outputs))
                       for source in inputs)


@lru_cache(None)
def elementary_rank(n, degree):
    if not 0 <= degree <= n:
        return 0
    return sum(inclusion_formula(n, low, low+degree) for low in range(n-degree+1))


def omega_rank(width, degree):
    return sum(comb(width, k)*(1 << (width-k))*elementary_rank(k, degree)
               for k in range(degree, width+1))


def binomial_tail(width, limit):
    return sum(comb(2*width, j) for j in range(limit+1))


def check_elementary_identification():
    combinations_checked = 0
    for width in range(1, 6):
        canonical = canonical_kernels(width)
        assert canonical == list(reversed(general_kernels(width, slope_family(width, 0), 0)))
        for flags in range(1 << (width+1)):
            tensor = combine(canonical, flags)
            assert quotient_flags(tensor, width) == flags
            assert quotient_flags(transpose(tensor, width), width) == reverse_flags(flags, width)
            expected = (exterior_rank(width, flags), exterior_rank(width, reverse_flags(flags, width)))
            assert rank_profile([tensor], width)[0] == expected
            combinations_checked += 1
        for degree in range(width+1):
            assert exterior_rank(width, 1 << degree) == omega_rank(width, degree)
    return combinations_checked


def check_frames(rng):
    profile_checks = 0
    mixed_checks = 0
    structural_coefficients = 0
    for width in range(1, 6):
        dimension = 1 << width
        canonical = canonical_kernels(width)
        frames = []
        for _ in range(32):
            free = rng.getrandbits(width*(dimension-1-width))
            slopes = slope_family(width, free)
            constants = rng.getrandbits(dimension-1)
            validate_slopes(slopes, width)
            frame = list(reversed(general_kernels(width, slopes, constants)))
            assert combine(frame, (1 << (width+1))-1) == [1 << s for s in range(dimension)]
            for degree, kernel in enumerate(frame):
                assert quotient_flags(kernel, width) == 1 << degree
                assert quotient_flags(transpose(kernel, width), width) == 1 << (width-degree)
                actual = rank_profile([kernel], width)[0]
                expected = (omega_rank(width, degree), omega_rank(width, width-degree))
                assert actual[0] >= expected[0] and actual[1] >= expected[1]
                structural_coefficients += 1
            frames.append(frame)
            profile_checks += 1
        for _ in range(300):
            tensor, flags = [0]*dimension, 0
            for _ in range(rng.randrange(1, 9)):
                selected = rng.getrandbits(width+1)
                local = combine(rng.choice(frames), selected)
                if rng.getrandbits(1):
                    local = transpose(local, width)
                    selected = reverse_flags(selected, width)
                flags ^= selected
                tensor = [a ^ b for a, b in zip(tensor, local)]
            assert quotient_flags(tensor, width) == flags
            assert quotient_flags(transpose(tensor, width), width) == reverse_flags(flags, width)
            actual = rank_profile([tensor], width)[0]
            expected = rank_profile([combine(canonical, flags)], width)[0]
            assert actual[0] >= expected[0] and actual[1] >= expected[1]
            mixed_checks += 1
        convolution_rank.cache_clear()
    return {'frame_profiles': profile_checks, 'structural_coefficient_checks': structural_coefficients,
            'mixed_frame_and_symmetry_rank_comparisons': mixed_checks}


def symmetry_check():
    partition = {frozenset((0, 3)), frozenset((1, 2))}
    kinds = []
    for p in permutations(range(4)):
        if {frozenset((p[0], p[3])), frozenset((p[1], p[2]))} != partition:
            continue
        transformations = []
        for s in range(2):
            for w in range(2):
                quartet = (0, s, w, s ^ w)
                changed = tuple(quartet[j] for j in p)
                transformations.append((changed[0] ^ changed[1], changed[0] ^ changed[2]))
        if transformations == [(s, w) for s in range(2) for w in range(2)]:
            kinds.append('identity')
        elif transformations == [(w, s) for s in range(2) for w in range(2)]:
            kinds.append('transpose')
        else:
            raise AssertionError('unexpected crossing symmetry on conserved kernel')
    assert kinds.count('identity') == kinds.count('transpose') == 4
    return {'identity_on_kernel': 4, 'transpose_kernel': 4}


def check_entropy():
    direct_cases = 0
    for n in range(9):
        for low in range(n+1):
            for high in range(low, n+1):
                assert inclusion_formula(n, low, high) == inclusion_direct(n, low, high)
                direct_cases += 1
    tail_checks = 0
    examples = []
    for width in range(1, 97):
        for degree in range(width+1):
            rank = omega_rank(width, degree)
            tail = binomial_tail(width, width-degree)
            assert tail <= (width+1)*rank <= 2*(width+1)*tail
            tail_checks += 1
        if width in (8, 16, 32, 64, 96):
            j = width//2
            examples.append({'width': width, 'forced_row_tail': binomial_tail(width, width-j),
                             'forced_column_tail': binomial_tail(width, j),
                             'denominator_per_rank': width+1,
                             'central_literal_ranks': [omega_rank(width, j), omega_rank(width, width-j)]})
    assert inclusion_direct(4, 1, 2) == 3 < choose(4, 1)
    return {'direct_inclusion_cases': direct_cases, 'exact_tail_bounds': tail_checks,
            'limiting_rho': 2-0.75*log2(3), 'selected_exact_bounds': examples,
            'characteristic_two_negative_control': {'n': 4, 'low': 1, 'high': 2, 'rank': 3}}


def nonlinear_negative_control():
    width = 3
    dimension = 1 << width
    kernels = [[0]*dimension for _ in range(width+1)]
    for s in range(dimension):
        selected = {1 << j for j in range(width) if s >> j & 1}
        if s & 1 and s & 4:
            selected.symmetric_difference_update((1, 2, 3))
        total = 0
        for p in selected:
            total ^= p
        assert total == s
        coefficients = [1]+[0]*width
        for p in selected:
            for i in range(width, 0, -1):
                coefficients[i] ^= coefficients[i-1] ^ xor_shift(coefficients[i-1], p, width)
        for i, value in enumerate(coefficients):
            kernels[i][s] = value
    assert combine(kernels, (1 << (width+1))-1) == [1 << s for s in range(dimension)]
    coordinates = polynomial_coordinates(kernels[1], width)
    # The h(S)=S1*S3 toggle on the relation {e1,e2,e1+e2}
    # contributes S1*S3*N1*N2: a nonsymmetric weight-zero term.
    assert coordinates[5] >> 3 & 1
    try:
        quotient_flags(kernels[1], width)
    except AssertionError:
        pass
    else:
        raise AssertionError('nonlinear selector unexpectedly admitted by affine quotient')
    return {'width': width, 'toggle_relation': [1, 2, 3],
            'quadratic_selector_monomial': 5,
            'offdiagonal_equal_degree_witness': {'S_monomial': 5, 'N_monomial': 3},
            'crossing_identity': True, 'rank_profile': rank_profile(kernels, width)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    start = perf_counter()
    result = {'elementary_combination_identifications': check_elementary_identification(),
              'crossing_symmetries': symmetry_check(),
              'frames': check_frames(random.Random(261035)),
              'entropy': check_entropy(),
              'nonlinear_selector_negative_control': nonlinear_negative_control(),
              'source_pins': {'binary_entropy': '147652574e7ce1ff68b2c378712846cbf438da34',
                              'binary_entropy_path': 'research/binary-coefficient-entropy.md'}}
    result['sha256'] = {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                        for name in ('check_frame_filtration.py', 'code_kernel_crossing.py',
                                     'general_section_crossing.py')}
    result['seconds'] = perf_counter()-start
    text = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text)


if __name__ == '__main__':
    main()
