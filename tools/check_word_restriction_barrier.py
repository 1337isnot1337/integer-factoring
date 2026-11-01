#!/usr/bin/env python3
"""Exact replay for the arbitrary-word rectangular crossing restriction."""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path
import random
import sys
from time import perf_counter

sys.path.insert(0, str(Path(__file__).with_name('reference')))

from check_quadratic_crossing import local_support_checks
from hamming_restrictions import histogram_profile
from hamming_tensor_power import local_terms, tensor_profile
from quadratic_crossing import (alternating_rank, binary_rank,
                                binomial_profile_from_histogram,
                                expanded_profile, truth_word)


def crossing(bits):
    x, y, u, v = bits
    return x == v and y == u


def bits(pattern):
    return tuple((pattern >> (3-j)) & 1 for j in range(4))


def source_support_table():
    supports = [sum(bool(value) << p for p, value in enumerate(local))
                for local, _ in local_terms()]
    relations = local_support_checks()
    contained = []
    for normal in range(1, 16):
        for offset in range(2):
            flat = sum(((p & normal).bit_count() & 1 == offset) << p for p in range(16))
            indices = [j for j, support in enumerate(supports) if not flat & ~support]
            if indices:
                contained.append({'normal': normal, 'offset': offset, 'terms': indices})
    assert contained == [{'normal': 3, 'offset': 0, 'terms': [13, 14]},
                         {'normal': 5, 'offset': 0, 'terms': [0, 5]}]
    assert all(support != 65535 for support in supports)
    symmetries = [p for p in permutations(range(4))
                  if all(crossing(tuple(bits(v)[j] for j in p)) == crossing(bits(v))
                         for v in range(16))]
    assert len(symmetries) == 8
    orbit = {tuple(sorted((p[a], p[b]))) for p in symmetries
             for a, b in ((1, 3), (2, 3))}
    assert orbit == {(0, 1), (0, 2), (1, 3), (2, 3)}
    cross_pairs = ((0, 1), (0, 2), (3, 1), (3, 2))
    witnesses = []
    specified = {
        (0, 5): ((3, 1), 1),
        (1, 6): ((3, 2), 3),
        (2, 7): ((3, 1), 2),
        (3, 8): ((3, 1), 3),
        (4, 10, 19, 20): ((3, 2), 3),
        (9, 21): ((3, 2), 1),
        (11, 12, 15, 16, 17, 18): ((3, 2), 2),
        (13, 14): ((3, 2), 1),
    }
    specified = {j: value for terms, value in specified.items() for j in terms}
    for j, (local, scalar) in enumerate(local_terms()):
        missing = []
        for a, b in cross_pairs:
            projection = {(bits(p)[a] << 1) | bits(p)[b]
                          for p, entry in enumerate(local) if entry}
            missing.extend((a, b, pattern) for pattern in range(4) if pattern not in projection)
        assert bool(missing) == (j != 22)
        if j != 22:
            (a, b), pattern = specified[j]
            assert (a, b, pattern) in missing
            witnesses.append({'term': j, 'pair': (a, b), 'missing_pattern': pattern})
        else:
            assert scalar == 16
            assert local == tuple(4 if crossing(bits(p)) else 0 for p in range(16))
    for permutation in symmetries:
        for witness in witnesses:
            a, b = witness['pair']
            transformed = {permutation[a], permutation[b]}
            assert len(transformed & {0, 3}) == 1
            assert len(transformed & {1, 2}) == 1
    return symmetries, {'local_product_equation_sets': len(relations),
                        'affine_hyperplanes': 30, 'full_spaces': 1,
                        'contained_hyperplanes': contained,
                        'crossing_symmetries': [list(p) for p in symmetries],
                        'exceptional_equality_orbit': sorted(orbit),
                        'noncrossing_term_projection_witnesses': witnesses,
                        'crossing_term': {'index': 22, 'local_polynomial': 'z^2 E',
                                          'low_scalar': 'z^4', 'leading_entry': 'z^2'}}


def span_polynomials(generators):
    result = [frozenset()]
    for generator in generators:
        result += [p.symmetric_difference(generator) for p in result]
    assert len(set(result)) == len(result)
    return result


def balanced_example():
    variables = 6
    rng = random.Random(610312026)
    pairs = tuple((1 << i) | (1 << j) for i, j in combinations(range(variables), 2))
    generators = [frozenset((3, 12, 48))]
    for _ in range(10000):
        if len(generators) == 3:
            break
        candidate = frozenset(p for p in pairs if rng.getrandbits(1))
        existing = span_polynomials(generators)
        if min(alternating_rank(variables, candidate.symmetric_difference(p))
               for p in existing) >= 4:
            generators.append(candidate)
    assert len(generators) == 3
    polynomial_labels = span_polynomials([frozenset([0]), *generators])
    words = [truth_word(variables, p) for p in polynomial_labels]
    length = 1 << variables
    weights = Counter(word.bit_count() for word in words)
    assert all(word in (0, (1 << length)-1)
               or length < 3*word.bit_count() < 2*length for word in words)
    disjoint_pairs = 0
    for f, g in product(words, repeat=2):
        if f & g == 0:
            assert f == 0 or g == 0 or f ^ g == (1 << length)-1
            disjoint_pairs += 1
    complete_pairs = 0
    for i, f in enumerate(words):
        for j, g in enumerate(words):
            if i >> 1 and j >> 1 and i >> 1 != j >> 1:
                projection = {(((f >> t) & 1) << 1) | ((g >> t) & 1)
                              for t in range(length)}
                assert projection == {0, 1, 2, 3}
                complete_pairs += 1
    full_affine = span_polynomials([frozenset([0]),
                                   *(frozenset([1 << j]) for j in range(variables)),
                                   *generators])
    affine_weights = Counter(truth_word(variables, p).bit_count() for p in full_affine)
    assert all(weight in (0, length) or length < 3*weight < 2*length
               for weight in affine_weights)
    return variables, words, {'variables': variables, 'label_dimension': 4,
                              'quadratic_generators': [sorted(p) for p in generators],
                              'nonzero_quadratic_ranks': [alternating_rank(variables, p)
                                                         for p in span_polynomials(generators)[1:]],
                              'weight_distribution': dict(sorted(weights.items())),
                              'disjoint_pairs': disjoint_pairs,
                              'complete_independent_pair_projections': complete_pairs,
                              'all_affine_extension_words': len(full_affine),
                              'all_affine_extension_weight_distribution': dict(sorted(affine_weights.items()))}


def fast_histogram(length, quartet):
    cells = [(1 << length)-1]
    for word in quartet:
        cells = [part for cell in cells for part in (cell & ~word, cell & word)]
    return tuple(cell.bit_count() for cell in cells)


def pair_span(a, b):
    a >>= 1
    b >>= 1
    assert a and b and a != b
    return tuple(sorted((a, b, a ^ b)))


def restricted_matrices(tensor, xs, ys):
    row_count = len(xs)*len(ys)
    rows, cols = [0]*row_count, [0]*row_count
    nx, ny = len(xs), len(ys)
    while tensor:
        low = tensor & -tensor
        index = low.bit_length()-1
        index, vi = divmod(index, nx)
        index, ui = divmod(index, ny)
        xi, yi = divmod(index, ny)
        rows[xi*ny+yi] |= 1 << (ui*nx+vi)
        cols[xi*ny+ui] |= 1 << (yi*nx+vi)
        tensor ^= low
    return rows, cols


def replay_rectangle(variables, words, symmetries):
    dimension = len(words)
    length = 1 << variables
    h = lambda label: (label >> 1) & 1
    xs = [label for label in range(dimension) if not h(label) and label >> 1]
    ys = [label for label in range(dimension) if h(label)]
    retained = len(xs)*len(ys)
    assert retained == dimension**2//4-dimension
    row_blocks = Counter(pair_span(x, y) for x in xs for y in ys)
    col_blocks = Counter(pair_span(u, v) for u in ys for v in xs)
    assert row_blocks == col_blocks and set(row_blocks.values()) == {8}
    tensors = [defaultdict(int) for _ in symmetries]
    crossing_tensor = 0
    histograms = set()
    nonuniform = 0
    larger_span_quartets = 0
    profile_evaluations = 0
    cache = {}
    index = 0
    for x, y, u, v in product(xs, ys, ys, xs):
        quartet = x, y, u, v
        joint_rank = binary_rank(label >> 1 for label in quartet)
        larger_span_quartets += joint_rank > 2
        if x == v and y == u:
            crossing_tensor |= 1 << index
        for which, perm in enumerate(symmetries):
            histogram = fast_histogram(length, tuple(words[quartet[j]] for j in perm))
            histograms.add(histogram)
            if which == 0:
                nonuniform += len({value for value in histogram if value}) > 1
            if histogram not in cache:
                # The source expansion and the compact oracle are separate
                # algorithms and must agree on every distinct histogram.
                profile = histogram_profile(histogram)
                assert profile == expanded_profile(binomial_profile_from_histogram(histogram))
                cache[histogram] = profile
            profile = cache[histogram]
            expected = (1 << (2*length+6)) if x == v and y == u else 0
            assert profile == expected
            assert profile.bit_count() & 1 == int(x == v and y == u)
            if joint_rank > 2:
                assert profile == 0
            if profile:
                assert pair_span(x, y) == pair_span(u, v)
                assert pair_span(x, u) == pair_span(y, v)
            while profile:
                bit = profile & -profile
                tensors[which][bit.bit_length()-1] |= 1 << index
                profile ^= bit
            profile_evaluations += 1
        index += 1
    assert nonuniform and larger_span_quartets
    assert crossing_tensor.bit_count() == retained
    for images in tensors:
        assert set(images) == {2*length+6}
        assert images[2*length+6] == crossing_tensor
        total = 0
        for tensor in images.values():
            total ^= tensor
            rows, cols = restricted_matrices(tensor, xs, ys)
            for r, row in enumerate(rows):
                x, y = xs[r//len(ys)], ys[r % len(ys)]
                for c in range(retained):
                    if row >> c & 1:
                        u, v = ys[c//len(xs)], xs[c % len(xs)]
                        assert pair_span(x, y) == pair_span(u, v)
            for r, row in enumerate(cols):
                x, u = xs[r//len(ys)], ys[r % len(ys)]
                for c in range(retained):
                    if row >> c & 1:
                        y, v = ys[c//len(xs)], xs[c % len(xs)]
                        assert pair_span(x, u) == pair_span(y, v)
        assert total == crossing_tensor
        assert tuple(binary_rank(matrix) for matrix in restricted_matrices(total, xs, ys)) == (retained, retained)
    # Exact decompositions by collective recombinations of source images.
    all_images = [tensor for family in tensors for tensor in family.values()]
    rng = random.Random(261031)
    combination_checks = 0
    for _ in range(32):
        pieces = [0, 0]
        for tensor in all_images:
            for i in range(2):
                if rng.getrandbits(1):
                    pieces[i] ^= tensor
        pieces.append(crossing_tensor ^ pieces[0] ^ pieces[1])
        detects = [(tensor & crossing_tensor).bit_count() for tensor in pieces]
        chosen = max(range(3), key=detects.__getitem__)
        assert 3*detects[chosen] >= retained
        ranks = tuple(binary_rank(matrix) for matrix in restricted_matrices(pieces[chosen], xs, ys))
        assert pieces[chosen] == crossing_tensor
        assert ranks == (retained, retained)
        assert min(ranks)*8 >= detects[chosen]
        assert min(ranks)*24 >= retained
        combination_checks += 1
    return {'dimension': dimension, 'rectangular_quartets': index,
            'retained_crossing_entries': retained, 'blocks_per_cut': len(row_blocks),
            'block_size': 8, 'quartets_of_joint_quotient_rank_above_two': larger_span_quartets,
            'source_profile_evaluations_including_symmetries': profile_evaluations,
            'unique_coefficient_degree': 2*length+6,
            'both_restricted_crossing_ranks': [retained, retained],
            'distinct_histograms_checked_by_both_algorithms': len(histograms),
            'nonuniform_original_histograms': nonuniform,
            'nonzero_coefficient_images_per_symmetry': [len(family) for family in tensors],
            'collective_three_piece_decompositions': combination_checks,
            'crossing_sha256': sha256(crossing_tensor.to_bytes((index+7)//8, 'little')).hexdigest()}


def negative_control():
    polynomials = [frozenset([1]), frozenset([5]), frozenset([2, 3]), frozenset([1])]
    words = [truth_word(3, p) for p in polynomials]
    histogram = fast_histogram(8, words)
    compact = binomial_profile_from_histogram(histogram)
    assert compact == ((11, 13), (15, 13))
    profile = histogram_profile(histogram)
    assert profile == expanded_profile(compact) and profile != 0
    # Basis order is 1,t1,t1t3,t2+t1t2.  h sees the last two basis bits.
    quartet = (2, 4, 8, 2)
    h = lambda a: (a & 12).bit_count() & 1
    assert tuple(h(a) for a in quartet) == (0, 1, 1, 0)
    assert binary_rank(a >> 1 for a in quartet) == 3
    assert pair_span(quartet[0], quartet[1]) != pair_span(quartet[2], quartet[3])
    return {'variables': 3, 'quartet_labels': quartet,
            'word_weights': [word.bit_count() for word in words],
            'joint_quotient_rank': 3, 'binomial_profile': compact,
            'expanded_polynomial_integer': profile}


def ceiling(value):
    return -(-value.numerator // value.denominator)


def fiber_masks(length):
    result = []
    for weight in range(1, length):
        for a in range(length):
            for b in range(length):
                if a == b:
                    continue
                mask = sum(1 << word for word in range(1 << length)
                           if word.bit_count() == weight and word >> a & 1
                           and not word >> b & 1)
                result.append(((weight, a, b), mask))
    return result


def exhaustive_fiber_bound():
    results = []
    for length in range(2, 5):
        harmonic = sum((Fraction(1, j) for j in range(1, length)), Fraction(0))
        denominator = 2*(length-1)*harmonic
        fibers = fiber_masks(length)
        constant_mask = 1 | (1 << ((1 << length)-1))
        checked, nontrivial = 0, 0
        for subset in range(1 << (1 << length)):
            nonconstants = (subset & ~constant_mask).bit_count()
            bound = ceiling(Fraction(nonconstants, 1)/denominator)
            selected = max((subset & mask for _, mask in fibers), key=int.bit_count)
            assert selected.bit_count() >= bound
            words = [word for word in range(1 << length) if selected >> word & 1]
            for f, g in combinations(words, 2):
                assert {(((f >> j) & 1) << 1) | ((g >> j) & 1)
                        for j in range(length)} == {0, 1, 2, 3}
            nontrivial += bound >= 2
            checked += 1
        results.append({'word_length': length, 'all_label_subsets': checked,
                        'harmonic_denominator': [denominator.numerator, denominator.denominator],
                        'nontrivial_lower_bounds': nontrivial})
    return results


def arbitrary_word_rectangle(symmetries):
    length = 8
    descriptor, selected = max(fiber_masks(length), key=lambda item: item[1].bit_count())
    words = [word for word in range(1 << length) if selected >> word & 1]
    xs, ys = words[:len(words)//2], words[len(words)//2:]
    assert len(words) == 20 and len(xs) == len(ys) == 10
    cache = {}
    cases = 0
    crossing_tensor = 0
    for index, quartet in enumerate(product(xs, ys, ys, xs)):
        x, y, u, v = quartet
        expected = (1 << (2*length+6)) if x == v and y == u else 0
        if expected:
            crossing_tensor |= 1 << index
        for permutation in symmetries:
            histogram = fast_histogram(length, tuple(quartet[j] for j in permutation))
            if histogram not in cache:
                profile = histogram_profile(histogram)
                assert profile == expanded_profile(binomial_profile_from_histogram(histogram))
                cache[histogram] = profile
            assert cache[histogram] == expected
            cases += 1
    ranks = tuple(binary_rank(matrix) for matrix in restricted_matrices(crossing_tensor, xs, ys))
    assert ranks == (100, 100)
    harmonic = sum((Fraction(1, j) for j in range(1, length)), Fraction(0))
    lower = ceiling(Fraction(254, 1)/(2*(length-1)*harmonic))
    assert len(words) >= lower
    weighted_cases = 0
    rng = random.Random(261032)
    for _ in range(4):
        weights = (1, *(rng.randrange(1, 5) for _ in range(length+1)))
        degree = 4+2*sum(weights[1:])
        for trial in range(32):
            x, y = rng.choice(xs), rng.choice(ys)
            quartet = ((x, y, y, x) if trial < 16 else
                       (x, y, rng.choice(ys), rng.choice(xs)))
            expected = (1 << degree) if crossing(quartet) else 0
            for permutation in symmetries:
                assert tensor_profile(length, tuple(quartet[j] for j in permutation), weights) == expected
                weighted_cases += 1
    return {'word_length': length, 'input_words': 256,
            'input_is_balanced': False, 'fiber_weight_and_coordinates': descriptor,
            'fiber_size': len(words), 'harmonic_guaranteed_fiber_size': lower,
            'rectangle_sizes': [len(xs), len(ys)],
            'source_profile_evaluations_including_symmetries': cases,
            'distinct_histograms_checked_by_both_algorithms': len(cache),
            'weighted_source_profile_checks': weighted_cases,
            'unique_coefficient_degree': 2*length+6,
            'both_restricted_crossing_ranks': ranks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    start = perf_counter()
    symmetries, table = source_support_table()
    variables, words, code = balanced_example()
    result = {'source_support': table, 'balanced_example': code,
              'rectangle': replay_rectangle(variables, words, symmetries),
              'harmonic_fiber_bound': exhaustive_fiber_bound(),
              'arbitrary_word_rectangle': arbitrary_word_rectangle(symmetries),
              'unbalanced_negative_control': negative_control(),
              'source_commit': '05d9c1b8aed4f8494ce230f0869f9048ef879474'}
    result['sha256'] = {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                        for name in ('check_word_restriction_barrier.py',
                                     'check_quadratic_crossing.py', 'quadratic_crossing.py')}
    result['seconds'] = perf_counter()-start
    text = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text)


if __name__ == '__main__':
    main()
