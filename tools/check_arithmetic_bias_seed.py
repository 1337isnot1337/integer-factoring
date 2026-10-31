#!/usr/bin/env python3
"""Exact angle-escape and character controls; no expander implementation claim."""
from fractions import Fraction as Q
import json

from arithmetic_bias_seed import powering_seed


def main():
    angles_checked = 0
    for denominator in range(2, 513):
        bound = (denominator-1).bit_length()
        for numerator in range(1, denominator):
            distances = [min(value, denominator-value)
                         for j in range(bound+1)
                         for value in [(numerator << j) % denominator]]
            assert any(4*distance >= denominator for distance in distances)
            angles_checked += 1
    character_controls = 0
    for dimension in range(1, 9):
        # An oversized order bound and redundant generators exercise the input
        # model; constructing the seed never uses these checker coordinates.
        generators = [1 << i for i in range(dimension)] + [(1 << dimension)-1, 0]
        bound = dimension+3
        words = powering_seed(len(generators), bound)
        assert len(words) == 4*len(generators)*(bound+1)
        assert max(sum(map(abs, word)) for word in words) == 1 << bound
        for character in range(1, 1 << dimension):
            signs = [1 if sum(value*(character & generator).bit_count()
                              for value, generator in zip(word, generators)) % 2 == 0 else -1
                     for word in words]
            eigenvalue = Q(sum(signs), len(signs))
            assert 0 <= eigenvalue <= 1-Q(2, len(words))
            character_controls += 1
    for numerator in range(1, 250):
        gap = Q(numerator, 1000)
        improved_gap = (1-Q(1, 100))*(2*gap-gap*gap)
        assert improved_gap >= Q(3, 2)*gap
    print(json.dumps({
        'scope': 'Exact powering-seed controls; published sparse expander amplification is not implemented here',
        'cyclic_character_angles_checked': angles_checked,
        'maximum_cyclic_denominator': 512,
        'elementary_two_group_characters_checked': character_controls,
        'maximum_elementary_dimension': 8,
        'scalar_gap_growth_controls': 249,
        'known_group_decomposition_used_by_constructor': False,
        'factoring_bound_improved': False,
    }, indent=2))


if __name__ == '__main__':
    main()
