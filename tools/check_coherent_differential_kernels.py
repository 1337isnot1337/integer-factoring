#!/usr/bin/env python3
"""Exhaustive small differential-kernel controls and literal q10 cut replay."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
import json
from pathlib import Path

from check_coherent_hot_compounds import reverse_map
from coherent_differential_kernels import differential_kernel
from coherent_hot_compounds import height, support
from coherent_hot_compounds_fast import profiles
from eight_letter_reflection_retraction import homogeneous_restriction


def literal_action(vector, operator):
    answer = 0
    for word in range(16):
        if vector >> word & 1:
            for removed in range(16):
                if operator >> removed & 1 and word & removed == removed:
                    answer ^= 1 << (word ^ removed)
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cells = defaultdict(list)
    for word in range(1, 16):
        cells[height(word), word.bit_count()].append(word)
    operators = [sum(1 << word for i, word in enumerate(words) if mask >> i & 1)
                 for words in cells.values() for mask in range(1, 1 << len(words))]
    counts = Counter()
    for operator in operators:
        for cap in range(17):
            code = differential_kernel([operator], 4, cap)
            words = [word for word in range(16) if height(word) <= cap]
            columns = [literal_action(1 << word, operator) for word in words]
            valid = 0
            # Enumerate the complete ambient vector space using Gray updates.
            previous, value, image = 0, 0, 0
            for index in range(1 << len(words)):
                gray = index ^ (index >> 1)
                if index:
                    changed = (gray ^ previous).bit_length() - 1
                    value ^= 1 << words[changed]
                    image ^= columns[changed]
                previous = gray
                counts['ambient_vectors_enumerated'] += 1
                if image == 0:
                    valid += 1
                    rebuilt = 0
                    for p, l in zip(code['embedding'], code['retraction']):
                        if (value & l).bit_count() % 2:
                            rebuilt ^= p
                    assert rebuilt == value
                    counts['kernel_vectors_reconstructed'] += 1
            assert valid == 1 << code['q']
            counts['complete_kernel_spaces'] += 1
    operator = (1 << 6) ^ (1 << 9)
    code = differential_kernel([operator], 4, 16)
    assert code['q'] == 10
    records = []
    for transposed in (False, True):
        p, l = code['embedding'], code['retraction']
        if transposed:
            p, l = l, p
        result = profiles(p, l, 4)
        literal = homogeneous_restriction(tuple(map(reverse_map, p)), tuple(map(reverse_map, l)),
                                           [height(support(value)[0]) for value in p])
        assert result['zero_orders_t'] == list(literal['zero_orders'])
        assert result['infinity_orders_t'] == list(literal['infinity_orders'])
        counts['complete_literal_profile_comparisons'] += 2
        counts['literal_endpoint_entries'] += literal['literal_endpoint_entries']
        counts['literal_infinity_restoration_entries'] += literal['literal_infinity_restoration_entries']
        records.append({'transposed': transposed, 'profile': result,
                        'literal_raw_degree': literal['actual_polynomial_degree']})
    names = (Path(__file__).name, 'coherent_differential_kernels.py', 'coherent_inverse_systems.py',
             'coherent_hot_compounds.py', 'coherent_hot_compounds_fast.py',
             'coherent_hot_compounds_fast.cpp', 'eight_letter_reflection_retraction.py',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp')
    result = {'status': 'pass', 'scope': __doc__, 'counts': dict(counts), 'records': records,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}, 'factoring_exponent_improved': False}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': 'pass', 'counts': dict(counts)}, indent=2))


if __name__ == '__main__':
    main()
