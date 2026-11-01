#!/usr/bin/env python3
"""Check graded quotient actions and the dyadic zero-cut factorization.

The one-particle polynomial identity is checked without Smith elimination.
The small complete quotient family is replayed and compared with its saved
complete profiles. The dyadic theorem applies to arbitrary larger inputs;
the finite checks do not substitute for its algebraic proof.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from check_cancellation_module_automata import pmul
from coherent_graded_quotients import parameter_images, quotient_maps
from coherent_hot_compounds import support
from coherent_hot_mixed import profiles
from search_coherent_hot_retractions import _atomic_write


def cross_entry(a, b, side_a, side_b):
    if side_a == side_b:
        return 3 << (2 * (a | b) + 1)
    return (1 << (2 * a + 1)) if a == b else 0


def column_coordinates(k, b, side):
    if b < k:
        return [int(index == side * k + b) for index in range(2 * k)]
    c, high = b % k, b - b % k
    shift = 2 * (high - k)
    result = []
    for target_side in (0, 1):
        for i in range(k):
            if target_side == side:
                result.append(((1 | (1 << (2 * k))) << shift) if i == c else 0)
            else:
                result.append(3 << (shift + 2 * (c & ~i)))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'evidence/coherent-graded-quotients-4-8.json')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    counts = Counter()
    for k in (1, 2, 4, 8):
        for b in range(32):
            for side in (0, 1):
                coordinates = column_coordinates(k, b, side)
                for row_side in (0, 1):
                    for a in range(k):
                        value = 0
                        for index, coefficient in enumerate(coordinates):
                            target_side, c = divmod(index, k)
                            value ^= pmul(cross_entry(a, c, row_side, target_side), coefficient)
                        assert value == cross_entry(a, b, row_side, side)
                        counts['one_particle_polynomial_identities'] += 1
    saved = json.loads(args.catalogue.read_text())
    assert saved['completed'] and saved['finite_list_exhausted']
    assert saved['virtual_modes'] == 4 and saved['physical_modes'] == 8
    assert len(saved['records']) == 16
    baseline = saved['records'][0]['profile']['zero_orders_t']
    records = []
    for record in saved['records']:
        mask = record['mask']
        images = parameter_images(4, 8, mask)
        code = quotient_maps(images, 4)
        assert code == record['code']
        embedding = code['embedding']
        for a in range(8):
            for virtual_word, physical_vector in enumerate(embedding):
                deleted = 0
                for word in support(physical_vector):
                    if word >> a & 1:
                        deleted ^= 1 << (word ^ (1 << a))
                represented = 0
                for monomial in support(images[a]):
                    if monomial & virtual_word == monomial:
                        represented ^= embedding[virtual_word ^ monomial]
                assert deleted == represented
                counts['induced_differential_action_identities'] += 1
        current = profiles(embedding, code['retraction'], 8)
        for key in ('zero_orders_t', 'infinity_orders_t'):
            assert current[key] == record['profile'][key]
            counts['complete_profiles_replayed'] += 1
            counts['complete_smith_orders_replayed'] += len(current[key])
        assert current['zero_orders_t'] == baseline
        records.append({'mask': mask, 'zero_profile_equals_bare_full_fock': True,
                        'infinity_profile_equals_bare_full_fock':
                            current['infinity_orders_t'] == saved['records'][0]['profile']['infinity_orders_t']})
    names = (Path(__file__).name, 'coherent_graded_quotients.py', 'coherent_hot_mixed.py',
             'coherent_hot_mixed.cpp', 'coherent_hot_compounds_fast.cpp',
             'coherent_hot_compounds.py', 'check_cancellation_module_automata.py',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp')
    report = {'status': 'pass', 'scope': __doc__, 'counts': dict(counts), 'records': records,
              'catalogue_sha256': sha256(args.catalogue.read_bytes()).hexdigest(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}, 'factoring_exponent_improved': False}
    _atomic_write(args.output, report)
    print(json.dumps({'status': 'pass', 'counts': dict(counts)}, indent=2))


if __name__ == '__main__':
    main()
