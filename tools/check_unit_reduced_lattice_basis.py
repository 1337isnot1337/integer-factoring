#!/usr/bin/env python3
"""Exact unit-removal identities and full constraint-rowspace comparisons."""
import json
from pathlib import Path
import random

import weyl_lattice_filter_synthesis as original
from unit_reduced_lattice_basis import lattice_basis as reduced, original_lattice_basis


def canonical(pivots):
    return tuple((pivot, (1 << pivot) ^ original.reduce_row(row ^ (1 << pivot), pivots))
                 for pivot, row in sorted(pivots.items()))


def main():
    rng = random.Random(0x01C0DE)
    records = []
    comparisons = identities = witnesses = 0
    for n, P, geometry, modes in ((3, 1, 'initial', 2), (5, 1, 'initial', 3),
                                   (31, 1, 'initial', 3), (32, 1, 'initial', 3),
                                   (4, 3, 'mixed', 2)):
        old = original_lattice_basis(n, P, geometry, modes)
        new = reduced(n, P, geometry, modes)
        row_units = [1] * len(old['states'])
        column_units = [1] * len(old['states'])
        for kind, index, unit in new['removed_zero_units']:
            assert unit & 1
            factors = row_units if kind == 'row' else column_units
            factors[index] = original.polynomial_multiply(factors[index], unit)
        row_block = {block['row_offset'] + i: block['row_offset']
                     for block in old['blocks'] for i in range(len(block['indices']))}
        for first, second in zip(old['coefficients'], new['coefficients']):
            assert set(first) == set(second)
            for (i, j), value in first.items():
                restored = original.polynomial_multiply(second[i, j], row_units[i])
                restored = original.polynomial_multiply(restored, column_units[row_block[i] + j])
                assert restored == value
                identities += 1
        minimum = tuple(min(original.polynomial_order(value) for matrix in old['coefficients']
                            for (i, j), value in matrix.items() if i == row)
                        for row in range(len(old['states'])))
        for trial in range(8):
            degrees = tuple(rng.choice((0, 1, n * P, n * P + 1)) for state in old['states'])
            targets = tuple(rng.randrange(a, b + 2) for a, b in zip(minimum, old['gaussian_orders']))
            original.lattice_basis = original_lattice_basis
            original.coefficient_tables.cache_clear()
            first = original.homogeneous_system(n, P, geometry, degrees, targets, modes)
            first_witness = original.synthesize(n, P, geometry, degrees, targets, modes=modes)
            original.lattice_basis = reduced
            original.coefficient_tables.cache_clear()
            second = original.homogeneous_system(n, P, geometry, degrees, targets, modes)
            second_witness = original.synthesize(n, P, geometry, degrees, targets, modes=modes)
            assert first[0] == second[0] and first[2] == second[2]
            assert canonical(first[1]) == canonical(second[1])
            assert bool(first_witness) == bool(second_witness)
            if first_witness:
                assert first_witness['alpha'] == second_witness['alpha']
                witnesses += 1
            comparisons += 1
        records.append({'n': n, 'ramification': P, 'geometry': geometry, 'modes': modes,
                        'removed_zero_units': len(new['removed_zero_units']),
                        'maximum_source_degree_before': max(value.bit_length() - 1
                                                            for matrix in old['coefficients'] for value in matrix.values()),
                        'maximum_source_degree_after': max(value.bit_length() - 1
                                                           for matrix in new['coefficients'] for value in matrix.values()),
                        'source_terms_before': sum(value.bit_count() for matrix in old['coefficients'] for value in matrix.values()),
                        'source_terms_after': sum(value.bit_count() for matrix in new['coefficients'] for value in matrix.values())})
    original.lattice_basis = original_lattice_basis
    original.coefficient_tables.cache_clear()
    result = {'status': 'Exact common-unit removal preserves every tested constraint rowspace and literal witness.',
              'complete_source_polynomial_identities': identities,
              'full_constraint_rowspace_comparisons': comparisons,
              'identical_literal_witnesses': witnesses, 'cases': records}
    path = Path(__file__).resolve().parents[1] / 'evidence' / 'unit-reduced-lattice-basis.json'
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
