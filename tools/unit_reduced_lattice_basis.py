#!/usr/bin/env python3
"""Remove exact common zero units from a complete polynomial source basis.

Each division is a regular invertible diagonal row/column change at zero.
It preserves every row-divisibility target, as well as all Smith profiles.
The original full coefficient polynomials are divided exactly; no series
coefficient is truncated.
"""
from functools import lru_cache

from weyl_lattice_filter_synthesis import lattice_basis as original_lattice_basis


def divide(first, second):
    quotient = 0
    while first and first.bit_length() >= second.bit_length():
        shift = first.bit_length() - second.bit_length()
        quotient ^= 1 << shift
        first ^= second << shift
    return quotient, first


def gcd(first, second):
    while second:
        first, second = second, divide(first, second)[1]
    return first


def common_unit(values):
    values = tuple(values)
    if not values:
        return 1
    result = values[0]
    for value in values[1:]:
        result = gcd(result, value)
        if result & (result - 1) == 0:
            return 1
    order = (result & -result).bit_length() - 1
    return result >> order


@lru_cache(maxsize=32)
def lattice_basis(n, ramification=1, geometry='initial', modes=2):
    source = original_lattice_basis(n, ramification, geometry, modes)
    coefficients = tuple(dict(value) for value in source['coefficients'])
    changes = []
    while True:
        before = len(changes)
        for row in range(len(source['states'])):
            unit = common_unit(value for matrix in coefficients
                               for (i, j), value in matrix.items() if i == row)
            if unit == 1:
                continue
            assert unit & 1
            for matrix in coefficients:
                for entry in tuple(matrix):
                    if entry[0] == row:
                        quotient, remainder = divide(matrix[entry], unit)
                        assert not remainder and quotient
                        matrix[entry] = quotient
            changes.append(('row', row, unit))
        for block in source['blocks']:
            offset, size = block['row_offset'], len(block['indices'])
            for column in range(size):
                unit = common_unit(value for matrix in coefficients for (i, j), value in matrix.items()
                                   if offset <= i < offset + size and j == column)
                if unit == 1:
                    continue
                assert unit & 1
                for matrix in coefficients:
                    for entry in tuple(matrix):
                        if offset <= entry[0] < offset + size and entry[1] == column:
                            quotient, remainder = divide(matrix[entry], unit)
                            assert not remainder and quotient
                            matrix[entry] = quotient
                changes.append(('column', offset + column, unit))
        if len(changes) == before:
            break
    result = dict(source)
    result['coefficients'] = coefficients
    result['removed_zero_units'] = tuple(changes)
    return result


def install():
    """Use the equivalent smaller source within the current Python process."""
    import weyl_lattice_filter_synthesis as lattice
    lattice.lattice_basis = lattice_basis
    lattice.coefficient_tables.cache_clear()


def annotate_output():
    import json
    from pathlib import Path
    import sys
    if '--output' not in sys.argv:
        return
    path = Path(sys.argv[sys.argv.index('--output') + 1])
    record = json.loads(path.read_text())
    data = lattice_basis(record['n'], record['ramification'], record['geometry'], record['modes'])
    record['unit_reduced_reference_basis'] = True
    record['removed_reference_zero_units'] = len(data['removed_zero_units'])
    path.write_text(json.dumps(record, separators=(',', ':')) + '\n')


if __name__ == '__main__':
    install()
    from search_weyl_lattice_filters import main
    main()
    annotate_output()
