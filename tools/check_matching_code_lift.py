#!/usr/bin/env python3
"""Exact checks for the uniform-set-partition multigraded crossing lift.

The rank profiles are finite certificates. In particular, the a=4 gain does
not continue at a=6; this checker does not certify a scalable contraction.
"""

import argparse
import json
from itertools import product

from check_modular_crossing import determinant_two, normalized_piece, rank_two
from fermionic_piece_polynomial import monomial_piece_polynomial
from matching_code_lift import (
    decode, labels, lifted_polynomial, modes, orbit_tables, specialize,
)


def translated_dual(label):
    """Explicit Lucas expansion, reduced in the valid operator-label quotient."""
    possibilities = []
    for block in label:
        choices = []
        subset = block
        while True:
            choices.append(1 + 2 * subset)
            if subset == 0:
                break
            subset = (subset - 1) & block
        possibilities.append(choices)
    terms = set()
    for choice in product(*possibilities):
        if len(set(choice)) < len(choice):
            continue
        monomial = tuple(sorted(choice))
        if monomial in terms:
            terms.remove(monomial)
        else:
            terms.add(monomial)
    return tuple(sorted(terms))


def unfolded_polynomial(x, v, u, y):
    """Independent folding check: explicitly sum monomial Wick polynomials."""
    weight = sum(modes(x))
    answer = 0
    for term in translated_dual(x):
        polynomial, offset = monomial_piece_polynomial(term, modes(v), modes(u), modes(y))
        answer ^= polynomial << (offset + 2 * weight)
    return answer


def direct_canonical_tables(vertices, block_size):
    labs = labels(vertices, block_size)
    dimension = len(labs)
    matrices = {}
    for xi, x in enumerate(labs):
        for yi, y in enumerate(labs):
            for vi, v in enumerate(labs):
                for ui, u in enumerate(labs):
                    polynomial = lifted_polynomial(vertices, x, v, u, y)
                    while polynomial:
                        bit = polynomial & -polynomial
                        polynomial ^= bit
                        low, high = decode(vertices, bit.bit_length() - 1)
                        if tuple(sorted(high)) == high:
                            matrix = matrices.setdefault((low, high), [0] * dimension ** 2)
                            matrix[xi * dimension + yi] |= 1 << (vi * dimension + ui)
    return matrices


def matching_partner(label, vertices):
    partner = [0] * vertices
    for block in label:
        endpoints = [i for i in range(vertices) if (block >> i) & 1]
        assert len(endpoints) == 2
        first, second = endpoints
        partner[first], partner[second] = second, first
    return tuple(partner)


def cycle_type(first, second):
    seen = 0
    lengths = []
    for start in range(len(first)):
        if (seen >> start) & 1:
            continue
        vertex, length = start, 0
        while True:
            seen |= (1 << vertex) | (1 << first[vertex])
            vertex = second[first[vertex]]
            length += 1
            if vertex == start:
                break
        lengths.append(length)
    return tuple(sorted(lengths))


def gram_cycle_value(first, second):
    nontrivial = [length for length in cycle_type(first, second) if length > 1]
    if not nontrivial:
        return int(len(first) // 2 % 2 == 0)
    return int(len(nontrivial) == 1 and nontrivial[0] % 2 == 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output')
    args = parser.parse_args()
    checks = {'unfolded_specializations': 0, 'literal_source_entries_both_orders': 0,
              'direct_orbit_matrix_entries': 0, 'gram_determinant_entries': 0}

    def verify(vertices, x, v, u, y, polynomial):
        assert specialize(vertices, polynomial) == unfolded_polynomial(x, v, u, y)
        checks['unfolded_specializations'] += 1

    profiles = []
    for vertices, block_size in ((2, 2), (4, 2), (6, 2), (6, 3)):
        report, tables = orbit_tables(vertices, block_size, verify_entry=verify)
        if vertices <= 4:
            direct = direct_canonical_tables(vertices, block_size)
            assert tables == direct
            checks['direct_orbit_matrix_entries'] += len(tables) * report['dimension'] ** 4
        profiles.append(report)
        print(json.dumps({key: report[key] for key in
                          ('vertices', 'block_size', 'dimension', 'degree_orbits',
                           'maximum_rank_product', 'worst')}), flush=True)

    one_label = (modes(labels(2)[0]),)
    weight = sum(one_label[0])
    polynomial = specialize(2, lifted_polynomial(2, *([labels(2)[0]] * 4)))
    for shift in range(-2 * weight, 2 * weight + 1):
        coefficient = (polynomial >> (shift + 2 * weight)) & 1
        for exchanged in (False, True):
            assert normalized_piece(2, weight, shift, one_label, exchanged)[0][0] == coefficient
            checks['literal_source_entries_both_orders'] += 1

    gram_reports = []
    for vertices in (2, 4, 6, 8):
        labs = labels(vertices)
        partners = [matching_partner(label, vertices) for label in labs]
        rows = []
        for x, first in zip(labs, partners):
            row = 0
            for j, (y, second) in enumerate(zip(labs, partners)):
                matrix = [sum(int(not (e & f)) << column for column, f in enumerate(y))
                          for e in x]
                value = determinant_two(matrix)
                assert value == gram_cycle_value(first, second)
                row |= value << j
                checks['gram_determinant_entries'] += 1
            rows.append(row)
        gram_reports.append({'vertices': vertices, 'dimension': len(labs),
                             'rank': rank_two(rows)})

    assert profiles[1]['maximum_rank_product'] == 18
    assert profiles[2]['maximum_rank_product'] == 32574
    assert profiles[3]['maximum_rank_product'] == 2116
    evidence = {
        'scope': 'exact finite profiles and algebraic-lift checks; no asymptotic contraction claim',
        'source_revision': '16a5ce41304fcdaea4e0be7c132e1c361dc26da5',
        'folded_formula_prior_revision': 'c12da897fcc6081685e8c663217ce9083db87b4a',
        'checks': checks, 'profiles': profiles, 'homogeneous_gram': gram_reports,
    }
    if args.output:
        with open(args.output, 'w') as handle:
            json.dump(evidence, handle, indent=2)
            handle.write('\n')
    print(json.dumps({'checks': checks, 'homogeneous_gram': gram_reports}), flush=True)


if __name__ == '__main__':
    main()
