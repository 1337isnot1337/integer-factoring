#!/usr/bin/env python3
"""Complete six-bit quadratic-feedback-plus-zero-indicator crossing screen."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

from check_controlled_permutation_smith import binary_rank, phi_rows
from check_cyclic_algebra_crossings import jet_orders
from check_shift_register_crossings import encode, full_cycle, literal_jets
from cyclic_algebra_crossing import coefficient_rows, construct


MONOMIALS = [()] + [(i,) for i in range(5)] + list(combinations(range(5), 2))


def truth_tables():
    features = [sum(int(all(state >> j & 1 for j in term)) << state
                    for state in range(32)) for term in MONOMIALS]
    tables = [0] * (1 << len(features))
    for code in range(1, len(tables)):
        bit = code & -code
        tables[code] = tables[code ^ bit] ^ features[bit.bit_length()-1]
    return tables


def register(table):
    feedback = table ^ 1  # Add the indicator of the all-zero lower tail.
    return [((state << 1) & 63) | ((state >> 5) ^ ((feedback >> (state & 31)) & 1))
            for state in range(64)]


def partition_cycle(permutation, encoding):
    result = [0] * 64
    for state, target in enumerate(permutation):
        result[encoding[state]] = encoding[target]
    return result


def check():
    q, n = 8, 64
    tables = truth_tables()
    # Independent monomial evaluation verifies the entire truth-table map.
    for code in range(len(tables)):
        direct = 0
        for state in range(32):
            value = 0
            for j, term in enumerate(MONOMIALS):
                if code >> j & 1:
                    value ^= int(all(state >> bit & 1 for bit in term))
            direct |= value << state
        assert direct == tables[code]
    partitions = []
    for first in combinations(range(6), 3):
        if 0 not in first:
            continue
        second = tuple(i for i in range(6) if i not in first)
        encoding = [encode(s, first, second) for s in range(n)]
        partitions.append((first, second, encoding))

    accepted, histogram, minimum = [], Counter(), None
    for code, table in enumerate(tables):
        permutation = register(table)
        assert sorted(permutation) == list(range(n))
        if not full_cycle(permutation):
            continue
        accepted.append(code)
        for first, second, encoding in partitions:
            cycle = partition_cycle(permutation, encoding)
            a0 = phi_rows(cycle, q)
            r, k = binary_rank(a0), None
            floor = n*r
            if floor**6 < q**17:
                r_again, k, _, _ = literal_jets(cycle, q)
                assert r_again == r
                floor += (n-1)*k
            assert floor**6 >= q**17
            histogram[r, k] += 1
            if minimum is None or floor < minimum['exact_floor']:
                minimum = dict(code=code, first=list(first), second=list(second),
                               rank=r, order_one_pivots=k, exact_floor=floor)

    assert (len(accepted), len(partitions)) == (1952, 10)
    witness = minimum
    first, second = witness['first'], witness['second']
    encoding = [encode(s, first, second) for s in range(n)]
    cycle = partition_cycle(register(tables[witness['code']]), encoding)
    row, column, infinity = construct(q, 2, cycle)
    r, k, a0, a1 = literal_jets(cycle, q)
    assert coefficient_rows(infinity, 0) == a0
    assert coefficient_rows(infinity, 1) == a1
    orders = jet_orders(infinity, 16)
    assert (r, k) == (orders.count(0), orders.count(1))
    for x in range(q):
        for y in range(q):
            for u in range(q):
                for v in range(q):
                    value = row[q*x+y][q*u+v]
                    assert value.bit_count() % 2 == int(x == v and y == u)
                    assert column[q*x+u][q*y+v] == value
    return dict(functions=len(tables), full_cycles=len(accepted),
                balanced_partitions=len(accepted)*len(partitions),
                accepted_codes_sha256=sha256(json.dumps(accepted).encode()).hexdigest(),
                monomial_order=[list(term) for term in MONOMIALS],
                minimum_certificate=minimum,
                minimum_full_row_orders=orders, minimum_order_one_pivots=k,
                histogram=[dict(rank=r, order_one_pivots=k, count=c)
                           for (r, k), c in sorted(histogram.items(), key=lambda item: (item[0][0], -1 if item[0][1] is None else item[0][1]))],
                certified_rho_lower_bound='17/24',
                scope='All six-bit quadratic lower-tail feedbacks plus the zero-tail indicator, all balanced coordinate partitions modulo complement; bare two-place criterion.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ('functions', 'full_cycles', 'balanced_partitions', 'factoring_exponent')}))


if __name__ == '__main__':
    main()
