#!/usr/bin/env python3
"""Exact full-cycle shift-register crossings and bounded two-jet exclusions.

The optional larger catalogue enumerates every maximal binary feedback
rule and every balanced physical bit partition, modulo complementary cuts.
It does not infer a factoring bound from finite spectra.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

from check_controlled_permutation_smith import binary_rank, phi_rows, toeplitz_rank
from check_cyclic_algebra_crossings import jet_orders
from cyclic_algebra_crossing import coefficient_rows, construct


def maximal_registers(bits):
    size = 1 << bits
    for taps in range(size//2, size):
        permutation = [((s << 1) & (size-1)) | ((s & taps).bit_count() & 1)
                       for s in range(size)]
        seen, current = set(), 1
        while current not in seen:
            seen.add(current)
            current = permutation[current]
        if len(seen) == size-1:
            yield taps, permutation


def encode(state, first, second):
    a = sum(((state >> bit) & 1) << j for j, bit in enumerate(first))
    b = sum(((state >> bit) & 1) << j for j, bit in enumerate(second))
    return (a << len(second)) | b


def joined_cycle(permutation, first, second):
    size = len(permutation)
    joined = list(permutation)
    joined[0], joined[size//2] = joined[size//2], joined[0]
    result = [0] * size
    for source, target in enumerate(joined):
        result[encode(source, first, second)] = encode(target, first, second)
    return result


def full_cycle(permutation):
    seen, current = set(), 0
    while current not in seen:
        seen.add(current)
        current = permutation[current]
    return current == 0 and len(seen) == len(permutation)


def linear_exponent(permutation, first, second):
    """Sum of the two diagonal block ranks of the binary linear map."""
    first_mask = sum(1 << i for i in first)
    second_mask = sum(1 << i for i in second)
    return (binary_rank([permutation[1 << i] & first_mask for i in first])
            + binary_rank([permutation[1 << i] & second_mask for i in second]))


def literal_jets(cycle, alphabet):
    a0 = phi_rows(cycle, alphabet)
    a2 = phi_rows([cycle[j] for j in cycle], alphabet)
    a1 = [a ^ b for a, b in zip(a0, a2)]
    r = binary_rank(a0)
    k = binary_rank(a0 + [a ^ (b << len(cycle)) for a, b in zip(a1, a0)]) - 2*r
    assert k >= 0
    return r, k, a0, a1


def catalogue(bits):
    q, size = 1 << (bits//2), 1 << bits
    count, registers, histogram, minimum = 0, 0, Counter(), None
    for taps, permutation in maximal_registers(bits):
        registers += 1
        for first in combinations(range(bits), bits//2):
            if 0 not in first:
                continue
            second = tuple(i for i in range(bits) if i not in first)
            cycle = joined_cycle(permutation, first, second)
            assert full_cycle(cycle)
            a0 = phi_rows(cycle, q)
            r, k = binary_rank(a0), None
            # Every case uses an exact rational exponent comparison.
            floor = size*r
            if floor**6 < q**17:
                r_again, k, _, _ = literal_jets(cycle, q)
                assert r_again == r
                floor += (size-1)*k
            assert floor**6 >= q**17
            exponent = linear_exponent(permutation, first, second)
            linear_cycle = [0] * size
            for source, target in enumerate(permutation):
                linear_cycle[encode(source, first, second)] = encode(target, first, second)
            assert binary_rank(phi_rows(linear_cycle, q)) == 1 << exponent
            assert abs(r-(1 << exponent)) <= 1
            histogram[r, k] += 1
            candidate = dict(taps=taps, first=list(first), second=list(second),
                             rank=r, order_one_pivots=k, exact_floor=floor)
            if minimum is None or floor < minimum['exact_floor']:
                minimum = candidate
            count += 1
    return dict(bits=bits, alphabet=q, maximal_registers=registers, partitions=count,
                minimum_certificate=minimum,
                histogram=[dict(rank=r, order_one_pivots=k, count=c)
                           for (r, k), c in sorted(histogram.items(), key=lambda item: (item[0][0], -1 if item[0][1] is None else item[0][1]))],
                certified_rho_lower_bound='17/24')


def controls():
    records = []
    for bits in (4, 6):
        taps, permutation = next(maximal_registers(bits))
        q, n = 1 << (bits//2), 1 << bits
        first = tuple(range(bits//2)) if bits == 4 else (0, 1, 3)
        second = tuple(i for i in range(bits) if i not in first)
        cycle = joined_cycle(permutation, first, second)
        assert full_cycle(cycle)
        row, _, infinity = construct(q, 2, cycle)
        for x in range(q):
            for y in range(q):
                for u in range(q):
                    for v in range(q):
                        assert row[q*x+y][q*u+v].bit_count() % 2 == int(x == v and y == u)
        r, k, a0, a1 = literal_jets(cycle, q)
        assert coefficient_rows(infinity, 0) == a0
        assert coefficient_rows(infinity, 1) == a1
        assert toeplitz_rank([a0, a1], 2) == 2*r+k
        # A representative complete row profile independently checks all later jets.
        orders = jet_orders(infinity, 16)
        assert (r, k) == (orders.count(0), orders.count(1))
        records.append(dict(bits=bits, taps=taps, first=list(first),
                            row_infinity_orders=orders, endpoint_entries=q**4))
    return records


def focused_twelve():
    bits, q, size = 12, 64, 4096
    scanned, registers, records = 0, 0, []
    for taps, permutation in maximal_registers(bits):
        registers += 1
        squared = [permutation[s] for s in permutation]
        for first in combinations(range(bits), bits//2):
            if 0 not in first:
                continue
            second = tuple(i for i in range(bits) if i not in first)
            scanned += 1
            a = linear_exponent(permutation, first, second)
            if a > 4:
                continue
            b = linear_exponent(squared, first, second)
            if b > 4:
                continue
            cycle = joined_cycle(permutation, first, second)
            assert full_cycle(cycle)
            r, k, _, _ = literal_jets(cycle, q)
            floor = size*r+(size-1)*k
            assert floor**6 >= q**17
            records.append(dict(taps=taps, first=list(first), second=list(second),
                                linear_exponents=[a, b], rank=r,
                                order_one_pivots=k, exact_floor=floor))
    assert (registers, scanned, len(records)) == (144, 66528, 138)
    return dict(maximal_registers=registers, partitions_scanned=scanned,
                selection='Both diagonal-block rank sums for L and L^2 are at most four.',
                selected_partitions=len(records), records=records,
                certified_rho_lower_bound='17/24',
                scope='Only the selected 138 systems; the other 66390 partitions are not excluded here.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-bits', type=int, choices=(4, 6, 8, 10), default=6)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--focused-twelve', action='store_true')
    args = parser.parse_args()
    result = dict(controls=controls(), catalogues=[catalogue(b) for b in range(4, args.max_bits+1, 2)],
                  source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  scope='All maximal feedbacks and balanced coordinate partitions at the listed sizes; bare two-place certificate only.',
                  factoring_exponent='1/2; unchanged')
    if args.focused_twelve:
        result['focused_twelve'] = focused_twelve()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(partitions=sum(c['partitions'] for c in result['catalogues']),
                          maximal_registers=sum(c['maximal_registers'] for c in result['catalogues']),
                          factoring_exponent=result['factoring_exponent'])))


if __name__ == '__main__':
    main()
