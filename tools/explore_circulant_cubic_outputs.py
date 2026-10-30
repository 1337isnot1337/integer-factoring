#!/usr/bin/env python3
"""Search circulant linear output maps for the cubic radical input code."""

import argparse
from functools import lru_cache
import hashlib
import json
from math import isqrt, sqrt
from pathlib import Path
import random

from check_modular_crossing import rank_two
from explore_complement_edge_pencils import minimum_decomposition
from explore_cyclic_radical_cubic import canonical_entry, circulant_rank


def general_entry(pairs, x, y, u, v, first_offset, second_offset):
    raw = (x, y, u, v, (x + first_offset) % pairs, (u + second_offset) % pairs)
    names = {}
    pattern = tuple(names.setdefault(value, len(names)) for value in raw)
    return canonical_entry(pattern, pairs % 2)


@lru_cache(None)
def forms(pairs):
    cases = [(0, 0, j, j) for j in range(pairs)] + [(0, j, 0, j) for j in range(pairs)]
    for delta in range(1, (pairs + 1) // 2):
        indices = ((0, delta), (delta, 0))
        cases += [(x, y, u, v) for x, y in indices for v, u in indices]
        cases += [(x, y, u, v) for x, u in indices for v, y in indices]
    result = []
    for x, y, u, v in cases:
        terms = []
        for a in range(1, pairs):
            value = general_entry(pairs, x, y, u, v, a, a)
            if value:
                terms.append((1 << (a - 1), value))
            for b in range(a + 1, pairs):
                value = (general_entry(pairs, x, y, u, v, a, b)
                         ^ general_entry(pairs, x, y, u, v, b, a))
                if value:
                    terms.append(((1 << (a - 1)) | (1 << (b - 1)), value))
        result.append(terms)
    return result


def profile(pairs, output_mask):
    assert output_mask.bit_count() % 2 == 1 and output_mask < 1 << (pairs - 1)
    values = []
    for form in forms(pairs):
        value = 0
        for support, polynomial in form:
            if output_mask & support == support:
                value ^= polynomial
        values.append(value)
    diagonal = [values[:pairs], values[pairs:2 * pairs]]
    kinds = []
    for i in range((pairs - 1) // 2):
        block = values[2 * pairs + 8 * i:2 * pairs + 8 * i + 8]
        kinds.append([[block[:2], block[2:4]], [block[4:6], block[6:8]]])
    for i, value in enumerate(diagonal[0]):
        assert value & 1 == int(i == 0)
    assert all(not (value >> 3) & 1 and not (value >> 12) & 1 for value in values)
    source = diagonal[0] + [v for cuts in kinds for row in cuts[0] for v in row]
    pivots, chosen = {}, []
    for degree in range(13):
        value = sum(((v >> degree) & 1) << i for i, v in enumerate(source))
        while value:
            lead = value.bit_length() - 1
            if lead in pivots:
                value ^= pivots[lead]
            else:
                pivots[lead] = value
                chosen.append(degree)
                break
    assert chosen[0] == 0
    both = []
    for mask in range(1 << len(chosen)):
        functional = sum(1 << degree for i, degree in enumerate(chosen) if mask >> i & 1)
        cut_ranks = []
        for cut in range(2):
            first_row = sum(((v & functional).bit_count() % 2) << i for i, v in enumerate(diagonal[cut]))
            rank = circulant_rank(first_row, pairs)
            for blocks in kinds:
                rows = [sum(((v & functional).bit_count() % 2) << j for j, v in enumerate(row))
                        for row in blocks[cut]]
                rank += pairs * rank_two(rows)
            cut_ranks.append(rank)
        both.append(cut_ranks)
    assert both[1] == [pairs * pairs, pairs * pairs]
    cost, witness = minimum_decomposition([sqrt(r * c) for r, c in both], target=1)
    lower, _ = minimum_decomposition([isqrt(r * c) for r, c in both], target=1)
    return {'pairs': pairs, 'output_mask': output_mask, 'basis_degrees': chosen,
            'cost': cost, 'integer_lower': lower, 'witness': witness,
            'witness_ranks': [both[i] for i in witness],
            'rank_table_sha256': hashlib.sha256(json.dumps(both, separators=(',', ':')).encode()).hexdigest(),
            'exclusion': lower ** 70 > pairs ** 99}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pairs', default='7,9,11')
    parser.add_argument('--limit', type=int, default=0)
    parser.add_argument('--seed', type=int, default=20261030)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = {'scope': 'Listed binary circulant output selectors and complete coefficient spaces',
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'cases': []}
    for pairs in map(int, args.pairs.split(',')):
        count = 1 << (pairs - 2)
        # Each lower (P-2)-bit word has one unique final bit making odd parity.
        indices = list(range(count)) if not args.limit or count <= args.limit else sorted(random.Random(args.seed).sample(range(count), args.limit))
        records = []
        for index in indices:
            output = index | ((1 ^ (index.bit_count() % 2)) << (pairs - 2))
            records.append(profile(pairs, output))
        case = {'pairs': pairs, 'selectors_available': count, 'selectors_checked': len(records),
                'minimum_integer_lower': min(r['integer_lower'] for r in records),
                'best': min(records, key=lambda r: r['cost']), 'records': records}
        result['cases'].append(case)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
        print(json.dumps({k: v for k, v in case.items() if k != 'records'}), flush=True)


if __name__ == '__main__':
    main()
