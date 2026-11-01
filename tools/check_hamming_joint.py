#!/usr/bin/env python3
"""Replay a joint-sector crossing witness; no optimum or rate claim."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_hamming_pfaffian import independent_ranks
from check_modular_recombination import candidates, crossing, encode, insert, permute, ranks
from hamming_sectors import binary_sector_labels, sector_matrices

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = json.loads((ROOT/'evidence/hamming-joint.json').read_text())
    width = result['binary_width']
    d = 1 << width
    sources, basis, dimensions = {}, {}, []
    for count in result['pair_counts']:
        free, domain = binary_sector_labels(width, count)
        source = sector_matrices(free, count, d)
        sources[count] = source
        for vector in candidates(source, d, False):
            insert(basis, vector, 0)
        dimensions.append([count, len(basis)])
    assert dimensions == result['successive_span_dimensions']
    assert len(basis) == result['span_dimension'] == 33
    generators = list(basis.values())
    total = 0
    costs = []
    for piece in result['witness']:
        rows = [0]*(d*d)
        for term in piece['terms']:
            image = permute(sources[term['pair_count']][term['degree']],
                            d, term['permutation'])
            rows = [a ^ b for a, b in zip(rows, image)]
        vector = encode(rows, d)
        mask = int(piece['basis_coefficient_mask'], 16)
        expected = 0
        for i, (other, description) in enumerate(generators):
            if (mask >> i) & 1:
                expected ^= other
        assert vector == expected
        actual = independent_ranks(vector, d)
        assert actual == ranks(rows, d) == (piece['row_rank'], piece['column_rank'])
        total ^= vector
        costs.append(actual[0]*actual[1])
    assert total == encode(crossing(d), d)
    assert max(costs) == result['maximum_rank_product'] == 25
    assert len(costs) == result['witness_piece_count'] == 6
    result['maximum_binary_width_sum'] = max(
        (p['row_rank']-1).bit_length()+(p['column_rank']-1).bit_length()
        for p in result['witness'])
    result['crossing_entries_checked'] = d**4
    paths = [Path(__file__), ROOT/'tools/hamming_sectors.py', ROOT/'tools/hamming_pfaffian.py']
    result['sha256'] = {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('witness', 'sha256')}))


if __name__ == '__main__':
    main()
