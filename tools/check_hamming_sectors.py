#!/usr/bin/env python3
"""Exact finite checks of the antipodal-subset crossing construction."""

import argparse
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path

from check_modular_recombination import (
    candidates, insert, ranks, replay, select,
)
from check_partition_channels import rank_mod
from hamming_pfaffian import matrices
from hamming_sectors import binary_sector_labels, sector_labels, sector_matrices

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    # AB alone must supply exact equality on all subsets, not just the prefix.
    gram_entries = 0
    for free_width in range(1, 5):
        pairs = 1 << free_width
        for count in range(1, min(4, pairs)+1):
            sets = list(combinations(range(pairs), count))
            if len(sets) > 100:
                sets = sets[:100]
            for left in sets:
                for right in sets:
                    forward = [[int(a & ~b == 0) for b in right] for a in left]
                    reverse = [[int(b & ~a == 0) for b in right] for a in left]
                    assert (rank_mod(forward, 2) == count and
                            rank_mod(reverse, 2) == count) == (left == right)
                    gram_entries += 1
    schedules = []
    for count in range(1, 7):
        for width in range(1, 9):
            free, domain = binary_sector_labels(width, count)
            assert len(domain) == 1 << width
            assert len(set(domain)) == len(domain)
            assert all(len(a) == 2*count for a in domain)
            assert comb(1 << free, count) >= 1 << width
            if free > 1:
                assert comb(1 << (free-1), count) < 1 << width
            assert count*free <= width+count*((count-1).bit_length()+1)
            schedules.append((width, count, free))
    records = []
    cases = [(2, 2, 4), (2, 2, None), (2, 3, None), (3, 2, 8)]
    if args.full:
        cases.append((3, 2, 16))
    for free, count, limit in cases:
        domain = sector_labels(free, count, limit)
        d = len(domain)
        source = sector_matrices(free, count, limit)
        profile = [(degree, *ranks(a, d)) for degree, a in sorted(source.items())]
        basis = {}
        for v in candidates(source, d, False):
            insert(basis, v, 0)
        record = {'free_width': free, 'pair_count': count, 'labels': domain,
                  'dimension': d, 'crossing_entries': d**4,
                  'span_dimension': len(basis), 'coefficient_pieces': len(source),
                  'ranks': profile, 'raw_max_product': max(r*c for k, r, c in profile)}
        if d <= 6:
            witness = select(candidates(source, d, True), d)
            replay(witness, source, d)
            record['pair_candidate_witness'] = witness
        records.append(record)
        print(json.dumps({k: v for k, v in record.items()
                          if k not in ('ranks', 'labels', 'pair_candidate_witness')}), flush=True)
    # Changing the label cardinality must not alter the original sector.
    for free in (1, 2, 3):
        assert sector_matrices(free, 1) == matrices(free)
    weighted = sector_matrices(2, 2, 4, (1, 2, 1, 2))
    assert weighted
    paths = [Path(__file__), ROOT/'tools/hamming_sectors.py', ROOT/'tools/hamming_pfaffian.py']
    result = {'status': 'effective crossing sectors; no improved asymptotic rank bound',
              'full': args.full, 'subset_gram_pairs': gram_entries,
              'binary_label_schedules': schedules, 'cases': records,
              'weighted_crossing_entries': 256,
              'sha256': {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
