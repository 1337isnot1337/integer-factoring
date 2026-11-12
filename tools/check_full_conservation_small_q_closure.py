#!/usr/bin/env python3
"""Close the small-q full-conservation lanes by combining existing bounds."""
from __future__ import annotations

from hashlib import sha256
import json
from math import floor, sqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'evidence' / 'full-conservation-small-q-closure.json'
PRIMES = (5, 7, 11, 13, 17, 19)


def main() -> None:
    rows = []
    for q in PRIMES:
        target = q ** sqrt(2)
        floor_cost = q
        one_piece = q * q
        max_pieces = floor(target / floor_cost)
        two_piece = (sqrt(q * q - q + 1)
                     + sqrt((q * q - 1) * (q - 1)))
        if max_pieces <= 1:
            closed = one_piece > target and 2 * floor_cost > target
            mechanism = 'one-piece cost and two-piece floor'
        elif max_pieces == 2:
            closed = two_piece > target
            mechanism = 'floor caps at two pieces; all-field two-piece bound'
        else:
            closed = False
            mechanism = 'not closed by this composition: three pieces possible'
        closed_by_minimum_cut = q <= 21
        assert closed == (q in (5, 7, 11, 13))
        assert two_piece > q ** 1.5
        rows.append({
            'q': q,
            'target_q_sqrt2': target,
            'piece_floor_q': floor_cost,
            'one_piece_cost_q_squared': one_piece,
            'max_pieces_from_floor': max_pieces,
            'two_piece_lower_bound': two_piece,
            'closed_by_floor_and_two_piece': closed,
            'closed_by_minimum_cut': closed_by_minimum_cut,
            'mechanism': mechanism,
        })
    receipt = {
        'status': 'pass',
        'prime_labels': list(PRIMES),
        'closure_by_floor_and_two_piece': [row for row in rows
                                          if row['closed_by_floor_and_two_piece']],
        'not_closed_by_floor_composition': [row for row in rows
                                           if not row['closed_by_floor_and_two_piece']],
        'minimum_cut_cutoff': 21,
        'all_rows_closed_by_minimum_cut': all(row['closed_by_minimum_cut']
                                              for row in rows),
        'scope': ('Combines the full-conservation piece floor q (worker 4) with the '
                  'all-field two-piece crossing bound. For q=5 the floor excludes two '
                  'pieces and q^2 exceeds the target; for q=7,11,13 the floor caps a '
                  'crossing at two pieces and the two-piece bound exceeds the target. '
                  'The q=17 and q=19 rows are open for this floor composition only; '
                  'the existing minimum-cut theorem independently closes every '
                  'full-conservation seed with q<=21. All fields; no score change.'),
        'source_sha256': {
            'two_piece_note': sha256((ROOT / 'research' /
                                      'two-piece-rank-sum-obstruction.md').read_bytes()).hexdigest(),
        },
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True,
                      'closed_by_floor_q': [row['q'] for row in
                                            receipt['closure_by_floor_and_two_piece']],
                      'not_closed_by_floor_q': [row['q'] for row in
                                                receipt['not_closed_by_floor_composition']],
                      'q13_two_piece_lower': next(row['two_piece_lower_bound']
                                                  for row in rows if row['q'] == 13),
                      'q13_target': next(row['target_q_sqrt2']
                                         for row in rows if row['q'] == 13)}))


if __name__ == '__main__':
    main()
