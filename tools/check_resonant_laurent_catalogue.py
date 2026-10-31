#!/usr/bin/env python3
"""Exhaustive offdiagonal and full-entry controls for the Laurent grouping."""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path
import random

from f4_repaired_identity_filters import evaluate_one
from resonant_identity_operators import determinant, signature
from resonant_laurent_catalogue import (
    complete_catalogue, diagonal_signature_groups, offdiagonal_classes, polynomial,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = random.Random(31446080)
    decisions, admissible, records = 0, 0, []
    for lower in (-2, -1):
        pairs = {(0, 0), (0, 63), (63, 0), (63, 63),
                 (1 << (-2*lower), 1 << (-2*lower))}
        while len(pairs) < 16:
            pairs.add((rng.randrange(64), rng.randrange(64)))
        for first, second in sorted(pairs):
            direct = Counter()
            for a, b in product(range(64), repeat=2):
                matrix = ((polynomial(first, lower), polynomial(a, lower)),
                          (polynomial(b, lower), polynomial(second, lower)))
                decisions += 1
                if evaluate_one(determinant(matrix)):
                    direct[signature(4, 1, matrix)] += 1
                    admissible += 1
            grouped = Counter()
            for values, count, offdiagonal in diagonal_signature_groups(4, 1, lower, first, second):
                grouped[values] += count
                matrix = ((polynomial(first, lower), polynomial(offdiagonal[0], lower)),
                          (polynomial(offdiagonal[1], lower), polynomial(second, lower)))
                assert signature(4, 1, matrix) == values
            assert direct == grouped
        classes = complete_catalogue(32, 1, lower)
        for row in classes:
            a, b, c, d = row['entry_codes']
            matrix = ((polynomial(a, lower), polynomial(b, lower)),
                      (polynomial(c, lower), polynomial(d, lower)))
            assert signature(32, 1, matrix) == row['signature']
        records.append({'window_lower': lower, 'complete_offdiagonal_pairs': 4096,
                        'exact_offdiagonal_product_classes': len(offdiagonal_classes(lower)),
                        'complete_diagonal_pair_controls': len(pairs),
                        'full_catalogue_kernels': sum(row['count'] for row in classes),
                        'full_catalogue_signatures': len(classes),
                        'all_representative_signatures_replayed': True})
    result = {'status': 'Exact complete Laurent catalogue grouping with independent literal polynomial controls; no pressure claim.',
              'literal_matrix_endpoint_decisions': decisions,
              'literal_admissible_profile_signatures': admissible,
              'windows': records}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
