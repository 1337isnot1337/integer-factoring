#!/usr/bin/env python3
"""Exact support witnesses excluding all F4 hyperoval orthogonal scalings."""
import argparse
from collections import Counter
from functools import reduce
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

from coherent_hot_compounds_f4 import mul
from orthogonal_differential_pairs_f4 import rank


def gram(A):
    return [[reduce(int.__xor__, (mul(x, y) for x, y in zip(row, other)), 0)
             for other in A] for row in A]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, default=Path('evidence/tensor-f4-hyperoval-obstruction.json'))
    args = ap.parse_args()
    columns = [(1, t, mul(t, t)) for t in range(4)] + [(0, 1, 0), (0, 0, 1)]
    C = [[column[i] for column in columns] for i in range(3)]
    for selected in combinations(range(6), 3):
        assert rank([[row[j] for j in selected] for row in C]) == 3
    H = [(a ^ b ^ c, a, b, c, a ^ mul(2, b) ^ mul(3, c),
          a ^ mul(3, b) ^ mul(2, c)) for a, b, c in product(range(4), repeat=3)]
    assert len(set(H)) == 64
    for h in H:
        assert all(reduce(int.__xor__, (mul(a, b) for a, b in zip(row, h)), 0) == 0 for row in C)
    weight_counts = Counter(sum(x != 0 for x in h) for h in H)
    assert weight_counts == {0: 1, 4: 45, 6: 18}
    zero_coordinate_covers = []
    for j in range(6):
        ambient = {h for h in H if h[j] == 0 and any(h)}
        lines = [{h for h in ambient if h[k] == 0} for k in range(6) if k != j]
        assert len(ambient) == 15 and all(len(line) == 3 for line in lines)
        assert all(not (a & b) for a, b in combinations(lines, 2))
        assert set().union(*lines) == ambient
        zero_coordinate_covers.append({'coordinate': j, 'nonzero_words': 15,
                                       'further_zero_lines': 5, 'words_per_line': 3})
    records, bad_support_histogram = [], Counter()
    for tail in product((1, 2, 3), repeat=5):
        scales = (1,) + tail
        A = [[mul(a, scales[j]) for j, a in enumerate(row)] for row in C]
        bad = {}
        for h in H:
            word = tuple(mul(a, mul(scale, scale)) for a, scale in zip(h, scales))
            support = tuple(i for i, value in enumerate(word) if value)
            if len(support) == 4 and reduce(int.__xor__, word, 0) == 0:
                bad.setdefault(support, word)
        assert len(bad) in (3, 5, 15)
        bad_support_histogram[len(bad)] += 1
        support, witness = min(bad.items())
        restricted = [[row[i] for i in support] for row in A]
        assert rank(restricted) == 3 and rank(gram(restricted)) == 2
        assert reduce(int.__xor__, (mul(x, x) for x in witness), 0) == 0
        assert all(reduce(int.__xor__, (mul(x, y) for x, y in zip(row, witness)), 0) == 0 for row in A)
        records.append({'scaling': scales, 'degenerate_four_subsets': len(bad),
                        'witness_support': support, 'isotropic_kernel_word': witness,
                        'column_rank': 3, 'gram_rank': 2})
    assert bad_support_histogram == {3: 165, 5: 72, 15: 6}
    report = {'scope': __doc__, 'passed': True, 'coefficient_matrix': C,
              'independent_triple_checks': 20, 'kernel_words': 64,
              'kernel_weight_counts': dict(weight_counts), 'zero_coordinate_covers': zero_coordinate_covers,
              'normalized_scalings': len(records), 'all_scalings_before_normalization': 3 * len(records),
              'degenerate_four_subset_histogram': dict(bad_support_histogram), 'witnesses': records,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest()
                                for n in (Path(__file__).name, 'coherent_hot_compounds_f4.py',
                                          'orthogonal_differential_pairs_f4.py')}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('witnesses', 'source_sha256', 'zero_coordinate_covers')}))


if __name__ == '__main__':
    main()
