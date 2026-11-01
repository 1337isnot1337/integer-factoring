#!/usr/bin/env python3
"""Exact controls for moving unequal-height cells and stable Gibbs roots."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_coherent_hot_boundary import logical
from check_moving_boundary_crossings import reference_compress
from moving_boundary_crossings import compress, reorder
from restrict_implicit_reflection import restricted_tensor
from search_moving_interior_crossings import pressure

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raw = restricted_tensor(16, list(map(logical, [0, 1, 2, 3, 4, 5, 8, 6, 9])))
    reordered_entries = contracted_entries = 0
    for pair in ((0, 1), (1, 5), (3, 7)):
        order = [i for i in range(9) if i not in pair] + list(pair)
        changed = reorder(raw, order)
        for index, value in enumerate(changed['values']):
            labels = [(index // 9**axis) % 9 for axis in range(4)]
            source = sum(order[label] * 9**axis for axis, label in enumerate(labels))
            assert value == raw['values'][source]
            reordered_entries += 1
        a, b = (3, 2), (2, 1)
        actual = compress(changed, 7, a, b)
        assert actual['values'] == reference_compress(changed, 7, a, b)
        contracted_entries += len(actual['values'])
    # The unnormalized logsumexp would lose exp(-150) in this example.
    base = pressure([[0, 100], [0, 100]], [0, 100])
    cut_shift = pressure([[-100000, -99900], [100000, 100100]], [0, 100])
    physical_shift = pressure([[300, 400], [300, 400]], [300, 400])
    assert abs(base['rho'] - .5) < 1e-12
    assert base == cut_shift == physical_shift
    sources = ['tools/check_moving_interior_crossings.py', 'tools/moving_boundary_crossings.py',
               'tools/search_moving_interior_crossings.py', 'tools/check_moving_boundary_crossings.py']
    report = {'scope': 'Independent full polynomial contractions after unequal-height relabeling and exact scalar-shift invariance of the numerical root; no improved bound.',
              'literal_relabeling_entries': reordered_entries,
              'independent_polynomial_entries': contracted_entries,
              'stable_root_control': base, 'cut_shift': 100000, 'physical_shift': 300,
              'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in sources}}
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: report[key] for key in ('literal_relabeling_entries', 'independent_polynomial_entries', 'stable_root_control')}))


if __name__ == '__main__':
    main()
