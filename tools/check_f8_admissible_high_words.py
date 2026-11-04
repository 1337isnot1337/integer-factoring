#!/usr/bin/env python3
"""Enumerate the single-word extensions of the q96 F8 hyperoval map."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coherent_f8_maps import basis, in_span, validate_maps
from coherent_hot_compounds import height

INPUT = Path('evidence/f8-admissible-high-words.json')
BASE = Path('evidence/f8-hyperoval-q96-full-profile.json')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=INPUT)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()

    saved = json.loads(args.input.read_text())
    assert saved['completed'] is True
    assert saved['base_profile_sha256'] == sha256(BASE.read_bytes()).hexdigest()
    base = json.loads(BASE.read_text())
    embedding = [{int(k): v for k, v in row.items()} for row in base['embedding']]
    retraction = [{int(k): v for k, v in row.items()} for row in base['retraction']]
    left = basis(embedding)
    right = basis(retraction)

    rows = []
    for word in range(1, 1 << 12):
        vector = {word: 1}
        if in_span(vector, left):
            continue
        if any(row.get(word, 0) for row in retraction):
            continue
        if any(row.get(word, 0) for row in embedding):
            continue
        if any(not in_span({word ^ (1 << bit): 1}, left)
               or not in_span({word ^ (1 << bit): 1}, right)
               for bit in range(12) if word >> bit & 1):
            continue
        try:
            validate_maps(embedding + [vector], retraction + [vector], 12)
        except Exception:
            continue
        rows.append({'word': word, 'height': height(word), 'charge': word.bit_count()})
    rows.sort(key=lambda row: (row['height'], row['word']))

    assert rows == saved['admissible_single_words']
    assert len(rows) == saved['count'] == 45
    assert min(row['height'] for row in rows) == saved['min_height'] == 26
    assert max(row['height'] for row in rows) == saved['max_height'] == 44
    assert all(row['height'] >= 26 for row in rows)
    words = {row['word'] for row in rows}
    assert {160, 192} <= words
    assert 224 not in words
    assert saved['q99_extension_words'] == [160, 192, 224]

    support = set()
    for vector in embedding + retraction:
        support.update(vector)
    low = [word for word in range(1, 1 << 12) if height(word) <= 25]
    free = [word for word in low if word not in support]
    assert saved['low_height_words_total'] == len(low)
    assert saved['low_height_words_in_support'] == len(low) - len(free)
    assert saved['free_low_words'] == free == []
    assert len(low) == saved['low_height_words_total']

    receipt = {
        'passed': True,
        'count': len(rows),
        'min_height': 26,
        'max_height': 44,
        'no_low_single_word_extension': True,
        'low_height_support_saturated': True,
        'q99_pair_words_admissible': True,
        'q99_triple_requires_pairs': True,
        'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
        'scope': (
            'Census of single-word vectors that extend the q96 map while '
            'preserving LP=I, deletion closure, and the graded-letter axioms. '
            'All 45 have height at least 26. The q99 triple 224 is admissible only '
            'after the two pairs 160 and 192 are added, so it is not in this '
            'single-word census. Multi-word additions and different low blocks '
            'are outside this census. No score change.'
        ),
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True, 'count': len(rows), 'min_height': 26,
                      'max_height': 44}))


if __name__ == '__main__':
    main()
