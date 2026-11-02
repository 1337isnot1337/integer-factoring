#!/usr/bin/env python3
"""Exact global quartic-gate expansion controls; no cut-rank calculation."""
import argparse
from collections import Counter
from functools import reduce
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from operator import xor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/charge-creation-square.json'
PARTITIONS = tuple((left, tuple(i for i in range(8) if i not in left))
                   for left in combinations(range(8), 4) if 0 in left)


def zero_quartet(colors):
    return len(set(colors)) == 4 and 0 not in colors and reduce(xor, colors) == 0


def literal_partitions(colors):
    """Count each unordered complementary spatial pair once, as integers."""
    return [(left, right) for left, right in PARTITIONS
            if zero_quartet(tuple(colors[i] for i in left))
            and zero_quartet(tuple(colors[i] for i in right))]


def binary_rank(columns):
    pivots = {}
    for column in columns:
        while column:
            lead = column.bit_length()-1
            if lead not in pivots:
                pivots[lead] = column
                break
            column ^= pivots[lead]
    return len(pivots)


def code_coefficient(colors):
    """The proved closed formula, independent of quartet enumeration."""
    if 0 in colors or len(set(colors)) < 8 or reduce(xor, colors):
        return 0
    return int(binary_rank([(color << 1) | 1 for color in colors]) <= 6)


def run():
    quartets = tuple(word for subset in combinations(range(1, 8), 4)
                     if reduce(xor, subset) == 0 for word in permutations(subset))
    assert len(quartets) == 168 and len(PARTITIONS) == 35
    # Directly expand every disjoint pair of q8 creation gates. This does
    # not use the repeated-color involution or the relation-code formula.
    coefficients = Counter()
    for left, right in PARTITIONS:
        for a, b in product(quartets, repeat=2):
            word = [0]*8
            for position, color in zip(left+right, a+b):
                word[position] = color
            coefficients[tuple(word)] += 1
    histogram = Counter(coefficients.values())
    assert sum(coefficients.values()) == 35*168**2
    assert all(count % 2 == 0 for count in coefficients.values())
    digest = sha256()
    for word, count in sorted(coefficients.items()):
        digest.update(bytes(word)+count.to_bytes(2, 'little'))

    # Check the actual free involution, with the SAME repeated pair for
    # every partition of a given word. Include an impossible multiplicity.
    involutions = []
    for colors in ((1, 2, 4, 7, 1, 2, 4, 7),
                   (1, 2, 4, 7, 1, 3, 5, 7),
                   (1, 1, 1, 2, 3, 4, 5, 6)):
        repeated = next(c for c in colors if colors.count(c) > 1)
        i, j = [k for k, c in enumerate(colors) if c == repeated][:2]
        valid = {left for left, _ in literal_partitions(colors)}
        images = {}
        for left in valid:
            assert (i in left) != (j in left)
            image = set(left) ^ {i, j}
            if 0 not in image:
                image = set(range(8))-image
            image = tuple(sorted(image))
            assert image in valid and image != left
            images[left] = image
        assert all(images[images[left]] == left for left in valid)
        involutions.append(dict(colors=list(colors), valid_partitions=len(valid)))

    panel, code_digest = Counter(), sha256()
    for colors in combinations(range(1, 16), 8):
        total = reduce(xor, colors)
        rank = binary_rank([(color << 1) | 1 for color in colors])
        count = len(literal_partitions(colors))
        assert count % 2 == code_coefficient(colors)
        if total == 0:
            assert count == (2**(8-rank)-2)//2
        else:
            assert count == 0
        panel[(total == 0, rank, count)] += 1
        code_digest.update(bytes(colors)+bytes((rank, count)))
    affine_plane = tuple(range(8, 16))
    assert len(literal_partitions(affine_plane)) == 7
    # A global XOR-zero condition alone is insufficient at larger q.
    only_full_relation = (1, 2, 4, 8, 16, 32, 64, 127)
    assert reduce(xor, only_full_relation) == 0
    assert binary_rank([(c << 1) | 1 for c in only_full_relation]) == 7
    assert literal_partitions(only_full_relation) == []
    assert code_coefficient(only_full_relation) == 0
    return dict(
        scope='Exact quartic creation-product coefficients over characteristic two. '
              'No Smith profile, pressure, global contraction or factoring bound.',
        q8_full_disjoint_expansion=dict(ordered_quartets=len(quartets),
            spatial_partitions=len(PARTITIONS),
            integer_terms=sum(coefficients.values()),
            supported_words=len(coefficients), odd_coefficients=0,
            multiplicity_histogram={str(k): v for k, v in sorted(histogram.items())},
            word_counts_sha256=digest.hexdigest()),
        repeated_color_involutions=involutions,
        q16_distinct_subsets=dict(subsets=sum(panel.values()),
            panels=[dict(total_xor_zero=zero, augmented_rank=rank,
                         integer_partitions=count, subsets=n)
                    for (zero, rank, count), n in sorted(panel.items())],
            records_sha256=code_digest.hexdigest()),
        nonzero_divided_square=dict(colors=list(affine_plane), integer_partitions=7),
        total_xor_alone_insufficient=dict(colors=list(only_full_relation),
            augmented_rank=7, integer_partitions=0),
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    else:
        assert result == json.loads(RECEIPT.read_text()), 'receipt differs; inspect and regenerate'
    print('Quartic creation square: complete q8 cancellation and relation-code controls PASS')
