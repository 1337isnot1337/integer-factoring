#!/usr/bin/env python3
"""Exact type-channel controls for the quantitative raw rank-sum transfer."""
import argparse
from collections import defaultdict
from hashlib import sha256
from itertools import product
import json
from math import comb, factorial, prod
from pathlib import Path


def multinomial(counts):
    value = 1
    remaining = sum(counts)
    for count in counts:
        value *= comb(remaining, count)
        remaining -= count
    return value


def unrank_word(counts, ordinal):
    counts = list(counts)
    remaining = sum(counts)
    total = multinomial(counts)
    assert 0 <= ordinal < total
    answer = []
    while remaining:
        for letter, count in enumerate(counts):
            first = total * count // remaining
            if ordinal >= first:
                ordinal -= first
            else:
                answer.append(letter)
                counts[letter] -= 1
                remaining -= 1
                total = first
                break
    assert ordinal == 0
    return tuple(answer)


def unrank_channel(counts, dimensions, ordinal):
    per_word = prod(d**c for d, c in zip(dimensions, counts))
    word_ordinal, local_ordinal = divmod(ordinal, per_word)
    word = unrank_word(counts, word_ordinal)
    channels = [0] * len(word)
    for position in range(len(word)-1, -1, -1):
        local_ordinal, channels[position] = divmod(local_ordinal, dimensions[word[position]])
    assert local_ordinal == 0
    return word, tuple(channels)


def channel_controls():
    checked = 0
    for dimensions in ((2, 3, 1), (3, 2, 4), (1, 3, 2, 5)):
        for length in range(1, 5):
            explicit = defaultdict(list)
            for word in product(range(len(dimensions)), repeat=length):
                counts = tuple(word.count(i) for i in range(len(dimensions)))
                for channels in product(*(range(dimensions[i]) for i in word)):
                    explicit[counts].append((word, channels))
            for counts, expected in explicit.items():
                count = multinomial(counts) * prod(d**c for d, c in zip(dimensions, counts))
                assert count == len(expected)
                actual = [unrank_channel(counts, dimensions, i) for i in range(count)]
                assert actual == expected
                checked += count
    return checked


def tensor_controls():
    # Partition the sixteen nonzero crossing entries into pieces of actual
    # row and column ranks 3, 5, and 8. Both cuts are partial permutations.
    groups = (tuple(range(3)), tuple(range(3, 8)), tuple(range(8, 16)))
    dimensions = tuple(map(len, groups))
    owner = {entry: i for i, group in enumerate(groups) for entry in group}
    length, label_count = 2, 16
    types = sorted({tuple(word.count(i) for i in range(3))
                    for word in product(range(3), repeat=length)})
    totals = [[0] * (label_count**2) for _ in range(label_count**2)]
    entries = 0
    exact_rank_sum = 0
    for counts in types:
        channel_count = multinomial(counts) * prod(d**c for d, c in zip(dimensions, counts))
        channels = [unrank_channel(counts, dimensions, i) for i in range(channel_count)]
        left, right = [], []
        for a, b in product(range(label_count), repeat=2):
            aa, bb = divmod(a, 4), divmod(b, 4)
            lvalue = rvalue = 0
            for index, (word, local) in enumerate(channels):
                selected = [groups[i][j] for i, j in zip(word, local)]
                if all(aa[p]*4+bb[p] == selected[p] for p in range(length)):
                    lvalue |= 1 << index
                if all(bb[p]*4+aa[p] == selected[p] for p in range(length)):
                    rvalue |= 1 << index
            left.append(lvalue)
            right.append(rvalue)
        # Every channel is used once on each side, proving both factor
        # matrices have full channel rank independently of a rank routine.
        assert sorted(v for v in left if v) == [1 << i for i in range(channel_count)]
        assert sorted(v for v in right if v) == [1 << i for i in range(channel_count)]
        exact_rank_sum += channel_count
        for x, y, u, v in product(range(label_count), repeat=4):
            row_value = (left[x*label_count+y] & right[u*label_count+v]).bit_count() & 1
            col_value = (left[x*label_count+u] & right[y*label_count+v]).bit_count() & 1
            pair_owners = (owner[(x//4)*4+y//4], owner[(x%4)*4+y%4])
            expected_counts = tuple(pair_owners.count(i) for i in range(3))
            expected = int(x == v and y == u and expected_counts == counts)
            assert row_value == col_value == expected
            totals[x*label_count+y][u*label_count+v] ^= row_value
            entries += 2
    for x, y, u, v in product(range(label_count), repeat=4):
        assert totals[x*label_count+y][u*label_count+v] == int(x == v and y == u)
    assert exact_rank_sum == sum(dimensions)**length == 256
    prematurely_padded_sum = sum(1 << (d-1).bit_length() for d in dimensions)
    assert prematurely_padded_sum == 20 > sum(dimensions) == 16
    return dict(type_tensors=len(types), individual_factorization_entries=entries,
                crossing_entries=label_count**4, raw_ranks=list(dimensions),
                exact_amplified_rank_sum=exact_rank_sum,
                prematurely_padded_raw_rank_sum=prematurely_padded_sum)


def budget_controls():
    # Integer square roots here are known exactly, so these are integer
    # certificates for unequal row/column dimensions without real arithmetic.
    profiles = ((1, 9, 3), (3, 12, 6), (2, 8, 4), (5, 20, 10))
    total_weight = sum(w for _, _, w in profiles)
    assert all(r*c == w*w for r, c, w in profiles)
    classes = 0
    for length in range(1, 9):
        weight_sum = 0
        for counts in product(range(length+1), repeat=len(profiles)):
            if sum(counts) != length:
                continue
            mult = multinomial(counts)
            assert mult == factorial(length) // prod(factorial(c) for c in counts)
            row = mult * prod(r**c for (r, _, _), c in zip(profiles, counts))
            col = mult * prod(s**c for (_, s, _), c in zip(profiles, counts))
            term = mult * prod(w**c for (_, _, w), c in zip(profiles, counts))
            assert row*col == term*term <= total_weight**(2*length)
            weight_sum += term
            classes += 1
        assert weight_sum == total_weight**length
    global_classes = 0
    for power in range(1, 5):
        for logarithm in (256, 512, 1024, 2048, 4096):
            m = 1
            while (m+1)**power <= logarithm//256:
                m += 1
            raw_budget = 1 << (m**power+1)
            assert 128*(m**power+1) <= logarithm+128
            entropy = raw_budget*(3*logarithm+20)
            # Stronger than the active-length bound with delta=4.
            assert 8*entropy*(logarithm+1)**4 <= 1 << (logarithm//4)
            global_classes += 1
    return dict(multinomial_rank_sum_classes=classes,
                integer_global_entropy_classes=global_classes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = dict(status='PASS', scope='Quantitative conditional transfer; no new exponent',
                  independently_enumerated_channels=channel_controls(),
                  **tensor_controls(), **budget_controls(),
                  checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    text = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text)


if __name__ == '__main__':
    main()
