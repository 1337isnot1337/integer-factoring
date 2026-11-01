#!/usr/bin/env python3
"""Replay the complete binary normalized four-label rank-two channel search."""

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile

from check_nilpotent_flag_crossings import apply_columns, compose_columns, phi, rank
from check_nilpotent_flag_crossings import nilpotent_family, reverse_at_degree, smith_series

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / 'evidence' / 'binary-rank-two-channels.json'
IDENTITY = 0x8421


def multiply(a, b):
    result = 0
    for i in range(4):
        row = 0
        for j in range(4):
            if (a >> (4 * i + j)) & 1:
                row ^= (b >> (4 * j)) & 15
        result |= row << (4 * i)
    return result


def inverse_mask(a):
    rows = [((a >> (4 * i)) & 15) | (1 << (4 + i)) for i in range(4)]
    for j in range(4):
        pivot = next((i for i in range(j, 4) if rows[i] >> j & 1), None)
        if pivot is None:
            return None
        rows[j], rows[pivot] = rows[pivot], rows[j]
        for i in range(4):
            if i != j and rows[i] >> j & 1:
                rows[i] ^= rows[j]
    return sum((row >> 4) << (4 * i) for i, row in enumerate(rows))


def centralizer(a):
    images = [multiply(a, 1 << i) ^ multiply(1 << i, a) for i in range(16)]
    rows = [sum(((image >> j) & 1) << i for i, image in enumerate(images)) for j in range(16)]
    pivots, used = [], 0
    for column in range(16):
        pivot = next((i for i in range(used, 16) if rows[i] >> column & 1), None)
        if pivot is None:
            continue
        rows[used], rows[pivot] = rows[pivot], rows[used]
        for i in range(16):
            if i != used and rows[i] >> column & 1:
                rows[i] ^= rows[used]
        pivots.append(column)
        used += 1
    basis = []
    for free in range(16):
        if free not in pivots:
            basis.append((1 << free) | sum(((rows[i] >> free) & 1) << p for i, p in enumerate(pivots)))
    result = [0]
    for vector in basis:
        result += [value ^ vector for value in result]
    return result


def literal_channel(a, b, p, p_inverse):
    # Build each matrix entry directly from the two transpose-channel factors.
    q0 = multiply(IDENTITY ^ multiply(a, b), p_inverse)
    p1, q1 = multiply(p, a), multiply(b, p_inverse)
    return [sum(((((p >> (4 * k + j)) & 1) & ((q0 >> (4 * i + ell)) & 1))
                 ^ (((p1 >> (4 * k + j)) & 1) & ((q1 >> (4 * i + ell)) & 1)))
                << (4 * k + ell) for k in range(4) for ell in range(4))
            for i in range(4) for j in range(4)]


def row_masks(operator):
    return [sum(value << j for j, value in enumerate(row)) for row in phi(operator, 4)]


def import_search(path):
    records = [json.loads(line) for line in path.read_text().splitlines() if line]
    summaries = [r for r in records if r.get('summary')]
    assert len(summaries) == 1 and summaries[0]['all_commuting']
    summary = dict(summaries[0]); summary.pop('seconds', None)
    return {'classes': [r for r in records if r.get('class')],
            'regular_cases': [r for r in records if r.get('regular') is True],
            'witnesses': [r for r in records if r.get('candidate')], 'summary': summary}


def verify(record):
    gl = {a: inverse_mask(a) for a in range(65536)}
    gl = {a: b for a, b in gl.items() if b is not None}
    assert len(gl) == 20160
    assert all(multiply(a, b) == multiply(b, a) == IDENTITY for a, b in gl.items())
    covered, canonical, pairs, rank_one, distinct = set(), set(), 0, 0, set()
    for item in record['classes']:
        a = item['A']
        orbit = {multiply(multiply(p, a), p_inverse) for p, p_inverse in gl.items()}
        assert len(orbit) == item['orbit_size'] and a == min(orbit)
        assert covered.isdisjoint(orbit)
        covered.update(orbit); canonical.add(a)
        commutants = centralizer(a)
        assert len(commutants) == item['commuting_matrices']
        assert all(multiply(a, b) == multiply(b, a) for b in commutants)
        pairs += len(commutants)
        if a in (0, IDENTITY):
            continue
        for b in commutants:
            base = tuple(literal_channel(a, b, IDENTITY, IDENTITY))
            if rank(row_masks(base)) < 2:
                rank_one += 1
            else:
                distinct.add(base)
    assert covered == set(range(65536))
    summary = record['summary']
    assert len(canonical) == summary['conjugacy_classes'] == 34
    assert pairs == summary['canonical_commuting_pairs'] == 136960
    assert rank_one == summary['rank_one_bases_skipped'] == 58
    assert len(distinct) == summary['distinct_rank_two_bases'] == 5830
    assert summary['tested'] == len(distinct) * len(gl) == 117532800
    histogram, seen, best = Counter(), set(), None
    for item in record['regular_cases']:
        a, b, p = item['A'], item['B'], item['P']
        assert a in canonical and p in gl and multiply(a, b) == multiply(b, a)
        assert (a, b, p) not in seen
        seen.add((a, b, p))
        u = literal_channel(a, b, p, gl[p])
        assert apply_columns(u, IDENTITY) == IDENTITY
        assert all(sum((column >> (5 * j)) & 1 for j in range(4)) % 2 == int(i % 5 == 0)
                   for i, column in enumerate(u))
        s = [column ^ (1 << i) for i, column in enumerate(u)]
        power = [1 << i for i in range(16)]
        for _ in range(15):
            power = compose_columns(s, power)
        assert any(power) and not any(compose_columns(s, power))
        h0, h2 = row_masks(u), row_masks(compose_columns(u, u))
        zeros = rank(h0)
        toeplitz = h0 + [(h0[i] ^ h2[i]) | (h0[i] << 16) for i in range(16)]
        ones = rank(toeplitz) - 2 * zeros
        assert zeros == item['zero_orders'] == 2
        assert ones == item['order_one_count'] >= 11
        histogram[ones] += 1
        if best is None or ones < best[0]:
            best = ones, item, s
    assert len(seen) == summary['regular'] == 5888
    assert [histogram[i] for i in range(17)] == summary['order_one_histogram']
    crossing = nilpotent_family(best[2], 4, 16)
    complete_orders = smith_series(reverse_at_degree(crossing, 15), 257)
    assert len(complete_orders) == 16 and max(complete_orders) < 257
    assert complete_orders.count(0) == 2 and complete_orders.count(1) == 11
    assert 197 == 2 * 16 + 11 * 15 and 197 > 64 and 197 ** 5 > 4 ** 19
    return {'matrix_similarity_classes': len(canonical), 'matrices_covered': len(covered),
            'canonical_commuting_pairs': pairs, 'distinct_rank_two_base_channels': len(distinct),
            'rank_one_nonscalar_bases': rank_one, 'regular_cases_replayed': len(seen),
            'order_one_histogram': {str(k): value for k, value in sorted(histogram.items())},
            'minimum_order_one_count': 11, 'profile_cost_lower_bound': 197,
            'best_retained_witness': best[1], 'best_complete_row_orders': complete_orders}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, default=DEFAULT)
    parser.add_argument('--search-log', type=Path)
    parser.add_argument('--regenerate', action='store_true')
    parser.add_argument('--replay-search', action='store_true')
    args = parser.parse_args()
    sources = (Path(__file__).name, 'search_cyclic_rank_two_channels.cpp', 'check_nilpotent_flag_crossings.py')
    hashes = {name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest() for name in sources}
    if args.regenerate:
        if args.search_log is None:
            parser.error('--regenerate requires --search-log')
        record = import_search(args.search_log)
        record['source_sha256'] = hashes
        record['controls'] = verify(record)
        record['scope'] = 'All binary normalized four-label commuting two-channel forms with binary invertible P; not arbitrary fields or nonunital forms.'
        record['factoring_exponent_improved'] = False
        args.evidence.write_text(json.dumps(record, separators=(',', ':')) + '\n')
    else:
        record = json.loads(args.evidence.read_text())
        assert record['source_sha256'] == hashes
        assert record['controls'] == verify(record)
    if args.replay_search:
        with tempfile.TemporaryDirectory(prefix='binary-rank-two-') as temporary:
            binary = Path(temporary) / 'search'
            subprocess.run(['g++', '-O3', '-std=c++17', str(ROOT / 'tools' / sources[1]), '-o', str(binary)], check=True)
            log = Path(temporary) / 'search.jsonl'
            with log.open('w') as output:
                subprocess.run([str(binary), '--all-commuting'], stdout=output, check=True)
            rerun = import_search(log)
            assert all(record[key] == value for key, value in rerun.items())
    print(json.dumps({'status': 'pass', 'full_search_replayed': args.replay_search, **record['controls']}, indent=2))


if __name__ == '__main__':
    main()
