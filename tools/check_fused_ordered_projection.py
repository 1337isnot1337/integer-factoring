#!/usr/bin/env python3
"""Validate the complete ordered binary envelope audit and its rank-sum argument."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
from tempfile import TemporaryDirectory

from check_fused_scalar_projection import (DUAL, TARGET, coefficient_space, combine,
                                         dense_binary_rank, partitions, write_blocks)

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT/'evidence/fused-ordered-scalar-audit.json'
SOURCE = ROOT/'evidence/fused-ordered-scalar-audit-source.cpp'
REPORT = re.compile(r'dimension (\d+) visited (\d+) cheap (\d+) seconds ([0-9.]+)\n')


def zero_source(source):
    assert source.count('offset=1ull<<p;') == source.count('rr*rc>501') == 1
    return source.replace('offset=1ull<<p;', 'offset=0;').replace('rr*rc>501', 'rr*rc>55')


def full_replay(audit, data, directory, jobs):
    compiler = shlex.split(os.environ.get('CXX', 'c++'))
    source = SOURCE.read_text()
    for side in ('one', 'zero'):
        path = directory/f'{side}.cpp'
        path.write_text(source if side == 'one' else zero_source(source))
        native = directory/side
        subprocess.run(compiler+['-O3', '-std=c++17', str(path), '-o', str(native)], check=True)

        def shard(index):
            output = directory/f'{side}-{index}.csv'
            run = subprocess.run(['nice', '-n', '10', str(native), str(data), str(output),
                                  'dual', '2', str(index)], check=True, capture_output=True, text=True)
            report = REPORT.fullmatch(run.stdout)
            assert report is not None, run.stdout+run.stderr
            expected = audit['shards'][side][index]
            assert [int(report[i]) for i in (1, 2, 3)] == [33, 1 << 30, expected['eligible']]
            with output.open() as stream:
                records = [[int(row[key]) for key in ('coordinates', 'row_rank', 'column_rank')]
                           for row in csv.DictReader(stream)]
            assert len(records) == expected['eligible']
            print(f'{side} shard {index}: {1 << 30} states, {len(records)} eligible', flush=True)
            return records

        with ThreadPoolExecutor(max_workers=min(jobs, 4)) as executor:
            records = sorted(row for rows in executor.map(shard, range(4)) for row in rows)
        assert records == (audit['affine_one_catalogue'] if side == 'one' else [])


def certificate(replay=False, jobs=4):
    audit = json.loads(AUDIT.read_text())
    source = SOURCE.read_bytes()
    assert sha256(source).hexdigest() == audit['source_one_sha256']
    assert sha256(zero_source(source.decode()).encode()).hexdigest() == audit['source_zero_sha256']
    assert audit['dimension'] == 33 and audit['dual'] == DUAL and audit['target'] == TARGET
    assert audit['thresholds'] == {'one': 501, 'zero': 55} and audit['shard_bits'] == 2
    assert (DUAL & TARGET).bit_count() % 2 == 1
    total = 0
    for side in ('one', 'zero'):
        assert len(audit['shards'][side]) == 4
        for index, report in enumerate(audit['shards'][side]):
            assert report['index'] == index and report['dimension'] == 33
            assert report['states'] == 1 << 30 and report['exit_code'] == 0
            if side == 'zero':
                assert report['eligible'] == 0
            total += report['states']
    assert total == audit['audited_states'] == 1 << 33
    cuts, labels = coefficient_space(ordered=True)
    assert len(labels) == 33
    with TemporaryDirectory(prefix='factoring-ordered-envelope-') as temporary:
        directory = Path(temporary)
        data = directory/'blocks.bin'
        write_blocks(data, cuts, partitions())
        assert sha256(data.read_bytes()).hexdigest() == audit['data_sha256']
        if replay:
            full_replay(audit, data, directory, jobs)
    records = audit['affine_one_catalogue']
    assert len(records) == 21 and records == sorted(records)
    assert len({row[0] for row in records}) == 21
    pairs = []
    for coordinate, rr, rc in records:
        assert 0 < coordinate < 1 << 33 and (DUAL & coordinate).bit_count() % 2 == 1
        assert rr*rc <= 501 and rr*rc >= 56
        assert [dense_binary_rank(combine(basis, coordinate)) for basis in cuts] == [rr, rc]
        complement = [dense_binary_rank(combine(basis, TARGET ^ coordinate)) for basis in cuts]
        assert complement[0]*complement[1] >= 224
        pairs.append(dict(coordinates=coordinate, ranks=[rr, rc], complement_ranks=complement))
    for index, report in enumerate(audit['shards']['one']):
        # The two fixed free coordinates are the original bits 31 and 32.
        assert sum((coordinate >> 31) == index for coordinate, _, _ in records) == report['eligible']
    assert min(row['complement_ranks'][0]*row['complement_ranks'][1] for row in pairs) == 3969
    assert min(rr*rc for _, rr, rc in records) == 162
    # At least three pieces cost >=3 sqrt(56). For two pieces, either the
    # affine-one piece has product >=502, or its listed complement has >=224.
    assert 502 >= 4*56 and 3969 >= 4*56
    assert 81**2 >= 9*56  # the only one-piece decomposition is the crossing
    assert 99**2 > 2*70**2 and 504**35 > 9**99
    names = ('check_fused_ordered_projection.py', 'check_fused_scalar_projection.py',
             'f4_lattice_linear.py', 'joint_energy_conjugations.py')
    return dict(dimension=33, physical_alphabet=9, basis_labels=labels,
                audited_states=total, affine_one_eligible=21, affine_zero_small_products=0,
                minimum_nonzero_product_lower_bound=56, smallest_listed_product=162,
                minimum_listed_complement_product=3969, independently_replayed_full_ranks=84,
                rank_sum_lower_bound='3*sqrt(56)=sqrt(504)', complement_pairs=pairs,
                data_sha256=audit['data_sha256'], audit_sha256=sha256(AUDIT.read_bytes()).hexdigest(),
                source_one_sha256=audit['source_one_sha256'], source_zero_sha256=audit['source_zero_sha256'],
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in names},
                scope='Every binary linear decomposition of the crossing in the full ordered two-copy scalar-projection envelope of the literal q3 theta4 fused seed.',
                conclusion='Every decomposition has rank-sum cost at least sqrt(504), strictly above 9^sqrt(2), so this entire envelope fails the finite-block threshold.',
                quick_check_scope='Reconstructs the complete space and exact audited input, validates frozen-source provenance and all 84 retained matrix ranks; does not repeat the exhaustive native screen unless --full-replay is supplied.',
                limit='Other fields, higher tensor powers and physical type restrictions remain open; this is local author computation without independent full audit replay.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--full-replay', action='store_true')
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/fused-ordered-scalar-projection.json')
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    result = certificate(args.full_replay, args.jobs)
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('basis_labels', 'complement_pairs', 'source_sha256')}))


if __name__ == '__main__':
    main()
