#!/usr/bin/env python3
"""Input identities, retained-source replays, and independent exact scan controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import lzma
from pathlib import Path
import random
import subprocess

from run_two_place_f4_census import q_planes, read_input, replay_class
from two_place_f4_census import build_input

ROOT = Path(__file__).resolve().parents[1]


def code_for_q(planes):
    r = []
    for q in planes:
        coefficient = value = 0
        for i in range(11):
            coefficient ^= (q >> i) & 1
            value |= coefficient << i
        assert q == value ^ (value << 1)
        r.append(value)
    return sum(((r[plane] >> i) & 1) << (2*i + plane)
               for i in range(11) for plane in range(2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    text, construction = build_input()
    assert text == (ROOT / 'evidence/two-place-f4-census-input.txt').read_bytes()
    source = json.loads(lzma.decompress((ROOT / 'evidence/two-place-pencil-source-8a3ea3e.json.xz').read_bytes()))
    records = source['receipt']['records']
    assert len(records) == 198
    retained = 0
    for record in records:
        code = code_for_q(record['q_planes'])
        assert q_planes(code) == record['q_planes']
        result = subprocess.run([str(args.binary), str(code), str(code+1), '1', '0'],
                                input=text, capture_output=True, check=True)
        response = json.loads(result.stdout)
        assert response['visited'] == 1 and not response['unresolved'] and not response['endpoint_failures']
        assert response['classes'][0]['orders'] == sum(record['profiles'], [])
        retained += 1
    entries, _ = read_input(ROOT / 'evidence/two-place-f4-census-input.txt')
    randomizer = random.Random(2026110104)
    independent = []
    for code in [0, (1 << 22)-1] + [randomizer.randrange(1 << 22) for _ in range(10)]:
        result = subprocess.run([str(args.binary), str(code), str(code+1), '1', '0'],
                                input=text, capture_output=True, check=True)
        response = json.loads(result.stdout)
        assert not response['unresolved'] and not response['endpoint_failures']
        record = response['classes'][0]
        assert record['representative'] == code and record['q_planes'] == q_planes(code)
        independent.append({'code': code, 'q_planes': record['q_planes'],
                            'orders': record['orders'], 'replay': replay_class(record, entries)})
    # A deliberately nonidentity affine candidate order validates both the
    # C++ enumeration mapping and the entire chunk outcome ledger.
    multiplier, offset, start, stop = 2800681, 1926619, 47, 559
    result = subprocess.run([str(args.binary), str(start), str(stop), str(multiplier), str(offset)],
                            input=text, capture_output=True, check=True)
    response = json.loads(result.stdout)
    assert response['visited'] == stop-start
    assert sum(record['count'] for record in response['classes']) == stop-start
    assert not response['unresolved'] and not response['endpoint_failures']
    inverse = pow(multiplier, -1, 1 << 22)
    for record in response['classes']:
        assert start <= inverse * (record['representative']-offset) % (1 << 22) < stop
        assert q_planes(record['representative']) == record['q_planes']
    sources = ['tools/check_two_place_f4_census.py', 'tools/two_place_f4_census.py',
               'tools/scan_two_place_f4_pencil.cpp', 'tools/run_two_place_f4_census.py',
               'tools/native_weighted_smith_f4.cpp', 'tools/weighted_smith_oracle_f4.py',
               'tools/check_joint_energy_conjugations.py', 'tools/joint_energy_conjugations.py',
               'tools/search_f4_recentered_crossings.py']
    report = {'status': 'pass', 'construction': construction,
              'retained_source_cases_replayed': retained,
              'retained_source_complete_orders_replayed': retained * 50,
              'independent_fraction_free_cases': independent,
              'independent_complete_orders_replayed': len(independent) * 50,
              'nonidentity_affine_chunk': response,
              'source_sha256': {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in sources},
              'scope': 'Exact construction identities, retained-source consistency, independent '
              'packed-F4 fraction-free profile controls, and candidate coverage controls. '
              'The full coefficient census is separate; no stronger factoring exponent is claimed.'}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'pass', 'retained_cases': retained,
                      'independent_complete_orders': len(independent) * 50,
                      'nonidentity_candidate_visits': stop-start}))


if __name__ == '__main__':
    main()
