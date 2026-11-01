#!/usr/bin/env python3
"""Bounded exact probe of the next Ward coefficient layer of the84-label fan."""

import argparse
import hashlib
import json
from pathlib import Path
import random
import time

from check_ward_fan_leading_block import build_code
from check_ward_pfaffian_queries import direct_coefficient
from check_ward_fan_entry_queries import columns
from explore_complement_edge_pencils import evaluate_local
from ward_fan_entry_queries import fan_entry_query
from ward_fan_jet_queries import FanJetOracle, STRIDE, original_jet


def controls(code, oracle):
    rng = random.Random(781231)
    weights = columns(3)
    scalar_cases = nonzero_scalars = 0
    for _ in range(120):
        values = [(rng.randrange(1 << 24) & oracle.mask,
                   rng.randrange(1 << 24) & oracle.mask) for _ in weights]
        sites = [(a, m, b) for a, (m, b) in zip(weights, values)]
        expected = direct_coefficient(3, sites, oracle)
        assert fan_entry_query(3, iter(values), oracle.mul) == expected
        scalar_cases += 1
        nonzero_scalars += bool(expected)
    assert nonzero_scalars
    endpoint = FanJetOracle(code, 1, endpoint=True)
    local = evaluate_local((0, 1 << STRIDE, 1, 1, 1, 1))
    records = []
    for _ in range(8):
        x, y = rng.randrange(84), rng.randrange(84)
        labels = x, y, y, x
        value = endpoint.polynomial(*labels)
        assert value & 255 == 1
        expected = endpoint.mul(value, 1 ^ (1 << STRIDE))  # (1+y)^3 mod y^2.
        actual = original_jet(code, labels, local, 1)
        assert expected == actual
        records.append({'labels': labels, 'quotient_jet': hex(value), 'original_jet': hex(actual)})
    return {'ring_cases': scalar_cases, 'nonzero_ring_cases': nonzero_scalars,
            'positive_crossing_jet_controls': records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', type=int, default=2048)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    code = build_code()
    oracle = FanJetOracle(code, 1)
    checked = controls(code, oracle)
    groups = {tau: [(a, b) for a in range(84) for b in range(84)
                    if code[a]['type_syndrome'] ^ code[b]['type_syndrome'] == tau]
              for tau in range(32)}
    groups = {tau: pairs for tau, pairs in groups.items() if pairs}
    rng = random.Random(191723)
    candidates = []
    residual = [i for i, label in enumerate(code) if label['residual'] == 1]
    # Prefer the two unicyclic input directions predicted by the local argument.
    for y in residual:
        for v in residual:
            for x in residual:
                tau = code[x]['type_syndrome'] ^ code[y]['type_syndrome']
                possible = [u for u in residual
                            if code[v]['type_syndrome'] ^ code[u]['type_syndrome'] == tau]
                if possible:
                    candidates.append((x, y, rng.choice(possible), v))
    rng.shuffle(candidates)
    while len(candidates) < args.samples:
        pairs = rng.choice(list(groups.values()))
        x, y = rng.choice(pairs)
        v, u = rng.choice(pairs)
        candidates.append((x, y, u, v))
    candidates = candidates[:args.samples]
    original = evaluate_local((3, 1 ^ (1 << STRIDE), 1, 1, 1, 1))
    found, original_controls = [], []
    for position, labels in enumerate(candidates):
        value = oracle.polynomial(*labels)
        coefficient = value >> STRIDE & 127
        if coefficient:
            found.append({'labels': labels, 'Gamma_y1_univariate': hex(coefficient)})
        if position < 16 or (coefficient and len(original_controls) < 48):
            actual = original_jet(code, labels, original, 4)
            assert actual == value << (3 * STRIDE)
            original_controls.append({'labels': labels, 'quotient_jet': hex(value)})
    names = (Path(__file__).name, 'ward_fan_jet_queries.py', 'check_ward_fan_leading_block.py',
             'ward_fan_entry_queries.py', 'check_ward_fan_entry_queries.py',
             'check_ward_pfaffian_queries.py', 'complement_shell_tensor.py',
             'explore_complement_edge_pencils.py', 'check_ward_stable_operators.py')
    result = {'source_sha256': {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                               for name in names},
              'controls': checked, 'sampled_logical_entries': len(candidates),
              'nonzero_next_boundary_entries': found,
              'independent_original_kernel_controls': original_controls,
              'query_counts': dict(oracle.counts),
              'scope': 'A bounded exact next-boundary probe; each displayed nonzero coefficient is exact, but no complete rank is inferred from sampling.'}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(json.dumps({'samples': len(candidates), 'nonzero': len(found), 'first': found[:3],
                      'query_counts': dict(oracle.counts),
                      'elapsed_seconds': round(time.monotonic() - started, 3)}), flush=True)


if __name__ == '__main__':
    main()
