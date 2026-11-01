#!/usr/bin/env python3
"""Verify global block covers without rerunning gcd or factor discovery."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
import lzma
from pathlib import Path

from compressed_fused_six_vertex import construct_compressed
from f16_identity_operator_profiles import MUL, add, multiply
from check_fused_closed_profiles import verify_irreducible
from fused_q12_endpoint_block_cover import (Q, M, PLACES, CHARGES, build_case,
                                            digest, normalise_at)

ROOT = Path(__file__).resolve().parents[1]
PACKET = ROOT/'evidence/fused-q12-endpoint-block-cover.json.xz'
SOURCES = ['tools/fused_q12_endpoint_block_cover.py', 'tools/check_fused_q12_endpoint_block_cover.py',
           'tools/compressed_fused_six_vertex.py', 'tools/f16_identity_operator_profiles.py',
           'tools/fused_six_vertex_crossings.py', 'tools/f4_lattice_linear.py',
           'tools/check_fused_multipoint_profiles.py', 'tools/check_fused_closed_profiles.py']


def polynomial(data):
    result = {int(k): v for k, v in data.items()}
    assert all(k >= 0 and type(v) is int and 0 < v < 16 for k, v in result.items())
    return result


def fibre(p, point, degree):
    if point is None:
        return p.get(degree, 0)
    result = 0
    for k in range(max(p, default=-1), -1, -1):
        result = MUL[result][point] ^ p.get(k, 0)
    return result


def determinant(matrix):
    """Division-free subset determinant in characteristic two."""
    n = len(matrix)
    values = [0]*(1 << n)
    values[0] = 1
    for mask in range((1 << n)-1):
        if not values[mask]:
            continue
        row = mask.bit_count()
        for column in range(n):
            if not (mask >> column & 1):
                values[mask | 1 << column] ^= MUL[values[mask]][matrix[row][column]]
    return values[-1]


def verify(packet):
    assert [case['theta'] for case in packet['cases']] == list(range(1, 16))
    counters = Counter()
    summaries = []
    all_factors = set()
    for case in packet['cases']:
        source = construct_compressed(M, case['theta'], normalize_joint=False)
        matrix, degree = source['matrix'], source['degree']
        assert degree == case['degree'] and digest(matrix) == case['source_sha256']
        assert source['rectangular_symmetry_checks'] == case['rectangular_symmetry_checks']
        assert source['output_symmetry_checks'] == case['output_symmetry_checks']
        assert all(not p or i//Q+i % Q == j//Q+j % Q
                   for i, row in enumerate(matrix) for j, p in enumerate(row))
        regular = []
        for point in PLACES:
            endpoint = [[fibre(p, point, degree) for p in row] for row in matrix]
            determinants = [determinant([[endpoint[i][j] for j in indices] for i in indices])
                            for indices in CHARGES]
            regular.append(all(determinants))
            counters['endpoint_block_determinants'] += len(determinants)
        assert regular == case['regular_fibres']
        expected = [point for point, ok in zip(PLACES, regular) if ok]
        assert [record['endpoint'] for record in case['endpoints']] == expected
        finite_bounds, infinity_bounds = [], []
        for record in case['endpoints']:
            normalised, _ = normalise_at(matrix, record['endpoint'], degree)
            assert digest(normalised) == record['normalised_sha256']
            endpoint = [[fibre(p, record['endpoint'], degree) for p in row] for row in matrix]
            for indices in CHARGES:
                for i, j in product(indices, repeat=2):
                    value = {}
                    for k in indices:
                        value = add(value, normalised[i][k], endpoint[k][j])
                    assert value == matrix[i][j]
            assert all(not p or i//Q+i % Q == j//Q+j % Q
                       for i, row in enumerate(normalised) for j, p in enumerate(row))
            assert all(fibre(p, record['endpoint'], degree) == int(i == j)
                       for i, row in enumerate(normalised) for j, p in enumerate(row))
            assert max(k for row in normalised for p in row for k in p) == degree
            counters['strict_endpoint_entries'] += Q**4
            assert [block['difference'] for block in record['blocks']] == list(range(1-Q, Q))
            supports = Counter()
            for block in record['blocks']:
                difference, value = block['difference'], {}
                for term in block['identity']['terms']:
                    x, y, u, v = term['coordinates']
                    assert all(type(a) is int and 0 <= a < Q for a in (x, y, u, v))
                    assert u-x == y-v == difference
                    value = add(value, multiply(polynomial(term['coefficient']), normalised[Q*y+x][Q*u+v]))
                g = polynomial(block['identity']['polynomial'])
                assert value == g
                reconstructed, distinct = {0: 1}, set()
                for factor in block['factors']:
                    assert len(factor) >= 2 and factor[-1] == 1
                    assert verify_irreducible(tuple(factor))
                    reconstructed = multiply(reconstructed, {k: v for k, v in enumerate(factor) if v})
                    distinct.add(tuple(factor))
                    all_factors.add(tuple(factor))
                assert reconstructed == g
                for factor in distinct:
                    supports[factor] += 1
                counters['polynomial_bezout_identities'] += 1
            lower = 2*Q-1-max(supports.values(), default=0)
            assert lower == record['finite_block_lower'] and lower >= Q
            finite_bounds.append(lower)
            covered, seen = 0, set()
            for witness in record['infinity_minors']:
                difference = witness['difference']
                assert difference not in seen
                seen.add(difference)
                rows = [(x, u) for x, u in product(range(Q), repeat=2) if u-x == difference]
                cols = [(y, v) for y, v in product(range(Q), repeat=2) if y-v == difference]
                ri, ci = witness['rows'], witness['columns']
                assert len(ri) == len(ci) > 0 and len(set(ri)) == len(ri) and len(set(ci)) == len(ci)
                assert all(0 <= a < len(rows) for a in ri) and all(0 <= a < len(cols) for a in ci)
                minor = [[normalised[Q*cols[j][0]+rows[i][0]][Q*rows[i][1]+cols[j][1]].get(degree, 0)
                          for j in ci] for i in ri]
                assert determinant(minor)
                covered += len(ri)
                counters['independent_infinity_minors'] += 1
            assert covered == record['infinity_rank_lower'] and covered >= Q
            infinity_bounds.append(covered)
            counters['normalised_endpoint_instances'] += 1
        summaries.append(dict(theta=case['theta'], degree=degree, regular_endpoints=len(expected),
                              finite_block_lower=min(finite_bounds), infinity_rank_lower=min(infinity_bounds)))
    assert counters['normalised_endpoint_instances'] == 123
    assert counters['polynomial_bezout_identities'] == 2829
    assert 2*3**2 > 4**2
    return dict(scope='All algebraic spectral places and scales for these fixed q12 F16 seeds and regular base-field endpoint normalizations; bare rho>=3/4, not a factoring lower bound.',
                alphabet=Q, fusion_parameter=M, parameter_cases=15, counters=dict(counters),
                distinct_irreducible_factors=len(all_factors), cases=summaries,
                opposite_rank_lower=Q, restored_bare_product_lower=Q**3, bare_rho_lower='3/4')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.generate:
        packet = dict(cases=[build_case(theta) for theta in range(1, 16)],
                      source_sha256={name: sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES})
    else:
        packet = json.loads(lzma.decompress(PACKET.read_bytes()))
    for name, expected in packet['source_sha256'].items():
        assert sha256((ROOT/name).read_bytes()).hexdigest() == expected, name
    report = verify(packet)
    if args.generate:
        PACKET.write_bytes(lzma.compress((json.dumps(packet, sort_keys=True)+'\n').encode()))
    report['packet_sha256'] = sha256(PACKET.read_bytes()).hexdigest()
    report['source_sha256'] = packet['source_sha256']
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    else:
        assert report == json.loads((ROOT/'evidence/fused-q12-endpoint-block-cover-controls.json').read_text())
    print(json.dumps({key: report[key] for key in ('scope', 'counters', 'opposite_rank_lower', 'bare_rho_lower')}))


if __name__ == '__main__':
    main()
