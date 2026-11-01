#!/usr/bin/env python3
"""Division-free polynomial checks of the spin-11 annihilator packet."""
import argparse
from hashlib import sha256
from itertools import product
import json
import lzma
from pathlib import Path

from integral_spin_annihilators import integer_source, primitive_source, cut_blocks, generate
from check_modular_spin_crossings import rational_casimir_crossing
from check_all_field_pascal_crossings import multiply, trim

ROOT = Path(__file__).resolve().parents[1]
PACKET = ROOT/'evidence/integral-spin-annihilators.json.xz'
RECEIPT = ROOT/'evidence/integral-spin-annihilator-controls.json'


def check(packet):
    assert packet['spin'] == 11 and packet['alphabet'] == 12
    assert packet['denominator'] == 39916800
    assert [case['prime'] for case in packet['cases']] == [3, 5, 7, 11]
    integer = integer_source()
    counts = dict(polynomial_inverse_entries=0, spectator_entries=0, blocks=0,
                  primitive_entries=0, strict_endpoint_entries=0)
    summaries = []
    for case in packet['cases']:
        p, degree = case['prime'], case['degree']
        common, primitive = primitive_source(integer, p)
        assert list(common) == case['common_polynomial']
        assert degree == {3: 2, 5: 4, 7: 6, 11: 10}[p]
        assert max(len(v)-1 for v in primitive.values()) == degree
        for key, value in primitive.items():
            assert multiply(common, value, p) == trim(integer[key], p)
            counts['primitive_entries'] += 1
        for x, y, u, v in product(range(12), repeat=4):
            polynomial = primitive.get((x, y, v, u), ())
            leading = polynomial[degree] if len(polynomial) > degree else 0
            assert leading == (x == v and y == u)
            counts['strict_endpoint_entries'] += 1
        if p == 3:
            small, denominator = rational_casimir_crossing(2)
            assert denominator == 2
            small = {k: trim([int(2*x) for x in value], 3) for k, value in small.items()}
            for x, y, u, v in product(range(12), repeat=4):
                expected = small.get((x % 3, y % 3, u % 3, v % 3), ()) \
                    if x//3 == u//3 and y//3 == v//3 else ()
                assert primitive.get((x, y, u, v), ()) == expected
                counts['spectator_entries'] += 1
            assert 'witnesses' not in case
            summaries.append(dict(prime=p, degree=degree, factorization='spin2 tensor spectator4'))
            continue
        witnesses = case['witnesses']
        assert [(w['cut'], w['block']) for w in witnesses] == [
            (cut, key) for cut in ('row', 'opposite') for key, _ in cut_blocks(primitive, cut)]
        target = [0]*(2*p+1)
        target[2], target[p+1], target[2*p] = 1, -2 % p, 1
        for cut in ('row', 'opposite'):
            for key, matrix in cut_blocks(primitive, cut):
                witness = next(w for w in witnesses if w['cut'] == cut and w['block'] == key)
                inverse, size = witness['inverse'], len(matrix)
                assert len(inverse) == size and all(len(row) == size for row in inverse)
                assert all(len(v) <= 2*p-degree+1 for row in inverse for v in row)
                for i, j in product(range(size), repeat=2):
                    result = [0]*(2*p+1)
                    for k in range(size):
                        value = multiply(matrix[i][k], inverse[k][j], p)
                        assert len(value) <= len(result)
                        for h, coefficient in enumerate(value):
                            result[h] = (result[h]+coefficient) % p
                    assert result == (target if i == j else [0]*len(result))
                    counts['polynomial_inverse_entries'] += 1
                counts['blocks'] += 1
        summaries.append(dict(prime=p, degree=degree, inverse_degree_at_most=2*p-degree,
                              all_place_orders_at_most=2, bare_rho_lower='1'))
    pins = {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
            for path in (ROOT/'tools/check_modular_spin_crossings.py',
                         ROOT/'tools/check_all_field_pascal_crossings.py',
                         ROOT/'tools/integral_spin_annihilators.py', Path(__file__))}
    return dict(scope='Fixed integral n=11 seeds only; polynomial inverse witnesses, no jet-window or runtime claim.',
                counts=counts, cases=summaries, source_sha256=pins,
                packet_sha256=sha256(PACKET.read_bytes()).hexdigest(), factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    args = parser.parse_args()
    if args.generate:
        PACKET.write_bytes(lzma.compress(json.dumps(generate(), separators=(',', ':')).encode()))
    packet = json.loads(lzma.decompress(PACKET.read_bytes()))
    receipt = check(packet)
    if args.generate:
        RECEIPT.write_text(json.dumps(receipt, indent=2)+'\n')
    else:
        assert receipt == json.loads(RECEIPT.read_text())
    print(json.dumps({key: receipt[key] for key in ('scope', 'counts', 'cases')}))


if __name__ == '__main__':
    main()
