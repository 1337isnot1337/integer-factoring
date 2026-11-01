#!/usr/bin/env python3
"""Direct proper-subfield system counts and acquired carry-budget controls."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from unittest.mock import patch

from teichmuller_cyclic_carries import UnramifiedRing, polynomial_sum
from check_teichmuller_cyclic_carries import POLYNOMIALS
import precision_rational_carries as backend

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/teichmuller-packed-system-controls.json'


def evaluate(ring, terms, point):
    answer = 0
    for coefficient, exponents in terms:
        value = coefficient
        for x, exponent in zip(point, exponents):
            value = ring.mul(value, ring.power(x, exponent))
        answer = ring.add(answer, value)
    return answer


def equations_for(name, ring, domain):
    a = next((x for x in domain if ring.mod_two(x) not in (0, 1)), 1)
    if name == 'empty':
        return []
    if name == 'binary-intersection':
        return [[(1, (1, 0)), (1, (0, 1))], [(1, (1, 1))]]
    if name == 'redundant-binary':
        return [[(1, (1, 0)), (1, (0, 1))],
                [(1, (2, 0)), (1, (1, 0))], [(1, (0, 2)), (1, (0, 1))]]
    if name == 'subfield-two-roots':
        # The second coefficient has deliberately non-subfield higher digits.
        lifted = ring.add(a, ring.mul(ring.integer(2), ring.generator))
        return [[(1, (1, 0)), (1, (0, 1))], [(1, (1, 1)), (lifted, (1, 0))]]
    if name == 'four-tags':
        return [[(1, (1,))], [(1, (2,))], [(ring.integer(2), (0,))], [(1, (3,))]]
    if name == 'full-field-quadratic':
        return [[(1, (2,)), (1, (1,)), (1, (0,))]]
    if name == 'inconsistent':
        return [[(1, (1,)), (a, (0,))], [(1, (1,)), (ring.add(a, 1), (0,))]]
    if name == 'nonzero-even-constant':
        return [[(ring.integer(2), (0, 0))]]
    raise ValueError(name)


def run():
    panel = [(1, 1, 2, 'empty'), (2, 1, 2, 'binary-intersection'),
             (3, 1, 2, 'redundant-binary'), (4, 2, 2, 'subfield-two-roots'),
             (4, 1, 1, 'four-tags'), (2, 2, 1, 'full-field-quadratic'),
             (4, 2, 1, 'inconsistent'), (4, 2, 2, 'nonzero-even-constant')]
    counts = dict(system_queries=0, direct_assignments=0, carry_queries=0,
                  tag_combinations=0, exact_integer_queries=0)
    records, digest = [], sha256()
    for precision in (1, 2, 3, 5):
        for degree, width, variables, name in panel:
            if precision == 5 and variables > 1 and name not in ('empty', 'nonzero-even-constant'):
                continue
            ring = UnramifiedRing(precision, POLYNOMIALS[degree])
            domain = [ring.teichmuller(x) for x in range(ring.field_size)]
            domain = [x for x in domain if ring.power(x, 1 << width) == x]
            assert len(domain) == 1 << width
            equations = equations_for(name, ring, domain)
            tags, tag = [], 1
            for _ in equations:
                tags.append(tag)
                tag = ring.mul(tag, ring.generator)
            seen = set()
            for coefficients in product(domain, repeat=len(tags)):
                combined = 0
                for c, tag in zip(coefficients, tags):
                    combined = ring.add(combined, ring.mul(c, tag))
                seen.add(ring.mod_two(combined))
                counts['tag_combinations'] += 1
            assert len(seen) == (1 << width)**len(tags)
            delta = max((sum(h) for eq in equations for c, h in eq), default=0)
            bound = (2*precision-1)*(degree//width)*delta
            for mask in product((False, True), repeat=variables):
                expected = 0
                for point in product(*([x for x in domain if x or not nz] for nz in mask)):
                    expected += all(not ring.mod_two(evaluate(ring, eq, point)) for eq in equations)
                    counts['direct_assignments'] += 1
                queries = []
                def acquire(R, w, t, v, factors, nonzero):
                    assert R is ring and w == t == width and v == variables
                    assert all(len(h) == variables for row, terms in factors for c, h in terms)
                    row_degrees = [sum(max((sum(h) for c, h in terms), default=0)
                                       for row, terms in factors if row == rho) for rho in range(width)]
                    carry_bounds = [sum(max((h[i] << row for c, h in terms), default=0)
                                        for row, terms in factors)//((1 << width)-1)
                                    for i in range(variables)]
                    assert max(row_degrees+carry_bounds, default=0) <= bound
                    queries.append(dict(row_degrees=row_degrees, carry_bounds=carry_bounds))
                    return polynomial_sum(R, w, t, v, factors, nonzero)
                with patch.object(backend, 'polynomial_sum', side_effect=acquire):
                    actual, invoice = backend.packed_system_count(ring, width, variables, equations, mask)
                assert actual == expected % (1 << precision), (name, precision, mask, actual, expected)
                assert len(queries) <= precision+1
                assert invoice['queries'] == len(queries) and invoice['max_initial_carry'] <= bound
                assert invoice['variables'] == variables and invoice['precision'] == precision
                if precision > width*variables:
                    assert actual == expected
                    counts['exact_integer_queries'] += 1
                counts['system_queries'] += 1
                counts['carry_queries'] += len(queries)
                record = dict(name=name, precision=precision, degree=degree, width=width,
                              variables=variables, equations=len(equations), delta=delta, mask=mask,
                              exact_count=expected, residue=actual, carry_bound=bound,
                              acquired_queries=len(queries),
                              max_row_degree=max((max(q['row_degrees']) for q in queries), default=0),
                              max_carry_bound=max((max(q['carry_bounds']) for q in queries), default=0))
                digest.update(json.dumps(record).encode())
                # Lists keep the stored control data identical after JSON round trip.
                records.append(json.loads(json.dumps(record)))
    ring = UnramifiedRing(3, POLYNOMIALS[2])
    invalid = [(1, 1, [[(ring.generator, (1,))]], 'coefficient-subfield'),
               (2, 1, [[], []], 'extension-capacity'), (3, 1, [], 'domain-width'),
               (1, 1, [[(1, (-1,))]], 'negative-exponent'),
               (1, 1, [[(1, (0, 0))]], 'exponent-dimension')]
    rejected = []
    for width, variables, equations, name in invalid:
        try:
            backend.packed_system_count(ring, width, variables, equations)
        except ValueError:
            rejected.append(name)
        else:
            raise AssertionError(name)
    names = ['tools/precision_rational_carries.py', 'tools/check_teichmuller_packed_systems.py',
             'tools/teichmuller_cyclic_carries.py']
    return dict(scope='Acquired coefficient-extension packing for residue systems; no extra summation variables or division precision. Supplied ring degree and carry budget charged; no integer-word compiler.',
                counts=counts, rejected=rejected, controls=records, query_sha256=digest.hexdigest(),
                source_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in names},
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    args = parser.parse_args()
    receipt = run()
    if args.generate:
        RECEIPT.write_text(json.dumps(receipt, indent=2)+'\n')
    else:
        assert receipt == json.loads(RECEIPT.read_text())
    print(json.dumps(dict(scope=receipt['scope'], counts=receipt['counts'], rejected=receipt['rejected'])))


if __name__ == '__main__':
    main()
