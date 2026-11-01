#!/usr/bin/env python3
"""Complete small-ring inverse checks and direct rational-sum controls."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random
from unittest.mock import patch

from teichmuller_cyclic_carries import UnramifiedRing, polynomial_sum
import teichmuller_unit_inverses as inverse_backend
from teichmuller_unit_inverses import unit_inverse_coefficients, conjugate_factor, rational_sum
from check_teichmuller_cyclic_carries import evaluate, POLYNOMIALS
from joint_frobenius_carries import joint_trace
from check_carry_fourier_elimination import BinaryField

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/teichmuller-unit-inverse-controls.json'


def run():
    complete = []
    for precision, degree in ((1, 1), (2, 2), (3, 2), (4, 2), (2, 3)):
        ring = UnramifiedRing(precision, POLYNOMIALS[degree])
        coefficients = unit_inverse_coefficients(precision)
        units = 0
        for a in range(1 << (precision*degree)):
            norm, partial = 1, 1
            for r in range(degree):
                conjugate = ring.frobenius(a, r)
                norm = ring.mul(norm, conjugate)
                if r:
                    partial = ring.mul(partial, conjugate)
            assert not any(ring.coefficients(norm)[1:])
            actual = ring.mul(partial, ring.evaluate(coefficients, norm))
            if ring.mod_two(a):
                assert ring.mul(a, actual) == 1 and actual == ring.unit_inverse(a)
                units += 1
            else:
                assert actual == 0
        complete.append(dict(precision=precision, degree=degree, elements=1 << (precision*degree),
                             units=units, inverse_polynomial_degree=len(coefficients)-1))
    rng = random.Random(202611019)
    counts = dict(rational_queries=0, direct_assignments=0, nonunit_values=0,
                  field_reductions=0, conjugate_evaluations=0, carry_queries=0)
    digest = sha256()
    for precision in (1, 2, 3):
        for degree, width, t, variables in ((1, 1, 1, 2), (2, 2, 1, 1),
                                          (2, 2, 2, 1), (3, 3, 3, 1), (4, 2, 2, 1)):
            ring = UnramifiedRing(precision, POLYNOMIALS[degree])
            field = BinaryField(sum(c << i for i, c in enumerate(POLYNOMIALS[degree])))
            dimension = variables*(width//t)
            domain = [ring.teichmuller(x) for x in range(ring.field_size)]
            domain = [x for x in domain if ring.power(x, 1 << width) == x]
            assert len(domain) == 1 << width
            for trial in range(5):
                def make_factor():
                    monomial = [0]*dimension
                    monomial[rng.randrange(dimension)] = 1
                    terms = [(rng.randrange(1 << (precision*degree)), (0,)*dimension),
                             (1, tuple(monomial))]
                    return (rng.randrange(t), terms)
                numerators = [make_factor()]
                if trial == 0:
                    denominators = []
                elif trial == 1:
                    denominators = [(0, [(ring.integer(2), (0,)*dimension)])]
                elif trial == 2:
                    denominators = [(0, [(ring.integer(3), (0,)*dimension)]), make_factor()]
                elif trial == 3:
                    denominators = [make_factor(), make_factor()]
                else:
                    denominators = [make_factor() for _ in range(8)]
                for factor in numerators+denominators:
                    for r in range(degree):
                        transformed = conjugate_factor(ring, width, t, variables, factor, r)
                        for point in product(domain, repeat=variables):
                            assert evaluate(ring, variables, t, transformed, point) == ring.frobenius(
                                evaluate(ring, variables, t, factor, point), r)
                            counts['conjugate_evaluations'] += 1
                for mask in product((False, True), repeat=variables):
                    expected = 0
                    for point in product(*([x for x in domain if x or not nz] for nz in mask)):
                        value = 1
                        for f in numerators:
                            value = ring.mul(value, evaluate(ring, variables, t, f, point))
                        for f in denominators:
                            a = evaluate(ring, variables, t, f, point)
                            if not ring.mod_two(a):
                                value = 0
                                counts['nonunit_values'] += 1
                            else:
                                value = ring.mul(value, ring.unit_inverse(a))
                        expected = ring.add(expected, value)
                        counts['direct_assignments'] += 1
                    def degree_budget(factors):
                        return sum(max((sum(h) for c, h in terms), default=0) for row, terms in factors)
                    bound = degree_budget(numerators)+((2*precision+1)*degree-1)*degree_budget(denominators)
                    queries = []
                    def acquire(ring, width, t, variables, factors, nonzero):
                        budget = degree_budget(factors)
                        assert budget <= bound
                        queries.append(budget)
                        return polynomial_sum(ring, width, t, variables, factors, nonzero)
                    with patch.object(inverse_backend, 'polynomial_sum', side_effect=acquire):
                        actual = rational_sum(ring, width, t, variables, numerators, denominators, mask)
                    assert len(queries) <= precision+1
                    counts['carry_queries'] += len(queries)
                    assert actual == expected, (precision, degree, width, t, variables, trial, mask)
                    if precision == 1:
                        assert joint_trace(field, degree, width, t, variables, numerators, denominators, mask) == actual
                        counts['field_reductions'] += 1
                    counts['rational_queries'] += 1
                    digest.update(json.dumps([precision, degree, width, t, variables, numerators, denominators, mask, actual]).encode())
    branches = []
    for precision in range(1, 9):
        coefficients = unit_inverse_coefficients(precision)
        positions = [j for j, c in enumerate(coefficients) if c % (1 << precision)]
        assert min(positions) >= precision and max(positions) <= 2*precision
        assert len(positions) <= precision+1
        branches.append(dict(precision=precision, nonzero_positions=positions,
                             one_factor_branches=len(positions)))
    names = ['tools/teichmuller_cyclic_carries.py', 'tools/teichmuller_unit_inverses.py',
             'tools/check_teichmuller_unit_inverses.py', 'tools/check_teichmuller_cyclic_carries.py',
             'tools/joint_frobenius_carries.py']
    return dict(scope='Collected unit inverses with zero on nonunits, acquired by at most h+1 norm-polynomial branches; denominator degree and count remain charged.',
                seed=202611019, complete_ring_controls=complete, counts=counts, branches=branches,
                query_sha256=digest.hexdigest(),
                source_sha256={name: sha256((ROOT/name).read_bytes()).hexdigest() for name in names},
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
    print(json.dumps(dict(scope=receipt['scope'], complete_elements=sum(x['elements'] for x in receipt['complete_ring_controls']),
                         counts=receipt['counts'], branch_counts=[x['one_factor_branches'] for x in receipt['branches']])))


if __name__ == '__main__':
    main()
