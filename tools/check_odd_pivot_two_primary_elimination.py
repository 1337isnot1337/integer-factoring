#!/usr/bin/env python3
"""Independent exact controls for odd-pivot presentation compression."""

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import gcd
from pathlib import Path
import random

from odd_relation_presentation import OddRelationPresentation
from relation_two_saturation import first_useful_relation_top


ROOT = Path(__file__).resolve().parents[1]


def rank(rows):
    matrix = [[Fraction(v) for v in row] for row in rows]
    current = 0
    for j in range(len(matrix[0]) if matrix else 0):
        pivot = next((i for i in range(current, len(matrix)) if matrix[i][j]), None)
        if pivot is None:
            continue
        matrix[current], matrix[pivot] = matrix[pivot], matrix[current]
        scale = matrix[current][j]
        matrix[current] = [v / scale for v in matrix[current]]
        for i in range(current + 1, len(matrix)):
            scale = matrix[i][j]
            matrix[i] = [a - scale * b for a, b in zip(matrix[i], matrix[current])]
        current += 1
        if current == len(matrix):
            break
    return current


def bezout(a, b):
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        q = old_r // r
        old_r, r, old_s, s, old_t, t = r, old_r-q*r, s, old_s-q*s, t, old_t-q*t
    sign = 1 if old_r >= 0 else -1
    return sign*old_r, sign*old_s, sign*old_t


def saturated_basis(rows, dimension):
    """Exact integer basis of Q-rowspan intersect Z^d, for d<=3."""
    r = rank(rows)
    assert r == len(rows) and dimension <= 3
    if r == 0:
        return []
    if r == dimension:
        return [[int(i == j) for j in range(dimension)] for i in range(dimension)]
    if r == 1:
        divisor = gcd(*rows[0])
        return [[v // divisor for v in rows[0]]]
    assert dimension == 3 and r == 2
    a, b = rows
    normal = [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
    divisor = gcd(*normal)
    normal = [v // divisor for v in normal]
    g, x, y = bezout(normal[0], normal[1])
    if not g:
        assert abs(normal[2]) == 1
        return [[1, 0, 0], [0, 1, 0]]
    assert gcd(g, normal[2]) == 1
    basis = [[normal[1]//g, -normal[0]//g, 0],
             [-normal[2]*x, -normal[2]*y, g]]
    assert all(sum(u*v for u, v in zip(row, normal)) == 0 for row in basis)
    return basis


def group_operations(moduli):
    identity = tuple(0 for _ in moduli)
    def add(a, b):
        return tuple((x+y) % m for x, y, m in zip(a, b, moduli))
    def negate(a):
        return tuple(-x % m for x, m in zip(a, moduli))
    return identity, add, negate


def evaluate(word, values, moduli):
    return tuple(sum(coefficient * values[label][i] for label, coefficient in word.items()) % m
                 for i, m in enumerate(moduli))


def check_invariants(engine, values, moduli):
    identity = tuple(0 for _ in moduli)
    assert engine.scale > 0 and engine.scale % 2
    assert len(set(engine.active_labels)) == len(engine.active_labels)
    for label, word in engine.expressions.items():
        assert set(word).issubset(engine.active_labels)
        assert evaluate(word, values, moduli) == tuple(engine.scale*x % m
                                                     for x, m in zip(values[label], moduli))
    assert all(evaluate(row, values, moduli) == identity for row in engine.rows)
    assert rank(engine.retained_matrix()) == len(engine.rows)


def torsion_tops(rows, labels, values, moduli):
    identity, add, _ = group_operations(moduli)
    images = [evaluate(dict(zip(labels, row)), values, moduli)
              for row in saturated_basis(rows, len(labels))]
    subgroup, frontier = {identity}, [identity]
    while frontier:
        value = frontier.pop()
        for generator in images:
            candidate = add(value, generator)
            if candidate not in subgroup:
                subgroup.add(candidate)
                frontier.append(candidate)
    return {v for v in subgroup if v != identity and add(v, v) == identity}


def test_product_groups():
    rng = random.Random(903172)
    controls = extractor_runs = negative = large_odd = rectangular = 0
    for moduli in product((2, 4, 6, 8, 10, 12), repeat=2):
        identity, add, inverse = group_operations(moduli)
        ambient_tops = [v for v in product(*(range(m) for m in moduli))
                        if v != identity and add(v, v) == identity]
        for repetition in range(8):
            coefficient = rng.choice((-7, -5, -3, -1, 1, 3, 5, 7))
            q = tuple(rng.randrange(m) for m in moduli)
            values = {2: (1, 0), 3: (0, 1), 5: q}
            first = {2: -(coefficient*q[0] % moduli[0]),
                     3: -(coefficient*q[1] % moduli[1]), 5: coefficient}
            rows = [[first.get(p, 0) for p in (2, 3, 5)]]
            engine = OddRelationPresentation((2, 3))
            engine.add_relation(first)
            check_invariants(engine, values, moduli)
            desired = 1 + repetition % 3
            while len(rows) < desired:
                c = rng.randrange(-9, 10)
                candidate = [(-c*q[i] % m) + rng.randrange(-2, 3)*m
                             for i, m in enumerate(moduli)] + [c]
                if rank(rows + [candidate]) != len(rows)+1:
                    continue
                assert evaluate(dict(zip((2, 3, 5), candidate)), values, moduli) == identity
                rows.append(candidate)
                engine.add_relation(dict(zip((2, 3, 5), candidate)))
                check_invariants(engine, values, moduli)
            old_tops = torsion_tops(rows, (2, 3, 5), values, moduli)
            new_tops = torsion_tops(engine.retained_matrix(), engine.active_labels, values, moduli)
            assert old_tops == new_tops
            for bad in [None] + ambient_tops:
                for labels, matrix in (((2, 3, 5), rows),
                                       (engine.active_labels, engine.retained_matrix())):
                    top, certificate, stats = first_useful_relation_top(
                        [values[p] for p in labels], matrix, identity, add, inverse,
                        lambda v: None if v == bad else 'useful', verify_relations=True)
                    assert (top is not None) == bool(old_tops - {bad})
                    if top is not None:
                        assert top in old_tops - {bad} and certificate == 'useful'
                    extractor_runs += 1
            controls += 1
            negative += coefficient < 0
            large_odd += abs(coefficient) > 1
            rectangular += desired < 3
    return {'product_group_presentations': controls, 'extractor_replays': extractor_runs,
            'negative_odd_pivots': negative, 'nonunit_odd_pivots': large_odd,
            'rectangular_presentations': rectangular}


def counterexample_and_deep_controls():
    values, moduli = {2: (1,), 5: (4,)}, (6,)
    engine = OddRelationPresentation((2,))
    engine.add_relation({5: 3})
    assert engine.scale == 3 and engine.expressions[5] == {}
    engine.add_relation({2: 2, 5: 1})
    assert engine.retained_matrix() == [[6]]
    check_invariants(engine, values, moduli)
    assert evaluate({2: 2}, values, moduli) != (0,)
    top, _, _ = first_useful_relation_top([(1,)], [[6]], (0,),
                                         lambda a,b: ((a[0]+b[0])%6,),
                                         lambda a: (-a[0]%6,), lambda _: 'useful')
    assert top == (3,)
    deep = OddRelationPresentation((2,))
    values = {2: (17,)}
    expected = 1
    for i in range(128):
        label = 100+i
        c = (-1 if i % 2 else 1) * (3, 5, 7, 9)[i % 4]
        values[label] = (0,)
        deep.add_relation({label: c})
        expected *= abs(c)
        assert deep.scale == expected and deep.expressions[label] == {}
        check_invariants(deep, values, (2310,))
    deep.add_relation({2: 2310})
    check_invariants(deep, values, (2310,))
    assert deep.retained_matrix() == [[2310*expected]]
    # Several fresh labels: largest fresh odd coefficient is the pivot;
    # a larger even coefficient must remain active.
    many = OddRelationPresentation((2,))
    many_values = {2: (1,), 5: (2,), 7: (3,), 11: (4,)}
    record = many.add_relation({2: -27, 5: 5, 7: 3, 11: 2})
    assert record['eliminated_label'] == 7 and many.active_labels == [2, 5, 11]
    check_invariants(many, many_values, (30,))
    return {'order_six_counterexample': 'genuine row 6g=0; naive row 2g=0 rejected',
            'deep_odd_pivots': 128, 'deep_scale_bits': deep.scale.bit_length(),
            'multiple_fresh_label_control': 1}


def dynamic_history_controls():
    source = ROOT / 'evidence/dynamic-leading-norm-relations.json'
    fixture = json.loads(source.read_text())
    from arithmetic_two_group import power
    from binary_quadratic_forms import inverse_form
    from quadratic_ideal_hnf import compose_forms_hnf
    records_out = []
    total_rows = total_pivots = 0
    for record in fixture['cases']:
        records = {int(p): value for p, value in record['prime_records'].items()}
        initial = sorted(p for p in records if p <= record['cutoff'])
        assert len(initial) == record['initial_generators']
        engine = OddRelationPresentation(initial)
        for row in record['history']:
            raw = {int(p): -v for p, v in row['word'].items()}
            p = row['pivot_prime']
            raw[p] = raw.get(p, 0) + row['exponent']
            result = engine.add_relation(raw)
            assert result['eliminated_label'] == row['eliminated_prime']
            if result['odd_pivot'] is not None:
                assert abs(result['odd_pivot']) == 1
                total_pivots += 1
            total_rows += 1
        assert engine.scale == 1
        assert engine.active_labels == record['active_primes']
        assert engine.rows == [{int(p): v for p,v in row.items()} for row in record['relation_rows']]
        assert engine.expressions == {int(p): {int(q): v for q,v in word.items()}
                                      for p,word in record['expressions'].items()}
        assert rank(engine.retained_matrix()) == len(engine.rows)
        delta = record['discriminant']
        identity = (1, delta % 2, (delta % 2-delta)//4)
        def class_word(word):
            value = identity
            for p,e in word.items():
                generator = tuple(records[p][1])
                if e < 0:
                    generator, e = inverse_form(generator), -e
                value = compose_forms_hnf(value, power(generator, e, identity, compose_forms_hnf))
            return value
        for p, word in engine.expressions.items():
            assert class_word(word) == tuple(records[p][1])
        assert all(class_word(row) == identity for row in engine.rows)
        records_out.append({'integer': record['integer'], 'raw_rows': len(record['history']),
                            'active_labels': len(engine.active_labels), 'retained_rows': len(engine.rows)})
    return {'dynamic_histories': len(records_out), 'dynamic_raw_rows': total_rows,
            'dynamic_unit_pivots': total_pivots, 'records': records_out,
            'fixture_sha256': sha256(source.read_bytes()).hexdigest()}


def cubic_prime_class_control():
    from arithmetic_two_group import power
    from binary_quadratic_forms import ambiguous_factor, inverse_form
    from quadratic_ideal_hnf import compose_forms_hnf
    integer = 6743
    identity = (1, 1, 1686)
    g2, g3 = (2, 1, 843), (3, 1, 562)
    cube_target = power(g2, 6, identity, compose_forms_hnf)
    second_target = power(g3, 7, identity, compose_forms_hnf)
    assert cube_target == (27, -13, 64)
    assert second_target == (16, 13, 108)
    assert power(inverse_form(g3), 3, identity, compose_forms_hnf) == cube_target
    assert power(g2, 4, identity, compose_forms_hnf) == second_target
    # The discarded generator really has a nonzero odd-order error if one
    # incorrectly treats the rational expression as a full-group identity.
    odd_error = compose_forms_hnf(g3, power(g2, 2, identity, compose_forms_hnf))
    assert odd_error == (12, -11, 143) and odd_error != identity
    assert power(odd_error, 3, identity, compose_forms_hnf) == identity
    assert power(g2, 18, identity, compose_forms_hnf) != identity
    engine = OddRelationPresentation((2,))
    engine.add_relation({2: 6, 3: 3})
    engine.add_relation({2: -4, 3: 7})
    assert engine.scale == 3 and engine.retained_matrix() == [[-54]]
    top, factor, stats = first_useful_relation_top(
        [g2], engine.retained_matrix(), identity, compose_forms_hnf,
        inverse_form, lambda form: ambiguous_factor(form, integer))
    assert factor == 11 and integer == 11 * 613
    assert top != identity and compose_forms_hnf(top, top) == identity
    return {'cubic_prime_class_control': {
        'integer': integer, 'discriminant': -integer,
        'raw_relations': [[6, 3], [-4, 7]], 'fresh_odd_pivot': 3,
        'selected_leading_norms': [27, 16], 'common_odd_scale': 3,
        'compressed_relation': [-54], 'odd_error_form': odd_error,
        'odd_error_order': 3, 'factor': factor, 'top': top,
        'extractor_stats': stats}}


def rejection_controls():
    rejected = 0
    for callback in (lambda: OddRelationPresentation((2,2)),
                     lambda: OddRelationPresentation(('2',)),
                     lambda: OddRelationPresentation().add_relation({}),
                     lambda: OddRelationPresentation().add_relation({2:0}),
                     lambda: OddRelationPresentation().add_relation({2:Fraction(1,2)}),
                     lambda: OddRelationPresentation().add_relation({'2':1})):
        try:
            callback()
        except ValueError:
            rejected += 1
        else:
            raise AssertionError('Malformed input accepted.')
    return rejected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'status': 'pass',
              'scope': 'Exact odd-pivot compression; no fast relation-selection claim',
              **test_product_groups(), **counterexample_and_deep_controls(),
              **dynamic_history_controls(), **cubic_prime_class_control(),
              'malformed_inputs_rejected': rejection_controls(),
              'engine_sha256': sha256((ROOT/'tools/odd_relation_presentation.py').read_bytes()).hexdigest(),
              'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'records'}, indent=2))


if __name__ == '__main__':
    main()
