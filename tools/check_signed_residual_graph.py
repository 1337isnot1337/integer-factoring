#!/usr/bin/env python3
"""Exact presentation maps, finite character counts, and ideal-factor controls."""
import argparse
from copy import deepcopy
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

from binary_quadratic_forms import ambiguous_factor, inverse_form
from check_oriented_coprime_basis import strict_gain
from cofactor_collision_relations import RationalRowBasis
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top
from signed_residual_graph import eliminate_signed_residual_graph, verify_signed_graph_elimination

ROOT = Path(__file__).resolve().parents[1]


def character_count(rows, dimension, modulus):
    return sum(all(sum(a*b for a, b in zip(row, vector)) % modulus == 0 for row in rows)
               for vector in product(range(modulus), repeat=dimension))


def exhaustive_graph_controls():
    cases = character_tests = 0
    digest = sha256()
    edges = ((0, 1), (0, 2), (1, 2))
    for choices in product((-1, 0, 1), repeat=3):
        for loops in range(8):
            for anchors in range(8):
                rows = []
                for (u, v), sign in zip(edges, choices):
                    if sign:
                        row = [0]*3
                        row[u], row[v] = 1, sign
                        rows.append(row)
                for vertex in range(3):
                    for mask, coefficient in ((loops, 2), (anchors, 1)):
                        if mask >> vertex & 1:
                            row = [0]*3
                            row[vertex] = coefficient
                            rows.append(row)
                result = eliminate_signed_residual_graph(rows, 0, generator_count=3)
                free = result['counts']['free_components']
                unbalanced = result['counts']['unbalanced_components']
                for modulus in (3, 4):
                    observed = character_count(rows, 3, modulus)
                    expected = modulus**free * (2 if modulus == 4 else 1)**unbalanced
                    assert observed == expected
                    character_tests += 1
                digest.update(json.dumps([rows, result['torsion_rows'], free, unbalanced],
                                         separators=(',', ':')).encode())
                cases += 1
    return {'complete_three_vertex_graphs': cases, 'independent_character_counts': character_tests,
            'canonical_result_sha256': digest.hexdigest()}


def mixed_row_controls():
    rng = Random(20261101)
    cases = character_tests = 0
    digest = sha256()
    for trial in range(420):
        base, residuals = rng.randrange(4), rng.randrange(7)
        dimension = base + residuals
        rows = []
        for _ in range(rng.randrange(15)):
            row = [rng.randrange(-5, 6) for _ in range(base)] + [0]*residuals
            if residuals:
                if rng.randrange(4):
                    for _ in range(rng.randrange(3)):
                        row[base+rng.randrange(residuals)] += rng.choice((-1, 1))
                else:
                    for i in range(residuals):
                        row[base+i] = rng.randrange(-3, 4)
            rows.append(row)
        result = eliminate_signed_residual_graph(rows, base, generator_count=dimension,
                                                 allow_extra_rows=True)
        verify_signed_graph_elimination(rows, json.loads(json.dumps(result)))
        if trial < 100 and dimension <= 5:
            for modulus in (3, 4):
                left = character_count(rows, dimension, modulus)
                right = character_count(result['torsion_rows'], len(result['torsion_columns']), modulus)
                right *= modulus**len(result['free_root_columns'])
                assert left == right
                character_tests += 1
        digest.update(json.dumps([rows, result['torsion_rows']], separators=(',', ':')).encode())
        cases += 1
    return {'mixed_base_and_extra_row_cases': cases, 'independent_character_counts': character_tests,
            'canonical_result_sha256': digest.hexdigest()}


def extract(integer, generators, rows):
    delta = generators[0][1]**2-4*generators[0][0]*generators[0][2]
    parity = delta % 2
    identity = (1, parity, (parity-delta)//4)
    basis = RationalRowBasis(len(generators))
    for row in rows:
        basis.add(row)
    top, factor, stats = first_useful_relation_top(
        generators, basis.rows, identity, compose_forms_hnf, inverse_form,
        lambda form: ambiguous_factor(form, integer))
    return {'factor': factor, 'ambiguous_form': top, 'independent_rows': basis.rows, 'stats': stats}


def ideal_controls():
    cycle = strict_gain()
    old = cycle['transported_rows']
    result = eliminate_signed_residual_graph(old, len(cycle['active_primes']))
    generators = [tuple(cycle['prime_records'][p][1]) for p in cycle['active_primes']]
    generators += [tuple(atom['form']) for atom in cycle['refinement']['atoms']]
    selected = [generators[i] for i in result['torsion_generator_indices']]
    extraction = extract(cycle['integer'], selected, result['torsion_rows'])
    assert result['torsion_rows'] == [[-2, 22, 64]] and extraction['factor'] == 1009

    records = {r.prime: r for r in build_prime_form_factor_base(-84, 41)['records']}
    primes = (5, 17, 41)
    triangle_generators = [records[p].positive for p in primes]
    triangle_rows = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
    assert all(form == (5, 4, 5) for form in triangle_generators)
    assert compose_forms_hnf((5, 4, 5), (5, 4, 5)) == (1, 0, 21)
    triangle = eliminate_signed_residual_graph(triangle_rows, 0)
    assert triangle['torsion_rows'] == [[2]] and triangle['counts']['unbalanced_components'] == 1
    triangle_extraction = extract(21, [triangle_generators[i] for i in triangle['torsion_generator_indices']],
                                  triangle['torsion_rows'])
    assert triangle_extraction['factor'] == 3
    return {
        'balanced_four_cycle': {'integer': cycle['integer'], 'rows': old,
                               'reduction': result, 'extraction': extraction},
        'necessary_unbalanced_root': {'integer': 21, 'residual_norms': primes,
                                     'roots': [records[p].root for p in primes],
                                     'generators': triangle_generators, 'rows': triangle_rows,
                                     'reduction': triangle, 'extraction': triangle_extraction,
                                     'base_only_relation_count': 0}}


def prefix_control(path):
    data = json.loads(path.read_text())
    if 'cases' in data:
        matches = [case['result'] for case in data['cases']
                   if case['result']['integer'] == 1000036000099]
        assert len(matches) == 1
        data = matches[0]
    base = len(data['active_primes'])
    result = eliminate_signed_residual_graph(data['transported_rows'], base, allow_extra_rows=True)
    generators = [tuple(data['prime_records'][str(p)][1]) for p in data['active_primes']]
    generators += [tuple(atom['form']) for atom in data['refinement']['atoms']]
    selected = [generators[i] for i in result['torsion_generator_indices']]
    extraction = extract(data['integer'], selected, result['torsion_rows'])
    assert extraction['factor'] is not None and data['integer'] % extraction['factor'] == 0
    return {'input_path': str(path), 'input_sha256': sha256(path.read_bytes()).hexdigest(),
            'integer': data['integer'], 'norm_queries': data['norm_queries'],
            'reduction': result, 'extraction': extraction}


def rejection_controls():
    bad = [lambda: eliminate_signed_residual_graph([[1, 1, 1]], 0),
           lambda: eliminate_signed_residual_graph([[3]], 0),
           lambda: eliminate_signed_residual_graph([[1], [1, 0]], 0),
           lambda: eliminate_signed_residual_graph([[1.5]], 0),
           lambda: eliminate_signed_residual_graph([[1]], 2),
           lambda: eliminate_signed_residual_graph([], -1)]
    for call in bad:
        try:
            call()
        except ValueError:
            continue
        raise AssertionError('An invalid graph contract was accepted.')
    rows = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
    good = eliminate_signed_residual_graph(rows, 0)
    tampered = []
    item = deepcopy(good); item['output_sources'][0][0] += 1; tampered.append(item)
    item = deepcopy(good); item['generator_words'][1][0] *= -1; tampered.append(item)
    item = deepcopy(good)
    changed = next(i for i, combination in enumerate(item['old_row_combinations']) if combination)
    item['old_row_combinations'][changed] = {}
    tampered.append(item)
    item = deepcopy(good); item['new_generator_indices'][0] = 2; tampered.append(item)
    for item in tampered:
        try:
            verify_signed_graph_elimination(rows, item)
        except AssertionError:
            continue
        raise AssertionError('A corrupted presentation certificate was accepted.')
    return {'invalid_contracts_rejected': len(bad), 'corrupted_certificates_rejected': len(tampered)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = {'claim': 'Integer signed-graph presentation reduction with torsion preservation',
              'scope': 'Supplied relation rows; no bound on witness search and no exponent improvement.',
              'exhaustive': exhaustive_graph_controls(), 'mixed': mixed_row_controls(),
              'ideals': ideal_controls(), 'rejections': rejection_controls()}
    if args.prefix:
        report['prefix'] = prefix_control(args.prefix)
    sources = ('signed_residual_graph.py', 'check_signed_residual_graph.py',
               'check_oriented_coprime_basis.py', 'oriented_coprime_basis.py',
               'relation_two_saturation.py', 'quadratic_ideal_hnf.py')
    report['source_sha256'] = {'tools/'+name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in sources}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: report[key] for key in ('exhaustive', 'mixed', 'rejections')}, indent=2))
    if args.prefix:
        print(json.dumps({'prefix': report['prefix']['reduction']['counts'],
                          'factor': report['prefix']['extraction']['factor']}))


if __name__ == '__main__':
    main()
