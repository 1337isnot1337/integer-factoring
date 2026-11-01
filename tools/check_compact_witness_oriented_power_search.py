#!/usr/bin/env python3
"""Check compact witnesses against complete reference searches and literal powers."""

import argparse
import ast
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
from time import process_time

from check_fresh_prime_cofactor_control import generate_input
from check_projected_oriented_power_search import (
    SOURCES as PROJECTED_SOURCES, apply_patch as apply_projected_patch,
    make_patch as make_projected_patch,
)
from check_sparse_oriented_power_search import canonical, read_data
from compact_witness_oriented_power_search import factor_with_compact_witnesses
from compact_witness_oriented_verifier import verify_compact_witness_search
from projected_oriented_power_search import factor_with_projected_oriented_refinement

ROOT = Path(__file__).resolve().parents[1]
SOURCES = PROJECTED_SOURCES + ('compact_witness_oriented_power_search.py',
    'compact_witness_oriented_verifier.py', 'check_compact_witness_oriented_power_search.py')
ENCODING = 'base-word-key-orientation-v1'


def query_logic_controls():
    """Compare the complete query syntax after three declared storage edits."""
    path = ROOT / 'tools/sparse_oriented_power_search.py'
    assert sha256(path.read_bytes()).hexdigest() == (
        '57a5178c9afb628ef102776e4353c68763366e80cb104a70287cb306c0000c07')
    old_tree = ast.parse(path.read_text())
    new_tree = ast.parse((ROOT / 'tools/compact_witness_oriented_power_search.py').read_text())
    old = next(n for n in old_tree.body if isinstance(n, ast.FunctionDef)
               and n.name == 'factor_with_sparse_oriented_refinement')
    new = next(n for n in new_tree.body if isinstance(n, ast.FunctionDef)
               and n.name == 'factor_with_compact_witnesses')
    counts = {'witness_encodings': 0, 'checkpoint_releases': 0, 'transport_calls': 0}
    append = ast.parse("witnesses.append((sorted(word.items()), "
        "key_indices[key] if key is not None else None, split['orientation']))").body[0]
    release = ast.parse('final_refinement = final_rows = final_graph = '
                        'final_units = final_independent = None').body[0]

    class StorageEdits(ast.NodeTransformer):
        def visit_FunctionDef(self, node):
            if (node.body and isinstance(node.body[0], ast.Expr)
                    and isinstance(node.body[0].value, ast.Constant)
                    and isinstance(node.body[0].value.value, str)):
                node.body = node.body[1:]
            return self.generic_visit(node)

        def visit_Expr(self, node):
            value = node.value
            if (isinstance(value, ast.Call) and isinstance(value.func, ast.Attribute)
                    and isinstance(value.func.value, ast.Name)
                    and value.func.value.id == 'witnesses' and value.func.attr == 'append'):
                counts['witness_encodings'] += 1
                return deepcopy(append)
            return self.generic_visit(node)

        def visit_Assign(self, node):
            value = node.value
            if (isinstance(value, ast.Call) and isinstance(value.func, ast.Attribute)
                    and isinstance(value.func.value, ast.Name)
                    and value.func.value.id == 'state' and value.func.attr == 'result'):
                counts['checkpoint_releases'] += 1
                return [deepcopy(release), node]
            return self.generic_visit(node)

        def visit_Name(self, node):
            if node.id == 'transport_sparse_witnesses':
                counts['transport_calls'] += 1
                node.id = 'transport_compact_witnesses'
            return node

    expected = StorageEdits().visit(old)
    expected.name = new.name

    class StripMetadata(ast.NodeTransformer):
        def visit_FunctionDef(self, node):
            if (node.body and isinstance(node.body[0], ast.Expr)
                    and isinstance(node.body[0].value, ast.Constant)
                    and isinstance(node.body[0].value.value, str)):
                node.body = node.body[1:]
            return self.generic_visit(node)

        def visit_Dict(self, node):
            kept = [(k, v) for k, v in zip(node.keys, node.values)
                    if not (isinstance(k, ast.Constant) and k.value in
                            ('witness_encoding', 'component_replacements'))]
            node.keys = [k for k, _ in kept]
            node.values = [v for _, v in kept]
            return self.generic_visit(node)

    actual = StripMetadata().visit(new)
    assert counts == {'witness_encodings': 1, 'checkpoint_releases': 1, 'transport_calls': 1}
    assert ast.dump(expected, include_attributes=False) == ast.dump(actual, include_attributes=False)
    return {'complete_query_syntax_unchanged_after_declared_storage_edits': True,
            'declared_storage_edits': counts,
            'reference_sha256': sha256(path.read_bytes()).hexdigest()}


def compact_witnesses(witnesses):
    return [[w['base_word'], w['key_index'], w['split']['orientation']] for w in witnesses]


def apply_patch(reference, patch):
    assert set(patch) == {'projected_patch', 'witness_encoding'}
    assert patch['witness_encoding'] == ENCODING
    result = apply_projected_patch(reference, patch['projected_patch'])
    result['witnesses'] = compact_witnesses(result['witnesses'])
    result['witness_encoding'] = ENCODING
    return result


def make_patch(reference, new):
    assert new['witness_encoding'] == ENCODING
    assert compact_witnesses(reference['witnesses']) == new['witnesses']
    expanded = deepcopy(new)
    expanded['witnesses'] = reference['witnesses']
    del expanded['witness_encoding']
    patch = {'projected_patch': make_projected_patch(reference, expanded),
             'witness_encoding': ENCODING}
    assert apply_patch(reference, patch) == new
    return patch


def mathematical_record(record):
    """Ignore only runtime metadata and normalize the redundant witness fields."""
    result = deepcopy(record)
    if 'witness_encoding' in result:
        assert result.pop('witness_encoding') == ENCODING
    else:
        result['witnesses'] = compact_witnesses(result['witnesses'])
    for name in ('elapsed_seconds', 'phase_seconds', 'component_replacements'):
        result.pop(name)
    for checkpoint in result['checkpoints']:
        checkpoint.pop('elapsed_seconds')
    return result


def verify_case(case):
    path = ROOT / case['reference_file']
    assert sha256(path.read_bytes()).hexdigest() == case['reference_sha256']
    reference = read_data(path)
    data = {**reference, 'result': apply_patch(reference['result'], case['certificate_patch'])}
    verification = verify_compact_witness_search(data)
    assert verification['units']['complete_original_graph_checked'] is True
    return verification


def storage_measure(value):
    """Measure literal canonical JSON bytes and scalar/container cells."""
    if isinstance(value, dict):
        cells = 1 + sum(1 + storage_measure(v)[1] for v in value.values())
    elif isinstance(value, (list, tuple)):
        cells = 1 + sum(storage_measure(v)[1] for v in value)
    else:
        cells = 1
    return len(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()), cells


def finite_controls():
    checks = literal_powers = 0
    methods = set()
    saved = None
    for bits, cutoff in ((32, 5), (34, 5), (36, 13), (38, 13), (40, 31), (42, 31)):
        source = generate_input(bits)
        policies = [(1, 1), (9, 4), (64, 8), (160, 16)]
        if bits == 32:
            policies.append((4096, 4096))
        for cap, interval in policies:
            options = {'checkpoint_interval': interval, 'max_extraction_columns': 64}
            old = canonical(factor_with_projected_oriented_refinement(
                source['integer'], cutoff, cap, **options))
            new = canonical(factor_with_compact_witnesses(
                source['integer'], cutoff, cap, **options))
            assert mathematical_record(old) == mathematical_record(new)
            data = {'input': source, 'result': new}
            replay = verify_compact_witness_search(data)
            literal_powers += replay['literal_power_checks']
            checks += 1
            methods.add(new['method'])
            if saved is None and new['witnesses'] and new['witnesses'][0][1] is not None:
                saved = data
    assert saved is not None
    assert methods == {'direct target factor', 'oriented residual presentation', 'query cap exhausted'}
    rejected = 0
    for mutate in (
        lambda r: r.__setitem__('witness_encoding', 'wrong'),
        lambda r: r['witnesses'][0].append(0),
        lambda r: r['witnesses'][0].__setitem__(2, -r['witnesses'][0][2]),
        lambda r: r['witnesses'][0].__setitem__(1, None),
        lambda r: r['witnesses'][0][0][0].__setitem__(1, r['witnesses'][0][0][0][1] + 1),
        lambda r: r['witnesses'][0][0].append([0, 0]),
    ):
        bad = deepcopy(saved)
        mutate(bad['result'])
        try:
            verify_compact_witness_search(bad)
        except (AssertionError, TypeError, ValueError, IndexError, KeyError):
            rejected += 1
        else:
            raise AssertionError('A corrupted compact witness was accepted.')
    return {'complete_small_search_equalities': checks,
            'independent_literal_power_checks': literal_powers,
            'observed_stopping_methods': sorted(methods),
            'corrupted_witness_rejections': rejected}


def source_hashes():
    return {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
            for name in SOURCES}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    hashes = source_hashes()
    controls = {'query_logic': query_logic_controls(), **finite_controls()}
    if args.replay:
        report = read_data(args.replay)
        assert report['source_sha256'] == hashes
        assert report['finite_controls'] == controls
        for case in report['cases']:
            verify_case(case)
    else:
        cases = []
        for bits in (64, 80):
            path = ROOT / f'evidence/sparse-oriented-power-{bits}.json.xz'
            reference = read_data(path)
            old = reference['result']
            started = process_time()
            result = factor_with_compact_witnesses(old['integer'], old['cutoff'],
                old['max_queries'], checkpoint_interval=old['checkpoint_interval'],
                max_extraction_columns=old['max_extraction_columns'])
            search_cpu = process_time() - started
            new = canonical(result)
            old_bytes, old_cells = storage_measure(old['witnesses'])
            new_bytes, new_cells = storage_measure(new['witnesses'])
            case = {'bits': bits, 'integer': old['integer'],
                    'reference_file': str(path.relative_to(ROOT)),
                    'reference_sha256': sha256(path.read_bytes()).hexdigest(),
                    'certificate_patch': make_patch(old, new),
                    'search_cpu_seconds': search_cpu,
                    'search_elapsed_seconds': new['elapsed_seconds'],
                    'same_queries': new['norm_queries'], 'same_factor': new['factor'],
                    'same_complete_mathematical_result': True,
                    'witness_json_bytes_before': old_bytes,
                    'witness_json_bytes_after': new_bytes,
                    'witness_cells_before': old_cells, 'witness_cells_after': new_cells}
            case['independent_verification'] = verify_case(case)
            cases.append(case)
        report = {'claim': 'Compact witnesses reconstruct every omitted field exactly',
                  'finite_controls': controls, 'cases': cases, 'source_sha256': hashes,
                  'scope': 'Finite complete search equality and literal replay; no uniform exponent gain.',
                  'timing_scope': 'Search CPU excludes JSON conversion and independent replay; concurrent research was active.',
                  'storage_scope': 'Canonical JSON bytes and logical cells are exact; they are not process RSS.'}
        if args.output:
            args.output.write_text(json.dumps(report, indent=2) + '\n')
    assert source_hashes() == hashes
    print(json.dumps({'finite_controls': controls,
        'cases': [{k: v for k, v in case.items() if k not in
                   ('certificate_patch', 'independent_verification')}
                  for case in report['cases']]}, indent=2))


if __name__ == '__main__':
    main()
