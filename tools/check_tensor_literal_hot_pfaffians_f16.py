#!/usr/bin/env python3
"""Independent original-kernel controls for literal F16 hot Pfaffians.

The reference kernel is extracted from the existing interval reflection
formula, not from the new hot closed form. Noncoherent physical maps are
included; literal matrices undergo independent binary Toeplitz ranks.
"""
import argparse
import ast
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random
import time

import check_tensor_hot_compounds_f16 as controls
import tensor_hot_compounds_f16 as field
import tensor_literal_hot_pfaffians_f16 as literal
from tensor_hot_boundary_retractions import boundary_code, validate_height_maps


def reference_kernel(modes):
    # Only these four definitions are needed. Their unrelated Schur-runtime
    # imports are not part of the reflected kernel calculation.
    source = Path(__file__).with_name('reflection_dual_rail_gaussian.py')
    tree = ast.parse(source.read_text())
    names = {'validate', 'logical_heights', 'laurent_kernel', 'gauged_kernel'}
    selected = ast.Module(body=[node for node in tree.body
                               if isinstance(node, ast.FunctionDef) and node.name in names], type_ignores=[])
    namespace = {}
    exec(compile(selected, str(source), 'exec'), namespace)
    interval = 4
    while interval < 4 * modes:
        interval *= 2
    ambient = interval // 2
    indices = [ambient * leg + ambient - 1 - hot for leg in range(4) for hot in range(modes)]
    zero = namespace['gauged_kernel'](interval, 'zero')
    infinity = namespace['gauged_kernel'](interval, 'infinity')

    def halve(polynomial):
        assert all(not (polynomial >> exponent) & 1 for exponent in range(1, polynomial.bit_length(), 2))
        return sum(((polynomial >> exponent) & 1) << (exponent // 2)
                   for exponent in range(0, polynomial.bit_length(), 2))

    kernel = [[halve(zero[i][j]) for j in indices] for i in indices]
    other = [[halve(infinity[i][j]) for j in indices] for i in indices]
    swap = lambda i: (0, 2, 1, 3)[i // modes] * modes + i % modes
    assert all(other[i][j] == kernel[swap(i)][swap(j)]
               for i in range(4 * modes) for j in range(4 * modes))
    return kernel


def reference_pfaffian(modes):
    kernel = reference_kernel(modes)

    @lru_cache(None)
    def pfaffian(mask):
        if not mask:
            return 1
        if mask.bit_count() % 2:
            return 0
        first = mask & -mask
        i = first.bit_length() - 1
        rest, result = mask ^ first, 0
        remaining = rest
        while remaining:
            chosen = remaining & -remaining
            remaining ^= chosen
            j = chosen.bit_length() - 1
            result ^= controls.binary_product(kernel[i][j], pfaffian(rest ^ chosen))
        return result

    return pfaffian


def reverse(polynomial, degree):
    assert polynomial.bit_length() - 1 <= degree
    return sum(((polynomial >> i) & 1) << (degree - i) for i in range(degree + 1))


def literal_entry(modes, row, column, pfaffian):
    planes = [0] * 4
    for rmask, a in row.items():
        for cmask, b in column.items():
            value = pfaffian(rmask | (cmask << (2 * modes)))
            coefficient = field.mul(a, b)
            for bit in range(4):
                if coefficient >> bit & 1:
                    planes[bit] ^= value
    return tuple(planes)


def check_map(code, label):
    p, left, modes = code['embedding'], code['retraction'], code['modes']
    heights, diagnostics = validate_height_maps(p, left, modes)
    q = len(p)
    pfaffian = reference_pfaffian(modes)
    records, entries, jets, cap_orders = [], 0, 0, 0
    matrices = []
    orders_by_cut = []
    for cut in ('zero', 'infinity'):
        dense = [[(0, 0, 0, 0) for _ in range(q*q)] for _ in range(q*q)]
        cut_orders = []
        for parity in (0, 1):
            pairs = [(i, j) for i in range(q) for j in range(q) if (heights[i] + heights[j]) % 2 == parity]
            rows = [field.tensor_vector(left[i], left[j] if cut == 'zero' else p[j], modes) for i, j in pairs]
            columns = [field.tensor_vector(p[j] if cut == 'zero' else left[j], p[i], modes) for i, j in pairs]
            native = literal.compound_block(modes, rows, columns, coefficients=True)
            reference = [[literal_entry(modes, row, column, pfaffian) for column in columns] for row in rows]
            assert native['matrix'] == reference
            checks = controls.toeplitz_controls(reference, native['orders_t'])
            for cap in (8, 24, 64):
                capped = literal.compound_block(modes, rows, columns, order_cap=cap)
                assert capped['orders_t'] == [min(cap, x) for x in native['orders_t']]
                cap_orders += len(pairs)
            for a, (i, j) in enumerate(pairs):
                for b, (k, l) in enumerate(pairs):
                    dense[q*i+j][q*k+l] = reference[a][b]
            entries += len(pairs)**2
            jets += len(checks)
            cut_orders.extend(native['orders_t'])
            records.append({'cut': cut, 'parity': parity, 'orders': native['orders_t'],
                            'toeplitz_checks': checks, 'cache': native['binary_pfaffian_cache']})
        matrices.append(dense)
        orders_by_cut.append(sorted(cut_orders))
    endpoint_checks, reciprocal_checks = 0, 0
    for x, y, u, v in product(range(q), repeat=4):
        zero = matrices[0][q*x+y][q*v+u]
        infinity = matrices[1][q*x+u][q*v+y]
        endpoint = sum((plane.bit_count() % 2) << bit for bit, plane in enumerate(zero))
        assert endpoint == int(x == v and y == u)
        total_height = heights[x] + heights[y] + heights[u] + heights[v]
        assert infinity == tuple(reverse(plane, total_height) for plane in zero)
        endpoint_checks += 1
        reciprocal_checks += 1
    if code.get('retraction_kind') == 'orthogonal':
        assert orders_by_cut[0] == orders_by_cut[1]
    if all(value['deletion_closed'] for value in diagnostics.values()):
        reference = field.profiles(p, left, modes)
        assert orders_by_cut[0] == reference['zero_orders_t']
        assert orders_by_cut[1] == reference['infinity_orders_t']
    return {'label': label, 'code': code, 'deletion_diagnostics': diagnostics,
            'literal_polynomial_entries': entries, 'independent_toeplitz_rank_checks': jets,
            'capped_order_comparisons': cap_orders, 'endpoint_checks': endpoint_checks,
            'infinity_reciprocity_checks': reciprocal_checks, 'cut_orders': orders_by_cut,
            'profiles': records, 'reference_pfaffian_cache_size': pfaffian.cache_info().currsize}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/tensor-literal-hot-pfaffians-f16-controls.json'))
    args = parser.parse_args()
    started = time.monotonic()
    edge_checks = 0
    for modes in (1, 4, 5, 8, 10, 16):
        kernel = reference_kernel(modes)
        for i in range(4*modes):
            for j in range(i+1, 4*modes):
                assert literal.literal_pfaffian(modes, (1 << i) | (1 << j)) == kernel[i][j]
                edge_checks += 1
    pfaffian = reference_pfaffian(4)
    for mask in range(1 << 16):
        assert literal.literal_pfaffian(4, mask) == pfaffian(mask)
    random_checks = 0
    rng = random.Random(20261101)
    for modes in (5, 8, 10):
        pfaffian = reference_pfaffian(modes)
        for _ in range(64):
            selected = rng.sample(range(4*modes), 2*rng.randrange(6))
            mask = sum(1 << i for i in selected)
            actual = literal.literal_pfaffian(modes, mask)
            assert actual == pfaffian(mask)
            jet = literal.literal_pfaffian(modes, mask, precision_words=1)
            assert jet == actual & ((1 << 64) - 1)
            random_checks += 1
    p = [{0: 1}] + [{1 << i: 1} for i in range(4)] + [{9: 1, 6: 4}]
    scale = field.INV[field.dot(p[-1], p[-1])]
    left = p[:-1] + [{w: field.mul(scale, c) for w, c in p[-1].items()}]
    coherent = {'modes': 4, 'embedding': p, 'retraction': left, 'retraction_kind': 'orthogonal'}
    noncoherent = boundary_code(10, 9, [2], modes=5)
    assert all(not value['deletion_closed'] for value in noncoherent['deletion_diagnostics'].values())
    cases = [check_map(coherent, 'coherent_control'), check_map(noncoherent, 'noncoherent_boundary')]
    names = [Path(__file__).name, 'tensor_literal_hot_pfaffians_f16.py',
             'tensor_literal_hot_pfaffians_f16.cpp', 'tensor_hot_boundary_retractions.py',
             'reflection_dual_rail_gaussian.py', 'check_tensor_hot_compounds_f16.py',
             'tensor_hot_compounds_f16_fast.cpp', 'tensor_hot_compounds_f16.py',
             'tensor_hot_compounds_f16.cpp', 'tensor_orthogonal_scaling.py',
             'coherent_hot_compounds_f4.cpp', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp']
    report = {'passed': True, 'scope': __doc__, 'original_kernel_edge_checks': edge_checks,
              'complete_four_mode_pfaffian_checks': 1 << 16,
              'larger_random_pfaffian_and_jet_checks': random_checks, 'cases': cases,
              'literal_polynomial_entries': sum(c['literal_polynomial_entries'] for c in cases),
              'independent_toeplitz_rank_checks': sum(c['independent_toeplitz_rank_checks'] for c in cases),
              'capped_order_comparisons': sum(c['capped_order_comparisons'] for c in cases),
              'endpoint_checks': sum(c['endpoint_checks'] for c in cases),
              'infinity_reciprocity_checks': sum(c['infinity_reciprocity_checks'] for c in cases),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'seconds': time.monotonic() - started}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps({key: value for key, value in report.items() if key not in ('cases', 'source_sha256')}))


if __name__ == '__main__':
    main()
