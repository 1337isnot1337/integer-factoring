#!/usr/bin/env python3
"""Exact bounded leading-coefficient blocks for the P10 duplicated fan code."""

import argparse
from collections import Counter
import hashlib
from itertools import product
import json
from math import isqrt, log2
from pathlib import Path
import random

from check_modular_crossing import rank_two
from check_ward_fan_entry_queries import restricted_interface, wedge
from check_ward_fan_family import fan_types
from check_ward_stable_operators import bits, local_operators
from complement_shell_tensor import multiply
from explore_complement_edge_pencils import evaluate_local
from pfaffian_particle_hole import particle_hole_mask
from ward_fan_entry_queries import fan_entry_query


FORMS, PAIRS, LENGTH = 3, 10, 10
ROOT = Path(__file__).resolve().parents[1]


def build_code():
    types, _ = fan_types(FORMS)
    interface = restricted_interface()
    labels = []
    for sector in interface['sectors']:
        selected = sector['selected_type_indices']
        even = [2 * i + copy for i in selected for copy in range(2)]
        singles = [i for i in range(PAIRS) if i not in even]
        weights = [sum(((types[site // 2] >> a) & 1) << j for j, site in enumerate(even))
                   for a in range(FORMS)]
        if len(selected) == FORMS:
            preimages = [{0}]
            images = [wedge(weights, len(even))]
        else:
            values = [types[i] for i in selected]
            cycle = next(mask for mask in range(1, 1 << len(selected))
                         if all(sum(((values[i] >> a) & 1) for i in range(len(values))
                                    if mask >> i & 1) % 2 == 0 for a in range(FORMS)))
            edge = next(i for i in range(len(selected)) if cycle >> i & 1)
            residuals = [3 << (2 * edge),
                         sum(1 << (2 * i + 1) for i in range(len(selected)) if cycle >> i & 1)]
            preimages = [{1 << j for j in bits(residual)} for residual in residuals]
            images = [wedge(weights + [residual], len(even)) for residual in residuals]
        assert list(map(len, preimages)) == sector['preimage_support_sizes']
        assert list(map(len, images)) == sector['input_support_sizes']
        for orientation in range(sector['orientation_count']):
            base = 1 << (2 * singles[0])
            for j, site in enumerate(singles[1:]):
                base |= 1 << (2 * site + ((orientation >> j) & 1))

            def embed(mask):
                return base | sum(3 << (2 * even[j]) for j in bits(mask))

            for residual, (before, after, pivot) in enumerate(
                    zip(preimages, images, sector['coordinate_output_pivots'])):
                source = {embed(mask) for mask in before}
                image = {embed(mask) for mask in after}
                # Literal product Q_1 Q_2 Q_3, including the local diagonal terms.
                actual = set(source)
                for a in range(FORMS):
                    following = set()
                    for word in actual:
                        for site in range(PAIRS):
                            local = word >> (2 * site) & 3
                            if not (types[site // 2] >> a & 1) or local in (1, 2):
                                continue
                            for term in (word, word ^ (3 << (2 * site))):
                                if term in following:
                                    following.remove(term)
                                else:
                                    following.add(term)
                    actual = following
                assert actual == image
                assert all(word.bit_count() == LENGTH - 2 * FORMS for word in source)
                assert all(word.bit_count() == LENGTH for word in image)
                labels.append({'type_syndrome': sector['type_syndrome'],
                               'orientation': orientation, 'residual': residual,
                               'preimage': sorted(source), 'input': sorted(image),
                               'output': embed(pivot)})
    assert len(labels) == 84
    for i, first in enumerate(labels):
        for j, second in enumerate(labels):
            assert (first['output'] in second['input']) == (i == j)
            assert len(set(first['input']) & set(second['input'])) % 2 == 0
    return labels


def local_query_table(p, q, r, eta):
    table = evaluate_local((p, q, r, 1, 1, 1))
    ordinary = [[0] * 16 for _ in range(16)]
    for x, y, u, v in product(range(4), repeat=4):
        ordinary[particle_hole_mask(x, y)][particle_hole_mask(v, u)] = table[((x * 4 + y) * 4 + u) * 4 + v]
    incidence, creation = local_operators()[:2]
    acted = [[0] * 16 for _ in range(16)]
    for source in range(16):
        for target in bits(incidence[source]):
            for col in range(16):
                acted[target][col] ^= ordinary[source][col]
        for target in bits(creation[source]):
            for col in range(16):
                acted[target][col] ^= multiply(eta, ordinary[source][col])
    return [(table[((x * 4 + y) * 4 + u) * 4 + v],
             acted[particle_hole_mask(x, y)][particle_hole_mask(v, u)])
            for x, y, u, v in product(range(4), repeat=4)]


class RestrictedOracle:
    def __init__(self, code, local=None, boundary=True):
        self.code = code
        self.local = local_query_table(3, 1, 1, 0) if local is None else local
        self.boundary = boundary
        self.pair_cache = {}
        self.physical_cache = {}
        self.logical_cache = {}
        self.counts = Counter()

    def physical(self, x, y, u, v):
        self.counts['physical_query_requests'] += 1
        keys = []
        for edge in range(2 * FORMS - 1):
            entries = []
            for copy in range(2):
                site = 2 * edge + copy
                index = (((x >> (2 * site) & 3) * 4 + (y >> (2 * site) & 3)) * 4
                         + (u >> (2 * site) & 3)) * 4 + (v >> (2 * site) & 3)
                entries.append(self.local[index])
            key = tuple(sorted(entries))
            if key not in self.pair_cache:
                (m1, b1), (m2, b2) = key
                self.pair_cache[key] = (multiply(m1, m2),
                                        multiply(b1, m2) ^ multiply(m1, b2))
            value = self.pair_cache[key]
            if value == (0, 0):
                self.counts['zero_local_pair_queries'] += 1
                return 0
            keys.append(value)
        key = tuple(keys)
        if key not in self.physical_cache:
            # (M,B),(1,0) has the same collapsed pair as the original sites.
            sites = [value for pair in key for value in (pair, (1, 0))]
            answer = fan_entry_query(FORMS, iter(sites), multiply)
            assert answer < 1 << 8
            self.physical_cache[key] = answer
            self.counts['distinct_physical_polynomials_computed'] += 1
        return self.physical_cache[key]

    def polynomial(self, x, y, u, v):
        key = x, y, u, v
        if key not in self.logical_cache:
            result = 0
            for yy in self.code[y]['preimage']:
                for vv in self.code[v]['input']:
                    result ^= self.physical(self.code[x]['output'], yy,
                                            self.code[u]['output'], vv)
            # The code has joint order six and the displayed Ward factor is y^3.
            if self.boundary:
                assert result & 7 == 0
            assert result < 1 << 8
            self.logical_cache[key] = result
        return self.logical_cache[key]


def independent_original_entry(code, x, y, u, v, local):
    # Work directly in F2[x,y]/(x^8,y^4), with stride sixteen to avoid carries.
    mask = sum(1 << (a + 16 * b) for a in range(8) for b in range(4))
    xx, uu = code[x]['output'], code[u]['output']
    result = 0
    for yy in code[y]['input']:
        for vv in code[v]['input']:
            value = 1
            for site in range(PAIRS):
                index = (((xx >> (2 * site) & 3) * 4 + (yy >> (2 * site) & 3)) * 4
                         + (uu >> (2 * site) & 3)) * 4 + (vv >> (2 * site) & 3)
                value = multiply(value, local[index]) & mask
                if not value:
                    break
            result ^= value
    return result


def high_pivot_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    return len(pivots)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--syndrome', type=int, default=7)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    code = build_code()
    positions = [(a, b) for a in range(len(code)) for b in range(len(code))
                 if code[a]['type_syndrome'] ^ code[b]['type_syndrome'] == args.syndrome]
    assert len(positions) <= 592
    oracle = RestrictedOracle(code)
    reports = []
    rng = random.Random(829117)
    original_local = evaluate_local((3, 1 ^ (1 << 16), 1, 1, 1, 1))
    truncation = sum(1 << (a + 16 * b) for a in range(8) for b in range(4))
    original_local = [value & truncation for value in original_local]
    controls = []
    endpoint = RestrictedOracle(code, local_query_table(0, 0, 1, 1), boundary=False)
    endpoint_local = evaluate_local((0, 0, 1, 1, 1, 1))
    endpoint_controls = []
    for x, y in rng.sample(positions, min(16, len(positions))):
        assert endpoint.polynomial(x, y, y, x) == 1
        assert independent_original_entry(code, x, y, y, x, endpoint_local) == 1
        endpoint_controls.append({'labels': [x, y, y, x], 'value': 1})
        v, u = rng.choice([pair for pair in positions if pair != (x, y)])
        assert endpoint.polynomial(x, y, u, v) == 0
        assert independent_original_entry(code, x, y, u, v, endpoint_local) == 0
        endpoint_controls.append({'labels': [x, y, u, v], 'value': 0})
    for cut in ('row', 'column'):
        rows = []
        polynomials = []
        for a, b in positions:
            row, values = 0, []
            for j, (c, d) in enumerate(positions):
                labels = (a, b, d, c) if cut == 'row' else (a, d, b, c)
                value = oracle.polynomial(*labels)
                values.append(value)
                row |= ((value >> 3) & 1) << j
            rows.append(row)
            polynomials.append(values)
        rank = rank_two(rows)
        assert rank == high_pivot_rank(rows)
        nonzero_positions = [(i, j) for i, row in enumerate(rows) for j in bits(row)]
        chosen = [(rng.randrange(len(rows)), rng.randrange(len(rows))) for _ in range(16)]
        if nonzero_positions:
            chosen.extend(rng.sample(nonzero_positions, min(16, len(nonzero_positions))))
        for i, j in chosen:
            a, b = positions[i]
            c, d = positions[j]
            labels = (a, b, d, c) if cut == 'row' else (a, d, b, c)
            actual = independent_original_entry(code, *labels, original_local)
            expected = polynomials[i][j] << 48
            assert actual == expected, (cut, i, j, actual, expected)
            controls.append({'cut': cut, 'row': i, 'column': j,
                             'univariate_quotient': polynomials[i][j]})
        boundary_ranks = []
        for degree in range(3, 8):
            coefficient = [sum(((value >> degree) & 1) << j for j, value in enumerate(values))
                           for values in polynomials]
            actual_rank = rank_two(coefficient)
            assert actual_rank == high_pivot_rank(coefficient)
            boundary_ranks.append({'centered_coefficient': [degree, 3], 'rank': actual_rank,
                                   'complete_coefficient_rank_upper_bound_using_only_this_block':
                                       len(code) ** 2 - len(positions) + actual_rank})
        reports.append({'cut': cut, 'block_dimension': len(positions), 'rank': rank,
                        'nonzero_entries': sum(row.bit_count() for row in rows),
                        'boundary_coefficient_ranks': boundary_ranks,
                        'whole_univariate_quotient_block_zero': not any(value for values in polynomials for value in values),
                        'coefficient_rows': [hex(row) for row in rows],
                        'univariate_quotient_rows': polynomials,
                        'complete_coefficient_rank_upper_bound_using_only_this_block': len(code) ** 2 - len(positions) + rank})
        print(json.dumps({key: value for key, value in reports[-1].items()
                          if key not in ('coefficient_rows', 'univariate_quotient_rows')}), flush=True)
    old_profile_path = ROOT / 'evidence/ward-fan-first-profile.json'
    old = next(record for record in json.loads(old_profile_path.read_text()) if record['orientation_selected'])
    assert (old['pairs'], old['length'], old['forms'], old['labels']) == (PAIRS, LENGTH, FORMS, len(code))
    boundary = {record['cut']: {item['centered_coefficient'][0]:
                                item['complete_coefficient_rank_upper_bound_using_only_this_block']
                                for item in record['boundary_coefficient_ranks']} for record in reports}
    updated = []
    for a, b, row, column in old['profile']:
        if b == 3 and a in boundary['row']:
            row, column = min(row, boundary['row'][a]), min(column, boundary['column'][a])
        updated.append([a, b, row, column])
    def rank_sum(profile):
        total = 0
        for _, _, row, column in profile:
            product_value = row * column
            root = isqrt(product_value)
            total += root + int(root * root < product_value)
        return total
    assert rank_sum(old['profile']) == old['integer_rank_sum_upper_bound'] == 192340
    new_cost = rank_sum(updated)
    assert new_cost <= old['integer_rank_sum_upper_bound']
    if args.syndrome == 7:
        assert new_cost == 191700
    names = (Path(__file__).name, 'ward_fan_entry_queries.py', 'check_ward_fan_entry_queries.py',
             'check_ward_fan_family.py', 'check_ward_stable_operators.py',
             'complement_shell_tensor.py', 'explore_complement_edge_pencils.py',
             'pfaffian_particle_hole.py')
    result = {'source_sha256': {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                               for name in names},
              'pairs': PAIRS, 'forms': FORMS, 'length': LENGTH, 'labels': len(code),
              'centered_coefficient': [3, 3], 'type_pair_syndrome': args.syndrome,
              'basis': code, 'logical_pair_positions': positions, 'cuts': reports,
              'independent_original_kernel_controls': controls,
              'independent_crossing_endpoint_controls': endpoint_controls,
              'crossing_endpoint_query_counts': dict(endpoint.counts),
              'query_counts': dict(oracle.counts),
              'complete_upper_profile': {
                  'prior_evidence_sha256': hashlib.sha256(old_profile_path.read_bytes()).hexdigest(),
                  'previous_integer_rank_sum_upper_bound': old['integer_rank_sum_upper_bound'],
                  'integer_rank_sum_upper_bound': new_cost,
                  'rho_of_upper_bound': log2(new_cost) / (2 * log2(len(code))),
                  'profile': updated},
              'scope': 'Exact ranks of the specified restricted coefficient blocks and the resulting complete upper profile; no factoring-score improvement asserted.'}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(json.dumps({'query_counts': result['query_counts'], 'original_kernel_controls': len(controls),
                      'crossing_endpoint_controls': len(endpoint_controls)}), flush=True)


if __name__ == '__main__':
    main()
