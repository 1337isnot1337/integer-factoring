#!/usr/bin/env python3
"""Close the finite per-target pair-lattice index by allowed row types.

Loads the committed window-4 sweep, takes every target whose HNF lattice
index in the relation kernel is above one, and tests two repairs:

* one order-derived row: an odd-order diagonal ``ord(p)e_p`` or an even-order
  half-order difference ``(ord(p)/2)e_p - (ord(q)/2)e_q``;
* for the remaining targets, pair rows at the first later exponent window.

The cycle order is ``198850 = 2*5^2*41*97``, so order-derived rows can only
remove 2-primary (and other divisor-primary) defects.  A 3-primary defect is
not order-derived and is tested separately.  Every added row is replayed.
Finite diagnostic only; no uniform supply, cost, or complete-bound claim.
"""
import argparse
import itertools
import json
from math import gcd

from arithmetic_two_group import power, project_to_two_part, two_group_basis
from binary_quadratic_forms import ambiguous_factor, inverse_form, reduce_form
from check_mixed_pair_gap_cycle import (
    BASE_PRIME, CUTOFF, INTEGER, cycle_forms, determinant, replay,
    row_hnf_basis, split_factorization,
)
from dynamic_leading_norm_relations import leading_norm_word
from leading_power_smooth_relations import compose_with_prime_form
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf


W4_PATH = 'evidence/mixed-pair-rank-sweep-w4.json'
MAX_REPAIR_WINDOW = 7


def build_rows(records, base_primes, forms, index, period, records_all, logs,
               target, max_window):
    targets = tuple(sorted(base_primes + (target,)))
    positions = {prime: position for position, prime in enumerate(targets)}
    target_logs = dict(logs)
    target_logs[target] = index[reduce_form(records_all[target][1])]
    rows_by_window = {}
    seen = set()
    for left_index, left in enumerate(targets):
        for right in targets[left_index:]:
            for exponent in range(1, max_window+1):
                form = forms[(exponent*(target_logs[left]+target_logs[right]))
                             % period]
                leading = form[0]
                if leading == 1:
                    continue
                parts = split_factorization(leading, records)
                if parts is None:
                    continue
                terms, new = leading_norm_word(
                    INTEGER, form, parts, records_all)
                assert not new
                row = [0]*len(targets)
                row[positions[left]] += exponent
                row[positions[right]] += exponent
                for prime, value in terms.items():
                    row[positions[prime]] -= value
                key = tuple(row)
                if any(row) and key not in seen:
                    seen.add(key)
                    rows_by_window.setdefault(exponent, []).append(row)
    return targets, target_logs, rows_by_window


def lattice_index(rows, width, period):
    basis, rank = row_hnf_basis(rows, width)
    assert rank == width
    covolume = abs(int(determinant(basis)))
    assert covolume % period == 0
    return covolume//period


def extract_factor(targets_set, records_all, identity, period):
    generators = [records_all[prime][1] for prime in targets_set]
    projected = project_to_two_part(
        generators, period, identity, compose_forms_hnf)
    basis = two_group_basis(
        projected, identity, compose_forms_hnf, inverse_form,
        period.bit_length())
    tops = [power(form, 1 << (power_order-1), identity, compose_forms_hnf)
            for form, power_order in basis]
    factors = [ambiguous_factor(form, INTEGER) for form in tops]
    factor = next((value for value in factors if value), None)
    assert factor in (1000003, 1000033) and INTEGER % factor == 0
    return factor


def run():
    with open(W4_PATH) as handle:
        w4 = json.load(handle)
    w4_by_target = {record['target']: record for record in w4['results']}
    targets = sorted(target for target, record in w4_by_target.items()
                     if record.get('lattice_multiple', 1) > 1)
    assert targets

    delta = -INTEGER
    base = build_prime_form_factor_base(delta, CUTOFF)
    records = {record.prime: (record.root, tuple(record.positive))
               for record in base['records']}
    base_primes = tuple(sorted(records))
    forms, index, identity = cycle_forms(records)
    period = len(forms)
    logs = {prime: index[reduce_form(records[prime][1])]
            for prime in base_primes}
    larger = build_prime_form_factor_base(delta, targets[-1])
    larger_records = {record.prime: (record.root, tuple(record.positive))
                      for record in larger['records']}

    order_records = []
    pair_records = []
    for target in targets:
        records_all = dict(records)
        records_all[target] = larger_records[target]
        targets_set, target_logs, rows_by_window = build_rows(
            records, base_primes, forms, index, period, records_all, logs,
            target, 4)
        width = len(targets_set)
        base_rows = [row for window in range(1, 5)
                     for row in rows_by_window.get(window, [])]
        original_index = lattice_index(base_rows, width, period)
        assert original_index == w4_by_target[target]['lattice_multiple']
        orders = {prime: period//gcd(period, target_logs[prime])
                  for prime in targets_set}
        candidates = []
        for position, prime in enumerate(targets_set):
            order = orders[prime]
            if order % 2:
                row = [0]*width
                row[position] = order
                candidates.append(('odd-diagonal', prime, None, order, row))
        even_positions = [
            (position, prime) for position, prime in enumerate(targets_set)
            if orders[prime] % 2 == 0]
        for (left_pos, left), (right_pos, right) in itertools.combinations(
                even_positions, 2):
            row = [0]*width
            row[left_pos] = orders[left]//2
            row[right_pos] = -orders[right]//2
            candidates.append(('half-difference', left, right,
                               orders[left]//2, row))
        chosen = None
        for kind, left, right, coefficient, row in candidates:
            if replay(row, targets_set, records_all, identity) != identity:
                continue
            candidate_index = lattice_index(
                base_rows + [row], width, period)
            if candidate_index == 1:
                chosen = {
                    'mechanism': kind,
                    'left': left,
                    'right': right,
                    'coefficient': coefficient,
                    'row': row,
                    'final_index': 1,
                    'extracted_factor': extract_factor(
                        targets_set, records_all, identity, period),
                }
                break
        if chosen is not None:
            order_records.append({
                'target': target,
                'w4_index': original_index,
                **chosen,
            })
            continue

        repair_rows_by_window = build_rows(
            records, base_primes, forms, index, period, records_all, logs,
            target, MAX_REPAIR_WINDOW)[2]
        repair_rows = []
        repair_window = None
        for window in range(1, MAX_REPAIR_WINDOW+1):
            repair_rows.extend(repair_rows_by_window.get(window, []))
            current = lattice_index(repair_rows, width, period)
            if current == 1:
                repair_window = window
                break
        assert repair_window is not None
        pair_records.append({
            'target': target,
            'w4_index': original_index,
            'mechanism': 'pair-window',
            'pair_window': repair_window,
            'final_index': 1,
            'extracted_factor': extract_factor(
                targets_set, records_all, identity, period),
        })

    final_indices = ([record['final_index'] for record in order_records]
                     + [record['final_index'] for record in pair_records])
    assert len(final_indices) == len(targets)
    assert set(final_indices) == {1}
    return {
        'scope': ('Finite per-target lattice-index closure on the 40-bit '
                  'order-198850 cycle. Every window-4 target with index above '
                  'one is repaired to index one by one order-derived row or '
                  'by pair rows at a later exponent window. No uniform '
                  'repair, cost, supply, or complete-bound claim.'),
        'cycle_order': period,
        'w4_targets': len(targets),
        'order_saturated': order_records,
        'order_saturated_count': len(order_records),
        'pair_saturated': pair_records,
        'pair_saturated_count': len(pair_records),
        'final_index_distribution': {str(value): final_indices.count(value)
                                     for value in sorted(set(final_indices))},
        'max_final_index': max(final_indices),
        'pair_window_distribution': {
            str(window): sum(1 for record in pair_records
                             if record['pair_window'] == window)
            for window in sorted({record['pair_window']
                                  for record in pair_records})},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output')
    arguments = parser.parse_args()
    receipt = run()
    text = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    if arguments.output:
        with open(arguments.output, 'w') as handle:
            handle.write(text)
    print(json.dumps({
        'w4_targets': receipt['w4_targets'],
        'order_saturated_count': receipt['order_saturated_count'],
        'pair_saturated_count': receipt['pair_saturated_count'],
        'pair_window_distribution': receipt['pair_window_distribution'],
        'final_index_distribution': receipt['final_index_distribution'],
        'max_final_index': receipt['max_final_index'],
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
