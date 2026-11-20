#!/usr/bin/env python3
"""Repair pair-lattice defects from an HNF annihilator.

Loads the committed window-4 sweep, takes the 77 targets with lattice index
above one, and computes the HNF covolume E of each pair-row lattice.  The
target order h is not enumerated.  For j=1,2,... the checker tests whether
E/2^j annihilates every target; the first j where it does not gives an odd
multiple a of h/2.  The half-power separates even from odd target logs, and
one candidate row a e_p or a(e_anchor-e_p) is added.

The first pass removes the 2-primary part of the index in 67 cases: 66 reach
index one and the index-six target drops to index three.  Eleven odd
residuals remain.  For them the checker divides a by the small odd divisor
three and retests; all eleven reach index one.  The finite sweep therefore
needs no cycle-length enumeration and no odd-unit data, but the odd divisor
bound is empirical.  Finite main-cycle diagnostic only; no uniform cost,
supply, or complete-bound claim.
"""
import argparse
import json

from arithmetic_two_group import power
from check_mixed_pair_gap_cycle import (
    CUTOFF, INTEGER, cycle_forms, determinant, replay, row_hnf_basis,
)
from check_mixed_pair_order_saturation import (
    W4_PATH, build_rows, extract_factor,
)
from binary_quadratic_forms import reduce_form
from leading_power_smooth_relations import compose_forms_hnf
from prime_form_factor_base import build_prime_form_factor_base


ODD_DIVISORS = (3,)


def covolume(rows, width):
    basis, rank = row_hnf_basis(rows, width)
    assert rank == width
    return abs(int(determinant(basis)))


def first_nonannihilator(annihilator, targets, records_all, identity):
    exponent_bits = (annihilator & -annihilator).bit_length()-1
    for level in range(1, exponent_bits+1):
        exponent = annihilator >> level
        values = {
            target: power(
                records_all[target][1], exponent, identity, compose_forms_hnf)
            for target in targets}
        if any(value != identity for value in values.values()):
            return level, exponent, values
    return None, None, None


def candidate_rows(targets, values, exponent, identity):
    even = {target: value == identity for target, value in values.items()}
    anchors = [target for target in targets if not even[target]]
    assert anchors
    anchor = anchors[0]
    positions = {target: position for position, target in enumerate(targets)}
    rows = []
    for position, target in enumerate(targets):
        if even[target]:
            row = [0]*len(targets)
            row[position] = exponent
            rows.append(('single', target, None, row))
        elif target != anchor:
            row = [0]*len(targets)
            row[positions[anchor]] = exponent
            row[position] = -exponent
            rows.append(('difference', anchor, target, row))
    return rows


def best_reduction(rows, annihilator, targets, records_all, identity,
                   exponent, width, period):
    candidates = candidate_rows(
        targets,
        {target: power(
            records_all[target][1], exponent, identity, compose_forms_hnf)
         for target in targets},
        exponent, identity)
    best = None
    for kind, left, right, row in candidates:
        if replay(row, targets, records_all, identity) != identity:
            continue
        new_annihilator = covolume(rows + [row], width)
        assert new_annihilator % period == 0
        if (new_annihilator < annihilator
                and annihilator % new_annihilator == 0):
            candidate = (
                new_annihilator//period, kind, left, right, row,
                new_annihilator)
            if best is None or candidate[0] < best[0]:
                best = candidate
    return best


def run():
    with open(W4_PATH) as handle:
        w4 = json.load(handle)
    w4_by_target = {record['target']: record for record in w4['results']}
    targets = sorted(
        target for target, record in w4_by_target.items()
        if record.get('lattice_multiple', 1) > 1)
    assert len(targets) == 77

    delta = -INTEGER
    base = build_prime_form_factor_base(delta, CUTOFF)
    records = {
        record.prime: (record.root, tuple(record.positive))
        for record in base['records']}
    base_primes = tuple(sorted(records))
    forms, index, identity = cycle_forms(records)
    period = len(forms)
    logs = {
        prime: index[reduce_form(records[prime][1])]
        for prime in base_primes}
    larger = build_prime_form_factor_base(delta, targets[-1])
    larger_records = {
        record.prime: (record.root, tuple(record.positive))
        for record in larger['records']}

    final_reductions = []
    two_primary = []
    two_primary_residual = []
    odd_repairs = []
    for target in targets:
        records_all = dict(records)
        records_all[target] = larger_records[target]
        targets_set, _, rows_by_window = build_rows(
            records, base_primes, forms, index, period, records_all, logs,
            target, 4)
        width = len(targets_set)
        rows = [row for window in range(1, 5)
                for row in rows_by_window.get(window, [])]
        annihilator = covolume(rows, width)
        original_index = annihilator//period
        assert original_index == w4_by_target[target]['lattice_multiple']

        level, exponent, values = first_nonannihilator(
            annihilator, targets_set, records_all, identity)
        assert level is not None
        best = best_reduction(
            rows, annihilator, targets_set, records_all, identity, exponent,
            width, period)
        two_primary_index = best[0] if best is not None else original_index
        two_primary.append({
            'target': target,
            'original_index': original_index,
            'level': level,
            'annihilator_bits': annihilator.bit_length(),
            'new_index': two_primary_index,
        })
        if two_primary_index != 1:
            two_primary_residual.append({
                'target': target,
                'original_index': original_index,
                'level': level,
                'annihilator_bits': annihilator.bit_length(),
                'new_index': two_primary_index,
            })

        row_divisor = 1
        if best is None or best[0] != 1:
            for divisor in ODD_DIVISORS:
                if exponent % divisor:
                    continue
                candidate = best_reduction(
                    rows, annihilator, targets_set, records_all, identity,
                    exponent//divisor, width, period)
                if candidate is not None and candidate[0] == 1:
                    best = candidate
                    row_divisor = divisor
                    break
        assert best is not None and best[0] == 1
        new_index, kind, left, right, row, new_annihilator = best
        factor = extract_factor(
            targets_set, records_all, identity, new_annihilator)
        assert factor in (1000003, 1000033)
        record = {
            'target': target,
            'original_index': original_index,
            'two_primary_index': two_primary_index,
            'new_index': new_index,
            'row_divisor': row_divisor,
            'row_kind': kind,
            'row_pair': [left, right] if right is not None else [left],
            'coefficient_bits': (exponent//row_divisor).bit_length(),
            'extracted_factor': factor,
        }
        final_reductions.append(record)
        if row_divisor != 1:
            odd_repairs.append(record)

    reached_one_two_primary = [
        record for record in two_primary if record['new_index'] == 1]
    assert len(reached_one_two_primary) == 66
    assert len(two_primary_residual) == 11
    assert sorted(record['target'] for record in two_primary_residual) == [
        229, 593, 617, 797, 1571, 2207, 2213, 2243, 4663, 5569, 6997]
    assert len(odd_repairs) == 11
    assert all(record['row_divisor'] == 3 for record in odd_repairs)
    assert all(record['new_index'] == 1 for record in final_reductions)
    assert all(record['extracted_factor'] in (1000003, 1000033)
               for record in final_reductions)

    return {
        'scope': (
            'Finite main-cycle annihilator repair: the window-4 HNF '
            'covolume supplies the half-cycle exponent without enumerating '
            'h. One row removes the 2-primary index for 67 of 77 targets '
            '(66 reach index one); dividing the exponent by the fixed small '
            'odd divisor three repairs the eleven odd residuals. The odd '
            'divisor bound is empirical, not uniform. No supply or complete-'
            'bound claim.'),
        'integer': INTEGER,
        'w4_targets': len(targets),
        'two_primary_reached_index_one': len(reached_one_two_primary),
        'two_primary_residual_targets': two_primary_residual,
        'two_primary_reductions': two_primary,
        'odd_divisors_tested': list(ODD_DIVISORS),
        'odd_divisor_repairs': odd_repairs,
        'targets_reached_index_one': len(final_reductions),
        'reductions': final_reductions,
        'h_enumeration_used_for_repair': False,
        'uniform_odd_bound_claim': False,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    arguments = parser.parse_args()
    result = run()
    with open(arguments.output, 'w') as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(json.dumps({
        'targets': result['w4_targets'],
        'two_primary_reached_index_one':
            result['two_primary_reached_index_one'],
        'two_primary_residual': [record['target']
                                 for record in
                                 result['two_primary_residual_targets']],
        'odd_divisor_repairs': [record['target']
                                for record in result['odd_divisor_repairs']],
        'final_reached_index_one': result['targets_reached_index_one'],
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
