#!/usr/bin/env python3
"""Derive the second-cycle layered repair from the HNF annihilator.

The second cycle has period ``h=20456=8*2557`` and four window-4 residual
targets of lattice index two.  This checker computes the HNF covolume
``E=covol(L)=h*2`` and walks the 2-adic levels ``a=E/2^j``.  At each level it
groups targets by the form value ``g_p^a`` and tests difference rows
``a(e_p-e_q)`` inside a nontrivial class, or single rows on the identity
class.  No cycle length and no target log is used for the construction.

Levels one through three do not reduce the index.  At level four,
``a=h/8``, the class rows are primitive layer-three differences and each
repairs its target to index one.  The repair row is selected from the base
pair rows and group powers; no cycle length or additional target-log data is
used.  All rows replay and the repaired covolume annihilators extract
``20323``.  Finite higher-2-adic control only; the class bound and uniform
cost remain open.
"""
import argparse
import json
from math import gcd

from arithmetic_two_group import power
from binary_quadratic_forms import reduce_form
from check_mixed_pair_gap_cycle import determinant, replay, row_hnf_basis
from check_mixed_pair_second_cycle import (
    BASE_PRIME, CUTOFF, EXPECTED_PERIOD, EXPECTED_RESIDUAL, INTEGER,
    TARGET_LIMIT, build_rows, cycle_forms, extract,
)
from leading_power_smooth_relations import compose_forms_hnf
from prime_form_factor_base import build_prime_form_factor_base, primes_through


def covolume(rows, width):
    basis, rank = row_hnf_basis(rows, width)
    assert rank == width
    return abs(int(determinant(basis)))


def class_rows(targets, records_all, exponent, identity):
    values = {
        target: power(
            records_all[target][1], exponent, identity, compose_forms_hnf)
        for target in targets}
    groups = {}
    for position, target in enumerate(targets):
        groups.setdefault(values[target], []).append(position)
    rows = []
    for value, positions in groups.items():
        if value == identity:
            for position in positions:
                row = [0]*len(targets)
                row[position] = exponent
                rows.append(('single', targets[position], None, row))
        else:
            for position in positions[1:]:
                row = [0]*len(targets)
                row[positions[0]] = exponent
                row[position] = -exponent
                rows.append((
                    'difference', targets[positions[0]],
                    targets[position], row))
    return rows


def best_reduction(rows, annihilator, targets, records_all, identity,
                   exponent, width, period):
    best = None
    for kind, left, right, row in class_rows(
            targets, records_all, exponent, identity):
        if replay(row, targets, records_all, identity) != identity:
            continue
        new_annihilator = covolume(rows + [row], width)
        assert new_annihilator % period == 0
        if (new_annihilator < annihilator
                and annihilator % new_annihilator == 0):
            candidate = (
                new_annihilator//period, kind, left, right, row,
                new_annihilator, exponent)
            if best is None or candidate[0] < best[0]:
                best = candidate
    return best


def run():
    delta = -INTEGER
    base = build_prime_form_factor_base(delta, CUTOFF)
    records = {
        record.prime: (record.root, tuple(record.positive))
        for record in base['records']}
    base_primes = tuple(sorted(records))
    primes = tuple(primes_through(CUTOFF))
    forms, index, identity = cycle_forms(delta, records, BASE_PRIME)
    period = len(forms)
    assert period == EXPECTED_PERIOD
    v2_period = (period & -period).bit_length()-1
    assert v2_period == 3
    logs = {
        prime: index[reduce_form(records[prime][1])]
        for prime in base_primes}
    extension = build_prime_form_factor_base(
        delta, max(EXPECTED_RESIDUAL))
    extension_records = {
        record.prime: (record.root, tuple(record.positive))
        for record in extension['records']}

    results = []
    for target in EXPECTED_RESIDUAL:
        targets = tuple(sorted(base_primes + (target,)))
        target_logs = dict(logs)
        target_logs[target] = index[
            reduce_form(extension_records[target][1])]
        records_all = dict(records)
        records_all[target] = extension_records[target]
        rows = build_rows(
            INTEGER, records, forms, index, period, records_all,
            target_logs, targets, primes)
        width = len(targets)
        annihilator = covolume(rows, width)
        original_index = annihilator//period
        assert original_index == 2

        level = None
        best = None
        for candidate_level in range(
                1, (annihilator & -annihilator).bit_length()):
            candidate = best_reduction(
                rows, annihilator, targets, records_all, identity,
                annihilator >> candidate_level, width, period)
            if candidate is not None:
                level = candidate_level
                best = candidate
                break
        assert level is not None and best is not None
        new_index, kind, left, right, row, new_annihilator, exponent = best
        assert new_index == 1
        assert level == 4
        basis, rank = row_hnf_basis(rows + [row], width)
        assert rank == width
        repaired_covolume = abs(int(determinant(basis)))
        assert repaired_covolume == period
        for basis_row in basis:
            assert replay(
                basis_row, targets, records_all, identity) == identity
        factor = extract(
            records_all, targets, basis, repaired_covolume, period, identity)
        assert factor in (20323, 22157)
        results.append({
            'target': target,
            'target_order': period//gcd(period, target_logs[target]),
            'original_index': original_index,
            'annihilator_bits': annihilator.bit_length(),
            'first_reduction_level': level,
            'coefficient_bits': exponent.bit_length(),
            'class_count': len({
                power(records_all[p][1], exponent, identity,
                      compose_forms_hnf)
                for p in targets}),
            'row_kind': kind,
            'row_pair': [left, right],
            'new_index': new_index,
            'repaired_covolume': repaired_covolume,
            'extracted_factor': factor,
        })

    assert [record['target'] for record in results] == list(EXPECTED_RESIDUAL)
    assert all(record['first_reduction_level'] == 4 for record in results)
    assert all(record['new_index'] == 1 for record in results)
    assert all(record['coefficient_bits'] == 12 for record in results)
    assert all(record['extracted_factor'] == 20323 for record in results)
    assert [record['row_pair'] for record in results] == [
        [5, 79], [31, 653], [2, 3911], [31, 17321]]

    return {
        'scope': (
            'Finite second-cycle annihilator repair: the HNF covolume '
            'E=h*2 drives a walk over 2-adic levels. Class rows grouped by '
            'g_p^{E/2^j} leave the four index-two targets unchanged at '
            'levels one through three and reach index one at level four '
            '(a=h/8). The repair selection uses the base pair rows and group '
            'powers, not the cycle length or additional target-log data. The '
            'class bound and uniform cost remain open.'),
        'integer': INTEGER,
        'period': period,
        'v2_period': v2_period,
        'cutoff': CUTOFF,
        'base_prime': BASE_PRIME,
        'residual_targets': list(EXPECTED_RESIDUAL),
        'results': results,
        'h_enumeration_used_for_repair': False,
        'target_logs_used_for_repair': False,
        'uniform_cost_claim': False,
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
        'targets': [record['target'] for record in result['results']],
        'first_reduction_levels': [
            record['first_reduction_level'] for record in result['results']],
        'new_indexes': [record['new_index'] for record in result['results']],
        'factors': [record['extracted_factor']
                    for record in result['results']],
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
