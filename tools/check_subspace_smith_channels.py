#!/usr/bin/env python3
"""Exact action-set and large-channel controls for the subspace oracle."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
from random import Random

import check_cancellation_module_automata as arithmetic
import check_exact_module_orbits as lattices
from dvr_smith_channels import DVRSmithLattices
from matroid_smith_channels import binary_plane_lattices
from subspace_smith_channels import SubspaceSmithLattices, regular_combination


def apply(matrix, vector):
    return tuple(sum(value * vector[target] for target, value in row.items())
                 for row in matrix)


def span(basis):
    result = {0}
    for vector in basis:
        result |= {old ^ vector for old in tuple(result)}
    return result


def shift_contained(source, target, shift):
    try:
        for row in source:
            lattices.membership(tuple(arithmetic.rshift(value, shift) for value in row), target)
    except AssertionError:
        return False
    return True


def compare_actions(graph, exhaustive):
    """Close independently enumerated exact targets under safe inclusions."""
    assert graph.states == exhaustive.states
    labels = residues = lift_memberships = 0
    inclusion_cache = {}
    expanded_families = []
    for state_index, family in enumerate(exhaustive.candidates):
        expected = {}
        for residue, actual_target, actual_order in family:
            for order in range(actual_order + 1):
                for target, basis in enumerate(graph.states):
                    key = actual_target, target, actual_order - order
                    if key not in inclusion_cache:
                        inclusion_cache[key] = shift_contained(
                            graph.states[actual_target], basis, actual_order-order,
                        )
                    if inclusion_cache[key]:
                        expected.setdefault((target, order), set()).add(residue)
        actual = {}
        for action in graph.actions[state_index]:
            label = action["target"], action["order"]
            actual[label] = span([residue for residue, _ in action["basis"]]) - {0}
            for residue, coefficients in action["basis"]:
                left = regular_combination(coefficients, graph.local_factors[state_index], action["order"])
                for start in range(0, len(left), graph.dimension):
                    lattices.membership(left[start:start+graph.dimension], graph.states[action["target"]])
                    lift_memberships += 1
        assert actual == expected
        labels += len(actual)
        residues += sum(len(values) for values in actual.values())
        expanded_families.append(actual)
    return {
        "matched_nonzero_label_spaces": labels,
        "matched_labeled_residues": residues,
        "polynomial_lift_memberships": lift_memberships,
        "exhaustive_jet_vectors": exhaustive.jet_vectors,
    }, expanded_families


def exhaustive_basis_cost(family, rank, potential, point, endpoint_degree):
    costs = {}
    for (target, order), residues in family.items():
        value = point ** (endpoint_degree-order) * potential[target]
        for residue in residues:
            costs[residue] = min(value, costs.get(residue, value))
    return min(sum(costs[residue] for residue in basis)
               for basis in combinations(costs, rank)
               if len(arithmetic.binary_basis(basis)) == rank)


def companion(alphabet):
    local = {}
    for x, y, u, v in product(range(alphabet), repeat=4):
        straight = int(x == u and y == v)
        crossed = int(x == v and y == u)
        local[x, y, u, v] = [
            [straight ^ ((straight ^ crossed) << 1), straight],
            [straight ^ (straight << 1), straight],
        ]
    return local


def check():
    assert companion(2) == lattices.companion_two()
    random = Random(20261031)
    action_controls, basis_controls = [], 0
    total_entries = total_memberships = 0
    for precision in (1, 2):
        states = binary_plane_lattices(precision)
        for cut, infinity in ((0, True), (1, False)):
            graph = SubspaceSmithLattices(companion(2), 1, seed_states=states,
                                        precision=precision, cut=cut, infinity=infinity)
            exhaustive = DVRSmithLattices(companion(2), 1, seed_states=states,
                                         precision=precision, cut=cut, infinity=infinity)
            comparison, families = compare_actions(graph, exhaustive)
            point = Fraction(5, 2)
            potential = tuple(random.randrange(1, 101) for _ in states)
            matrix, choices = graph.choose(potential, point)
            image = apply(matrix, potential)
            for state, family in enumerate(families):
                assert image[state] == exhaustive_basis_cost(
                    family, len(graph.local_factors[state]), potential, point,
                    graph.degree if infinity else 0,
                )
                basis_controls += 1
            certificate = graph.verify(choices)
            total_entries += certificate["local_entry_checks"]
            total_memberships += certificate["regular_membership_checks"]
            action_controls.append({"graph": graph.summary(), **comparison,
                                    "potential": potential, "point": str(point),
                                    "minimum_row_costs": list(map(str, image))})

    states, point = binary_plane_lattices(2), Fraction(5)
    row = SubspaceSmithLattices(companion(2), 1, seed_states=states, cut=0, infinity=True)
    column = SubspaceSmithLattices(companion(2), 1, seed_states=states, cut=1, infinity=False)
    row_potential = (10000, 5320, 10000, 10000, 5320, 10000, 5320, 10000, 10000, 10000)
    col_potential = (100000, 94634, 63719, 62303, 94634, 62828, 94634, 63719, 54764, 61230)
    rates = Fraction(266, 25), Fraction(241, 200)
    companion_certificates = []
    for graph, potential, rate in zip((row, column), (row_potential, col_potential), rates):
        matrix, choices = graph.choose(potential, point)
        image = apply(matrix, potential)
        assert all(value <= rate * height for value, height in zip(image, potential))
        verified = graph.verify(choices)
        total_entries += verified["local_entry_checks"]
        total_memberships += verified["regular_membership_checks"]
        companion_certificates.append({
            "graph": graph.summary(), "potential": potential, "rate": str(rate),
            "maximum_row_ratio": str(max(value/height for value, height in zip(image, potential))),
            "policy": choices,
            "transfer": [[str(values.get(j, Fraction(0))) for j in range(len(states))]
                         for values in matrix],
        })

    # Rank sixteen, precision two: an explicit expanded residue/jet family
    # is intentionally never created. Every chosen rational identity is
    # nevertheless reconstructed in full.
    large = SubspaceSmithLattices(companion(4), 1, seed_states=states,
                                 alphabet=4, precision=2, cut=1, infinity=False)
    assert all(len(factors) == 16 for factors in large.local_factors)
    matrix, choices = large.choose((1,) * len(states), point)
    verified = large.verify(choices)
    assert large.maximum_variables == 32
    assert verified["local_entry_checks"] == 10240
    assert verified["regular_membership_checks"] == 5120
    total_entries += verified["local_entry_checks"]
    total_memberships += verified["regular_membership_checks"]
    return {
        "status": "pass",
        "scope": "Polynomial-size exact oracle for the declared channel/lattice model; no stronger factoring exponent",
        "exhaustive_action_controls": action_controls,
        "independent_all_basis_minimum_controls": basis_controls,
        "companion_certificates": companion_certificates,
        "retained_companion_product_bound": str(rates[0] * rates[1]),
        "rank_sixteen_control": {
            "graph": large.summary(), "point": str(point),
            "local_ranks": [len(factors) for factors in large.local_factors],
            "maximum_row_sum": str(max(sum(values.values()) for values in matrix)),
            "verified_local_entries": verified["local_entry_checks"],
            "verified_memberships": verified["regular_membership_checks"],
            "policy": choices,
        },
        "total_local_identity_entries": total_entries,
        "total_selected_target_memberships": total_memberships,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    report = check()
    report["checker_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    report["engine_sha256"] = sha256(Path(__file__).with_name("subspace_smith_channels.py").read_bytes()).hexdigest()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
