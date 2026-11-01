#!/usr/bin/env python3
"""Exact block/precision hierarchy controls and companion improvements."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from random import Random

import check_cancellation_module_automata as arithmetic
from check_exact_module_orbits import companion_two
from block_subspace_channels import (
    block_local_crossing, conditional_right_product, regular_right_inverse,
    transition_product,
)
from channel_pressure_separation import minimum_transfer, verify_lower_potential
from matroid_smith_channels import binary_plane_lattices
from subspace_smith_channels import SubspaceSmithLattices


ROOT = Path(__file__).resolve().parents[1]


def apply(matrix, vector):
    return tuple(sum(value*vector[target] for target, value in row.items()) for row in matrix)


def check(maximum_block=2):
    point, local = Fraction(5), companion_two()
    records = []
    identity_entries = memberships = endpoint_entries = rank_controls = 0
    domination_controls = precision_controls = right_inverse_controls = 0
    random = Random(20261031)
    lower_evidence = json.loads((ROOT / "evidence" / "channel-pressure-separation.json").read_text())
    one_lower_product = Fraction(lower_evidence["companion_product_lower"])
    one_graphs = []
    for index, (cut, infinity) in enumerate(((0, True), (1, False))):
        graph = SubspaceSmithLattices(local, 1, seed_states=binary_plane_lattices(2),
                                     cut=cut, infinity=infinity)
        certificate = lower_evidence["companion_records"][index]["search_result"]["certificate"]
        potential = tuple(Fraction(value) for value in certificate["potential"])
        oracle = lambda h: minimum_transfer(graph, h, point)
        assert verify_lower_potential(oracle, potential, Fraction(certificate["rate"])) is not None
        one_graphs.append(graph)

        # Adding finer states is safe when enlarged targets and lower charges
        # are allowed. Extend any old positive potential by h_I on new states.
        coarse = SubspaceSmithLattices(local, 1, seed_states=binary_plane_lattices(1),
                                      precision=1, cut=cut, infinity=infinity)
        identity_index = 0
        for _ in range(5):
            coarse_potential = tuple(Fraction(random.randrange(1, 101)) for _ in coarse.states)
            old_matrix, _ = coarse.choose(coarse_potential, point)
            old_image = apply(old_matrix, coarse_potential)
            bound = max(value/height for value, height in zip(old_image, coarse_potential))
            old_positions = {state: i for i, state in enumerate(coarse.states)}
            finer_potential = tuple(coarse_potential[old_positions[state]] if state in old_positions
                                    else coarse_potential[identity_index] for state in graph.states)
            _, new_image = minimum_transfer(graph, finer_potential, point)
            assert all(value <= bound*height for value, height in zip(new_image, finer_potential))
            precision_controls += len(graph.states)

    # A rectangular varying-rank instance of the conditional right inverse.
    one, zero, tau = arithmetic.ONE, arithmetic.ZERO, arithmetic.rat(2)
    first = ((one, tau, zero), (zero, one, one))
    seconds = (((one, tau, one),), ((one, zero, tau), (zero, one, one)))
    for matrix in (first, *seconds, conditional_right_product(first, seconds)):
        regular_right_inverse(matrix)
        right_inverse_controls += 1

    fixtures = {
        2: {
            "precision": 2,
            "potentials": (
                (1000000, 512311, 1000000, 1000000, 512311, 1000000, 512311, 1000000, 1000000, 1000000),
                (1000000, 906557, 645834, 617226, 895729, 637956, 895729, 637956, 540672, 605830),
            ),
            "rates": (Fraction(13904, 125), Fraction(707, 500)),
            "per_site_upper": Fraction(6271, 500),
        },
        3: {
            "precision": 1,
            "potentials": ((1000000, 512311, 1000000, 1000000),
                           (1000000, 895993, 628166, 625925)),
            "rates": (Fraction(1140), Fraction(5, 3)),
            "per_site_upper": Fraction(6193, 500),
        },
    }
    for length in range(2, maximum_block+1):
        fixture = fixtures[length]
        precision = fixture["precision"]
        blocked = block_local_crossing(local, 2, 1, length)
        for (x, y, u, v), matrix in blocked.items():
            assert matrix[0][0].bit_count() % 2 == int(x == v and y == u)
            assert matrix[1][0].bit_count() % 2 == 0
            endpoint_entries += 2
        block_graphs, transfers, side_records = [], [], []
        for index, (cut, infinity) in enumerate(((0, True), (1, False))):
            graph = SubspaceSmithLattices(blocked, length,
                                         seed_states=binary_plane_lattices(precision),
                                         precision=precision, alphabet=2**length,
                                         cut=cut, infinity=infinity)
            potential, rate = fixture["potentials"][index], fixture["rates"][index]
            matrix, policy = graph.choose(potential, point)
            image = apply(matrix, potential)
            assert all(value <= rate*height for value, height in zip(image, potential))
            certificate = graph.verify(policy)
            identity_entries += certificate["local_entry_checks"]
            memberships += certificate["regular_membership_checks"]
            terminal = graph.terminal_weights(point)
            prefactor = potential[graph.start] * max(value/potential[i] for i, value in enumerate(terminal))
            side_records.append({
                "graph": graph.summary(), "potential": potential, "rate": str(rate),
                "maximum_row_ratio": str(max(value/height for value, height in zip(image, potential))),
                "prefactor": str(prefactor), "policy": policy,
                "transfer": [[str(row.get(j, Fraction(0))) for j in range(len(graph.states))] for row in matrix],
                "verified_identity_entries": certificate["local_entry_checks"],
            })
            transfers.append(matrix)
            block_graphs.append(graph)

            if length == 2:
                one_graph = one_graphs[index]
                for _ in range(5):
                    height = tuple(Fraction(random.randrange(1, 101)) for _ in one_graph.states)
                    one_matrix, _ = one_graph.choose(height, point)
                    inherited = apply(transition_product(one_matrix, one_matrix), height)
                    _, improved = minimum_transfer(graph, height, point)
                    assert all(value <= old for value, old in zip(improved, inherited))
                    domination_controls += len(graph.states)
                # Every selected local R has a regular right inverse. Their
                # path product remains surjective, as the proof requires.
                _, one_policy = one_graph.choose((Fraction(1),) * len(one_graph.states), point)
                one_certificate = one_graph.verify(one_policy)
                for state in one_certificate["states"]:
                    factors = state["factors"]
                    first_right = tuple(factor["right"] for factor in factors)
                    next_rights = [tuple(next_factor["right"] for next_factor in
                                        one_certificate["states"][factor["target"]]["factors"])
                                   for factor in factors]
                    regular_right_inverse(first_right)
                    regular_right_inverse(conditional_right_product(first_right, next_rights))
                    right_inverse_controls += 2

        product_bound = fixture["rates"][0] * fixture["rates"][1]
        assert product_bound < fixture["per_site_upper"]**length
        assert product_bound < one_lower_product**length
        finite = []
        for blocks in range(1, (2 if length == 2 else 1)+1):
            original_length = length*blocks
            ranks, entries = arithmetic.coefficient_ranks(local, 1, original_length)
            endpoint_entries += entries
            words = []
            for graph, matrix, side_record in zip(block_graphs, transfers, side_records):
                vector = graph.terminal_weights(point)
                for _ in range(blocks):
                    vector = apply(matrix, vector)
                value = vector[graph.start]
                assert value <= Fraction(side_record["prefactor"]) * Fraction(side_record["rate"])**blocks
                words.append(value)
            for degree in range(original_length+1):
                assert ranks[0][degree] <= (original_length-degree+1)*point**(-degree)*words[0]
                assert ranks[1][degree] <= (degree+1)*point**degree*words[1]
                rank_controls += 2
            finite.append({"blocks": blocks, "original_sites": original_length,
                           "row_ranks": ranks[0], "column_ranks": ranks[1],
                           "word_sums": list(map(str, words))})
        records.append({"block_length": length, "precision": precision,
                        "alphabet": 2**length, "point": str(point),
                        "product_bound": str(product_bound),
                        "per_original_site_upper": str(fixture["per_site_upper"]),
                        "sides": side_records, "finite_rank_controls": finite})
    return {
        "status": "pass", "maximum_block": maximum_block,
        "scope": "Exact improved companion certificates and complete-oracle hierarchy controls; no factoring-exponent improvement",
        "one_site_global_product_lower_at_s5": str(one_lower_product),
        "records": records,
        "total_local_identity_entries": identity_entries,
        "total_selected_memberships": memberships,
        "literal_endpoint_entries": endpoint_entries,
        "literal_coefficient_rank_inequalities": rank_controls,
        "block_domination_row_controls": domination_controls,
        "precision_domination_row_controls": precision_controls,
        "regular_right_inverse_controls": right_inverse_controls,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-block", type=int, choices=(2, 3), default=2)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    report = check(arguments.maximum_block)
    report["checker_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    report["engine_sha256"] = sha256(Path(__file__).with_name("block_subspace_channels.py").read_bytes()).hexdigest()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
