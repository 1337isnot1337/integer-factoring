#!/usr/bin/env python3
"""Replay exact lower potentials for all regular channel-basis policies.

The default checker needs no linear-programming package. --search uses
the floating cutting-plane driver to find new rational witnesses and
then applies the same exact verification gate.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from random import Random

from channel_pressure_separation import (
    find_lower_potential, minimum_transfer, verify_lower_potential,
)
from check_exact_module_orbits import companion_two
from check_subspace_smith_channels import exhaustive_basis_cost, span
from matroid_smith_channels import binary_plane_lattices
from subspace_smith_channels import SubspaceSmithLattices


ROOT = Path(__file__).resolve().parents[1]


def reducible_oracle(potential):
    return ({0: Fraction(1)}, {1: Fraction(2)}), (potential[0], 2*potential[1])


def check(search=False, evidence_path=None):
    evidence_path = evidence_path or ROOT / "evidence" / "channel-pressure-separation.json"
    if not search:
        saved = json.loads(evidence_path.read_text())
    states, point = binary_plane_lattices(2), Fraction(5)
    random = Random(20261031)
    controls = entries = 0
    records = []
    lower_rates = (Fraction(53197, 5000), Fraction(1204989, 1000000))
    upper_rates = (Fraction(266, 25), Fraction(241, 200))
    upper_potentials = (
        (10000, 5320, 10000, 10000, 5320, 10000, 5320, 10000, 10000, 10000),
        (100000, 94634, 63719, 62303, 94634, 62828, 94634, 63719, 54764, 61230),
    )
    for index, (cut, infinity) in enumerate(((0, True), (1, False))):
        graph = SubspaceSmithLattices(companion_two(), 1, seed_states=states,
                                     cut=cut, infinity=infinity)
        oracle = lambda potential: minimum_transfer(graph, potential, point)
        families = []
        for actions in graph.actions:
            families.append({(action["target"], action["order"]):
                             span([residue for residue, _ in action["basis"]])-{0}
                             for action in actions})
        for trial in range(3):
            potential = tuple(Fraction(random.randrange(1, 101)) if (trial == 0 or i % 3 == trial)
                              else Fraction(0) for i in range(len(states)))
            matrix, image = oracle(potential)
            for source, family in enumerate(families):
                assert image[source] == exhaustive_basis_cost(
                    family, len(graph.local_factors[source]), potential, point,
                    graph.degree if infinity else 0,
                )
                controls += 1
            if trial == 0:
                direct, _ = graph.choose(potential, point)
                assert image == tuple(sum(value*potential[target] for target, value in row.items())
                                      for row in direct)
        if search:
            found = find_lower_potential(oracle, len(states), lower_rates[index])
            assert found["status"] == "certified"
        else:
            found = saved["companion_records"][index]["search_result"]
        certificate = found["certificate"]
        potential = tuple(Fraction(value) for value in certificate["potential"])
        verified = verify_lower_potential(oracle, potential, lower_rates[index])
        assert verified is not None
        assert json.loads(json.dumps(verified)) == json.loads(json.dumps(certificate))
        upper_matrix, upper_choices = graph.choose(upper_potentials[index], point)
        upper_image = tuple(sum(value*upper_potentials[index][target] for target, value in row.items())
                            for row in upper_matrix)
        assert all(value <= upper_rates[index]*height
                   for value, height in zip(upper_image, upper_potentials[index]))
        entry_check = graph.verify(upper_choices)
        entries += entry_check["local_entry_checks"]
        records.append({
            "cut": cut, "infinity": infinity, "point": str(point),
            "lower_rate": str(lower_rates[index]), "upper_rate": str(upper_rates[index]),
            "upper_potential": upper_potentials[index],
            "search_result": found,
            "lower_witness_start_coordinate": str(potential[graph.start]),
        })
    reducible = verify_lower_potential(reducible_oracle, (0, 1), 2)
    assert reducible is not None and reducible["support"] == [1]
    assert verify_lower_potential(reducible_oracle, (1, 1), 2) is None
    assert verify_lower_potential(reducible_oracle, (0, 1), Fraction(201, 100)) is None
    if search:
        found = find_lower_potential(reducible_oracle, 2, 2)
        assert found["status"] == "certified"
        assert found["certificate"]["potential"] == ("0", "1")
        failed = find_lower_potential(reducible_oracle, 2, Fraction(201, 100))
        assert failed["status"] == "no_exact_certificate"
    lower_product = lower_rates[0] * lower_rates[1]
    upper_product = upper_rates[0] * upper_rates[1]
    assert lower_product == Fraction(64101799833, 5000000000)
    assert upper_product == Fraction(32053, 2500)
    return {
        "status": "pass",
        "scope": "Global all-state policy-rate certificates at the fixed scale s=5; no boundary-word, actual-rank, or factoring lower bound",
        "nonnegative_potential_all_basis_controls": controls,
        "exact_upper_policy_identity_entries": entries,
        "companion_records": records,
        "companion_product_lower": str(lower_product),
        "companion_product_upper": str(upper_product),
        "reducible_zero_coordinate_control": {
            "matrix": [[1, 0], [0, 2]], "certificate": reducible,
            "initial_state": 0, "initial_state_word_rate": "1",
            "global_spectral_rate": "2",
            "scope": "A global spectral lower certificate may have zero at the initial state and need not bound that state's word growth.",
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--evidence", type=Path)
    arguments = parser.parse_args()
    report = check(arguments.search, arguments.evidence)
    report["checker_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    report["engine_sha256"] = sha256(Path(__file__).with_name("channel_pressure_separation.py").read_bytes()).hexdigest()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
