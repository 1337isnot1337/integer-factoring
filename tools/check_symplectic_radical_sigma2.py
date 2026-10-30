#!/usr/bin/env python3
"""Exact inclusion and coefficient algebras of the 24-label two-even-site code."""

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

from check_symplectic_radical_seed_algebra import (
    BlockAlgebra, plus, rank_two, representation_filtration, verify_normalization,
)


def pair_groups(code):
    result = defaultdict(list)
    for first in range(code["labels"]):
        for second in range(code["labels"]):
            result[code["syndromes"][first] ^ code["syndromes"][second]].append((first, second))
    return result


def verify_inclusion(parent, child):
    original, restricted = parent["code"], child["code"]
    selected = [i for i, syndrome in enumerate(original["syndromes"]) if syndrome]
    assert selected == list(range(2, 26))
    assert restricted["labels"] == 24
    assert restricted["ambient_occupations"] == original["ambient_occupations"]
    for key in ("input_vectors", "output_indices", "syndromes"):
        assert restricted[key] == [original[key][i] for i in selected]
    assert child["basis_monomials_p_q_r"] == parent["basis_monomials_p_q_r"]
    assert child["basis_scalar_polynomials_p_q_r"] == parent["basis_scalar_polynomials_p_q_r"]
    old_groups, new_groups = pair_groups(original), pair_groups(restricted)
    assert sorted(old_groups) == sorted(new_groups)
    inclusions = []
    for key in sorted(old_groups):
        positions = {pair: i for i, pair in enumerate(old_groups[key])}
        inclusions.append([positions[selected[x], selected[y]] for x, y in new_groups[key]])
    checks = 0
    for cut in ("row", "column"):
        for first, second in zip(parent[f"basis_{cut}_blocks"], child[f"basis_{cut}_blocks"]):
            for old_block, new_block, indices in zip(first, second, inclusions):
                for row_index, encoded in enumerate(new_block):
                    row = int(old_block[indices[row_index]], 16)
                    value = sum(((row >> old_index) & 1) << i for i, old_index in enumerate(indices))
                    assert value == int(encoded, 16)
                    checks += 1
    assert checks == 6912
    return {"parent_labels": 26, "labels": 24, "label_inclusion": selected,
            "cut_block_inclusions": inclusions, "exact_bit_rows_checked": checks,
            "both_cuts_and_all_six_normalized_generators_inherited": True}


def analyze(data, cut):
    raw = [[[int(row, 16) for row in block] for block in matrix]
           for matrix in data[f"basis_{cut}_blocks"]]
    generators = raw[:4] + [plus(raw[1], raw[4]), plus(raw[1], raw[5])]
    algebra = BlockAlgebra(generators)
    commutators = [algebra.table[i][j] ^ algebra.table[j][i]
                   for i in range(1, 6) for j in range(i + 1, 6)]
    radical = algebra.ideal(commutators + [
        algebra.table[i][i] ^ (1 << i) for i in range(1, 4)] + [16, 32])
    powers = algebra.ideal_powers(radical)
    characters = algebra.characters()
    assert len(characters) == algebra.dimension - len(radical)
    assert rank_two([functional for values, functional in characters]) == len(characters)
    expected = {"row": (22, [20, 9, 0], [384, 192]),
                "column": (9, [8, 5, 0], [576])}[cut]
    assert algebra.dimension == expected[0]
    assert list(map(len, powers)) == expected[1]
    filtration, layers, adapted = representation_filtration(algebra, powers, characters)
    assert filtration[0]["character_multiplicities"] == expected[2]
    return {
        "cut": cut, "module_dimension": 576, "algebra_dimension": algebra.dimension,
        "basis_words": algebra.words, "multiplication_table": algebra.table,
        "jacobson_radical_basis": radical,
        "radical_power_dimensions": list(map(len, powers)),
        "characters_on_E_A_B_C_D_F": [list(values) for values, functional in characters],
        "representation_radical_filtration": filtration,
        "representation_Loewy_layers": layers,
        "adapted_original_blocks": adapted,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent", type=Path,
                        default=Path("evidence/symplectic-radical-P4-L4-basis.json"))
    parser.add_argument("--basis", type=Path,
                        default=Path("evidence/symplectic-radical-sigma2-basis.json"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--inclusion-output", type=Path, required=True)
    args = parser.parse_args()
    parent = json.loads(args.parent.read_text())
    data = json.loads(args.basis.read_text())
    verify_normalization(data)
    inclusion = verify_inclusion(parent, data)
    hashes = {"parent_basis_sha256": hashlib.sha256(args.parent.read_bytes()).hexdigest(),
              "basis_sha256": hashlib.sha256(args.basis.read_bytes()).hexdigest(),
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    inclusion.update(hashes)
    args.inclusion_output.write_text(json.dumps(inclusion, indent=2, sort_keys=True) + "\n")
    result = {"labels": 24, "normalization": "R=E+tA+t^2B+t^3C+zD+ztF",
              **hashes, "cuts": [analyze(data, cut) for cut in ("row", "column")]}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    for cut in result["cuts"]:
        print(json.dumps({key: cut[key] for key in (
            "cut", "algebra_dimension", "radical_power_dimensions",
            "characters_on_E_A_B_C_D_F", "representation_radical_filtration")}), flush=True)


if __name__ == "__main__":
    main()
