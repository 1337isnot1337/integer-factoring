#!/usr/bin/env python3
"""Audit a bounded-module, all-power certificate for the eight-site seed.

The row witness tensors verified seven-site restrictions with one full local
factor.  The column witness uses only four- and five-site cores with Jordan
factors.  The largest matrix has dimension 5120; no full eight-site matrix is
constructed.  --discover reproduces the fixed coordinate witnesses.
"""

import argparse
from itertools import combinations
import hashlib
import json
from math import comb
from pathlib import Path

from check_common_character_compressions import audit_block
from check_two_even_site_jordan_compressions import (
    build_cores, core_row_triangle, embedding_controls, local_similarity, module,
)
from find_tagged_common_character_compressions import search_block


ROOT = Path(__file__).resolve().parents[1]
P7_CERTIFICATE = ROOT / "evidence/two-even-site-P7-jordan-compression-certificate.json"
DEFAULT_CERTIFICATE = ROOT / "evidence/two-even-site-P8-compression-certificate.json"
PAIRS, LABELS = 8, 1792
CHARACTER = [comb(6, i) % 2 if i <= 6 else 0 for i in range(8)] + [0] * 6


def weighted(items):
    return sum(item["multiplicity"] * item["witness"]["compression"] for item in items)


def relative_syndrome(syndrome_sites, anchors, pairs=PAIRS):
    sites = sorted(syndrome_sites + anchors)
    assert len(sites) == len(set(sites)) == 5
    assert all(0 <= site < pairs for site in sites)
    return sum(1 << sites.index(site) for site in syndrome_sites)


def discover(cores, prior):
    old_rows = {item["core_syndrome"]: weighted(item["cuts"]["row"])
                for item in prior["types"] if item["core_pairs"] == 5}
    records, columns = [], {}
    types = [(4, 15)] + [(5, syndrome) for syndrome in cores[5]["syndromes"]
                         if syndrome.bit_count() in (0, 2)]
    for h, syndrome in types:
        free_copies = PAIRS - h
        items = []
        for j in range(free_copies + 1):
            size = len(cores[h]["groups"][syndrome]) * (1 << j)
            if j < 6 - (h - 1):
                witness = {"size": size, "compression": 0, "annihilated_span": 0,
                           "coordinate_indices": []}
            else:
                matrices = module(cores, PAIRS, h, syndrome, j, 1)
                witness = search_block(matrices, CHARACTER, seed=0)
            items.append({"j": j, "multiplicity": comb(free_copies, j) * (1 << (free_copies - j)),
                          "witness": witness})
        columns[h, syndrome] = weighted(items)
        records.append({"core_pairs": h, "core_syndrome": syndrome, "column": items,
                        "weighted_column_compression": columns[h, syndrome]})
    restrictions = []
    totals = [comb(PAIRS, 4) * 32 * 4 ** (PAIRS - 4), comb(PAIRS, 4) * columns[4, 15]]
    for syndrome in combinations(range(PAIRS), 2):
        choices = []
        for anchors in combinations([site for site in range(PAIRS) if site not in syndrome], 3):
            relative = relative_syndrome(list(syndrome), list(anchors))
            choices.append((relative, anchors))
        record = {"syndrome_sites": list(syndrome)}
        for cut, name in enumerate(("row", "column")):
            def value(choice):
                return 4 * old_rows[choice[0]] if cut == 0 else columns[5, choice[0]]
            relative, anchors = max(choices, key=value)
            amount = value((relative, anchors))
            totals[cut] += amount
            record[name] = {"anchors": list(anchors), "core_syndrome": relative, "compression": amount}
        restrictions.append(record)
    totals[0] += 4 * old_rows[0]
    totals[1] += columns[5, 0]
    return {"pairs": PAIRS, "labels": LABELS, "coefficient_directions": 14,
            "character": CHARACTER, "discovery": {"algorithm": "tagged coordinate-subspace search", "seed": 0},
            "row_source_certificate_sha256": hashlib.sha256(P7_CERTIFICATE.read_bytes()).hexdigest(),
            "types": records, "three_anchor_restrictions": restrictions,
            "zero_restriction_sites": list(range(5)), "rank_dimensions": totals,
            "rank_product": totals[0] * totals[1],
            "cost_exponent": {"numerator": 99, "denominator": 70}}


def audit_witness(matrices, character, witness):
    assert witness["size"] == len(matrices[0])
    indices = witness["coordinate_indices"]
    assert len(indices) == len(set(indices)) == witness["compression"]
    assert all(0 <= index < witness["size"] for index in indices)
    if not indices:
        assert witness["annihilated_span"] == 0
        return {"compression_dimension": 0, "all_generator_compression_entries_checked": 0}
    expanded = dict(witness, vectors=[hex(1 << index) for index in indices])
    return audit_block(matrices, character, expanded)


def audit(cores, certificate, prior):
    assert certificate["pairs"] == PAIRS and certificate["labels"] == LABELS
    assert certificate["coefficient_directions"] == 2 * PAIRS - 2
    assert certificate["character"] == CHARACTER
    assert certificate["row_source_certificate_sha256"] == hashlib.sha256(P7_CERTIFICATE.read_bytes()).hexdigest()
    assert prior["pairs"] == 7 and prior["labels"] == 672
    assert prior["character"] == [comb(5, i) % 2 if i <= 5 else 0 for i in range(7)] + [0] * 5
    similarity = local_similarity()
    triangle = core_row_triangle()
    literal = embedding_controls(cores, certificate)
    old_rows, records = {}, []
    for supplied in prior["types"]:
        if supplied["core_pairs"] != 5:
            continue
        syndrome = supplied["core_syndrome"]
        items = supplied["cuts"]["row"]
        assert len(items) == 3
        reports = []
        for j, item in enumerate(items):
            assert item["j"] == j and item["multiplicity"] == comb(2, j) * (1 << (2 - j))
            matrices = module(cores, 7, 5, syndrome, j, 0)
            report = audit_witness(matrices, prior["character"], item["witness"])
            reports.append(dict(report, jordan_copies=j, multiplicity=item["multiplicity"]))
        old_rows[syndrome] = weighted(items)
        records.append({"core_pairs": 5, "core_syndrome": syndrome, "cut": "row",
                        "source_pairs": 7, "local_triangular_lift_dimension": 4,
                        "weighted_rank_dimension": 4 * old_rows[syndrome], "modules": reports})
    columns = {}
    expected_types = {(4, 15)} | {(5, syndrome) for syndrome in cores[5]["syndromes"]
                                 if syndrome.bit_count() in (0, 2)}
    assert len(certificate["types"]) == len(expected_types)
    assert {(item["core_pairs"], item["core_syndrome"]) for item in certificate["types"]} == expected_types
    for supplied in certificate["types"]:
        h, syndrome = supplied["core_pairs"], supplied["core_syndrome"]
        free_copies = PAIRS - h
        core_size = len(cores[h]["groups"][syndrome])
        assert len(supplied["column"]) == free_copies + 1
        reports, dimension = [], 0
        for j, item in enumerate(supplied["column"]):
            multiplicity = comb(free_copies, j) * (1 << (free_copies - j))
            assert item["j"] == j and item["multiplicity"] == multiplicity
            assert item["witness"]["size"] == core_size * (1 << j)
            matrices = module(cores, PAIRS, h, syndrome, j, 1)
            report = audit_witness(matrices, CHARACTER, item["witness"])
            reports.append(dict(report, jordan_copies=j, multiplicity=multiplicity))
            dimension += multiplicity * item["witness"]["size"]
        assert dimension == core_size * 4 ** free_copies
        total = weighted(supplied["column"])
        assert total == supplied["weighted_column_compression"]
        columns[h, syndrome] = total
        records.append({"core_pairs": h, "core_syndrome": syndrome, "cut": "column",
                        "source_pairs": PAIRS, "weighted_rank_dimension": total, "modules": reports})
    weight_four = [comb(PAIRS, 4) * triangle["core_character_dimension"] * 4 ** (PAIRS - 4),
                   comb(PAIRS, 4) * columns[4, 15]]
    weight_two = [0, 0]
    restrictions = certificate["three_anchor_restrictions"]
    assert len(restrictions) == comb(PAIRS, 2)
    assert {tuple(item["syndrome_sites"]) for item in restrictions} == set(combinations(range(PAIRS), 2))
    for record in restrictions:
        for cut, name in enumerate(("row", "column")):
            supplied = record[name]
            relative = relative_syndrome(record["syndrome_sites"], supplied["anchors"])
            assert relative == supplied["core_syndrome"]
            value = 4 * old_rows[relative] if cut == 0 else columns[5, relative]
            assert value == supplied["compression"]
            weight_two[cut] += value
    assert certificate["zero_restriction_sites"] == list(range(5))
    zero = [4 * old_rows[0], columns[5, 0]]
    totals = [sum(values) for values in zip(weight_four, weight_two, zero)]
    assert weight_four == [573440, 1120]
    assert weight_two == [182732, 980] and zero == [12888, 6]
    assert totals == certificate["rank_dimensions"] == [769060, 2106]
    product = totals[0] * totals[1]
    assert product == certificate["rank_product"] == 1619640360
    assert certificate["cost_exponent"] == {"numerator": 99, "denominator": 70}
    assert product ** 35 > LABELS ** 99 and 99 ** 2 > 2 * 70 ** 2
    checks = sum(item["all_generator_compression_entries_checked"]
                 for record in records for item in record["modules"])
    return {"status": "pass", "labels": LABELS, "pairs": PAIRS, "coefficient_directions": 14,
            "character": CHARACTER, "character_scalar_polynomial": "(q+r)^2 (p+r)^6",
            "local_constant_similarities": similarity,
            "literal_eight_site_embedding_controls": literal,
            "weight_four_row_triangular_certificate": triangle,
            "modules": records, "weight_four_rank_dimensions": weight_four,
            "weight_two_restricted_rank_dimensions": weight_two, "weight_zero_restricted_rank_dimensions": zero,
            "final_rank_certificate_dimensions": totals, "rank_product": product,
            "all_generator_compression_entries_checked": checks,
            "largest_matrix_dimension": max(item["witness"]["size"]
                for supplied in certificate["types"] for item in supplied["column"]),
            "exact_margin": {"rank_product_power": 35, "label_power": 99,
                             "difference": product ** 35 - LABELS ** 99},
            "cost_exponent": {"numerator": 99, "denominator": 70},
            "scope": "All characteristic-two fields, every tensor power, arbitrary coefficient-space recombinations."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--basis5", type=Path)
    parser.add_argument("--discover", action="store_true", help="Reproduce and compare every stored coordinate witness.")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    prior = json.loads(P7_CERTIFICATE.read_text())
    cores = build_cores(args.basis5)
    certificate = json.loads(args.certificate.read_text())
    if args.discover:
        assert discover(cores, prior) == certificate
    result = audit(cores, certificate, prior)
    result["certificate_sha256"] = hashlib.sha256(args.certificate.read_bytes()).hexdigest()
    result["row_source_certificate_sha256"] = hashlib.sha256(P7_CERTIFICATE.read_bytes()).hexdigest()
    names = (Path(__file__).name, "check_two_even_site_jordan_compressions.py",
             "assemble_two_even_site_cores.py", "check_common_character_compressions.py",
             "find_tagged_common_character_compressions.py", "check_noncommutative_crossing_character_compression.py",
             "explore_symplectic_boundary_codes.py", "check_modular_crossing.py",
             "complement_shell_tensor.py", "explore_complement_edge_pencils.py")
    result["source_sha256"] = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                               for name in names}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
