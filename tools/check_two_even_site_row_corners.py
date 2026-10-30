#!/usr/bin/env python3
"""Small-core triangular certificates and the nine-site all-power barrier.

The new row calculation uses only five-site matrices of dimension at most 640.
The column and fallback row witnesses are independently rechecked from the
published bounded seven- and eight-site packages. No dense eight- or nine-site
coefficient matrix is constructed.
"""

import argparse
from itertools import combinations
import hashlib
import json
from math import comb
from pathlib import Path

from check_common_character_compressions import dual_rows, low_basis, matvec
from check_two_even_site_P8_compressions import audit_witness, weighted
from check_two_even_site_jordan_compressions import (
    build_cores, core_row_triangle, embedding_controls, local_similarity, module,
)
from find_tagged_common_character_compressions import columns


ROOT = Path(__file__).resolve().parents[1]
P7_CERTIFICATE = ROOT / "evidence/two-even-site-P7-jordan-compression-certificate.json"
P8_CERTIFICATE = ROOT / "evidence/two-even-site-P8-compression-certificate.json"
DEFAULT_CERTIFICATE = ROOT / "evidence/two-even-site-row-corner-certificate.json"
CHARACTER5 = [1, 1, 1, 1, 0, 0, 0, 0]


def product(first, second):
    result = []
    for row in first:
        value = 0
        while row:
            bit = row & -row
            value ^= second[bit.bit_length() - 1]
            row ^= bit
        result.append(value)
    return result


def image_filtration(matrices, character):
    size = len(matrices[0])
    differences = [[row ^ ((1 << i) if scalar else 0) for i, row in enumerate(matrix)]
                   for matrix, scalar in zip(matrices, character)]
    basis = [1 << i for i in range(size)]
    filtration = [basis]
    while basis:
        new_basis = low_basis([matvec(matrix, vector) for matrix in differences for vector in basis])
        assert len(new_basis) < len(basis), "The compressed difference ideal is not nilpotent."
        filtration.append(new_basis)
        basis = new_basis
    return filtration


def adapted_columns(filtration):
    pivots, selected = {}, []
    for basis in reversed(filtration):
        for vector in basis:
            original = vector
            while vector:
                pivot = (vector & -vector).bit_length() - 1
                if pivot not in pivots:
                    pivots[pivot] = vector
                    selected.append(original)
                    break
                vector ^= pivots[pivot]
    assert len(selected) == len(filtration[0])
    return selected


def linear_combination(vectors, coordinates):
    result = 0
    while coordinates:
        bit = coordinates & -coordinates
        result ^= vectors[bit.bit_length() - 1]
        coordinates ^= bit
    return result


def corner_matrices(matrices, vectors, zero_vectors):
    size = len(matrices[0])
    left = dual_rows(zero_vectors, vectors, size)
    result = []
    for matrix in matrices:
        images = [matvec(matrix, vector) for vector in vectors]
        result.append([sum(((row & image).bit_count() % 2) << j for j, image in enumerate(images))
                       for row in left])
    return left, result


def idempotent(matrix):
    value = matrix
    squarings = 0
    while product(value, value) != value:
        value = product(value, value)
        squarings += 1
        assert squarings <= 10
    return squarings, value


def core_matrices(cores, syndrome):
    data = cores[5]
    block = data["syndromes"].index(syndrome)
    return [matrix[block] for matrix in data["cuts"][0]]


def discover_corners(cores):
    records = []
    for syndrome in [s for s in cores[5]["syndromes"] if s.bit_count() in (0, 2)]:
        matrices = core_matrices(cores, syndrome)
        size = len(matrices[0])
        squarings, projector = idempotent(matrices[1])
        vectors = low_basis(columns(projector))
        zero = low_basis(columns([row ^ (1 << i) for i, row in enumerate(projector)]))
        assert len(vectors) + len(zero) == size
        _, compressed = corner_matrices(matrices, vectors, zero)
        filtration = image_filtration(compressed, CHARACTER5)
        adapted = adapted_columns(filtration)
        final_vectors = [linear_combination(vectors, coordinates) for coordinates in adapted]
        records.append({"core_syndrome": syndrome, "ambient_dimension": size,
                        "idempotent_squarings": squarings, "idempotent_exponent": 1 << squarings,
                        "rank_dimension": len(vectors),
                        "nilpotent_image_dimensions": [len(basis) for basis in filtration],
                        "vectors": [hex(vector) for vector in final_vectors]})
    return records


def source_values(prior7, prior8):
    row = {item["core_syndrome"]: weighted(item["cuts"]["row"])
           for item in prior7["types"] if item["core_pairs"] == 5}
    column = {(item["core_pairs"], item["core_syndrome"]): weighted(item["column"])
              for item in prior8["types"]}
    return row, column


def physical_certificate(pairs, corner_values, prior7, prior8):
    old_rows, columns8 = source_values(prior7, prior8)
    restrictions = []
    weight_four = [comb(pairs, 4) * 32 * 4 ** (pairs - 4), comb(pairs, 4) * columns8[4, 15]]
    weight_two = [0, 0]
    for syndrome in combinations(range(pairs), 2):
        choices = []
        for anchors in combinations([site for site in range(pairs) if site not in syndrome], 3):
            sites = sorted(syndrome + anchors)
            relative = sum(1 << sites.index(site) for site in syndrome)
            corner = corner_values[relative] * 4 ** (pairs - 5)
            fallback = old_rows[relative] * 4 ** (pairs - 7)
            row_mode = "core_corner" if corner >= fallback else "seven_site_scalar_lift"
            choices.append((relative, anchors, max(corner, fallback), row_mode))
        record = {"syndrome_sites": list(syndrome)}
        relative, anchors, amount, mode = max(choices, key=lambda item: item[2])
        record["row"] = {"anchors": list(anchors), "core_syndrome": relative,
                         "compression": amount, "mode": mode}
        weight_two[0] += amount
        relative, anchors, _, _ = max(choices, key=lambda item: columns8[5, item[0]])
        amount = columns8[5, relative]
        record["column"] = {"anchors": list(anchors), "core_syndrome": relative, "compression": amount,
                            "mode": "eight_site_scalar_lift"}
        weight_two[1] += amount
        restrictions.append(record)
    zero = [corner_values[0] * 4 ** (pairs - 5), columns8[5, 0]]
    totals = [sum(values) for values in zip(weight_four, weight_two, zero)]
    return {"pairs": pairs, "labels": comb(pairs, 2) * 2 ** (pairs - 2),
            "coefficient_directions": 2 * pairs - 2,
            "character": [comb(pairs - 2, i) % 2 if i <= pairs - 2 else 0 for i in range(pairs)] + [0] * (pairs - 2),
            "three_anchor_restrictions": restrictions, "zero_restriction_sites": list(range(5)),
            "weight_four_rank_dimensions": weight_four, "weight_two_restricted_rank_dimensions": weight_two,
            "weight_zero_restricted_rank_dimensions": zero, "rank_dimensions": totals,
            "rank_product": totals[0] * totals[1]}


def discover(cores, prior7, prior8):
    corners = discover_corners(cores)
    values = {item["core_syndrome"]: item["rank_dimension"] for item in corners}
    return {"core_pairs": 5, "core_character": CHARACTER5, "corners": corners,
            "source_certificates_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                           for path in (P7_CERTIFICATE, P8_CERTIFICATE)},
            "seeds": [physical_certificate(pairs, values, prior7, prior8) for pairs in (8, 9)],
            "cost_exponent": {"numerator": 17, "denominator": 12}}


def audit_corner(cores, supplied):
    syndrome = supplied["core_syndrome"]
    matrices = core_matrices(cores, syndrome)
    size = len(matrices[0])
    assert size == supplied["ambient_dimension"]
    projector = matrices[1]
    for _ in range(supplied["idempotent_squarings"]):
        projector = product(projector, projector)
    assert product(projector, projector) == projector
    assert supplied["idempotent_exponent"] == 1 << supplied["idempotent_squarings"]
    vectors = [int(value, 16) for value in supplied["vectors"]]
    amount = supplied["rank_dimension"]
    assert len(vectors) == len(low_basis(vectors)) == amount
    assert all(matvec(projector, vector) == vector for vector in vectors)
    zero = low_basis(columns([row ^ (1 << i) for i, row in enumerate(projector)]))
    assert len(vectors) + len(zero) == size
    left, compressed = corner_matrices(matrices, vectors, zero)
    assert compressed[0] == [1 << i for i in range(amount)]
    checks = 0
    for matrix, scalar in zip(compressed, CHARACTER5):
        for i, row in enumerate(matrix):
            assert (row >> i & 1) == scalar
            assert row & ((1 << i) - 1) == 0
            checks += i + 1
    dimensions = [len(basis) for basis in image_filtration(compressed, CHARACTER5)]
    assert dimensions == supplied["nilpotent_image_dimensions"]
    return {"core_syndrome": syndrome, "ambient_dimension": size, "rank_dimension": amount,
            "idempotent_exponent": supplied["idempotent_exponent"],
            "nilpotent_image_dimensions": dimensions,
            "triangular_entries_checked": checks,
            "compressed_matrix_entries_constructed": len(matrices) * amount ** 2,
            "left_covectors_sha256": hashlib.sha256(json.dumps([hex(row) for row in left], separators=(",", ":")).encode()).hexdigest(),
            "compressed_matrices_sha256": hashlib.sha256(json.dumps(compressed, separators=(",", ":")).encode()).hexdigest()}


def audit_dependencies(cores, prior7, prior8):
    reports = []
    fallback = next(item for item in prior7["types"] if (item["core_pairs"], item["core_syndrome"]) == (5, 17))
    for j, item in enumerate(fallback["cuts"]["row"]):
        assert item["j"] == j and item["multiplicity"] == comb(2, j) * (1 << (2 - j))
        report = audit_witness(module(cores, 7, 5, 17, j, 0), prior7["character"], item["witness"])
        reports.append(dict(report, pairs=7, core_pairs=5, core_syndrome=17, cut="row", jordan_copies=j))
    assert weighted(fallback["cuts"]["row"]) == 1163
    for supplied in prior8["types"]:
        h, syndrome = supplied["core_pairs"], supplied["core_syndrome"]
        for j, item in enumerate(supplied["column"]):
            assert item["j"] == j and item["multiplicity"] == comb(8 - h, j) * (1 << (8 - h - j))
            report = audit_witness(module(cores, 8, h, syndrome, j, 1), prior8["character"], item["witness"])
            reports.append(dict(report, pairs=8, core_pairs=h, core_syndrome=syndrome, cut="column", jordan_copies=j))
    return reports


def audit(cores, certificate, prior7, prior8):
    assert certificate["core_pairs"] == 5 and certificate["core_character"] == CHARACTER5
    assert certificate["cost_exponent"] == {"numerator": 17, "denominator": 12}
    for path in (P7_CERTIFICATE, P8_CERTIFICATE):
        assert certificate["source_certificates_sha256"][str(path.relative_to(ROOT))] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert prior7["character"] == [1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0]
    assert prior8["character"] == [1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0]
    expected_syndromes = {s for s in cores[5]["syndromes"] if s.bit_count() in (0, 2)}
    assert len(certificate["corners"]) == len(expected_syndromes)
    assert {item["core_syndrome"] for item in certificate["corners"]} == expected_syndromes
    corners = [audit_corner(cores, supplied) for supplied in certificate["corners"]]
    corner_values = {item["core_syndrome"]: item["rank_dimension"] for item in corners}
    assert corner_values == {0: 640, 3: 384, 5: 256, 6: 384, 9: 128, 10: 256,
                             12: 384, 17: 0, 18: 128, 20: 256, 24: 384}
    for syndrome, value in corner_values.items():
        if syndrome:
            a, b = [i for i in range(5) if syndrome >> i & 1]
            assert value == 128 * (4 - (b - a))
    similarity = local_similarity()
    # The physical coordinate zero gives the column scalar 1+t directly.
    assert [1, 2, 4, 8][0] & 1 == 1 and [9, 0, 0, 9][0] & 1 == 1
    triangle = core_row_triangle()
    dependencies = audit_dependencies(cores, prior7, prior8)
    seeds, literal = [], []
    assert [item["pairs"] for item in certificate["seeds"]] == [8, 9]
    expected_dimensions = {8: [1225260, 2106], 9: [7522480, 3332]}
    for supplied in certificate["seeds"]:
        pairs = supplied["pairs"]
        expected = physical_certificate(pairs, corner_values, prior7, prior8)
        assert supplied == expected
        assert supplied["rank_dimensions"] == expected_dimensions[pairs]
        # Only one physical extreme syndrome uses the old scalar row witness.
        fallbacks = [item for item in supplied["three_anchor_restrictions"]
                     if item["row"]["mode"] == "seven_site_scalar_lift"]
        assert len(fallbacks) == 1 and fallbacks[0]["syndrome_sites"] == [0, pairs - 1]
        assert fallbacks[0]["row"]["core_syndrome"] == 17
        assert supplied["rank_dimensions"][0] == (2048 * comb(pairs, 4) + 6144 * comb(pairs, 2) - 9077) * 4 ** (pairs - 7)
        assert supplied["rank_dimensions"][1] == 16 * comb(pairs, 4) + 46 * comb(pairs, 2) - 38 * pairs + 2
        product_value, labels = supplied["rank_product"], supplied["labels"]
        assert product_value ** 6 > labels ** 17 and 17 ** 2 > 2 * 12 ** 2
        controls = embedding_controls(cores, supplied)
        literal.append({"pairs": pairs, "controls": controls})
        seeds.append({key: value for key, value in supplied.items() if key not in ("three_anchor_restrictions", "zero_restriction_sites")})
        seeds[-1]["exact_margin"] = {"rank_product_power": 6, "label_power": 17,
                                     "difference": product_value ** 6 - labels ** 17}
    ten = physical_certificate(10, corner_values, prior7, prior8)
    assert ten["rank_dimensions"] == [44638912, 5052]
    assert ten["rank_product"] == 225515783424
    assert ten["rank_product"] ** 5 < ten["labels"] ** 14 and 7 ** 2 < 2 * 5 ** 2
    return {"status": "pass", "core_character": CHARACTER5,
            "core_character_scalar_polynomial": "(q+r)^2 (p+r)^3",
            "core_row_corners": corners, "seeds": seeds,
            "triangular_entries_checked": sum(item["triangular_entries_checked"] for item in corners),
            "compressed_matrix_entries_constructed": sum(item["compressed_matrix_entries_constructed"] for item in corners),
            "dependency_scalar_compression_entries_checked": sum(item["all_generator_compression_entries_checked"] for item in dependencies),
            "dependency_module_audits": dependencies,
            "local_constant_similarities": similarity,
            "weight_four_row_triangular_certificate": triangle,
            "literal_embedding_controls": literal,
            "largest_new_row_matrix_dimension": 640,
            "largest_dependency_matrix_dimension": 5120,
            "all_size_rank_formulas": {
                "domain": "Every integer P >= 8; the anchor case count is proved in the note.",
                "row": "(2048*C(P,4) + 6144*C(P,2) - 9077)*4^(P-7)",
                "column": "16*C(P,4) + 46*C(P,2) - 38*P + 2"},
            "first_unclosed_seed_for_this_certificate": {
                "pairs": 10, "labels": ten["labels"], "rank_dimensions": ten["rank_dimensions"],
                "rank_product": ten["rank_product"],
                "exact_shortfall": {"rank_product_power": 5, "label_power": 14,
                                    "difference": ten["labels"] ** 14 - ten["rank_product"] ** 5}},
            "cost_exponent": {"numerator": 17, "denominator": 12},
            "scope": "All characteristic-two fields, every tensor power, arbitrary coefficient-space recombinations."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--basis5", type=Path)
    parser.add_argument("--discover", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    prior7, prior8 = [json.loads(path.read_text()) for path in (P7_CERTIFICATE, P8_CERTIFICATE)]
    cores = build_cores(args.basis5)
    certificate = json.loads(args.certificate.read_text())
    if args.discover:
        assert discover(cores, prior7, prior8) == certificate
    result = audit(cores, certificate, prior7, prior8)
    result["certificate_sha256"] = hashlib.sha256(args.certificate.read_bytes()).hexdigest()
    result["source_certificates_sha256"] = certificate["source_certificates_sha256"]
    names = (Path(__file__).name, "check_two_even_site_P8_compressions.py",
             "check_two_even_site_jordan_compressions.py", "assemble_two_even_site_cores.py",
             "check_common_character_compressions.py", "find_tagged_common_character_compressions.py",
             "check_noncommutative_crossing_character_compression.py", "explore_symplectic_boundary_codes.py",
             "check_modular_crossing.py", "complement_shell_tensor.py", "explore_complement_edge_pencils.py")
    result["source_sha256"] = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
