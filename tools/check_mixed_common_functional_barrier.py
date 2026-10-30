#!/usr/bin/env python3
"""Check source dimensions and exact margins for the mixed-space corollary.

This small checker reads the already audited local certificates. The proof of
mixed tensor compression and the cyclic product-grid step is in the note; this
script does not replace the local matrix reconstruction checkers.
"""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    hashes = {}

    def read(name):
        path = ROOT / "evidence" / name
        hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        return json.loads(path.read_text())

    records = []

    def add(name, labels, row, column, mechanism="common triangular/scalar diagonal", scope=None):
        product = row * column
        assert product ** 6 > labels ** 17
        records.append({"family": name, "labels": labels, "row_rank_base": row,
                        "column_rank_base": column, "rank_product": product,
                        "column_mechanism": mechanism, "scope": scope,
                        "exact_margin": {"rank_product_power": 6, "label_power": 17,
                                         "difference": product ** 6 - labels ** 17}})

    data = read("noncommutative-crossing-character-compression.json")
    assert data["column_compression"]["dimension"] == 43
    for labels, row in ((24, 384), (26, 436)):
        assert data["row_triangular_certificates"][str(labels)]["multiplicity"] == row
        add("symplectic-" + str(labels), labels, row, 43)
    for filename, labels, row, column in (
        ("symplectic-two-even-P5-compression.json", 80, 2072, 286),
        ("two-even-site-P6-compression.json", 240, 15688, 718),
    ):
        data = read(filename)
        assert data["labels"] == labels
        assert data["cuts"]["row"]["compression_dimension"] == row
        assert data["cuts"]["column"]["compression_dimension"] == column
        add("two-even-" + str(labels), labels, row, column)
    data = read("two-even-site-P7-jordan-compression.json")
    assert data["labels"] == 672 and data["final_rank_certificate_dimensions"] == [108315, 1150]
    add("two-even-672", 672, 108315, 1150)
    data = read("two-even-site-row-corner.json")
    for labels, row, column in ((1792, 1225260, 2106), (4608, 7522480, 3332)):
        seed = next(item for item in data["seeds"] if item["labels"] == labels)
        assert seed["rank_dimensions"] == [row, column]
        add("two-even-" + str(labels), labels, row, column)
    data = read("cubic-cyclic-fourier-bound.json")
    assert data["row_scalar_compression_dimension"] == 14
    assert data["column_simultaneously_upper_triangular_restriction_dimension"] == 28
    assert data["tensor_power_grid_size_per_coordinate"] == 7
    assert data["tensor_power_polynomial_degree_per_coordinate"] == 1
    assert data["tensor_power_minimum_grid_support_base"] == 6
    assert data["tensor_power_column_rank_lower_base"] == 24
    add("cubic-cyclic", 7, 14, 24, "four diagonal copies per grid point; seven-point degree-one grid",
        "The full formal scalar-output space C=lambda*S for the fixed seven-cycle shift.")
    for filename, name, row, column, scope in (
        ("cubic-jordan-filtration.json", "cubic-jordan", 14, 28,
         "The full formal scalar-output space C=lambda*J for the fixed strict shift."),
        ("cubic-fixed-image-rank-one.json", "cubic-fixed-image-rank-one", 9, 28,
         "The image vector u is fixed within this local factor; covectors vary. Distinct tensor factors may use different fixed u."),
        ("cubic-coordinate-image-mixtures.json", "cubic-coordinate-images", 14, 28,
         "All coordinate image choices e_i and covectors may mix inside this local factor; general multiple-row cross terms are excluded."),
    ):
        data = read(filename)
        assert data["labels"] == 7 and data["tensor_power_rank_lower_bases"] == [row, column]
        add(name, 7, row, column, scope=scope)
    assert 17 ** 2 > 2 * 12 ** 2
    result = {"status": "pass", "families": records,
              "source_certificate_sha256": hashes,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "cost_exponent": {"numerator": 17, "denominator": 12},
              "scope": "Every characteristic-two field; every nonempty mixed tensor product of the listed fixed local spaces; arbitrary coefficient-space recombinations within that tensor product.",
              "local_space_limit": "This is a tensor-product closure, not an arbitrary union of distinct local coefficient spaces. Each rank-one image vector remains fixed within its own local factor."}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
