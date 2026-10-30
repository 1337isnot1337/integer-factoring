#!/usr/bin/env python3
"""Independent all-field scalar-compression certificate for the 240-label seed."""

import argparse
import hashlib
import json
from pathlib import Path

from assemble_two_even_site_cores import assemble, independent_coordinates, literal_controls
from check_common_character_compressions import audit_block


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "evidence/two-even-site-P6-compression-certificate.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=CERTIFICATE)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    certificate = json.loads(args.certificate.read_text())
    assert certificate["pairs"] == 6 and certificate["labels"] == 240
    character = certificate["character"]
    assert character == [1, 0, 0, 0, 1, 0, 0, 0, 0, 0]
    metadata, cuts = assemble(6)
    coordinates = independent_coordinates(cuts[0])
    controls = literal_controls(metadata, cuts)
    result_cuts = {}
    for name, matrices in zip(("row", "column"), cuts):
        supplied = certificate["cuts"][name]
        assert len(supplied) == len(metadata["block_sizes"])
        reports = []
        for block, witness in enumerate(supplied):
            assert len(set(witness["coordinate_indices"])) == witness["compression"]
            assert all(0 <= index < witness["size"] for index in witness["coordinate_indices"])
            expanded = dict(witness, vectors=[hex(1 << index) for index in witness["coordinate_indices"]])
            report = audit_block([matrix[block] for matrix in matrices], character, expanded)
            reports.append(report)
            print(name, block, report["compression_dimension"], flush=True)
        result_cuts[name] = {"blocks": reports,
                             "compression_dimension": sum(record["compression_dimension"] for record in reports),
                             "all_generator_compression_entries_checked": sum(
                                 record["all_generator_compression_entries_checked"] for record in reports)}
    first, second = (result_cuts[cut]["compression_dimension"] for cut in ("row", "column"))
    assert (first, second) == (15688, 718)
    rank_product = first * second
    assert rank_product ** 25 > 240 ** 74
    assert 37 ** 2 > 2 * 25 ** 2
    result = {"status": "pass", "labels": 240, "coefficient_dimension": 10,
              "character": character, "character_scalar_polynomial": "(q+r)^2 (p+r)^4",
              "cuts": result_cuts, "rank_product": rank_product,
              "certified_cost_exponent": {"numerator": 37, "denominator": 25},
              "exact_margin": {"rank_product_power": 25, "label_power": 74,
                               "difference": rank_product ** 25 - 240 ** 74},
              "independent_coordinate_witnesses": coordinates,
              "literal_polynomial_controls": controls,
              "scope": "Every characteristic-two field, every tensor power, arbitrary coefficient-space recombinations.",
              "certificate_sha256": hashlib.sha256(args.certificate.read_bytes()).hexdigest(),
              "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                  for name in (Path(__file__).name, "assemble_two_even_site_cores.py",
                               "check_common_character_compressions.py",
                               "find_tagged_common_character_compressions.py")}}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
