#!/usr/bin/env python3
"""Replay the d24 seed's constant module bases and export weighted powers."""

import argparse
import hashlib
import json
from pathlib import Path

from check_modular_crossing import rank_two
from check_symplectic_seed_tensor_modules import ROOT, apply_columns, transpose
from symplectic_seed_weighted_modules import (
    coefficient_block, monomials, safe_integer_gate, type_records,
)


CERTIFICATE = ROOT / "evidence/symplectic-radical-sigma2-modules-certificate.json"


def verify(certificate):
    path = ROOT / certificate["source"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == certificate["source_sha256"]
    source = json.loads(path.read_text())
    reports = {}
    for name, family in certificate["families"].items():
        checks = 0
        for cut in ("row", "column"):
            original = next(record for record in source["cuts"] if record["cut"] == cut)
            counts = [0] * len(family["types"][cut])
            for block in family["blocks"][cut]:
                old = original["adapted_original_blocks"][block["source_block"]]
                operators = [transpose([int(value, 16) for value in matrix])
                             for matrix in old["generator_rows"][:family["operator_count"]]]
                basis = []
                for summand in block["summands"]:
                    index = summand["type"]
                    counts[index] += 1
                    module = family["types"][cut][index]
                    columns = [int(value, 16) for value in summand["basis_columns"]]
                    assert len(columns) == module["size"]
                    basis.extend(columns)
                    for operator, target_rows in zip(operators, module["generator_rows"]):
                        target = transpose(target_rows)
                        assert all(apply_columns(operator, vector) == apply_columns(columns, image)
                                   for vector, image in zip(columns, target))
                        checks += len(columns)
                assert len(basis) == block["size"] and rank_two(basis) == block["size"]
            assert counts == [module["multiplicity"] for module in family["types"][cut]]
            assert sum(module["size"] * module["multiplicity"] for module in family["types"][cut]) == 576
        reports[name] = {"intertwining_column_checks": checks,
                         "type_dimensions_and_multiplicities": {
                             cut: [[module["size"], module["multiplicity"]] for module in modules]
                             for cut, modules in family["types"].items()}}
    return reports


def export(certificate, path, power, cap):
    labels = certificate["labels"] ** power
    coordinates = monomials(power)
    gate = safe_integer_gate(labels)
    result = {"labels": labels, "dimension": len(coordinates), "gate": gate,
              "rank_bound": "exact" if cap is None else "lower", "maximum_block_size": cap,
              "monomials_v_s": coordinates, "cuts": {}}
    with path.open("w") as stream:
        stream.write(f"{labels} {len(coordinates)} {gate}\n")
        for cut in ("row", "column"):
            types = certificate["families"]["full"]["types"][cut]
            records = type_records(types, power)
            kept = [record for record in records if cap is None or record["size"] <= cap]
            omitted = sum(record["size"] * record["multiplicity"] for record in records
                          if cap is not None and record["size"] > cap)
            stream.write(f"{len(kept) + int(bool(omitted))}\n")
            if omitted:
                stream.write(f"1 {omitted}\n" + "0\n" * len(coordinates))
            for record in kept:
                size = record["size"]
                stream.write(f"{size} {record['multiplicity']}\n")
                for matrix in coefficient_block(types, record["types"]):
                    for row in matrix:
                        stream.write(" ".join(f"{(row >> shift) & ((1 << 64) - 1):x}"
                                              for shift in range(0, size, 64)) + "\n")
            assert sum(record["size"] * record["multiplicity"] for record in kept) + omitted == labels ** 2
            result["cuts"][cut] = {"retained_blocks": kept, "omitted_weighted_dimension": omitted}
    result["input_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--native", type=Path)
    parser.add_argument("--power", type=int, default=2)
    parser.add_argument("--maximum-block-size", type=int)
    args = parser.parse_args()
    certificate = json.loads(CERTIFICATE.read_text())
    result = {"checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "certificate_sha256": hashlib.sha256(CERTIFICATE.read_bytes()).hexdigest(),
              "families": verify(certificate)}
    if args.native:
        result["native_export"] = export(certificate, args.native, args.power, args.maximum_block_size)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "native_export"}, indent=2))


if __name__ == "__main__":
    main()
