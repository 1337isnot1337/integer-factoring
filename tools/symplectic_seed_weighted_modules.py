#!/usr/bin/env python3
"""Exact small modules and weighted common-parameter powers of the radical seed."""

import argparse
import hashlib
from itertools import combinations_with_replacement
import json
from math import factorial, prod
from pathlib import Path

from check_modular_crossing import rank_two
from check_symplectic_seed_tensor_modules import (
    ROOT, add, apply_columns, bits, kronecker, transpose,
)


CERTIFICATE = ROOT / "evidence/symplectic-seed-weighted-modules-certificate.json"
SEED_MONOMIALS = ((0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1))


def load_certificate():
    return json.loads(CERTIFICATE.read_text())


def verify_certificate(certificate):
    source_path = ROOT / certificate["source"]
    assert hashlib.sha256(source_path.read_bytes()).hexdigest() == certificate["source_sha256"]
    source = json.loads(source_path.read_text())
    checks = 0
    for cut in ("row", "column"):
        original = next(record for record in source["cuts"] if record["cut"] == cut)
        counts = [0] * len(certificate["types"][cut])
        for block in certificate["blocks"][cut]:
            source_block = original["adapted_original_blocks"][block["source_block"]]
            operators = [transpose([int(value, 16) for value in matrix])
                         for matrix in source_block["generator_rows"]]
            basis = []
            for summand in block["summands"]:
                index = summand["type"]
                counts[index] += 1
                module = certificate["types"][cut][index]
                embedding = [int(value, 16) for value in summand["basis_columns"]]
                assert len(embedding) == module["size"]
                basis.extend(embedding)
                for operator, small_rows in zip(operators, module["generator_rows"]):
                    small = transpose(small_rows)
                    assert all(apply_columns(operator, vector) == apply_columns(embedding, image)
                               for vector, image in zip(embedding, small))
                    checks += len(embedding)
            assert len(basis) == block["size"] and rank_two(basis) == block["size"]
        assert counts == [module["multiplicity"] for module in certificate["types"][cut]]
        assert sum(module["size"] * module["multiplicity"]
                   for module in certificate["types"][cut]) == 676
    return {"intertwining_column_checks": checks,
            "type_dimensions_and_multiplicities": {
                cut: [[module["size"], module["multiplicity"]] for module in modules]
                for cut, modules in certificate["types"].items()}}


def polynomial_multiply(first, second):
    result = 0
    while second:
        if second & 1:
            result ^= first
        second >>= 1
        first <<= 1
    return result


def scalar_curve_polynomials(u, v, s):
    """Binary polynomial bitsets multiplying E,A,B,C,D,F in K(u,v,s)."""
    multiply = polynomial_multiply
    return (multiply(multiply(u, u), u), multiply(multiply(u, u), v),
            multiply(u, multiply(v, v)), multiply(multiply(v, v), v),
            multiply(u, s), multiply(v, s))


def specialize_types(types, u, v, s):
    scalars = scalar_curve_polynomials(u, v, s)
    degree = max(scalars).bit_length() - 1
    result = []
    for module in types:
        matrices = [[0] * module["size"] for _ in range(degree + 1)]
        for scalar, operator in zip(scalars, module["generator_rows"]):
            for index in bits(scalar):
                matrices[index] = add(matrices[index], operator)
        result.append({"size": module["size"], "multiplicity": module["multiplicity"],
                       "generator_rows": matrices})
    return result


def monomials(power):
    """(v-degree,s-degree); u-degree is 3*power-v_degree-2*s_degree."""
    return [(v, s) for s in range(power + 1) for v in range(3 * power - 2 * s + 1)]


def type_records(types, power):
    records = []
    for indices in combinations_with_replacement(range(len(types)), power):
        multiplicity = factorial(power)
        for index in set(indices):
            multiplicity //= factorial(indices.count(index))
        multiplicity *= prod(types[index]["multiplicity"] for index in indices)
        records.append({"types": indices, "size": prod(types[index]["size"] for index in indices),
                        "multiplicity": multiplicity})
    return sorted(records, key=lambda record: (record["size"], -record["multiplicity"]))


def coefficient_block(types, indices):
    coefficients = {(0, 0): [1]}
    for index in indices:
        module = types[index]
        size = len(next(iter(coefficients.values()))) * module["size"]
        successor = {}
        for (v, s), matrix in coefficients.items():
            for (extra_v, extra_s), operator in zip(SEED_MONOMIALS, module["generator_rows"]):
                key = (v + extra_v, s + extra_s)
                old = successor.setdefault(key, [0] * size)
                successor[key] = add(old, kronecker(matrix, operator))
        coefficients = successor
    return [coefficients.get(key, [0] * len(next(iter(coefficients.values()))))
            for key in monomials(len(indices))]


def safe_integer_gate(labels):
    """Ceiling(labels^(99/70)); 99/70 exceeds sqrt(2) by exact squaring."""
    assert 99 * 99 > 2 * 70 * 70
    target = labels ** 99
    low, high = 0, 1
    while high ** 70 < target:
        high *= 2
    while high - low > 1:
        middle = (low + high) // 2
        if middle ** 70 < target:
            low = middle
        else:
            high = middle
    return high


def export_native(certificate, path, power, maximum_block_size=None, gate=None):
    """Export complete or capped rank blocks for the weighted native screen.

    Omitted blocks are replaced by a zero summand of the same weighted
    dimension. Their computed ranks are therefore certified lower bounds.
    """
    labels = 26 ** power
    coordinates = monomials(power)
    gate = safe_integer_gate(labels) if gate is None else gate
    metadata = {"labels": labels, "dimension": len(coordinates), "gate": gate,
                "monomials_v_s": coordinates, "maximum_block_size": maximum_block_size,
                "rank_bound": "exact" if maximum_block_size is None else "lower",
                "cuts": {}}
    with Path(path).open("w") as stream:
        stream.write(f"{labels} {len(coordinates)} {gate}\n")
        for cut in ("row", "column"):
            types = certificate["types"][cut]
            records = type_records(types, power)
            retained = [record for record in records
                        if maximum_block_size is None or record["size"] <= maximum_block_size]
            omitted = sum(record["size"] * record["multiplicity"] for record in records
                          if maximum_block_size is not None and record["size"] > maximum_block_size)
            stream.write(f"{len(retained) + int(bool(omitted))}\n")
            if omitted:
                stream.write(f"1 {omitted}\n")
                stream.write("0\n" * len(coordinates))
            for record in retained:
                size, multiplicity = record["size"], record["multiplicity"]
                stream.write(f"{size} {multiplicity}\n")
                words = (size + 63) // 64
                for matrix in coefficient_block(types, record["types"]):
                    for row in matrix:
                        stream.write(" ".join(f"{(row >> (64 * word)) & ((1 << 64) - 1):x}"
                                              for word in range(words)) + "\n")
            assert sum(record["size"] * record["multiplicity"] for record in retained) + omitted == labels ** 2
            metadata["cuts"][cut] = {"retained_blocks": retained, "omitted_weighted_dimension": omitted}
    metadata["input_sha256"] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--native", type=Path)
    parser.add_argument("--power", type=int, default=2)
    parser.add_argument("--maximum-block-size", type=int)
    parser.add_argument("--gate", type=int)
    args = parser.parse_args()
    certificate = load_certificate()
    result = {"checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "certificate_sha256": hashlib.sha256(CERTIFICATE.read_bytes()).hexdigest(),
              "verification": verify_certificate(certificate)}
    if args.native:
        result["native_export"] = export_native(certificate, args.native, args.power,
                                                args.maximum_block_size, args.gate)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "native_export"},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
