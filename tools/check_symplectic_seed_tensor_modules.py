#!/usr/bin/env python3
"""Replay constant seed-module splittings and complete binary tensor-power costs.

The certificate contains explicit bases, so no random splitting search is part
of verification.  The two cuts may have different module types and sizes.
"""

import argparse
import hashlib
import heapq
from itertools import combinations_with_replacement
import json
from math import factorial, isqrt, log2, prod, sqrt
from pathlib import Path

from check_modular_crossing import rank_two


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "evidence/symplectic-seed-tensor-modules-certificate.json"


def bits(value):
    while value:
        bit = value & -value
        yield bit.bit_length() - 1
        value ^= bit


def apply_columns(matrix, vector):
    result = 0
    for i in bits(vector):
        result ^= matrix[i]
    return result


def transpose(rows):
    result = [0] * len(rows)
    for i, row in enumerate(rows):
        for j in bits(row):
            result[j] ^= 1 << i
    return result


def product_columns(first, second):
    return [apply_columns(first, column) for column in second]


def add(first, second):
    return [a ^ b for a, b in zip(first, second)]


def verify_certificate(certificate):
    for name, digest in certificate["sources"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest
    raw = json.loads((ROOT / "evidence/symplectic-radical-P4-L4-basis.json").read_text())
    algebra = json.loads((ROOT / "evidence/symplectic-radical-seed-algebra.json").read_text())
    row_source = next(record for record in algebra["cuts"] if record["cut"] == "row")
    checks = 0
    for cut in ("row", "column"):
        counts = [0] * len(certificate["types"][cut])
        for block in certificate["blocks"][cut]:
            index, size = block["source_block"], block["size"]
            if cut == "row":
                source_rows = row_source["adapted_original_blocks"][index]["generator_rows"][:4]
            else:
                source_rows = [matrix[index] for matrix in raw["basis_column_blocks"][:4]]
            source = [transpose([int(value, 16) for value in matrix]) for matrix in source_rows]
            basis = []
            for summand in block["summands"]:
                type_index = summand["type"]
                counts[type_index] += 1
                module = certificate["types"][cut][type_index]
                columns = [int(value, 16) for value in summand["basis_columns"]]
                assert len(columns) == module["size"]
                basis.extend(columns)
                for original, target_rows in zip(source, module["generator_rows"]):
                    target = transpose(target_rows)
                    assert all(apply_columns(original, vector) == apply_columns(columns, image)
                               for vector, image in zip(columns, target))
                    checks += len(columns)
            assert len(basis) == size and rank_two(basis) == size
        assert counts == [module["multiplicity"] for module in certificate["types"][cut]]
        assert sum(module["size"] * module["multiplicity"]
                   for module in certificate["types"][cut]) == 676

    # Exact relations in the column algebra, including both central summands.
    relation_ranks = {name: 0 for name in ("e", "A_complement", "A2_complement",
                                          "N", "AN", "eC")}
    for index, size in enumerate(raw["block_sizes"]):
        E, A, B, C = [transpose([int(value, 16) for value in matrix[index]])
                       for matrix in raw["basis_column_blocks"][:4]]
        N = add(A, B)
        A2 = product_columns(A, A)
        e = product_columns(A2, A)
        complement = add(E, e)
        zero = [0] * size
        assert product_columns(e, e) == e
        assert product_columns(e, A) == product_columns(e, B) == e
        assert product_columns(N, N) == product_columns(A2, N) == zero
        assert product_columns(C, C) == zero
        assert product_columns(complement, C) == product_columns(A, N)
        R1 = add(add(E, A), add(B, C))
        assert product_columns(R1, R1) == E
        assert product_columns(add(A, C), R1) == A
        assert product_columns(add(B, C), R1) == B
        assert product_columns(C, R1) == C
        matrices = (e, product_columns(complement, A), product_columns(complement, A2),
                    N, product_columns(A, N), product_columns(e, C))
        for name, matrix in zip(relation_ranks, matrices):
            relation_ranks[name] += rank_two(matrix)
    assert relation_ranks == {"e": 8, "A_complement": 179, "A2_complement": 16,
                              "N": 321, "AN": 22, "eC": 4}
    return {"intertwining_column_checks": checks, "column_relation_ranks": relation_ranks}


def kronecker(first, second):
    width = len(second)
    result = []
    for left in first:
        positions = list(bits(left))
        for right in second:
            row = 0
            for position in positions:
                row |= right << (position * width)
            result.append(row)
    return result


def polynomial_tensor(first, second):
    result = [[0] * (len(first[0]) * len(second[0]))
              for _ in range(len(first) + len(second) - 1)]
    for i, left in enumerate(first):
        for j, right in enumerate(second):
            result[i + j] = add(result[i + j], kronecker(left, right))
    return result


def rank_table(types, power):
    """All binary combinations of the coefficients of R(t,0)^tensor power."""
    coordinates = 3 * power + 1
    ranks = [0] * (1 << coordinates)
    maximum = count = 0
    for indices in combinations_with_replacement(range(len(types)), power):
        coefficient = [[1]]
        multiplicity = factorial(power)
        for index in set(indices):
            multiplicity //= factorial(indices.count(index))
        for index in indices:
            coefficient = polynomial_tensor(coefficient, types[index]["generator_rows"])
            multiplicity *= types[index]["multiplicity"]
        size = len(coefficient[0])
        maximum = max(maximum, size)
        count += 1
        rows, previous = [0] * size, 0
        for step in range(1, 1 << coordinates):
            mask = step ^ (step >> 1)
            changed = (mask ^ previous).bit_length() - 1
            rows = add(rows, coefficient[changed])
            previous = mask
            ranks[mask] += multiplicity * rank_two(rows)
    return ranks, {"type_blocks": count, "largest_type_block": maximum}


def shortest_cost(row, column, rounding):
    """Shortest XOR decomposition of coefficient zero, the crossing tensor.

    Integer floor weights give a rigorous lower bound on the real optimum;
    integer ceiling weights give a rigorous upper bound from an actual path.
    """
    costs = []
    for r, c in zip(row, column):
        value = r * c
        if rounding == "float":
            costs.append(sqrt(value))
        else:
            floor = isqrt(value)
            costs.append(floor + int(rounding == "ceiling" and floor * floor != value))
    distances = [None] * len(row)
    predecessors = [None] * len(row)
    distances[0] = 0
    queue = [(0, 0)]
    while queue:
        value, vertex = heapq.heappop(queue)
        if value != distances[vertex]:
            continue
        if vertex == 1:
            break
        for edge in range(1, len(row)):
            successor, candidate = vertex ^ edge, value + costs[edge]
            if distances[successor] is None or candidate < distances[successor]:
                distances[successor] = candidate
                predecessors[successor] = vertex, edge
                heapq.heappush(queue, (candidate, successor))
    path, vertex = [], 1
    while vertex:
        vertex, edge = predecessors[vertex]
        path.append(edge)
    return {"cost": distances[1], "masks": path,
            "rank_pairs": [[row[edge], column[edge]] for edge in path]}


def power_report(types, power):
    tables, metadata = {}, {}
    for cut in ("row", "column"):
        tables[cut], metadata[cut] = rank_table(types[cut], power)
    row, column = tables["row"], tables["column"]
    if power == 1:
        independent = [(0, 0), (676, 676), (404, 187), (560, 668), (656, 243),
                       (63, 668), (580, 321), (397, 676), (396, 26), (580, 676),
                       (62, 187), (656, 672), (560, 243), (405, 672), (676, 325), (1, 676)]
        assert list(zip(row, column)) == independent
    lower = shortest_cost(row, column, "floor")
    upper = shortest_cost(row, column, "ceiling")
    numerical = shortest_cost(row, column, "float")
    labels = 26 ** power
    # This stronger rational-exponent comparison implies failure at sqrt(2).
    assert lower["cost"] ** 2 > labels ** 3
    return {"power": power, "labels": labels, "coefficient_coordinates": 3 * power + 1,
            "rank_tables": tables, "module_evaluation": metadata,
            "integer_lower": lower, "integer_upper": upper, "numerical_optimum": numerical,
            "numerical_rho": log2(numerical["cost"]) / (2 * log2(labels)),
            "integer_lower_squared_exceeds_labels_cubed": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-power", type=int, default=2, choices=(1, 2))
    args = parser.parse_args()
    certificate = json.loads(CERTIFICATE.read_text())
    result = {"checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "certificate_sha256": hashlib.sha256(CERTIFICATE.read_bytes()).hexdigest(),
              "certificate_verification": verify_certificate(certificate),
              "type_dimensions_and_multiplicities": {
                  cut: [[module["size"], module["multiplicity"]] for module in modules]
                  for cut, modules in certificate["types"].items()},
              "powers": [power_report(certificate["types"], power)
                         for power in range(1, args.max_power + 1)]}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"certificate_verification": result["certificate_verification"],
                      "powers": [{key: value for key, value in record.items()
                                  if key not in ("rank_tables",)} for record in result["powers"]]},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
