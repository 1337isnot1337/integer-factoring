#!/usr/bin/env python3
"""Bounded Jordan-module certificates for restricted seven-site seed blocks.

No full seven-site coefficient matrix is constructed.  The largest audited
matrix has dimension2560.  The default rebuilds the five-site literal core;
--basis5 can reuse its independently certified JSON.
"""

import argparse
from collections import Counter, defaultdict
from itertools import combinations
import hashlib
import json
from math import comb
from pathlib import Path
import random

from assemble_two_even_site_cores import core, monomials, normalize
from check_common_character_compressions import audit_block, code80, coefficient_digest, matvec
from check_noncommutative_crossing_character_compression import reconstructed_cut
from check_modular_crossing import rank_two
from complement_shell_tensor import multiply
from explore_complement_edge_pencils import evaluate_local
from explore_symplectic_boundary_codes import coefficient_basis
from find_tagged_common_character_compressions import search_block


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERTIFICATE = ROOT / "evidence/two-even-site-P7-jordan-compression-certificate.json"


def add(first, second):
    return [a ^ b for a, b in zip(first, second)]


def tensor(first, second):
    size = len(second)
    result = []
    for row in first:
        for other in second:
            value, output = row, 0
            while value:
                bit = value & -value
                output ^= other << ((bit.bit_length() - 1) * size)
                value ^= bit
            result.append(output)
    return result


def free_jordan(copies, cut):
    result = [[1]]
    for _ in range(copies):
        old = result
        size = 2 * len(old[0])
        result = [[0] * size for _ in range(len(old) + 1)]
        linear = [3, 2] if cut == 0 else [2, 0]
        for degree, matrix in enumerate(old):
            result[degree] = add(result[degree], tensor(matrix, [1, 2]))
            result[degree + 1] = add(result[degree + 1], tensor(matrix, linear))
    return result


def groups(code):
    grouped = defaultdict(list)
    for i, first in enumerate(code["syndromes"]):
        for j, second in enumerate(code["syndromes"]):
            grouped[first ^ second].append((i, j))
    return grouped


def build_cores(path5):
    path4 = ROOT / "evidence/symplectic-radical-sigma2-basis.json"
    data4 = json.loads(path4.read_text())
    code = code80()
    if path5:
        data5 = json.loads(path5.read_text())
        for key in ("input_vectors", "output_indices", "syndromes", "ambient_occupations"):
            assert data5["code"][key] == json.loads(json.dumps(code[key]))
        rows, columns = [[[[int(value, 16) for value in block] for block in matrix]
                          for matrix in data5[name]]
                         for name in ("basis_row_blocks", "basis_column_blocks")]
    else:
        data5, rows, columns = coefficient_basis(code)
        data5["code"] = code
    certified = json.loads((ROOT / "evidence/symplectic-two-even-P5-compression-certificate.json").read_text())
    assert coefficient_digest(rows, columns) == certified["coefficient_matrices_sha256"]
    result = {4: dict(core(4), code=data4["code"]),
              5: {"code": data5["code"], "monomials": monomials(5),
                  "cuts": [normalize(data5, rows), normalize(data5, columns)]}}
    for value in result.values():
        value["groups"] = groups(value["code"])
        value["syndromes"] = sorted(value["groups"])
    return result


def core_row_triangle():
    basis_path = ROOT / "evidence/symplectic-radical-sigma2-basis.json"
    algebra_path = ROOT / "evidence/symplectic-radical-sigma2-algebra.json"
    basis = json.loads(basis_path.read_text())
    algebra = json.loads(algebra_path.read_text())
    assert hashlib.sha256(basis_path.read_bytes()).hexdigest() == algebra["basis_sha256"]
    row = next(cut for cut in algebra["cuts"] if cut["cut"] == "row")
    reconstructed = reconstructed_cut(basis, row)
    character = [1, 0, 1, 0, 0, 0]
    index = row["characters_on_E_A_B_C_D_F"].index(character)
    block = row["adapted_original_blocks"][-1]
    assert block["size"] == 96
    selected = [i for i, (_, value) in enumerate(block["basis_depth_and_character"]) if value == index]
    assert len(selected) == 32
    checks = 0
    for scalar, encoded in zip(character, block["generator_rows"]):
        rows = [int(value, 16) for value in encoded]
        for i, position in enumerate(selected):
            assert (rows[position] >> position & 1) == scalar
            checks += 1
            for old in selected[:i]:
                assert not rows[position] >> old & 1
                checks += 1
    return {"core_syndrome": 15, "core_character_dimension": 32,
            "triangular_matrix_entries_checked": checks,
            "independent_algebra_reconstruction": reconstructed,
            "basis_sha256": hashlib.sha256(basis_path.read_bytes()).hexdigest(),
            "algebra_sha256": hashlib.sha256(algebra_path.read_bytes()).hexdigest()}


def module(cores, pairs, core_pairs, syndrome, jordan_copies, cut):
    data = cores[core_pairs]
    block = data["syndromes"].index(syndrome)
    matrices = [matrix[block] for matrix in data["cuts"][cut]]
    free = free_jordan(jordan_copies, cut)
    size = len(matrices[0]) * (1 << jordan_copies)
    positions = {term: i for i, term in enumerate(monomials(pairs))}
    result = [[0] * size for _ in positions]
    scalar_copies = pairs - core_pairs - jordan_copies
    shifts = [i for i in range(scalar_copies + 1) if comb(scalar_copies, i) % 2] if cut == 0 else [0]
    for (degree, charge), matrix in zip(data["monomials"], matrices):
        for extra, other in enumerate(free):
            product = tensor(matrix, other)
            for shift in shifts:
                index = positions[degree + extra + shift, charge]
                result[index] = add(result[index], product)
    return result


def local_similarity():
    records = []
    for name, operator, basis, target in (
        ("row", [1, 4, 2, 8], [1, 8, 6, 2], [1, 2, 12, 8]),
        ("column", [9, 0, 0, 9], [2, 4, 9, 1], [0, 0, 8, 0]),
    ):
        assert rank_two(basis) == 4
        for j, vector in enumerate(basis):
            expected = 0
            for i in range(4):
                if target[i] >> j & 1:
                    expected ^= basis[i]
            assert matvec(operator, vector) == expected
        records.append({"cut": name, "basis_columns": basis, "linear_target_rows": target})
    return records


def embedding_controls(cores, certificate):
    pairs = certificate["pairs"]
    variable_p, variable_q = 2, 1 << (pairs + 1)
    local = evaluate_local((variable_p, variable_q, 1, 1, 1, 1))

    def power(value, exponent):
        result = 1
        for _ in range(exponent):
            result = multiply(result, value)
        return result

    def supports(word):
        edge = [i for i, value in enumerate(word) if value in (0, 3)]
        assert len(edge) == 2
        result = []
        for full in edge:
            value = list(word)
            for site in edge:
                value[site] = 3 * int(site == full)
            result.append(value)
        return result

    physical = []
    for sites in combinations(range(pairs), 4):
        for cut in range(2):
            physical.append((4, 15, sites, cut, "entire weight-four block"))
    for record in certificate["three_anchor_restrictions"]:
        for cut, name in enumerate(("row", "column")):
            sites = tuple(sorted(record["syndrome_sites"] + record[name]["anchors"]))
            physical.append((5, record[name]["core_syndrome"], sites, cut, "weight-two principal restriction"))
    for cut in range(2):
        physical.append((5, 0, tuple(certificate["zero_restriction_sites"]), cut, "weight-zero principal restriction"))
    rng = random.Random(709133)
    counts = Counter()
    for h, syndrome, sites, cut, kind in physical:
        data = cores[h]
        group = data["groups"][syndrome]
        block = data["syndromes"].index(syndrome)
        weights = [multiply(multiply(power(1 ^ variable_q, h - a - 2 * b),
                                     power(variable_p ^ variable_q, a)), power(variable_q, b))
                   for a, b in data["monomials"]]
        for _ in range(3):
            row, column = rng.randrange(len(group)), rng.randrange(len(group))
            first, second = group[row]
            third, fourth = group[column]
            if cut == 0:
                labels = (first, second, fourth, third)  # x,y,u,v
            else:
                labels = (first, fourth, second, third)
            words = [[0] * pairs for _ in range(4)]
            for role, label in enumerate(labels):
                source = data["code"]["ambient_occupations"][data["code"]["output_indices"][label]]
                for site, value in zip(sites, source):
                    words[role][site] = value
            expected = 0
            for matrix, weight in zip(data["cuts"][cut], weights):
                if matrix[block][row] >> column & 1:
                    expected ^= weight
            for site in range(pairs):
                if site in sites:
                    continue
                for word in words:
                    word[site] = 1 + rng.randrange(2)
                x, y, u, v = (word[site] - 1 for word in words)
                ri, ci = (2 * x + y, 2 * v + u) if cut == 0 else (2 * x + u, 2 * v + y)
                linear = ([1, 4, 2, 8] if cut == 0 else [9, 0, 0, 9])[ri]
                factor = ((1 ^ variable_q) if ri == ci else 0) ^ (
                    (variable_p ^ variable_q) if linear >> ci & 1 else 0)
                expected = multiply(expected, factor)
            actual = 0
            for yy in supports(words[1]):
                for vv in supports(words[3]):
                    value = 1
                    for x, y, u, v in zip(words[0], yy, words[2], vv):
                        value = multiply(value, local[((x * 4 + y) * 4 + u) * 4 + v])
                    actual ^= value
            assert actual == expected
            counts[kind, cut] += 1
    return [{"kind": kind, "cut": cut, "complete_polynomial_entries_checked": count}
            for (kind, cut), count in sorted(counts.items())]


def restriction_totals(certificate, values, key, core_pairs):
    result = [0, 0]
    records = certificate[key]
    assert len(records) == comb(7, 2)
    assert {tuple(record["syndrome_sites"]) for record in records} == set(combinations(range(7), 2))
    for record in records:
        syndrome = record["syndrome_sites"]
        for cut, name in enumerate(("row", "column")):
            supplied = record[name]
            sites = sorted(syndrome + supplied["anchors"])
            assert len(sites) == len(set(sites)) == core_pairs
            assert all(0 <= site < 7 for site in sites)
            relative = sum(1 << sites.index(site) for site in syndrome)
            assert relative == supplied["core_syndrome"]
            value = values[core_pairs, relative, name]
            assert value == supplied["compression"]
            result[cut] += value
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--basis5", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    certificate = json.loads(args.certificate.read_text())
    assert certificate["pairs"] == 7 and certificate["labels"] == 672
    character = certificate["character"]
    assert character == [comb(5, i) % 2 if i <= 5 else 0 for i in range(7)] + [0] * 5
    cores = build_cores(args.basis5)
    similarity = local_similarity()
    triangle = core_row_triangle()
    literal = embedding_controls(cores, certificate)
    values, records = {}, []
    for supplied in certificate["types"]:
        h, syndrome = supplied["core_pairs"], supplied["core_syndrome"]
        free_copies = 7 - h
        assert h in (4, 5) and syndrome in cores[h]["groups"]
        core_size = len(cores[h]["groups"][syndrome])
        for cut, name in enumerate(("row", "column")):
            reports, total, dimension = [], 0, 0
            assert len(supplied["cuts"][name]) == free_copies + 1
            for j, item in enumerate(supplied["cuts"][name]):
                assert item["j"] == j
                multiplicity = comb(free_copies, j) * (1 << (free_copies - j))
                assert item["multiplicity"] == multiplicity
                matrices = module(cores, 7, h, syndrome, j, cut)
                witness = item["witness"]
                assert witness["size"] == core_size * (1 << j)
                assert len(set(witness["coordinate_indices"])) == witness["compression"]
                expanded = dict(witness, vectors=[hex(1 << i) for i in witness["coordinate_indices"]])
                if witness["compression"]:
                    report = audit_block(matrices, character, expanded)
                else:
                    assert witness["annihilated_span"] == 0
                    report = {"compression_dimension": 0, "all_generator_compression_entries_checked": 0}
                total += multiplicity * witness["compression"]
                dimension += multiplicity * witness["size"]
                reports.append(dict(report, jordan_copies=j, multiplicity=multiplicity))
            assert dimension == core_size * 4 ** free_copies
            values[h, syndrome, name] = total
            records.append({"core_pairs": h, "core_syndrome": syndrome, "cut": name,
                            "weighted_compression_dimension": total, "modules": reports})
            print(h, syndrome, name, total, flush=True)
    weight_four = [comb(7, 4) * values[4, 15, name] for name in ("row", "column")]
    two = restriction_totals(certificate, values, "two_anchor_restrictions", 4)
    three = restriction_totals(certificate, values, "three_anchor_restrictions", 5)
    assert certificate["zero_restriction_sites"] == list(range(5))
    zero = [values[5, 0, name] for name in ("row", "column")]
    intermediate = [a + b for a, b in zip(weight_four, two)]
    scalar_totals = [a + b + c for a, b, c in zip(weight_four, three, zero)]
    assert weight_four == [50785, 560] and intermediate == [72483, 990]
    assert scalar_totals == [87420, 1150]
    stronger_row = comb(7, 4) * triangle["core_character_dimension"] * 4 ** 3
    totals = [scalar_totals[0] - weight_four[0] + stronger_row, scalar_totals[1]]
    assert totals == [108315, 1150]
    rank_product = totals[0] * totals[1]
    assert rank_product ** 7 > 672 ** 20 and 10 ** 2 > 2 * 7 ** 2
    checks = sum(item["all_generator_compression_entries_checked"]
                 for record in records for item in record["modules"])
    result = {"status": "pass", "labels": 672, "coefficient_directions": 12,
              "character": character, "character_scalar_polynomial": "(q+r)^2 (p+r)^5",
              "local_constant_similarities": similarity, "literal_seven_site_embedding_controls": literal,
              "weight_four_row_triangular_certificate": triangle,
              "modules": records, "weight_four_only_compressions": weight_four,
              "two_anchor_plus_weight_four_compressions": intermediate,
              "all_scalar_compression_dimensions": scalar_totals,
              "final_rank_certificate_dimensions": totals, "rank_product": rank_product,
              "all_generator_compression_entries_checked": checks,
              "largest_matrix_dimension": max(item["witness"]["size"]
                   for supplied in certificate["types"] for items in supplied["cuts"].values() for item in items),
              "exact_margin": {"rank_product_power": 7, "label_power": 20,
                               "difference": rank_product ** 7 - 672 ** 20},
              "cost_exponent": {"numerator": 10, "denominator": 7},
              "scope": "All characteristic-two fields, every tensor power, arbitrary coefficient-space recombinations.",
              "certificate_sha256": hashlib.sha256(args.certificate.read_bytes()).hexdigest(),
              "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                  for name in (Path(__file__).name, "assemble_two_even_site_cores.py",
                               "check_common_character_compressions.py", "find_tagged_common_character_compressions.py",
                               "check_noncommutative_crossing_character_compression.py",
                               "explore_symplectic_boundary_codes.py")}}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
