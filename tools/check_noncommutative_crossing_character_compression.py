#!/usr/bin/env python3
"""Independent common-character rank certificate for the 24/26-label seeds."""

from collections import Counter
from pathlib import Path
import argparse
import hashlib
import json

from check_noncommutative_crossing_seed_algebra_audit import (
    audit_cut, indices, rank, transposed, xor_all,
)


CHI = (1, 0, 1, 0, 0, 0)
WITNESSES = {0: (2, 3), 1: (16, 20), 5: (32, 36),
             6: (64, 80), 8: (32, 36), 9: (8, 10)}


def number(value):
    return int(value, 16) if isinstance(value, str) else value


def reconstructed_cut(data, supplied):
    # The d24 generator writes words as generator-index lists; the independent
    # d26 audit writes the same words as strings.  Normalize only this encoding.
    supplied = dict(supplied)
    supplied["basis_words"] = [
        "".join("EABCDF"[i] for i in word) or "E"
        if isinstance(word, list) else word
        for word in supplied["basis_words"]
    ]
    return audit_cut(data, supplied)


def matvec(rows, vector):
    return sum(((row & vector).bit_count() % 2) << i for i, row in enumerate(rows))


def inverse(rows):
    size = len(rows)
    work = [row | (1 << (size + i)) for i, row in enumerate(rows)]
    for col in range(size):
        pivot = next(i for i in range(col, size) if work[i] >> col & 1)
        work[col], work[pivot] = work[pivot], work[col]
        for i in range(size):
            if i != col and work[i] >> col & 1:
                work[i] ^= work[col]
    assert [row & ((1 << size) - 1) for row in work] == [1 << i for i in range(size)]
    return [row >> size for row in work]


def triangular_character(algebra):
    row = next(cut for cut in algebra["cuts"] if cut["cut"] == "row")
    character = row["characters_on_E_A_B_C_D_F"].index(list(CHI))
    count = 0
    checked = 0
    for block in row["adapted_original_blocks"]:
        selected = [i for i, (_, value) in enumerate(block["basis_depth_and_character"])
                    if value == character]
        for g, encoded in enumerate(block["generator_rows"]):
            rows = list(map(number, encoded))
            for i, position in enumerate(selected):
                assert (rows[position] >> position & 1) == CHI[g]
                assert not any(rows[position] >> old & 1 for old in selected[:i])
                checked += 1
        count += len(selected)
    assert count == row["representation_radical_filtration"][0]["character_multiplicities"][character]
    return {"multiplicity": count, "literal_triangular_diagonal_checks": checked}


def column_compression(algebra, certificate):
    column = next(cut for cut in algebra["cuts"] if cut["cut"] == "column")
    family = certificate["families"]["full"]
    types, blocks = family["types"]["column"], family["blocks"]["column"]
    multiplicities = Counter()
    left, right = [], []
    global_generators = [[] for _ in CHI]
    offset = 0
    similarities = 0
    records = []
    for block in blocks:
        size = block["size"]
        supplied = column["adapted_original_blocks"][block["source_block"]]
        assert supplied["size"] == size
        generators = [list(map(number, rows)) for rows in supplied["generator_rows"]]
        for g, rows in enumerate(generators):
            global_generators[g].extend(row << offset for row in rows)
        basis_columns = [number(value) for summand in block["summands"] for value in summand["basis_columns"]]
        assert len(basis_columns) == rank(basis_columns) == size
        inverse_rows = inverse(transposed(basis_columns))
        start = 0
        for summand in block["summands"]:
            kind = summand["type"]
            data = types[kind]
            columns = list(map(number, summand["basis_columns"]))
            dim = data["size"]
            assert len(columns) == dim
            multiplicities[kind] += 1
            for g, encoded in enumerate(data["generator_rows"]):
                source = list(map(number, encoded))
                for j, vector in enumerate(columns):
                    image = xor_all(columns[i] for i in range(dim) if source[i] >> j & 1)
                    assert matvec(generators[g], vector) == image
                    similarities += 1
            if kind in WITNESSES:
                v, w = WITNESSES[kind]
                values = [((matvec(list(map(number, rows)), v) & w).bit_count() % 2)
                          for rows in data["generator_rows"]]
                assert tuple(values) == CHI
                right.append(xor_all(columns[i] for i in indices(v)) << offset)
                left.append(xor_all(inverse_rows[start + i] for i in indices(w)) << offset)
            start += dim
        assert start == size
        offset += size
    assert offset == algebra["labels"] ** 2
    for kind, data in enumerate(types):
        assert multiplicities[kind] == data["multiplicity"]
        if kind in WITNESSES:
            records.append({"type": kind, "size": data["size"], "multiplicity": multiplicities[kind],
                            "right_vector": WITNESSES[kind][0], "left_covector": WITNESSES[kind][1]})
    equations = 0
    for g, generator in enumerate(global_generators):
        images = [matvec(generator, vector) for vector in right]
        for i, row in enumerate(left):
            for j, vector in enumerate(images):
                assert (row & vector).bit_count() % 2 == CHI[g] * int(i == j)
                equations += 1
    assert len(left) == len(right) == 43
    return {"dimension": 43, "type_witnesses": records,
            "all_full_column_intertwining_columns": similarities,
            "explicit_compression_matrix_entries_checked": equations,
            "left_covectors_in_adapted_direct_sum": [hex(row) for row in left],
            "right_vectors_in_adapted_direct_sum": [hex(vector) for vector in right]}


def audit_inclusion(smaller, larger, supplied):
    chosen = supplied["label_inclusion"]
    assert chosen == list(range(2, 26))
    assert smaller["basis_monomials_p_q_r"] == larger["basis_monomials_p_q_r"]
    assert smaller["basis_scalar_polynomials_p_q_r"] == larger["basis_scalar_polynomials_p_q_r"]
    for key in ("input_vectors", "output_indices", "syndromes"):
        assert smaller["code"][key] == [larger["code"][key][i] for i in chosen]
    count = 0
    for key in ("basis_row_blocks", "basis_column_blocks"):
        for low, high in zip(smaller[key], larger[key]):
            for first, second, positions in zip(low, high, supplied["cut_block_inclusions"]):
                assert len(first) == len(positions)
                for i, row in enumerate(first):
                    source = number(second[positions[i]])
                    restricted = sum(((source >> index) & 1) << j for j, index in enumerate(positions))
                    assert restricted == number(row)
                    count += 1
    assert count == 6912
    return {"labels": chosen, "all_six_generators_both_cuts_bit_rows_checked": count}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    paths = {
        "basis24": Path("evidence/symplectic-radical-sigma2-basis.json"),
        "algebra24": Path("evidence/symplectic-radical-sigma2-algebra.json"),
        "modules24": Path("evidence/symplectic-radical-sigma2-modules-certificate.json"),
        "inclusion": Path("evidence/symplectic-radical-sigma2-inclusion.json"),
        "basis26": Path("evidence/symplectic-radical-P4-L4-basis.json"),
        "algebra26": Path("evidence/symplectic-radical-seed-algebra.json"),
    }
    data = {key: json.loads(path.read_text()) for key, path in paths.items()}
    hashes = {key: hashlib.sha256(path.read_bytes()).hexdigest() for key, path in paths.items()}
    assert hashes["basis24"] == data["algebra24"]["basis_sha256"]
    assert hashes["basis26"] == data["algebra26"]["basis_sha256"]
    assert hashes["algebra24"] == data["modules24"]["source_sha256"]
    reconstructed = [reconstructed_cut(data["basis24"], cut) for cut in data["algebra24"]["cuts"]]
    reconstructed.append(reconstructed_cut(data["basis26"], next(
        cut for cut in data["algebra26"]["cuts"] if cut["cut"] == "row")))
    rows = {str(d): triangular_character(data[f"algebra{d}"]) for d in (24, 26)}
    assert rows["24"]["multiplicity"] == 384 and rows["26"]["multiplicity"] == 436
    compression = column_compression(data["algebra24"], data["modules24"])
    inclusion = audit_inclusion(data["basis24"], data["basis26"], data["inclusion"])
    gates = [{"labels": d, "row_character_multiplicity": rows[str(d)]["multiplicity"],
              "column_compression_dimension": 43, "rank_product": rows[str(d)]["multiplicity"] * 43,
              "label_cube": d ** 3} for d in (24, 26)]
    assert all(item["rank_product"] > item["label_cube"] for item in gates)
    result = {"status": "pass", "character_E_A_B_C_D_F": CHI,
              "scope": "Every characteristic-two field, every tensor power and mixed tensor product, arbitrary coefficient-space recombinations.",
              "row_triangular_certificates": rows, "column_compression": compression,
              "exact_label_restriction": inclusion, "integer_separation_certificates": gates,
              "independent_actual_matrix_reconstructions": reconstructed,
              "source_hashes": hashes,
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
