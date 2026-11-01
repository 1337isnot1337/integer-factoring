#!/usr/bin/env python3
"""Independent reconstruction audit of the 26-label radical coefficient algebras."""

from pathlib import Path
import argparse
import hashlib
import json


def indices(value):
    while value:
        bit = value & -value
        yield bit.bit_length() - 1
        value ^= bit


def xor_all(values):
    out = 0
    for value in values:
        out ^= value
    return out


def basis(values):
    pivots = {}
    for value in values:
        while value:
            lead = value.bit_length() - 1
            if lead not in pivots:
                pivots[lead] = value
                break
            value ^= pivots[lead]
    return list(pivots.values())


def rank(values):
    return len(basis(values))


def add(first, second):
    return [[a ^ b for a, b in zip(x, y)] for x, y in zip(first, second)]


def product(first, second):
    return [[xor_all(right[j] for j in indices(row)) for row in left]
            for left, right in zip(first, second)]


def zero_like(matrix):
    return [[0] * len(block) for block in matrix]


def combine(matrices, coordinate):
    result = zero_like(matrices[0])
    for i in indices(coordinate):
        result = add(result, matrices[i])
    return result


def flattened(matrix):
    out = offset = 0
    for block in matrix:
        for row in block:
            out ^= row << offset
            offset += len(block)
    return out


def transposed(rows):
    return [sum(((row >> col) & 1) << i for i, row in enumerate(rows))
            for col in range(len(rows))]


def apply_columns(columns, vector):
    return xor_all(columns[i] for i in indices(vector))


def coordinate_product(table, a, b):
    return xor_all(table[i][j] for i in indices(a) for j in indices(b))


def characterize(words, values):
    result = 0
    for index, word in enumerate(words):
        value = 1
        if word != "E":
            for letter in word:
                value &= values["EABCDF".index(letter)]
        if value:
            result |= 1 << index
    return result


def audit_cut(data, supplied):
    cut = supplied["cut"]
    key = "basis_row_blocks" if cut == "row" else "basis_column_blocks"
    raw = [[[int(row, 16) for row in block] for block in matrix] for matrix in data[key]]
    generators = raw[:4] + [add(raw[1], raw[4]), add(raw[1], raw[5])]
    assert rank([flattened(matrix) for matrix in generators]) == 6
    identity = generators[0]
    assert all(rows == [1 << i for i in range(len(rows))] for rows in identity)
    words = supplied["basis_words"]
    matrices = []
    for word in words:
        value = identity
        if word != "E":
            for letter in word:
                value = product(value, generators["EABCDF".index(letter)])
        matrices.append(value)
    dimension = len(matrices)
    assert rank([flattened(matrix) for matrix in matrices]) == dimension
    table = supplied["multiplication_table"]
    products_checked = 0
    for i in range(dimension):
        for j in range(dimension):
            assert product(matrices[i], matrices[j]) == combine(matrices, table[i][j])
            products_checked += 1

    characters = supplied["characters_on_E_A_B_C_D_F"]
    functionals = [characterize(words, values) for values in characters]
    assert rank(functionals) == len(functionals)
    for functional in functionals:
        for i in range(dimension):
            for j in range(dimension):
                assert ((table[i][j] & functional).bit_count() & 1) == (
                    ((functional >> i) & 1) * ((functional >> j) & 1))
    radical = supplied["jacobson_radical_basis"]
    assert rank(radical) == dimension - len(functionals)
    assert all((j & chi).bit_count() % 2 == 0 for j in radical for chi in functionals)
    for j in radical:
        for i in range(6):
            assert all((coordinate_product(table, j, 1 << i) & chi).bit_count() % 2 == 0
                       and (coordinate_product(table, 1 << i, j) & chi).bit_count() % 2 == 0
                       for chi in functionals)
    powers = [radical]
    while powers[-1]:
        powers.append(basis(coordinate_product(table, a, b)
                            for a in powers[-1] for b in radical))
        assert len(powers) <= dimension + 1
    assert list(map(len, powers)) == supplied["radical_power_dimensions"]

    sizes = [len(block) for block in identity]
    spaces = [[[1 << i for i in range(size)] for size in sizes]]
    for ideal in powers:
        images = [[] for _ in sizes]
        for element in ideal:
            for block, rows in enumerate(combine(matrices, element)):
                images[block].extend(transposed(rows))
        spaces.append([basis(vectors) for vectors in images])
    assert [sum(map(len, v)) for v in spaces] == [
        layer["dimension"] for layer in supplied["representation_radical_filtration"]]

    # Reconstruct the actual change-of-basis columns, but verify similarities
    # by direct multiplication rather than the original coordinate decoder.
    exponent = 1
    while exponent < len(powers):
        exponent *= 2
    projectors = []
    for values in characters:
        candidate = identity
        for i in range(1, 4):
            selector = generators[i] if values[i] else add(generators[i], identity)
            candidate = product(candidate, selector)
        power = exponent
        while power > 1:
            candidate = product(candidate, candidate)
            power //= 2
        assert product(candidate, candidate) == candidate
        projectors.append(candidate)
    projector_columns = [[transposed(block) for block in matrix] for matrix in projectors]
    counts = []
    for layer in spaces:
        counts.append([sum(rank(apply_columns(matrix[b], v) for v in vectors)
                           for b, vectors in enumerate(layer))
                       for matrix in projector_columns])
    assert counts == [x["character_multiplicities"]
                      for x in supplied["representation_radical_filtration"]]

    similarities = 0
    for block, size in enumerate(sizes):
        pivots, vectors, labels = {}, [], []
        for depth in range(len(spaces) - 2, -1, -1):
            for character, projector in enumerate(projector_columns):
                for source in spaces[depth][block]:
                    original = apply_columns(projector[block], source)
                    reduced = original
                    while reduced:
                        lead = reduced.bit_length() - 1
                        if lead not in pivots:
                            pivots[lead] = reduced
                            vectors.append(original)
                            labels.append([depth, character])
                            break
                        reduced ^= pivots[lead]
        assert len(vectors) == rank(vectors) == size
        expected = supplied["adapted_original_blocks"][block]
        assert labels == expected["basis_depth_and_character"]
        for i, generator in enumerate(generators):
            actual_columns = transposed(generator[block])
            target_rows = [int(row, 16) for row in expected["generator_rows"][i]]
            target_columns = transposed(target_rows)
            for j, vector in enumerate(vectors):
                assert apply_columns(actual_columns, vector) == xor_all(
                    vectors[k] for k in indices(target_columns[j]))
                similarities += 1
    return {"cut": cut, "algebra_dimension": dimension,
            "all_basis_products_reconstructed": products_checked,
            "quotient_is_split_field_product_of_size": len(characters),
            "radical_power_dimensions": list(map(len, powers)),
            "module_filtration_dimensions": [sum(map(len, v)) for v in spaces],
            "character_multiplicities_by_depth": counts,
            "direct_change_of_basis_column_identities": similarities}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basis", type=Path, default=Path("evidence/symplectic-radical-P4-L4-basis.json"))
    parser.add_argument("--algebra", type=Path, default=Path("evidence/symplectic-radical-seed-algebra.json"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.loads(args.basis.read_text())
    algebra = json.loads(args.algebra.read_text())
    assert hashlib.sha256(args.basis.read_bytes()).hexdigest() == algebra["basis_sha256"]
    result = {
        "status": "Independent actual-matrix products, split quotient, nilpotent radical, module images, and change-of-basis identities passed.",
        "cuts": [audit_cut(data, supplied) for supplied in algebra["cuts"]],
        "basis_sha256": hashlib.sha256(args.basis.read_bytes()).hexdigest(),
        "algebra_evidence_sha256": hashlib.sha256(args.algebra.read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()

