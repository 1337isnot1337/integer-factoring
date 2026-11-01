#!/usr/bin/env python3
"""Exact coefficient algebra and Loewy filtrations of the 26-label radical seed."""

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path

from check_modular_crossing import rank_two


NAMES = ("E", "A", "B", "C", "D", "F")


def polynomial_product(first, second):
    result = set()
    for a in first:
        for b in second:
            term = tuple(x + y for x, y in zip(a, b))
            result.symmetric_difference_update((term,))
    return result


def polynomial_power(value, exponent):
    result = {(0, 0, 0)}
    for _ in range(exponent):
        result = polynomial_product(result, value)
    return result


def verify_normalization(data):
    u = {(0, 1, 0), (0, 0, 1)}
    v = {(1, 0, 0), (0, 1, 0)}
    s = {(0, 1, 1)}
    normalized = [polynomial_power(u, 4),
                  polynomial_product(polynomial_power(u, 3), v),
                  polynomial_product(polynomial_power(u, 2), polynomial_power(v, 2)),
                  polynomial_product(u, polynomial_power(v, 3)),
                  polynomial_product(s, polynomial_power(u, 2)),
                  polynomial_product(s, polynomial_product(u, v))]
    # D=B1+B4 and F=B1+B5 in the exported independent coefficient basis.
    raw = normalized[:]
    raw[1] = raw[1] ^ raw[4] ^ raw[5]
    assert raw == [{tuple(term) for term in polynomial}
                   for polynomial in data["basis_scalar_polynomials_p_q_r"]]


def bits(value):
    while value:
        bit = value & -value
        yield bit.bit_length() - 1
        value ^= bit


def plus(first, second):
    return [[a ^ b for a, b in zip(x, y)] for x, y in zip(first, second)]


def matrix_product(first, second):
    result = []
    for left, right in zip(first, second):
        block = []
        for row in left:
            value = 0
            for i in bits(row):
                value ^= right[i]
            block.append(value)
        result.append(block)
    return result


def flatten(matrix):
    result = offset = 0
    for block in matrix:
        for row in block:
            result |= row << offset
            offset += len(block)
    return result


def matrix_rank(matrix):
    return sum(rank_two(block) for block in matrix)


def linear_basis(vectors):
    pivots = {}
    for value in vectors:
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return list(pivots.values())


def columns(rows):
    result = [0] * len(rows)
    for i, row in enumerate(rows):
        for j in bits(row):
            result[j] ^= 1 << i
    return result


def apply_columns(matrix_columns, vector):
    result = 0
    for i in bits(vector):
        result ^= matrix_columns[i]
    return result


class BlockAlgebra:
    def __init__(self, generators):
        self.generators = generators
        self.matrices = []
        self.words = []
        self.pivots = {}
        for i, matrix in enumerate(generators):
            self.insert(matrix, () if i == 0 else (i,))
        index = 0
        while index < len(self.matrices):
            for i, generator in enumerate(generators[1:], 1):
                self.insert(matrix_product(generator, self.matrices[index]),
                            (i,) + self.words[index])
            index += 1
            assert len(self.matrices) <= 128
        self.dimension = len(self.matrices)
        self.table = [[self.coordinates(matrix_product(first, second))
                       for second in self.matrices] for first in self.matrices]

    def insert(self, matrix, word):
        value = flatten(matrix)
        coordinate = 1 << len(self.matrices)
        while value:
            pivot = value.bit_length() - 1
            if pivot not in self.pivots:
                self.pivots[pivot] = value, coordinate
                self.matrices.append(matrix)
                self.words.append(word)
                return True
            old, combination = self.pivots[pivot]
            value ^= old
            coordinate ^= combination
        return False

    def coordinates(self, matrix):
        value, coordinate = flatten(matrix), 0
        while value:
            pivot = value.bit_length() - 1
            old, combination = self.pivots[pivot]
            value ^= old
            coordinate ^= combination
        return coordinate

    def multiply(self, first, second):
        result = 0
        for i in bits(first):
            for j in bits(second):
                result ^= self.table[i][j]
        return result

    def power(self, value, exponent):
        result = 1
        while exponent:
            if exponent & 1:
                result = self.multiply(result, value)
            value = self.multiply(value, value)
            exponent >>= 1
        return result

    def matrix(self, coordinate):
        result = [[0] * len(block) for block in self.generators[0]]
        for i in bits(coordinate):
            result = plus(result, self.matrices[i])
        return result

    def ideal(self, generators):
        vectors = linear_basis(generators)
        pivots = {value.bit_length() - 1: value for value in vectors}
        index = 0
        while index < len(vectors):
            for i in range(1, 6):
                for value in (self.multiply(1 << i, vectors[index]),
                              self.multiply(vectors[index], 1 << i)):
                    while value:
                        pivot = value.bit_length() - 1
                        if pivot not in pivots:
                            pivots[pivot] = value
                            vectors.append(value)
                            break
                        value ^= pivots[pivot]
            index += 1
        return vectors

    def ideal_powers(self, ideal):
        result = [ideal]
        while result[-1]:
            successor = linear_basis(self.multiply(a, b) for a in result[-1] for b in ideal)
            assert len(successor) < len(result[-1]), "the proposed ideal is not nilpotent"
            result.append(successor)
        return result

    def characters(self):
        result = []
        for triple in product((0, 1), repeat=3):
            generator_values = (1,) + triple + (0, 0)
            values = []
            for word in self.words:
                value = 1
                for i in word:
                    value &= generator_values[i]
                values.append(value)
            functional = sum(value << i for i, value in enumerate(values))
            if all((self.table[i][j] & functional).bit_count() % 2 == values[i] * values[j]
                   for i in range(self.dimension) for j in range(self.dimension)):
                result.append((generator_values, functional))
        return result


def representation_filtration(algebra, radical_powers, characters):
    sizes = [len(block) for block in algebra.generators[0]]
    image_spaces = [[[1 << i for i in range(size)] for size in sizes]]
    for ideal in radical_powers:
        spans = [[] for _ in sizes]
        for element in ideal:
            for block, rows in enumerate(algebra.matrix(element)):
                spans[block].extend(columns(rows))
        image_spaces.append([linear_basis(vectors) for vectors in spans])
    assert all(not block for block in image_spaces[-1])

    # A lift of each character selector becomes idempotent after a Frobenius
    # power at least the nilpotence index of the radical.
    exponent = 1
    while exponent < len(radical_powers):
        exponent *= 2
    projectors = []
    for values, _ in characters:
        value = 1
        for i in range(1, 4):
            value = algebra.multiply(value, (1 << i) ^ int(values[i] == 0))
        value = algebra.power(value, exponent)
        assert algebra.multiply(value, value) == value
        projectors.append(algebra.matrix(value))
    projector_columns = [[columns(block) for block in matrix] for matrix in projectors]

    filtration = []
    for spaces in image_spaces:
        multiplicities = [sum(rank_two([apply_columns(matrix[block], vector)
                                        for vector in vectors])
                              for block, vectors in enumerate(spaces))
                          for matrix in projector_columns]
        dimension = sum(map(len, spaces))
        assert sum(multiplicities) == dimension
        filtration.append({"dimension": dimension, "character_multiplicities": multiplicities})
    layers = [{"dimension": filtration[i]["dimension"] - filtration[i + 1]["dimension"],
               "character_multiplicities": [a - b for a, b in zip(
                   filtration[i]["character_multiplicities"],
                   filtration[i + 1]["character_multiplicities"])]}
              for i in range(len(filtration) - 1)]

    # Construct actual bases adapted to both radical depth and quotient character.
    adapted = []
    for block, size in enumerate(sizes):
        pivots = {}
        vectors = []
        labels = []
        for depth in range(len(image_spaces) - 2, -1, -1):
            for character, matrix in enumerate(projector_columns):
                for source in image_spaces[depth][block]:
                    original = apply_columns(matrix[block], source)
                    value = original
                    coordinates = 1 << len(vectors)
                    while value:
                        pivot = value.bit_length() - 1
                        if pivot not in pivots:
                            pivots[pivot] = value, coordinates
                            vectors.append(original)
                            labels.append((depth, character))
                            break
                        old, combination = pivots[pivot]
                        value ^= old
                        coordinates ^= combination
        assert len(vectors) == size
        transformed = []
        for generator_index, generator in enumerate(algebra.generators):
            original_columns = columns(generator[block])
            rows = [0] * size
            for source_index, source in enumerate(vectors):
                value = apply_columns(original_columns, source)
                coordinate = 0
                while value:
                    pivot = value.bit_length() - 1
                    old, combination = pivots[pivot]
                    value ^= old
                    coordinate ^= combination
                for target in bits(coordinate):
                    assert labels[target][0] >= labels[source_index][0]
                    if labels[target][0] == labels[source_index][0]:
                        assert target == source_index
                expected_diagonal = characters[labels[source_index][1]][0][generator_index]
                assert ((coordinate >> source_index) & 1) == expected_diagonal
                for target in bits(coordinate):
                    rows[target] ^= 1 << source_index
            transformed.append([hex(row) for row in rows])
        adapted.append({"size": size, "basis_depth_and_character": labels,
                        "generator_rows": transformed})
    return filtration, layers, adapted


def analyze(data, cut):
    supplied = data["basis_row_blocks" if cut == "row" else "basis_column_blocks"]
    raw = [[[int(row, 16) for row in block] for block in matrix] for matrix in supplied]
    generators = raw[:4] + [plus(raw[1], raw[4]), plus(raw[1], raw[5])]
    assert all(block == [1 << i for i in range(len(block))] for block in generators[0])
    algebra = BlockAlgebra(generators)
    commutators = [algebra.table[i][j] ^ algebra.table[j][i]
                   for i in range(1, 6) for j in range(i + 1, 6)]
    df_ideal = algebra.ideal((1 << 4, 1 << 5))
    df_powers = algebra.ideal_powers(df_ideal)
    radical = algebra.ideal(commutators + [
        algebra.table[i][i] ^ (1 << i) for i in range(1, 4)] + [1 << 4, 1 << 5])
    radical_powers = algebra.ideal_powers(radical)
    characters = algebra.characters()
    assert len(characters) == algebra.dimension - len(radical)
    # Distinct field characters are independent; this separately checks it.
    assert rank_two([functional for _, functional in characters]) == len(characters)
    expected = {"row": (45, [42, 26, 6, 0], [22, 0], 3),
                "column": (11, [9, 4, 0], [5, 0], 2)}[cut]
    assert algebra.dimension == expected[0]
    assert list(map(len, radical_powers)) == expected[1]
    assert list(map(len, df_powers)) == expected[2]
    assert len(characters) == expected[3]
    filtration, layers, adapted = representation_filtration(algebra, radical_powers, characters)
    minimum_polynomials = []
    for i in range(1, 6):
        value = 1
        pivots = {}
        ranks = []
        for degree in range(algebra.dimension + 1):
            ranks.append(matrix_rank(algebra.matrix(value)))
            vector, coefficients = value, 1 << degree
            while vector:
                pivot = vector.bit_length() - 1
                if pivot not in pivots:
                    pivots[pivot] = vector, coefficients
                    break
                old, combination = pivots[pivot]
                vector ^= old
                coefficients ^= combination
            if not vector:
                minimum_polynomials.append({"operator": NAMES[i],
                                            "nonzero_degrees": list(bits(coefficients)),
                                            "power_ranks": ranks})
                break
            value = algebra.multiply(1 << i, value)
    return {
        "cut": cut, "module_dimension": data["labels"] ** 2,
        "algebra_dimension": algebra.dimension,
        "basis_words": ["".join(NAMES[i] for i in word) or "E" for word in algebra.words],
        "multiplication_table": algebra.table,
        "jacobson_radical_basis": radical,
        "radical_power_dimensions": list(map(len, radical_powers)),
        "DF_ideal_power_dimensions": list(map(len, df_powers)),
        "characters_on_E_A_B_C_D_F": [list(values) for values, _ in characters],
        "minimal_polynomials": minimum_polynomials,
        "representation_radical_filtration": filtration,
        "representation_Loewy_layers": layers,
        "adapted_original_blocks": adapted,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basis", type=Path, default=Path("evidence/symplectic-radical-P4-L4-basis.json"))
    parser.add_argument("--cut", choices=("row", "column", "both"), default="both")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.loads(args.basis.read_text())
    assert data["labels"] == 26 and data["basis_dimension"] == 6
    assert data["basis_monomials_p_q_r"] == [[0, 0, 4], [1, 0, 3], [2, 0, 2],
                                            [3, 0, 1], [0, 1, 3], [1, 1, 2]]
    verify_normalization(data)
    cuts = ("row", "column") if args.cut == "both" else (args.cut,)
    result = {"labels": 26, "normalization": "G=u^4 R(v/u,qr/u^2); R=E+tA+t^2B+t^3C+zD+ztF",
              "scalar_normalization_verified": True,
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "basis_sha256": hashlib.sha256(args.basis.read_bytes()).hexdigest(), "cuts": []}
    for cut in cuts:
        record = analyze(data, cut)
        result["cuts"].append(record)
        print(json.dumps({key: record[key] for key in (
            "cut", "algebra_dimension", "radical_power_dimensions",
            "characters_on_E_A_B_C_D_F", "representation_radical_filtration",
            "representation_Loewy_layers")}), flush=True)
        if args.output:
            args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()

