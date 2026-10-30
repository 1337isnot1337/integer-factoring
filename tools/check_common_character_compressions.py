#!/usr/bin/env python3
"""Construct and audit constant common-functional compressions of two cuts.

The default reconstructs the complete 80-label, five-site, two-even-site
Pfaffian seed.  A supplied basis JSON avoids the literal reconstruction,
but receives the same exact matrix digest and compression audit.
"""

import argparse
import hashlib
import json
from pathlib import Path

from explore_symplectic_boundary_codes import boundary_code, coefficient_basis


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERTIFICATE = ROOT / "evidence/symplectic-two-even-P5-compression-certificate.json"


def bits(value):
    while value:
        bit = value & -value
        yield bit.bit_length() - 1
        value ^= bit


def matvec(rows, vector):
    return sum(((row & vector).bit_count() % 2) << i for i, row in enumerate(rows))


def low_basis(vectors):
    """Use low pivots, independently of the witness search's high pivots."""
    pivots = {}
    for vector in vectors:
        while vector:
            pivot = (vector & -vector).bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = vector
                break
            vector ^= pivots[pivot]
    return list(pivots.values())


def dual_rows(zero_vectors, vectors, ambient):
    """Solve all dual systems together, then check the actual scalar maps."""
    pivots = {}
    equations = [(vector, 0) for vector in zero_vectors]
    equations += [(vector, 1 << i) for i, vector in enumerate(vectors)]
    for vector, target in equations:
        while vector:
            pivot = (vector & -vector).bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = vector, target
                break
            old, image = pivots[pivot]
            vector ^= old
            target ^= image
        else:
            assert target == 0
    assert len(pivots) == len(zero_vectors) + len(vectors)
    solution = [0] * ambient
    for pivot in sorted(pivots, reverse=True):
        vector, target = pivots[pivot]
        for index in bits(vector ^ (1 << pivot)):
            target ^= solution[index]
        solution[pivot] = target
    rows = [0] * len(vectors)
    for position, value in enumerate(solution):
        for index in bits(value):
            rows[index] |= 1 << position
    for row in rows:
        assert all((row & vector).bit_count() % 2 == 0 for vector in zero_vectors)
    assert all((row & vector).bit_count() % 2 == int(i == j)
               for i, row in enumerate(rows) for j, vector in enumerate(vectors))
    return rows


def audit_block(matrices, character, witness):
    ambient = len(matrices[0])
    assert ambient == witness["size"]
    assert matrices[0] == [1 << i for i in range(ambient)]
    assert character[0] == 1 and len(character) == len(matrices)
    vectors = [int(value, 16) for value in witness["vectors"]]
    assert all(0 < value < 1 << ambient for value in vectors)
    assert len(vectors) == witness["compression"]
    images = [[matvec(matrix, vector) for vector in vectors] for matrix in matrices]
    differences = [image ^ (vector if scalar else 0)
                   for scalar, image_vectors in zip(character, images)
                   for image, vector in zip(image_vectors, vectors)]
    zero_vectors = low_basis(differences)
    assert len(zero_vectors) == witness["annihilated_span"]
    assert len(low_basis(zero_vectors + vectors)) == len(zero_vectors) + len(vectors)
    rows = dual_rows(zero_vectors, vectors, ambient)
    checks = 0
    for scalar, image_vectors in zip(character, images):
        for i, row in enumerate(rows):
            for j, image in enumerate(image_vectors):
                assert (row & image).bit_count() % 2 == scalar * int(i == j)
                checks += 1
    encoded = json.dumps([hex(row) for row in rows], separators=(",", ":")).encode()
    return {"ambient_dimension": ambient, "compression_dimension": len(vectors),
            "annihilated_span_dimension": len(zero_vectors),
            "all_generator_compression_entries_checked": checks,
            "constructed_left_covectors_sha256": hashlib.sha256(encoded).hexdigest()}


def code80():
    original = boundary_code(5, 5)
    selected = [i for i, syndrome in enumerate(original["syndromes"])
                if syndrome.bit_count() == 3]
    code = dict(original)
    code["labels"] = len(selected)
    for key in ("input_vectors", "output_indices", "syndromes"):
        code[key] = [original[key][i] for i in selected]
    for key in ("input_Gram_rows", "induced_mode_complement_rows"):
        code[key] = [sum(((original[key][i] >> j) & 1) << k
                        for k, j in enumerate(selected)) for i in selected]
    assert code["labels"] == 80
    assert all(vector.bit_count() == 2 for vector in code["input_vectors"])
    return code


def coefficient_digest(rows, columns):
    encoded = json.dumps([rows, columns], separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def polynomial_product(first, second):
    result = set()
    for a in first:
        for b in second:
            term = tuple(x + y for x, y in zip(a, b))
            if term in result:
                result.remove(term)
            else:
                result.add(term)
    return result


def audit_scalar_polynomial(data, character):
    actual = set()
    for scalar, polynomial in zip(character, data["basis_scalar_polynomials_p_q_r"]):
        if scalar:
            actual.symmetric_difference_update(tuple(term) for term in polynomial)
    expected = {(0, 0, 0)}
    for factor in [{(0, 1, 0), (0, 0, 1)}] * 2 + [{(1, 0, 0), (0, 0, 1)}] * 3:
        expected = polynomial_product(expected, factor)
    assert actual == expected
    return sorted(actual)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--basis", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    certificate = json.loads(args.certificate.read_text())
    character = certificate["character"]
    code = code80()
    if args.basis:
        data = json.loads(args.basis.read_text())
        for key in ("input_vectors", "output_indices", "syndromes", "ambient_occupations"):
            assert data["code"][key] == json.loads(json.dumps(code[key]))
        rows, columns = [[[list(map(lambda value: int(value, 16), block)) for block in matrix]
                          for matrix in data[name]]
                         for name in ("basis_row_blocks", "basis_column_blocks")]
    else:
        data, rows, columns = coefficient_basis(code)
    assert data["labels"] == certificate["labels"] == 80
    assert data["basis_dimension"] == len(character) == 8
    assert data["endpoint_coordinates"] == [1] + [0] * 7
    assert data["p_endpoint_zero"] and data["q_endpoint_ranks"] == [1, 6400]
    digest = coefficient_digest(rows, columns)
    assert digest == certificate["coefficient_matrices_sha256"]
    cuts = {}
    for name, matrices in (("row", rows), ("column", columns)):
        assert len(certificate["cuts"][name]) == len(data["block_sizes"])
        reports = [audit_block([matrix[block] for matrix in matrices], character, witness)
                   for block, witness in enumerate(certificate["cuts"][name])]
        cuts[name] = {"compression_dimension": sum(item["compression_dimension"] for item in reports),
                      "all_generator_compression_entries_checked": sum(
                          item["all_generator_compression_entries_checked"] for item in reports),
                      "blocks": reports}
    rank_product = cuts["row"]["compression_dimension"] * cuts["column"]["compression_dimension"]
    assert rank_product == 592592 and rank_product > 80 ** 3
    result = {"status": "pass", "labels": 80, "coefficient_space_dimension": 8,
              "character": character, "character_scalar_polynomial": "(q+r)^2 (p+r)^3",
              "character_scalar_monomials_p_q_r": audit_scalar_polynomial(data, character),
              "cuts": cuts, "rank_product": rank_product, "label_cube": 80 ** 3,
              "coefficient_matrices_sha256": digest,
              "literal_contracted_entries": data["contracted_entries_checked"],
              "literal_entries_summed": data["literal_entries_summed"],
              "scope": "Every characteristic-two field, every tensor power, arbitrary coefficient-space recombinations.",
              "certificate_sha256": hashlib.sha256(args.certificate.read_bytes()).hexdigest(),
              "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                 for name in (Path(__file__).name, "explore_symplectic_boundary_codes.py")}}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
