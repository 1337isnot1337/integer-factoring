#!/usr/bin/env python3
"""Faithful local spectral identities, negligible residuals and longer slides."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from bounded_spin_crossing import paths
from check_fixed_field_descent import Field
from paired_path_crossing import power
from quantum_block_braid import braid_word
from quantum_braid_path import local_action
from quantum_braid_spectral_fusion import (
    FIELD,
    apply_linear_gate,
    fused_string_coefficients,
    matrix_product as mm,
    normalized_string_scalar,
    normalized_three_one_coefficients,
    physical_generator,
    slide_scalars,
    string_scalar,
    symmetrizer,
    three_one_coefficients,
)

ROOT = Path(__file__).resolve().parents[1]
WORDS = (
    (),
    (0,),
    (1,),
    (2,),
    (0, 1),
    (0, 2),
    (1, 0),
    (1, 2),
    (2, 1),
    (0, 1, 2),
    (0, 2, 1),
    (1, 0, 2),
    (2, 1, 0),
    (1, 0, 2, 1),
)


def pivot_columns(matrix, field):
    a = matrix.tolist()
    columns, r = [], 0
    for j in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][j]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        inv = field.power(a[r][j], 14)
        a[r] = [field.mul(inv, v) for v in a[r]]
        for i in range(r + 1, len(a)):
            factor = a[i][j]
            a[i] = [x ^ field.mul(factor, y) for x, y in zip(a[i], a[r])]
        columns.append(j)
        r += 1
        if r == len(a):
            break
    return columns


def determinant(matrix, field):
    a, result = matrix.tolist(), 1
    for j in range(len(a)):
        p = next((i for i in range(j, len(a)) if a[i][j]), None)
        if p is None:
            return 0
        a[j], a[p] = a[p], a[j]  # Row-swap signs coincide in characteristic two.
        pivot = a[j][j]
        result = field.mul(result, pivot)
        inv = field.power(pivot, 14)
        a[j] = [field.mul(inv, v) for v in a[j]]
        for i in range(j + 1, len(a)):
            factor = a[i][j]
            a[i] = [x ^ field.mul(factor, y) for x, y in zip(a[i], a[j])]
    return result


def gram_projector(strands, q, field):
    frame = [[0] * (strands + 1) for _ in range(1 << strands)]
    for word in range(1 << strands):
        inversions = sum(
            j - (word & ((1 << j) - 1)).bit_count()
            for j in range(strands)
            if word >> j & 1
        )
        frame[word][word.bit_count()] = power(q, -inversions)
    transpose = list(map(list, zip(*frame)))
    gram = field.matmul(transpose, frame)
    inverse = [
        [field.power(gram[i][i], 14) if i == j else 0 for j in range(strands + 1)]
        for i in range(strands + 1)
    ]
    assert all(gram[i][i] for i in range(strands + 1))
    return np.array(
        field.matmul(field.matmul(frame, inverse), transpose), dtype=np.uint8
    )


def poly_product(a, b):
    result = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            key = i + k, j + l
            if key not in result:
                result[key] = np.zeros_like(x)
            result[key] ^= mm(x, y)
    return result


def scalar_eval(coefficients, value, field):
    answer = 0
    for c in coefficients[::-1]:
        answer = field.mul(answer, value) ^ c
    return answer


def literal_scalar_product(a, b, field):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] ^= field.mul(x, y)
    return result


def check():
    field = Field(0b10011)
    for a, b in product(range(16), repeat=2):
        assert FIELD[a, b] == field.mul(a, b)
    models, ambient_entries, long_entries, scalar_controls, yb_coefficients = (
        [],
        0,
        0,
        0,
        0,
    )
    for q in (8, 10, 12, 15):
        identity = np.eye(16, dtype=np.uint8)
        es = [physical_generator(4, i, q) for i in range(3)]
        bs = [FIELD[q, identity] ^ e for e in es]
        inverses = [FIELD[power(q, -1), identity] ^ e for e in es]
        images = []
        for word in WORDS:
            image = identity.copy()
            for site in word:
                image = mm(es[site], image)
            images.append(image.reshape(-1))
        images = np.array(images)
        columns = pivot_columns(images, field)
        assert len(columns) == 14
        minor = images[:, columns]
        det = determinant(minor, field)
        assert det and field.rank(minor.tolist()) == 14
        for i in range(3):
            assert np.array_equal(mm(bs[i], inverses[i]), identity)
            unitarity = (
                mm(inverses[i], bs[i]),
                mm(inverses[i], inverses[i]) ^ mm(bs[i], bs[i]),
                mm(bs[i], inverses[i]),
            )
            delta = q ^ power(q, -1)
            assert all(
                np.array_equal(a, FIELD[c, identity])
                for a, c in zip(unitarity, (1, field.mul(delta, delta), 1))
            )

        def r(site, monomial):
            return {(0, 0): inverses[site], monomial: bs[site]}

        lhs = poly_product(poly_product(r(0, (1, 0)), r(1, (1, 1))), r(0, (0, 1)))
        rhs = poly_product(poly_product(r(1, (0, 1)), r(0, (1, 1))), r(1, (1, 0)))
        keys = lhs.keys() | rhs.keys()
        for key in keys:
            assert np.array_equal(
                lhs.get(key, np.zeros_like(identity)),
                rhs.get(key, np.zeros_like(identity)),
            )
        yb_coefficients += len(keys)
        f3, f4 = symmetrizer(3, q), symmetrizer(4, q)
        assert np.array_equal(f3, gram_projector(3, q, field))
        assert np.array_equal(f4, gram_projector(4, q, field))
        left, right = np.kron(np.eye(2, dtype=np.uint8), f3), np.kron(
            f3, np.eye(2, dtype=np.uint8)
        )
        low = left ^ f4
        assert [field.rank(a.tolist()) for a in (left, f4, low)] == [8, 5, 3]
        assert np.array_equal(mm(left, f4), f4)
        assert not np.any(FIELD[q, f4[:8, :8]] ^ FIELD[power(q, -1), f4[8:, 8:]])
        coefficients, braid, projector = three_one_coefficients(q)
        assert np.array_equal(projector, f4)
        p, residual = slide_scalars(q)
        p_reference = [power(q, 2)]
        for exponent in (0, 3, 4):
            p_reference = literal_scalar_product(
                p_reference, [1, power(q, exponent)], field
            )
        assert p == p_reference
        common = literal_scalar_product([power(q, 2)], [1, 1], field)
        common = literal_scalar_product(common, [1, power(q, 3)], field)
        assert residual == [field.mul(1 ^ power(q, 3), v) for v in common] + [0]
        normalized = normalized_three_one_coefficients(q)
        expanded = np.zeros_like(coefficients)
        for i, c in enumerate(common):
            for j, a in enumerate(normalized):
                expanded[i + j] ^= FIELD[c, a]
        assert np.array_equal(expanded, coefficients)
        for c, x, y in zip(coefficients, p, residual):
            assert np.array_equal(c, FIELD[x, braid] ^ FIELD[y, f4])
            assert np.array_equal(mm(right, c), c)
            ambient_entries += c.size
        assert scalar_eval(p, q, field) == 0 and scalar_eval(residual, q, field)
        exceptional = np.zeros_like(identity)
        for c in coefficients[::-1]:
            exceptional = FIELD[q, exceptional] ^ c
        assert field.rank(exceptional.tolist()) == 5
        assert np.array_equal(exceptional, FIELD[scalar_eval(residual, q, field), f4])
        bad = left[None, ...]
        for site in (2, 1, 0):
            bad = apply_linear_gate(bad, inverses[site], bs[site])
        assert any(np.any(mm(identity ^ right, c)) for c in bad)
        for value in range(16):
            literal = left.tolist()
            for site, exponent in zip((2, 1, 0), (0, -2, -4)):
                local = (
                    inverses[site]
                    ^ FIELD[field.mul(value, power(q, exponent)), bs[site]]
                )
                literal = field.matmul(local.tolist(), literal)
            evaluated = np.zeros_like(identity)
            for c in coefficients[::-1]:
                evaluated = FIELD[value, evaluated] ^ c
            assert evaluated.tolist() == literal
            scalar_controls += evaluated.size
        low_braid = mm(braid, low)
        returned = low_braid.copy()
        for site in (0, 1, 2):
            returned = mm(bs[site], returned)
        assert np.array_equal(returned, FIELD[power(q, 3), low])

        longer = []
        for m, expected_rank in ((3, 1), (5, 3), (7, 8)):
            polynomial, constant = fused_string_coefficients(m, q)
            scalar = string_scalar(m, q)
            assert len(polynomial) == len(scalar)
            for a, c in zip(polynomial, scalar):
                assert np.array_equal(a, FIELD[c, constant])
                long_entries += a.size
            assert field.rank(constant.tolist()) == expected_rank
            longer.append(
                {
                    "free_strands": m,
                    "vacuum_states": len(constant),
                    "input_channel_rank": expected_rank,
                    "degree": len(polynomial) - 1,
                }
            )
        # At m=3 the vacuum input f3 necessarily pairs with the other f3.
        ps = paths(6)
        index = {v: i for i, v in enumerate(ps)}
        embed = (0, 1, delta, delta ^ 1)
        generators = []
        for site in range(5):
            e = np.zeros((len(ps), len(ps)), dtype=np.uint8)
            for j, path in enumerate(ps):
                for target, c in local_action(path, site):
                    e[index[target], j] = embed[c]
            generators.append(e)
        path_identity = np.eye(len(ps), dtype=np.uint8)

        def f3_path(site):
            a, b = generators[site : site + 2]
            return path_identity ^ a ^ b ^ FIELD[power(delta, -1), mm(a, b) ^ mm(b, a)]

        pair_projector = mm(f3_path(0), f3_path(3))
        assert field.rank(pair_projector.tolist()) == 1
        crossed = pair_projector.copy()
        for site in braid_word(3, 3):
            crossed = mm(FIELD[q, path_identity] ^ generators[site], crossed)
        assert np.array_equal(crossed, FIELD[power(q, 2), pair_projector])
        sparse = []
        for m in (5, 10, 20, 40):
            raw, normal = string_scalar(m, q), normalized_string_scalar(m, q)
            assert all(v == int(i in (0, m, 2 * m, 3 * m)) for i, v in enumerate(raw))
            assert all(v == int(i in (0, m)) for i, v in enumerate(normal))
            sparse.append(
                {
                    "free_strands": m,
                    "raw_support": [0, m, 2 * m, 3 * m],
                    "normalized_support": [0, m],
                }
            )
        models.append(
            {
                "q": q,
                "faithful_TL4_words_application_order": WORDS,
                "faithful_minor_columns": columns,
                "faithful_minor": minor.tolist(),
                "faithful_minor_determinant": det,
                "quotient_scalar": p,
                "negligible_residual": residual,
                "ordinary_channel_ranks": [8, 5, 3],
                "exceptional_ambient_rank_at_z_q": 5,
                "unshifted_negative_control": True,
                "raw_double_braid_phase": power(q, 3),
                "raw_self_braid_phase": power(q, 2),
                "longer_transport_controls": longer,
                "sparse_scalar_controls": sparse,
            }
        )
    sources = (
        "research/fused-spectral-transport.md",
        "tools/quantum_braid_spectral_fusion.py",
        "tools/check_quantum_braid_spectral_fusion.py",
        "tools/quantum_braid_path.py",
        "tools/check_fixed_field_descent.py",
        "tools/paired_path_crossing.py",
        "tools/quantum_block_braid.py",
        "tools/bounded_spin_crossing.py",
        "tools/bounded_spin_support.py",
        "tools/check_fourier_multiplication_gates.py",
    )
    return {
        "status": "pass",
        "field_multiplication_controls": 256,
        "ambient_slide_coefficient_entries": ambient_entries,
        "yang_baxter_coefficient_matrices": yb_coefficients,
        "independent_literal_scalar_matrix_entries": scalar_controls,
        "longer_slide_coefficient_entries": long_entries,
        "models": models,
        "scope": "Exact local TL polynomial fusion, explicit nonzero negligible residual, and uniform transport by composition; no complete crossing rank-sum or stronger factoring bound.",
        "source_sha256": {
            f: sha256((ROOT / f).read_bytes()).hexdigest() for f in sources
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
