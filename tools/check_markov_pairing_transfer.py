#!/usr/bin/env python3
"""Independent tensor-word verification of finite-state pairing constructions."""
from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import random
import subprocess
from tempfile import TemporaryDirectory


def bits(word):
    while word:
        bit = word & -word
        yield bit.bit_length() - 1
        word ^= bit


def poly_product(a, b):
    return xor_all(a << j for j in bits(b))


def xor_all(values):
    answer = 0
    for value in values:
        answer ^= value
    return answer


def matrix_product(a, b):
    return [
        xor_all(poly_product(a[2 * i + k], b[2 * k + j]) for k in range(2))
        for i in range(2) for j in range(2)
    ]


def rank(rows):
    basis = {}
    for word in rows:
        while word:
            pivot = word.bit_length() - 1
            if pivot in basis:
                word ^= basis[pivot]
            else:
                basis[pivot] = word
                break
    return len(basis)


def local(symbol, quad):
    x, y, u, v = quad
    return (
        int(x == y and u == v),
        int(x ^ y ^ u ^ v == 0),
        int(x == u and y == v),
    )[symbol]


def word_polynomial(word, vertex, parameters, mode, boundary):
    shift = sum(vertex[symbol] for symbol in word)
    if mode == "auxprofile":
        value = [1, 0, 0, 1]
        for symbol in word:
            mask = parameters[symbol]
            factor = [
                int(i == j) ^ (3 if mask >> (2 * i + j) & 1 else 0)
                for i in range(2) for j in range(2)
            ]
            value = matrix_product(value, factor)
        coefficient = value[0] if boundary == 0 else value[1] ^ value[2] ^ value[3]
        return coefficient << shift
    value = 1
    for previous, symbol in zip(word, word[1:]):
        parameter = parameters[3 * previous + symbol]
        if mode == "polyprofile":
            value = poly_product(value, parameter)
        else:
            value <<= parameter
    return value << shift


def direct_tensors(t, vertex, parameters, mode, boundary=0):
    d = 1 << t
    coefficients = {}
    supports = [
        [quad for quad in product(range(2), repeat=4) if local(symbol, quad)]
        for symbol in range(3)
    ]
    for word in product(range(3), repeat=t):
        polynomial = word_polynomial(word, vertex, parameters, mode, boundary)
        assert polynomial.bit_count() & 1 == 1
        tensor = 0
        for quads in product(*(supports[symbol] for symbol in word)):
            x, y, u, v = (
                sum(quads[i][leg] << i for i in range(t)) for leg in range(4)
            )
            tensor ^= 1 << (((x * d + y) * d + u) * d + v)
        for degree in bits(polynomial):
            coefficients[degree] = coefficients.get(degree, 0) ^ tensor
    endpoint = xor_all(coefficients.values())
    crossing = sum(
        1 << (((x * d + y) * d + y) * d + x)
        for x, y in product(range(d), repeat=2)
    )
    assert endpoint == crossing
    return coefficients


def ranks_of_coefficients(t, tensors, length):
    d = 1 << t
    mask = (1 << (d * d)) - 1
    rows, columns = [], []
    for degree in range(length):
        tensor = tensors.get(degree, 0)
        row_matrix = [(tensor >> (d * d * i)) & mask for i in range(d * d)]
        column_matrix = [0] * (d * d)
        for index in bits(tensor):
            v = index % d
            index //= d
            u = index % d
            index //= d
            y = index % d
            x = index // d
            column_matrix[x * d + u] ^= 1 << (y * d + v)
        rows.append(rank(row_matrix))
        columns.append(rank(column_matrix))
    return rows, columns


def native_profile(executable, directory, t, vertex, parameters, mode, boundary=0):
    path = directory / "current.json"
    command = [str(executable), mode, str(t), str(path),
               *map(str, vertex + parameters)]
    if mode == "auxprofile":
        command.append(str(boundary))
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    source_dir = Path(__file__).resolve().parent
    sources = [source_dir / "markov_pairing_search.cpp",
               source_dir / "markov_pairing_profile.cpp"]
    rng = random.Random(20261030)
    cases = [
        (2, [0, 1, 2], [0] * 9, "profile", 0),
        (3, [0, 0, 2], [0] * 8 + [1], "profile", 0),
    ]
    for _ in range(10):
        cases.append((3, [rng.randrange(3) for _ in range(3)],
                      [rng.randrange(2) for _ in range(9)], "profile", 0))
    for _ in range(6):
        cases.append((3, [rng.randrange(3) for _ in range(3)],
                      [rng.choice((1, 2, 7)) for _ in range(9)], "polyprofile", 0))
    for boundary in (0, 1):
        for matrices in ([2, 9, 4], [3, 12, 7], [rng.randrange(16) for _ in range(3)]):
            cases.append((3, [0, 1, 2], list(matrices), "auxprofile", boundary))
    with TemporaryDirectory(prefix="pairing-transfer-") as temporary:
        directory = Path(temporary)
        binaries = [directory / "small", directory / "large"]
        for source, binary in zip(sources, binaries):
            subprocess.run(["g++", "-O2", "-std=c++17", str(source), "-o", str(binary)],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        subprocess.run([str(binaries[0]), "selftest", "4", str(directory / "kernels.json")],
                       check=True, stdout=subprocess.DEVNULL)
        kernel_result = json.loads((directory / "kernels.json").read_text())
        for t, vertex, parameters, mode, boundary in cases:
            native = native_profile(binaries[0], directory, t, vertex, parameters, mode, boundary)
            tensors = direct_tensors(t, vertex, parameters, mode, boundary)
            rows, columns = ranks_of_coefficients(t, tensors, len(native["row_ranks"]))
            assert rows == native["row_ranks"] and columns == native["column_ranks"]
        wide_cases = [
            (6, [0, 0, 2], [0] * 8 + [1]),
            (6, [0, 1, 2], [1, 0, 1, 0, 1, 0, 1, 1, 0]),
        ]
        for t, vertex, edge in wide_cases:
            small = native_profile(binaries[0], directory, t, vertex, edge, "profile")
            large = native_profile(binaries[1], directory, t, vertex, edge, "profile")
            assert small["row_ranks"] == large["row_ranks"]
            assert small["column_ranks"] == large["column_ranks"]
        # Exercises polynomial exponents above 63 in the wide implementation.
        t, vertex, edge = 2, [0, 0, 40], [0] * 8 + [20]
        large = native_profile(binaries[1], directory, t, vertex, edge, "profile")
        tensors = direct_tensors(t, vertex, edge, "profile")
        rows, columns = ranks_of_coefficients(t, tensors, len(large["row_ranks"]))
        assert rows == large["row_ranks"] and columns == large["column_ranks"]

    # A fixed coupled model is genuinely not position invariant.
    domain_edges = [0, 0, 1, 0, 0, 1, 1, 1, 0]
    tensors = direct_tensors(3, [0, 0, 5], domain_edges, "profile")
    d = 8
    first = ((0 * d + 6) * d + 6) * d
    second = ((0 * d + 5) * d + 5) * d
    first_polynomial = sum(1 << e for e, tensor in tensors.items() if tensor >> first & 1)
    second_polynomial = sum(1 << e for e, tensor in tensors.items() if tensor >> second & 1)
    assert first_polynomial == 1 << 6 and second_polynomial == 1 << 7

    result = {
        "status": "passed",
        "scope": "Exact finite-state crossing identities and rank engine; no improved factoring exponent or asymptotic rank certificate.",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        "independent_tensor_word_cases": len(cases),
        "exact_circulant_shortcut_cases": kernel_result["all_kernels"],
        "wide_vs_single_word_cases": len(wide_cases),
        "wide_polynomial_degree_case": 100,
        "position_symmetry_counterexample_degrees": [6, 7],
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


if __name__ == "__main__":
    main()
