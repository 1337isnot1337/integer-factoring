#!/usr/bin/env python3
"""Test the finite source-shaped exchange system for weight-two Witt fields."""

import argparse
import json
from fractions import Fraction
from time import monotonic

from witt_artin_hasse import monomials, artin_hasse_coefficients, hall_gram, energy
from check_modular_crossing import inverse_two


def insert(basis, value, witness=0):
    """Insert a binary vector while recording its input combination."""
    while value:
        pivot = value.bit_length() - 1
        if pivot not in basis:
            basis[pivot] = (value, witness)
            return True
        vector, mask = basis[pivot]
        value ^= vector
        witness ^= mask
    return False


def represent(basis, value):
    """Return a binary combination for value, or None outside the span."""
    witness = 0
    while value:
        pivot = value.bit_length() - 1
        if pivot not in basis:
            return None
        vector, mask = basis[pivot]
        value ^= vector
        witness ^= mask
    return witness


def transpose(rows, n):
    result = [0] * n
    for i, row in enumerate(rows):
        while row:
            bit = row & -row
            result[bit.bit_length() - 1] |= 1 << i
            row ^= bit
    return result


def multiply(A, B):
    result = []
    for row in A:
        value = 0
        while row:
            bit = row & -row
            value ^= B[bit.bit_length() - 1]
            row ^= bit
        result.append(value)
    return result


def calculate(H, neutral=False):
    started = monotonic()
    variables = H.bit_length()
    states = [m for e in range(H + 1) for m in monomials(e, variables)]
    positions = {m: i for i, m in enumerate(states)}
    grades = [energy(m) for m in states]
    grade_indices = [[i for i, e in enumerate(grades) if e == k] for k in range(H + 1)]
    size = len(states)
    G, GI = [0] * size, [0] * size
    for e in range(H + 1):
        if e == 0:
            small, small_inverse = [1], [1]
        else:
            _, _, small = hall_gram(e)
            small_inverse = inverse_two(small, len(small))
        offset = grade_indices[e][0]
        for j, (a, b) in enumerate(zip(small, small_inverse)):
            G[offset + j], GI[offset + j] = a << offset, b << offset
    parity = [x.numerator & 1 for x in artin_hasse_coefficients(H)]
    coherent = {m for m in states if all(parity[a] for a in m)}

    def creation(polynomial):
        answer = [0] * size
        for j, state in enumerate(states):
            for monomial in polynomial:
                target = tuple(a + b for a, b in zip(state, monomial))
                if energy(target) <= H:
                    answer[positions[target]] ^= 1 << j
        return answer

    def coherent_pair(polynomial):
        plus = creation(polynomial)
        minus = multiply(multiply(GI, transpose(plus, size)), G)
        return multiply(plus, minus)

    ordinary = coherent_pair(coherent)
    squared = coherent_pair({tuple(2 * a for a in m) for m in coherent if 2 * energy(m) <= H})

    def derivative_mask(matrix):
        result = [0] * size
        for i, row in enumerate(matrix):
            while row:
                bit = row & -row
                j = bit.bit_length() - 1
                if ((grades[i] - grades[j]) >> 1) & 1:
                    result[i] |= bit
                row ^= bit
        return result

    Y11 = derivative_mask(squared)
    Y2 = derivative_mask(ordinary)
    if neutral:
        inverse_creation = inverse_two(creation(coherent), size)
        inverse_creation_squared = multiply(inverse_creation, inverse_creation)
        adjoint_inverse = multiply(multiply(GI, transpose(inverse_creation, size)), G)
        adjoint_inverse_squared = multiply(adjoint_inverse, adjoint_inverse)
        Y11 = multiply(multiply(inverse_creation_squared, Y11), adjoint_inverse_squared)
        Y2 = multiply(multiply(inverse_creation, Y2), adjoint_inverse)
    fields = [Y11, [a ^ b for a, b in zip(Y11, Y2)]]
    labels = grade_indices[2]
    assert len(labels) == 2
    for v in range(2):
        for u in range(2):
            assert ((fields[v][labels[u]] >> 0) & 1) == (v == u)
    gram = [sum(((fields[x][0] >> labels[v]) & 1) << v for v in range(2)) for x in range(2)]
    inverse = inverse_two(gram, 2)
    duals = []
    for mask in inverse:
        duals.append([sum(0 for _ in ())] * size)
        for x in range(2):
            if mask >> x & 1:
                duals[-1] = [a ^ b for a, b in zip(duals[-1], fields[x])]
    left, right = [], []
    for e in range(H + 1):
        L = R = 0
        for x in range(2):
            for y in range(2):
                for v in range(2):
                    for u in range(2):
                        a = b = 0
                        for k in grade_indices[e]:
                            a ^= ((duals[x][k] >> labels[y]) & 1) * ((fields[v][labels[u]] >> k) & 1)
                            b ^= ((duals[x][labels[u]] >> k) & 1) * ((fields[v][k] >> labels[y]) & 1)
                        index = (2 * x + y) * 4 + 2 * v + u
                        L |= a << index
                        R |= b << index
        left.append(L)
        right.append(R)
    target_crossing = sum(1 << (5 * j) for j in range(4))
    systems = []
    for degree in range(H - 3):
        shifts = range(-degree, 5)
        columns = []
        for i in range(degree + 1):
            constraint, total = 0, 0
            for position, shift in enumerate(shifts):
                e, f = shift + i, 4 - shift - i
                a = left[e] if 0 <= e <= H else 0
                b = right[f] if 0 <= f <= H else 0
                constraint |= (a ^ b) << (16 * position)
                total ^= a
            constraint |= total << (16 * len(shifts))
            columns.append(constraint)
        basis = {}
        for i, column in enumerate(columns):
            insert(basis, column, 1 << i)
        target = target_crossing << (16 * len(shifts))
        hit = represent(basis, target)
        systems.append({"angle_degree": degree, "rank": len(basis),
                        "solution": None if hit is None else hex(hit)})
    return {"h": 2, "neutral_ratio_fields": neutral, "energy_window": H, "total_states": size,
            "vacuum_gram_rows": gram, "left_energy_tensors": left,
            "right_energy_tensors": right, "finite_exchange_systems": systems,
            "seconds": round(monotonic() - started, 3)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--window", type=int, default=12)
    parser.add_argument("--output")
    parser.add_argument("--neutral",action="store_true")
    args = parser.parse_args()
    result = calculate(args.window,args.neutral)
    if args.output:
        with open(args.output, "w") as f:
            json.dump(result, f, indent=2)
            f.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
