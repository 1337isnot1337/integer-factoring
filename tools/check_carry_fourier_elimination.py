#!/usr/bin/env python3
"""Exact finite controls for the constant-state carry Fourier recurrence.

All arithmetic is in explicitly constructed binary finite fields. No
floating point roots, uncharged formal characters, or probabilistic
irreducibility tests are used.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from itertools import product
import json
from pathlib import Path
import random


def prime_factors(value: int) -> list[int]:
    result, candidate = [], 2
    while candidate * candidate <= value:
        if value % candidate == 0:
            result.append(candidate)
            while value % candidate == 0:
                value //= candidate
        candidate += 1
    if value > 1:
        result.append(value)
    return result


def polynomial_remainder(value: int, modulus: int) -> int:
    while value.bit_length() >= modulus.bit_length():
        value ^= modulus << (value.bit_length() - modulus.bit_length())
    return value


def polynomial_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, polynomial_remainder(a, b)
    return a


class BinaryField:
    def __init__(self, polynomial: int):
        self.polynomial = polynomial
        self.degree = polynomial.bit_length() - 1
        self.size = 1 << self.degree

    def mul(self, a: int, b: int) -> int:
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & self.size:
                a ^= self.polynomial
        return result

    def power(self, value: int, exponent: int) -> int:
        if exponent < 0:
            raise ValueError("nonnegative exponent required")
        result = 1
        while exponent:
            if exponent & 1:
                result = self.mul(result, value)
            value = self.mul(value, value)
            exponent >>= 1
        return result

    def inverse(self, value: int) -> int:
        if value == 0:
            raise ValueError("zero has no inverse")
        answer = self.power(value, self.size - 2)
        assert self.mul(value, answer) == 1
        return answer

    def product(self, *values: int) -> int:
        result = 1
        for value in values:
            result = self.mul(result, value)
        return result


def irreducible(polynomial: int) -> bool:
    field = BinaryField(polynomial)
    degree = field.degree
    powers = [2]
    for _ in range(degree):
        powers.append(field.mul(powers[-1], powers[-1]))
    return (powers[degree] == 2
            and all(polynomial_gcd(polynomial, powers[degree // p] ^ 2) == 1
                    for p in prime_factors(degree)))


def construct_field_and_root(m: int) -> tuple[BinaryField, int, dict]:
    degree, trials = 2 * m, 0
    for lower in range(1, 1 << degree, 2):
        polynomial = (1 << degree) | lower
        trials += 1
        if irreducible(polynomial):
            break
    else:
        raise AssertionError("no irreducible polynomial found")
    field = BinaryField(polynomial)
    modulus = (1 << m) + 1
    factors = prime_factors(modulus)
    for candidate in range(1, field.size):
        root = field.power(candidate, (1 << m) - 1)
        if all(field.power(root, modulus // p) != 1 for p in factors):
            break
    else:
        raise AssertionError("no primitive M-th root found")
    assert field.power(root, modulus) == 1
    assert all(field.power(root, modulus // p) != 1 for p in factors)
    return field, root, {
        "m": m, "M": modulus, "field_degree": degree,
        "irreducible_polynomial_binary": bin(polynomial),
        "irreducible_candidates_tested": trials,
        "prime_divisors_of_M": factors,
        "root_search_candidate": candidate,
        "primitive_root_binary": bin(root),
    }


def geometric(field: BinaryField, value: int, length: int) -> int:
    if length < 0:
        raise ValueError("negative geometric sum length")
    power, total = 1, 0
    step_power, step_total = value, 1
    while length:
        if length & 1:
            total ^= field.mul(power, step_total)
            power = field.mul(power, step_power)
        step_total ^= field.mul(step_power, step_total)
        step_power = field.mul(step_power, step_power)
        length >>= 1
    return total


def geometric_pair(field: BinaryField, z: int, beta: int, length: int) -> int:
    """Return sum_{r<length} z^r * sum_{a<r} beta^a without division."""
    def concatenate(left, right):
        pz, pb, sz, sb, joint = left
        qz, qb, tz, tb, other = right
        return (
            field.mul(pz, qz), field.mul(pb, qb),
            sz ^ field.mul(pz, tz), sb ^ field.mul(pb, tb),
            joint ^ field.mul(pz, field.mul(sb, tz) ^ field.mul(pb, other)),
        )
    answer, step = (1, 1, 0, 0, 0), (z, beta, 1, 1, 0)
    while length:
        if length & 1:
            answer = concatenate(answer, step)
        step = concatenate(step, step)
        length >>= 1
    return answer[-1]


@dataclass
class CarryState:
    zero: int = 1
    endpoint: int = 0
    u: int = 0
    v: int = 0
    tail: int = 0
    alpha: int = 1
    beta: int = 1


class CarryFourier:
    def __init__(self, digits: int, base: int):
        if digits < 2 or base <= digits or base & (base - 1):
            raise ValueError("requires D>=2 and a power-of-two base b>D")
        self.digits, self.base = digits, base
        self.bound = digits * (base - 1)
        self.m = (digits * base * base - 1).bit_length()
        self.modulus = (1 << self.m) + 1
        self.field, self.root, self.construction = construct_field_and_root(self.m)
        self.characters = [1]
        for _ in range(1, self.modulus):
            self.characters.append(self.field.mul(self.characters[-1], self.root))
        assert len(set(self.characters)) == self.modulus
        assert self.bound < self.modulus
        assert digits * base * base <= 1 << self.m

    def character(self, frequency: int, value: int) -> int:
        return self.characters[(frequency * value) % self.modulus]

    def state_value(self, state: CarryState, carry: int) -> int:
        f, d = self.field, self.digits
        if carry == 0:
            return state.zero
        if carry == d - 1:
            return state.endpoint
        if carry >= d:
            return f.mul(state.tail, f.power(state.alpha, carry))
        return f.mul(f.power(state.alpha, carry),
                     state.u ^ f.mul(state.v, geometric(f, state.beta, carry - 1)))

    def prefix(self, state: CarryState, frequency: int, upper: int) -> int:
        if upper < 0:
            return 0
        f, d = self.field, self.digits
        upper = min(upper, self.bound)
        z = f.mul(state.alpha, self.character(-frequency, 1))
        answer = state.zero
        core = min(upper, d - 2)
        if core >= 1:
            answer ^= f.product(state.u, z, geometric(f, z, core))
            answer ^= f.product(state.v, z, geometric_pair(f, z, state.beta, core))
        if upper >= d - 1:
            answer ^= f.mul(state.endpoint, self.character(-frequency, d - 1))
        if upper >= d:
            answer ^= f.product(state.tail, f.power(z, d),
                                geometric(f, z, upper - d + 1))
        return answer

    def step(self, state: CarryState, digit: int, frequency: int) -> CarryState:
        if not 0 <= digit < self.base:
            raise ValueError("invalid target digit")
        f, d, b, k = self.field, self.digits, self.base, self.bound
        z = f.mul(state.alpha, self.character(-frequency, 1))
        alpha = self.character(frequency, b)
        beta = f.power(z, b)
        phase = self.character(frequency, digit)
        small = self.prefix(state, frequency, d - 1)
        length = digit + b - d + 1
        assert length > 0
        return CarryState(
            zero=f.mul(phase, self.prefix(state, frequency, digit)),
            endpoint=f.product(phase, f.power(alpha, d - 1),
                               self.prefix(state, frequency, min(k, digit + b * (d - 1)))),
            u=f.mul(phase, small ^ f.product(state.tail, f.power(z, d), geometric(f, z, length))),
            v=f.product(phase, state.tail, f.power(z, d + length), geometric(f, z, b)),
            tail=f.mul(phase, self.prefix(state, frequency, k)),
            alpha=alpha, beta=beta,
        )

    def transform(self, digits: tuple[int, ...], frequencies: tuple[int, ...]) -> int:
        if len(digits) != len(frequencies):
            raise ValueError("frequency length differs from target length")
        state = CarryState()
        for digit, frequency in zip(digits, frequencies):
            state = self.step(state, digit, frequency)
        return state.zero

    def direct_step(self, vector: list[int], digit: int, frequency: int,
                    literal: bool) -> list[int]:
        f, k = self.field, self.bound
        if literal:
            result = []
            for outgoing in range(k + 1):
                value = 0
                for incoming, coefficient in enumerate(vector):
                    original = digit + self.base * outgoing - incoming
                    if 0 <= original < self.modulus:
                        value ^= f.mul(coefficient, self.character(frequency, original))
                result.append(value)
            return result
        # A full carry vector, independent of the seven-element state form.
        prefixes, value = [], 0
        for incoming, coefficient in enumerate(vector):
            value ^= f.mul(coefficient, self.character(-frequency, incoming))
            prefixes.append(value)
        return [f.mul(self.character(frequency, digit + self.base * outgoing),
                      prefixes[min(k, digit + self.base * outgoing)])
                for outgoing in range(k + 1)]


def field_rank(field: BinaryField, matrix: list[list[int]]) -> int:
    rows, rank = [row[:] for row in matrix], 0
    for column in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = field.inverse(rows[rank][column])
        rows[rank] = [field.mul(x, scale) for x in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][column]:
                multiple = rows[i][column]
                rows[i] = [x ^ field.mul(multiple, y) for x, y in zip(rows[i], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def geometric_checks(model: CarryFourier, totals: Counter) -> None:
    f = model.field
    values = (0, 1, model.root, f.power(model.root, 2), model.characters[-1])
    for z, beta, length in product(values, values, (0, 1, 2, 3, 7, 16, 31)):
        pz, pb, sz, sb, joint = 1, 1, 0, 0, 0
        for _ in range(length):
            joint ^= f.mul(pz, sb)
            sz ^= pz
            sb ^= pb
            pz, pb = f.mul(pz, z), f.mul(pb, beta)
        assert geometric(f, z, length) == sz
        assert geometric_pair(f, z, beta, length) == joint
        totals["geometric_pair_controls"] += 1


def state_checks(model: CarryFourier, totals: Counter, rng: random.Random) -> list[dict]:
    d, b, k, modulus = model.digits, model.base, model.bound, model.modulus
    boundary_digits = sorted({0, b - 1, b - d - 1, b - d, d - 2, d - 1})
    frequency_patterns = [
        (0,) * 6, (1,) * 6,
        tuple(pow(b, i, modulus) for i in range(6)),
        tuple(rng.randrange(modulus) for _ in range(6)),
    ]
    records = []
    for digits in [tuple([x] * 6) for x in boundary_digits] + [
            tuple(rng.randrange(b) for _ in range(6)) for _ in range(3)]:
        for frequencies in frequency_patterns:
            state = CarryState()
            vector = [1] + [0] * k
            for digit, frequency in zip(digits, frequencies):
                state = model.step(state, digit, frequency)
                literal = k <= 10
                vector = model.direct_step(vector, digit, frequency, literal)
                assert [model.state_value(state, carry) for carry in range(k + 1)] == vector
                for upper in sorted({0, d - 2, d - 1, d, k, rng.randrange(k + 1)}):
                    check_frequency = rng.randrange(modulus)
                    expected = 0
                    for r in range(min(upper, k) + 1):
                        expected ^= model.field.mul(vector[r], model.character(-check_frequency, r))
                    assert model.prefix(state, check_frequency, upper) == expected
                    totals["weighted_prefix_checks"] += 1
                totals["state_coordinates_checked"] += k + 1
                totals["literal_matrix_steps" if literal else "full_vector_steps"] += 1
            assert model.transform(digits, frequencies) == vector[0]
            totals["complete_frequency_queries"] += 1
            if len(records) < 4:
                records.append({"target_digits": list(digits), "frequencies": list(frequencies),
                                "answer_field_element": state.zero})
    return records


def rank_checks(model: CarryFourier, totals: Counter) -> dict:
    d, b, k, modulus = model.digits, model.base, model.bound, model.modulus
    records = []
    for digit in sorted({0, b - d - 1, b - d, b - 1}):
        for frequency in (0, 1):
            matrix = [[(model.character(frequency, digit + b * outgoing - incoming)
                        if incoming <= digit + b * outgoing else 0)
                       for outgoing in range(k + 1)] for incoming in range(k + 1)]
            expected = (k - digit + b - 1) // b + 1
            rank = field_rank(model.field, matrix)
            assert rank == expected <= d + 1
            records.append({"digit": digit, "frequency": frequency, "rank": rank})
            totals["fixed_frequency_matrix_ranks"] += 1
    matrix = [[model.character(frequency, b * k - incoming)
               for frequency in range(k + 1)] for incoming in range(k + 1)]
    assert field_rank(model.field, matrix) == k + 1
    totals["joint_frequency_full_rank_controls"] += 1
    return {"fixed_frequency_ranks": records, "joint_frequency_rank": k + 1}


def coefficient_and_dual_checks(model: CarryFourier, totals: Counter) -> dict:
    """Complete two-coefficient domains; no intermediate carry-state formula."""
    assert model.digits == 2 and model.base == 4
    b, k, modulus, f = model.base, model.bound, model.modulus, model.field
    selected = (2, 4, 9)
    result = []
    for target in selected:
        digits = (target % b, target // b)
        allowed = []
        for a0, a1 in product(range(modulus), repeat=2):
            first, d0 = divmod(a0, b)
            last, d1 = divmod(a1 + first, b)
            if first <= k and last == 0 and (d0, d1) == digits:
                allowed.append((a0, a1))
        for frequencies in ((0, 0), (1, 0), (0, 1), (1, 1), (3, 7), (modulus - 1, 2)):
            expected = 0
            for coefficients in allowed:
                expected ^= model.character(1, sum(s * a for s, a in zip(frequencies, coefficients)))
            assert model.transform(digits, frequencies) == expected
            totals["complete_coefficient_fourier_checks"] += 1

        # Weighted input polynomial: the coefficient at X^p must equal
        # the number modulo two of allowed q with p*q=target. Every p
        # here has one digit, so exponents remain explicitly distinguishable.
        polynomial = [0] * b
        for s0, s1 in product(range(modulus), repeat=2):
            carry = model.transform(digits, (s0, s1))
            for p in range(b):
                input_sum = 0
                for q in range(b):
                    input_sum ^= model.character(-s0, p * q)
                polynomial[p] ^= f.mul(carry, input_sum)
            totals["dual_frequency_tuples"] += 1
        expected = [sum(p * q == target for q in range(b)) % 2 for p in range(b)]
        assert polynomial == expected
        totals["complete_weighted_dual_identities"] += 1
        result.append({"target": target, "accepted_coefficient_tuples": len(allowed),
                       "weighted_factor_polynomial": polynomial})
    return {"models": result, "frequency_tuples_per_target": modulus ** 2}


def two_frequency_checks(model: CarryFourier, totals: Counter) -> list[dict]:
    records = []
    for start, length in ((0, 1), (0, model.base), (model.digits, model.bound - model.digits + 1),
                          (3, 2 * model.base)):
        interval = list(range(start, start + length))
        matrix = []
        for s in range(model.modulus):
            row = []
            for t in range(model.modulus):
                value = 0
                for r in interval:
                    value ^= model.character(model.base * s - t, r)
                row.append(value)
            matrix.append(row)
        assert field_rank(model.field, matrix) == length
        # Nonconsecutive frequency restrictions retain full rank when each
        # side has at least the interval length.
        rows = [(2 * i) % model.modulus for i in range(length)]
        columns = [(3 * i + 1) % model.modulus for i in range(length)]
        if len(set(columns)) != length:
            columns = [(2 * i + 1) % model.modulus for i in range(length)]
        assert len(set(rows)) == len(set(columns)) == length
        assert field_rank(model.field, [[matrix[s][t] for t in columns] for s in rows]) == length
        totals["two_frequency_interval_ranks"] += 1
        records.append({"interval_start": start, "interval_length": length,
                        "full_rank": length, "restricted_rows": rows,
                        "restricted_columns": columns, "restricted_rank": length})
    return records


def twisted_composition_checks(model: CarryFourier, totals: Counter) -> dict:
    d, b, k, modulus, f = model.digits, model.base, model.bound, model.modulus, model.field
    kernel = []
    for s in range(modulus):
        row = []
        for t in range(modulus):
            value = 0
            for r in range(k + 1):
                value ^= model.character(b * s - t, r)
            row.append(value)
        kernel.append(row)
    records = []
    for digit, coefficient in ((0, 0), (b - 1, 0), (0, b - 1), (b - 1, d * (b - 1) ** 2)):
        pairs = [(b * q + digit - coefficient, q) for q in range(k + 1)
                 if 0 <= b * q + digit - coefficient <= k]
        assert len(pairs) <= d
        composed = []
        for s in range(modulus):
            row = []
            for t in range(modulus):
                value = 0
                for u in range(modulus):
                    value ^= f.product(kernel[s][u], model.character(digit - coefficient, u), kernel[u][t])
                expected = 0
                for r, q in pairs:
                    expected ^= model.character(1, b * s * r - t * q)
                assert value == expected
                row.append(value)
            composed.append(row)
        assert field_rank(f, composed) == len(pairs)
        totals["twisted_composition_identities"] += 1
        records.append({"target_digit": digit, "coefficient": coefficient,
                        "surviving_integer_pairs": pairs, "rank": len(pairs)})
    # Keep the coefficient physical instead of fixing it: all carry
    # characters return to the row span.
    physical = []
    for s in range(modulus):
        row = []
        for coefficient in range(b):
            for t in range(modulus):
                value = 0
                for q in range(k + 1):
                    r = b * q - coefficient
                    if 0 <= r <= k:
                        value ^= model.character(1, b * s * r - t * q)
                row.append(value)
        physical.append(row)
    assert field_rank(f, physical) == k + 1
    totals["unresolved_coefficient_full_rank_controls"] += 1
    return {"fixed_coefficient_cases": records, "unresolved_coefficient_rank": k + 1}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    totals, records = Counter(), []
    rng = random.Random(20261101)
    smallest = None
    for digits, base in ((2, 4), (3, 4), (4, 8), (8, 16)):
        model = CarryFourier(digits, base)
        if smallest is None:
            smallest = model
        geometric_checks(model, totals)
        queries = state_checks(model, totals, rng)
        ranks = rank_checks(model, totals) if digits <= 4 else None
        records.append({"D": digits, "b": base, "K": model.bound,
                        "field_construction": model.construction,
                        "sample_queries": queries, "rank_controls": ranks})
        totals["explicit_field_root_constructions"] += 1
    report = {
        "result": "constant-state exact carry Fourier evaluation",
        "scope": "one specified frequency vector; no uniform small tensor bond or improved factoring exponent",
        "models": records,
        "complete_dual_checks": coefficient_and_dual_checks(smallest, totals),
        "two_frequency_kernel_checks": two_frequency_checks(smallest, totals),
        "twisted_composition_checks": twisted_composition_checks(smallest, totals),
        "summary": dict(sorted(totals.items())),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
