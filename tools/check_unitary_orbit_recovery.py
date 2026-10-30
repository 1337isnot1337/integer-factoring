#!/usr/bin/env python3
"""Exact finite checks for unitary orbit factoring recovery.

Ordinary factor tuples are enumerated independently of the orbit-count
formula. The script is a transcription check, not a fast factoring code.
Only the Python standard library is required.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import hashlib
from itertools import combinations, product
import json
from math import gcd, isqrt
from pathlib import Path


def factor_reference(value: int) -> list[int]:
    result: list[int] = []
    divisor = 2
    while divisor * divisor <= value:
        while value % divisor == 0:
            result.append(divisor)
            value //= divisor
        divisor += 1
    if value > 1:
        result.append(value)
    return result


@lru_cache(maxsize=None)
def divisors(value: int) -> tuple[int, ...]:
    low: list[int] = []
    high: list[int] = []
    for divisor in range(1, isqrt(value) + 1):
        if value % divisor == 0:
            low.append(divisor)
            if divisor * divisor != value:
                high.append(value // divisor)
    return tuple(low + high[::-1])


def ordinary_tuples(value: int, length: int):
    if length == 1:
        yield (value,)
        return
    for first in divisors(value):
        for tail in ordinary_tuples(value // first, length - 1):
            yield (first, *tail)


def rotations(values: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(values[shift:] + values[:shift] for shift in range(len(values)))


def pairwise_coprime(values: tuple[int, ...]) -> bool:
    return all(gcd(left, right) == 1 for left, right in combinations(values, 2))


def inverse_witness(a: int, b: int) -> tuple[int, int] | None:
    if gcd(a, b) != 1:
        return None
    if b == 1:
        return 0, 0
    u = pow(a, -1, b)
    v = (a * u - 1) // b
    assert 0 <= u < b and 0 <= v < a
    assert a * u == b * v + 1
    return u, v


def add_one_trace(value: int, delta: int, width: int) -> tuple[int, int, tuple[int, ...]]:
    assert delta in (0, 1)
    carry = delta
    result = 0
    trace = [carry]
    for position in range(width):
        total = ((value >> position) & 1) + carry
        result |= (total & 1) << position
        carry = total >> 1
        trace.append(carry)
    return result, carry, tuple(trace)


def compare_trace(left: int, right: int, width: int) -> tuple[int, tuple[int, ...]]:
    state = 0
    trace = [state]
    for position in range(width - 1, -1, -1):
        a = (left >> position) & 1
        b = (right >> position) & 1
        if state == 0:
            state = (a > b) - (a < b)
        trace.append(state)
    return state, tuple(trace)


def verify_prefix_inverse_network(values: tuple[int, ...]) -> tuple[bool, int]:
    total = 1
    for value in values:
        total *= value
    width = total.bit_length()
    prefix = values[0]
    equations = 0
    for value in values[1:]:
        witness = inverse_witness(prefix, value)
        if witness is None:
            return False, equations
        u, v = witness
        delta = int(value > 1)
        assert 0 <= u < value and 0 <= v < prefix
        au, bv = prefix * u, value * v
        product_value = prefix * value
        assert au < product_value <= total and bv < product_value <= total
        recovered, carry, trace = add_one_trace(bv, delta, width)
        assert recovered == au and carry == 0 and len(trace) == width + 1
        prefix = product_value
        equations += 1
    assert prefix == total
    return True, equations


def power_capped(base: int, exponent: int, cap: int) -> int:
    result = 1
    for _ in range(exponent):
        result *= base
        if result > cap:
            return cap + 1
    return result


def integer_root(value: int, exponent: int) -> int:
    low = 1
    high = 1 << ((value.bit_length() + exponent - 1) // exponent)
    while low < high:
        mid = (low + high + 1) // 2
        if power_capped(mid, exponent, value) <= value:
            low = mid
        else:
            high = mid - 1
    return low


@lru_cache(maxsize=None)
def perfect_power(value: int) -> tuple[int, int] | None:
    for exponent in range(2, value.bit_length()):
        base = integer_root(value, exponent)
        if base >= 2 and power_capped(base, exponent, value) == value:
            return base, exponent
    return None


def pack(values: tuple[int, ...], width: int) -> int:
    word = 0
    for value in values:
        assert 0 <= value < 1 << width
        word = (word << width) | value
    return word


def unpack(word: int, width: int, length: int) -> tuple[int, ...]:
    mask = (1 << width) - 1
    return tuple((word >> ((length - index - 1) * width)) & mask for index in range(length))


def residue(value: int, characteristic: int) -> int:
    return value if characteristic == 0 else value % characteristic


def recover_nonzero(
    witnesses: tuple[tuple[int, ...], ...], width: int, length: int, characteristic: int
) -> tuple[tuple[int, ...], int]:
    encoded = tuple(pack(values, width) for values in witnesses)
    assert len(set(encoded)) == len(encoded)
    assert residue(len(encoded), characteristic) != 0
    prefix = 0
    free = width * length
    calls = 0
    while free:
        free -= 1
        left_prefix = prefix * 2
        left = residue(sum((word >> free) == left_prefix for word in encoded), characteristic)
        calls += 1
        prefix = left_prefix if left else left_prefix + 1
        child = residue(sum((word >> free) == prefix for word in encoded), characteristic)
        assert child != 0
    result = unpack(prefix, width, length)
    assert result in witnesses
    return result, calls


def factor_by_counts(value: int, characteristic: int, representatives: dict[int, tuple[tuple[int, ...], ...]]):
    length = max(2, characteristic)
    pending = [(value, 1)]
    output: list[int] = []
    queries = 0
    power_nodes = 0
    while pending:
        current, multiplicity = pending.pop()
        root = perfect_power(current)
        if root is not None:
            base, exponent = root
            assert 1 < base < current and base**exponent == current
            pending.append((base, multiplicity * exponent))
            power_nodes += 1
            continue
        witnesses = representatives[current]
        queries += 1
        if residue(len(witnesses), characteristic) == 0:
            output.extend([current] * multiplicity)
            continue
        values, used = recover_nonzero(witnesses, current.bit_length(), length, characteristic)
        queries += used
        assert sum(component > 1 for component in values) >= 2
        divisor = next(component for component in values if component > 1)
        assert 1 < divisor < current and current % divisor == 0
        pending.extend(((divisor, multiplicity), (current // divisor, multiplicity)))
    return sorted(output), queries, power_nodes


def check_characteristic(characteristic: int, limit: int) -> dict[str, int]:
    length = max(2, characteristic)
    representatives: dict[int, tuple[tuple[int, ...], ...]] = {}
    stats: Counter[str] = Counter()
    for value in range(2, limit + 1):
        exponents = Counter(factor_reference(value))
        omega = len(exponents)
        unitary: list[tuple[int, ...]] = []
        canonical: set[tuple[int, ...]] = set()
        for values in ordinary_tuples(value, length):
            stats["ordinary_tuples"] += 1
            coprime = pairwise_coprime(values)
            valid, equations = verify_prefix_inverse_network(values)
            assert valid == coprime
            if not coprime:
                continue
            stats["inverse_equations"] += equations
            unitary.append(values)
            rotated = rotations(values)
            assert len(set(rotated)) == length
            canonical.add(min(rotated))
        assert len(unitary) == length**omega
        assert len(canonical) * length == len(unitary)
        trivial = (1,) * (length - 1) + (value,)
        assert trivial in canonical
        admitted = tuple(sorted(canonical - {trivial}))
        representatives[value] = admitted
        assert len(admitted) == length ** (omega - 1) - 1
        expected = 0 if omega == 1 else (-1 if characteristic else len(admitted))
        assert residue(len(admitted), characteristic) == residue(expected, characteristic)
        for values in canonical:
            word = pack(values, value.bit_length())
            for other in rotations(values)[1:]:
                state, _ = compare_trace(word, pack(other, value.bit_length()), length * value.bit_length())
                assert state == -1
        root = perfect_power(value)
        exponent_gcd = 0
        for exponent in exponents.values():
            exponent_gcd = gcd(exponent_gcd, exponent)
        assert (root is not None) == (exponent_gcd > 1)
        if root is not None:
            stats["perfect_power_inputs"] += 1
        elif residue(len(admitted), characteristic) == 0:
            assert factor_reference(value) == [value]
            stats["prime_inputs"] += 1
        else:
            stats["unitary_orbit_inputs"] += 1
        result, queries, power_nodes = factor_by_counts(value, characteristic, representatives)
        assert result == factor_reference(value)
        stats["maximum_factorization_queries"] = max(stats["maximum_factorization_queries"], queries)
        stats["maximum_power_nodes"] = max(stats["maximum_power_nodes"], power_nodes)
        stats["unitary_tuples"] += len(unitary)
        stats["canonical_representatives"] += len(canonical)
    return {"characteristic": characteristic, "tuple_length": length, "integers": limit - 1, **dict(stats)}


def check_local_witnesses() -> dict[str, int]:
    inverse_pairs = 0
    unit_pairs = 0
    for a, b in product(range(1, 49), repeat=2):
        delta = int(b > 1)
        literal = [(u, v) for u in range(b) for v in range(a) if a * u == b * v + delta]
        expected = inverse_witness(a, b)
        assert literal == ([] if expected is None else [expected])
        inverse_pairs += 1
        unit_pairs += int(a == 1 or b == 1)
    additions = 0
    overflows = 0
    for width in range(1, 10):
        for value, delta in product(range(1 << width), (0, 1)):
            output, carry, trace = add_one_trace(value, delta, width)
            assert output + (carry << width) == value + delta
            assert len(trace) == width + 1 and trace[0] == delta
            additions += 1
            overflows += carry
    comparisons = 0
    for width in range(1, 8):
        for a, b in product(range(1 << width), repeat=2):
            state, trace = compare_trace(a, b, width)
            assert state == (a > b) - (a < b)
            assert len(trace) == width + 1 and trace[0] == 0
            comparisons += 1
    return {"literal_inverse_pairs": inverse_pairs, "unit_argument_pairs": unit_pairs,
            "addition_traces": additions, "overflow_negative_cases": overflows,
            "comparison_traces": comparisons}


def check_negative_controls() -> dict[str, object]:
    ordinary12 = [values for values in ordinary_tuples(12, 2) if 1 < values[0] < values[1]]
    unitary12 = [values for values in ordinary12 if pairwise_coprime(values)]
    assert ordinary12 == [(2, 6), (3, 4)] and unitary12 == [(3, 4)]
    assert len(ordinary12) % 2 == 0 and len(unitary12) % 2 == 1
    prime_power8 = [values for values in ordinary_tuples(8, 2) if 1 < values[0] < values[1] and pairwise_coprime(values)]
    assert prime_power8 == [] and perfect_power(8) is not None
    assert 2 * 2 == 3 * 1 + 1 and 2 * 5 == 3 * 3 + 1
    assert inverse_witness(5, 1) == (0, 0)
    assert 5 * 0 != 1 * 0 + 1
    assert not pairwise_coprime((2, 3, 2, 3))
    assert len(set(rotations((2, 3, 2, 3)))) == 2
    assert (2 ** (3 - 1) - 1) % 3 == 0
    assert (3 ** (3 - 1) - 1) % 3 == 2
    return {
        "omitting_coprimality": {"N": 12, "ordinary_proper_pairs": ordinary12, "unitary_pairs": unitary12},
        "omitting_perfect_power_test": {"N": 8, "unitary_nontrivial_count": 0},
        "omitting_inverse_bound": {"a": 2, "b": 3, "distinct_witnesses": [[2, 1], [5, 3]]},
        "omitting_unit_guard": {"a": 5, "b": 1, "valid_witness": [0, 0]},
        "omitting_coprimality_with_composite_length": {"tuple": [2, 3, 2, 3], "orbit_size": 2},
        "wrong_tuple_length_for_characteristic": {"N": 30, "characteristic": 3, "length2_residue": 0, "length3_residue": 2},
        "trivial_orbit_must_be_removed": {"N": 7, "canonical_unitary_orbits": 1},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=384)
    parser.add_argument("--binary-limit", type=int, default=4096)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    reports = []
    for characteristic in (0, 2, 3, 4, 5, 6, 7, 8, 9):
        limit = args.binary_limit if characteristic in (0, 2) else args.limit
        report = check_characteristic(characteristic, limit)
        reports.append(report)
        print(json.dumps(report, sort_keys=True), flush=True)
    report = {
        "status": "PASS", "scope": "Exact finite transcription checks, not an asymptotic runtime measurement",
        "characteristics": reports, "local_witnesses": check_local_witnesses(),
        "negative_controls": check_negative_controls(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    if args.output is not None:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(f"PASS: wrote {args.output}")
    else:
        print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
