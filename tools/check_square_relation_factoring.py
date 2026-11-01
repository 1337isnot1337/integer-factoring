#!/usr/bin/env python3
"""Replay finite factoring certificates independently of the box search.

The replay recovers signed words from the integer coordinates, checks
every dominant row, computes the determinant with rational elimination
instead of the construction's Bareiss routine, verifies annihilation,
and checks the extracted factors.  It never reruns the MITM search.
"""

from fractions import Fraction
from math import gcd, isqrt, prod
from pathlib import Path
import argparse
import json

from early_stop_class_factoring import ambiguous_factorization, factor_from_annihilated_generators
from quadratic_form_arithmetic import (
    compose_primitive, inverse_primitive, power_primitive, reduce_primitive,
)
from universal_smooth_square import word_from_square_representation


def rational_determinant(matrix):
    data = [[Fraction(value) for value in row] for row in matrix]
    result = Fraction(1)
    for column in range(len(data)):
        pivot = next((index for index in range(column, len(data)) if data[index][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            data[column], data[pivot] = data[pivot], data[column]
            result = -result
        leading = data[column][column]
        result *= leading
        for index in range(column+1, len(data)):
            if not data[index][column]:
                continue
            multiple = data[index][column]/leading
            for position in range(column+1, len(data)):
                data[index][position] -= multiple*data[column][position]
            data[index][column] = 0
    assert result.denominator == 1
    return result.numerator


def trial_prime(value):
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    return all(value % divisor for divisor in range(3, isqrt(value)+1, 2))


def check(path):
    record = json.loads(path.read_text())
    odd_n = record["odd_n"]
    discriminant = -4*odd_n
    primes = record["factor_base_primes"]
    roots = {int(prime): root for prime, root in record["roots"].items()}
    cap = record["cap_fixed_before_targets"]
    count, identity = len(primes), (1, 0, odd_n)
    half_diagonal = cap*count+1
    assert odd_n > 1 and odd_n % 2 == 1 and cap > 0 and count > 0
    assert record["half_diagonal"] == half_diagonal
    assert len(record["targets"]) == count
    assert record["status"] == "factored" and record["missing_relation_indices"] == []
    assert not record["class_number_supplied"] and not record["annihilator_supplied"]
    assert len(set(primes)) == count and all(trial_prime(prime) for prime in primes)
    assert all(gcd(prime, discriminant) == 1 for prime in primes)
    generators = []
    for prime in primes:
        root = roots[prime]
        assert (root*root-discriminant) % (4*prime) == 0
        generators.append(reduce_primitive((prime, root, (root*root-discriminant)//(4*prime))))
    root_norm = prod(prime**cap for prime in primes)
    rows = []
    direct_factors = []
    for index, entry in enumerate(record["targets"]):
        assert entry["represented"] and entry["index"] == index
        target = tuple(entry["target"])
        assert target == power_primitive(generators[index], 2*half_diagonal)
        witness = entry["representation"]
        x, y = witness["vector"]
        assert target[0]*x*x+target[1]*x*y+target[2]*y*y == root_norm*root_norm
        assert witness["root_norm"] == root_norm and witness["fixed_norm"] == root_norm*root_norm
        content = gcd(x, y)
        assert content == witness["content"] and root_norm % content == 0
        assert list(witness["primitive_vector"]) == [x//content, y//content]
        assert witness["primitive_norm"] == (root_norm//content)**2
        recovered = word_from_square_representation(target, x, y, primes, cap)
        word = recovered["half_word"]
        assert word == entry["half_word"]
        row = [-2*value for value in word]
        row[index] += 2*half_diagonal
        assert row == entry["relation_row"]
        assert row[index] > sum(abs(value) for j, value in enumerate(row) if j != index)
        value, half_value = identity, identity
        for generator, exponent, half_exponent in zip(generators, row, word):
            value = compose_primitive(value, power_primitive(generator, exponent))
            half_value = compose_primitive(half_value, power_primitive(generator, half_exponent))
        assert value == identity
        difference = compose_primitive(power_primitive(generators[index], half_diagonal),
                                       inverse_primitive(half_value))
        assert difference == tuple(entry["root_difference"])
        factor = ambiguous_factorization(difference, odd_n)["proper_factor"]
        assert factor == entry["direct_factor_from_root_difference"]
        if factor is not None:
            direct_factors.append(factor)
        rows.append(row)
    assert rows == record["relation_matrix"]
    annihilator = abs(rational_determinant(rows))
    assert annihilator == record["annihilator"] > 0
    assert all(power_primitive(generator, annihilator) == identity for generator in generators)
    extraction = factor_from_annihilated_generators(odd_n, generators, annihilator, with_trace=True)
    assert json.loads(json.dumps(extraction)) == record["extraction"]
    factor = extraction["factor"]
    assert factor and 1 < factor < odd_n and odd_n % factor == 0
    factors = sorted((factor, odd_n//factor))
    assert factors == record["factorization"]
    # This finite replay's selected factors fit below 2^31.  Trial division
    # certifies their primality only after the search has already returned them.
    assert max(factors) < 1 << 31
    assert all(trial_prime(value) for value in factors)
    return {"certificate": path.name, "odd_n": odd_n, "input_bits": odd_n.bit_length(),
            "dominant_rows_checked": count, "fixed_norm_bits": (root_norm*root_norm).bit_length(),
            "rational_elimination_determinant": annihilator,
            "annihilator_bits": annihilator.bit_length(),
            "factorization": factors, "returned_factors_checked_prime": True,
            "direct_root_difference_factors": direct_factors,
            "complete_generation_claimed": False,
            "mitm_search_rerun": False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("certificates", nargs="*")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    paths = ([Path(path) for path in args.certificates] if args.certificates else
             [root/"evidence/square-relation-factoring-51bit.json",
              root/"evidence/square-relation-factoring-60bit.json"])
    print(json.dumps({"scope": "Independent arithmetic replay of finite square-relation factoring certificates",
                      "certificates": [check(path) for path in paths]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
