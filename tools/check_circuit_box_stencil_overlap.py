#!/usr/bin/env python3
"""Independent proper-root, prime-power, and literal stencil-gcd controls."""

import argparse
from hashlib import sha256
import json
from math import gcd
from pathlib import Path

from binary_quadratic_forms import reduce_form, reduced_forms
from quadratic_ideal_hnf import IdealBasis, extended_gcd, form_to_ideal, multiply_ideals


ROOT = Path(__file__).resolve().parents[1]
FOUR = [(1, 0), (0, 1), (1, 1), (1, -1)]


def _det(first, second):
    return first[0] * second[1] - first[1] * second[0]


def _norm(form, vector):
    a, b, c = form
    x, y = vector
    return a * x * x + b * x * y + c * y * y


def _polar(form, first, second):
    a, b, c = form
    x, y = first
    u, v = second
    return 2 * a * x * u + b * (x * v + y * u) + 2 * c * y * v


def _completion(vector):
    common, first, second = extended_gcd(*vector)
    assert common == 1
    result = (-second, first)
    assert _det(vector, result) == 1
    return result


def _key(form, vector):
    complement = _completion(vector)
    norm = _norm(form, vector)
    middle = _polar(form, vector, complement)
    trailing = _norm(form, complement)
    discriminant = form[1] ** 2 - 4 * form[0] * form[2]
    assert middle * middle - 4 * norm * trailing == discriminant
    assert reduce_form((norm, middle, trailing)) == reduce_form(form)
    root = middle % (2 * norm)
    canonical = min(root, (-root) % (2 * norm))
    sign = 1 if root == canonical else -1
    return {"norm": norm, "root": root, "middle": middle, "canonical": canonical, "sign": sign, "complement": complement}


def _ideal(discriminant, norm, root):
    root %= 2 * norm
    assert (root * root - discriminant) % (4 * norm) == 0
    return form_to_ideal((norm, root, (root * root - discriminant) // (4 * norm)))


def _remove_support(value, support):
    while True:
        common = gcd(value, support)
        if common == 1:
            return value
        value //= common


def _factor_small(value):
    result = {}
    prime = 2
    while prime * prime <= value:
        while value % prime == 0:
            result[prime] = result.get(prime, 0) + 1
            value //= prime
        prime += 1
    if value > 1:
        result[value] = result.get(value, 0) + 1
    return result


def _pair_identity(form, first_vector, second_vector, first, second, common):
    determinant = _det(first_vector, second_vector)
    alpha = _det(second_vector, first["complement"])
    beta = _det(first_vector, second["complement"])
    assert determinant * (first["middle"] + second["middle"]) == 2 * beta * second["norm"] - 2 * alpha * first["norm"]
    assert common > 0 and common % 2 == 1 and gcd(common, determinant) == 1
    assert first["norm"] % common == 0 and second["norm"] % common == 0
    assert (first["root"] + second["root"]) % (2 * common) == 0
    assert (first["sign"] * first["canonical"] + second["sign"] * second["canonical"]) % (2 * common) == 0
    discriminant = form[1] ** 2 - 4 * form[0] * form[2]
    assert gcd(common, discriminant) == 1
    product = multiply_ideals(_ideal(discriminant, common, first["root"]), _ideal(discriminant, common, second["root"]), discriminant)
    assert product == IdealBasis(common, 0, common)


def _graph_control(form, vectors, *, remove_determinants):
    discriminant = form[1] ** 2 - 4 * form[0] * form[2]
    keys = [_key(form, vector) for vector in vectors]
    determinant_product = 1
    for i, vector in enumerate(vectors):
        for other in vectors[i + 1:]:
            determinant = _det(vector, other)
            assert determinant
            determinant_product *= abs(determinant)
    support = abs(2 * discriminant) * (determinant_product if remove_determinants else 1)
    residuals = [_remove_support(key["norm"], support) for key in keys]
    for i, residual in enumerate(residuals):
        assert residual & 1 and gcd(residual, discriminant) == 1
        for j in range(i):
            assert gcd(gcd(residual, residuals[j]), _det(vectors[i], vectors[j])) == 1
    occurrences = {}
    for index, residual in enumerate(residuals):
        for prime in _factor_small(residual):
            occurrences.setdefault(prime, []).append(index)
    assert all(len(indices) <= 2 for indices in occurrences.values())
    incident = [1] * len(keys)
    edges, scalar = [], 1
    for i in range(len(keys)):
        for j in range(i):
            common = gcd(residuals[i], residuals[j])
            if common == 1:
                continue
            _pair_identity(form, vectors[i], vectors[j], keys[i], keys[j], common)
            edges.append((j, i, common))
            incident[i] *= common
            incident[j] *= common
            scalar *= common
    for index, (_, _, common) in enumerate(edges):
        assert all(gcd(common, previous[2]) == 1 for previous in edges[:index])
    leftovers = []
    for residual, divisor in zip(residuals, incident):
        quotient, remainder = divmod(residual, divisor)
        assert remainder == 0
        leftovers.append(quotient)
    for index, value in enumerate(leftovers):
        assert all(gcd(value, previous) == 1 for previous in leftovers[:index])
    left = right = IdealBasis(1, 0, 1)
    norm_left = norm_right = 1
    for key, residual, leftover in zip(keys, residuals, leftovers):
        left = multiply_ideals(left, _ideal(discriminant, residual, key["root"]), discriminant)
        right = multiply_ideals(right, _ideal(discriminant, leftover, key["root"]), discriminant)
        norm_left *= residual
        norm_right *= leftover
    assert norm_left == scalar * scalar * norm_right
    assert right.d == 1
    assert left == IdealBasis(scalar * right.a, scalar * right.r, scalar * right.d)
    assert left.d == scalar
    return {"form": form, "vertices": len(vectors), "residuals": residuals, "edges": edges, "leftovers": leftovers, "integer_content": scalar, "prime_occurrences": len(occurrences)}


def _prime_power_pair(prime, exponent, other_exponent):
    root = next(value for value in range(prime) if (value * value + 4) % prime == 0)
    modulus = prime
    for _ in range(1, exponent):
        root += ((-(root * root + 4) // modulus) * pow(2 * root, -1, prime) % prime) * modulus
        modulus *= prime
    if root & 1:
        root += modulus
    form = (modulus, root, (root * root + 4) // (4 * modulus))
    a, b, c = form
    coordinate = -c * pow(b, -1, prime) % prime
    power = prime
    for _ in range(1, other_exponent):
        value = a * coordinate * coordinate + b * coordinate + c
        coordinate += ((-value // power) * pow(2 * a * coordinate + b, -1, prime) % prime) * power
        power *= prime
    for adjustment in range(2 * prime):
        candidate = coordinate + adjustment * power
        norm = a * candidate * candidate + b * candidate + c
        if norm & 1 and norm // power % prime:
            coordinate = candidate
            break
    else:
        raise AssertionError("failed to choose an exact odd prime-power valuation")
    assert norm % power == 0 and norm // power % prime
    first_vector, second_vector = (1, 0), (coordinate, 1)
    first, second = _key(form, first_vector), _key(form, second_vector)
    common = gcd(first["norm"], second["norm"])
    assert common == prime ** min(exponent, other_exponent)
    _pair_identity(form, first_vector, second_vector, first, second, common)
    # If determinant coprimality is removed, an entire prime-power overlap
    # can instead retain the same orientation.
    bad_vector = (1, modulus)
    bad = _key(form, bad_vector)
    assert _det(first_vector, bad_vector) == modulus
    assert gcd(first["norm"], bad["norm"]) == modulus
    assert (first["root"] - bad["root"]) % (2 * modulus) == 0
    assert (first["root"] + bad["root"]) % (2 * modulus) != 0
    return {"prime": prime, "exponent": exponent, "other_exponent": other_exponent,
            "first_norm_bits": first["norm"].bit_length(), "second_norm_bits": second["norm"].bit_length(),
            "gcd_bits": common.bit_length(), "gcd": common,
            "proper_roots": [first["root"], second["root"]],
            "canonical_signs": [first["sign"], second["sign"]]}


def run_checks():
    sources = ["tools/check_circuit_box_stencil_overlap.py", "tools/binary_quadratic_forms.py", "tools/quadratic_ideal_hnf.py"]
    hashes = {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in sources}
    counts = {"primitive_forms": 0, "proper_keys": 0, "completion_shift_controls": 0,
              "polarization_identities": 0, "nontrivial_unramified_pair_controls": 0,
              "four_stencil_graph_identities": 0, "larger_stencil_graph_identities": 0,
              "independent_prime_incidence_controls": 0, "prime_power_pairs": 0,
              "determinant_counterexamples": 0}
    signs = {}
    stencil = [(x, y) for x in range(4) for y in range(-3, 4) if gcd(x, y) == 1 and (x > 0 or y > 0)]
    forms = [form for discriminant in (-4, -7, -15, -20, -23, -31, -39, -47, -84, -120, -260, -575, -1811) for form in reduced_forms(discriminant)]
    for form in forms:
        counts["primitive_forms"] += 1
        discriminant = form[1] ** 2 - 4 * form[0] * form[2]
        keys = [_key(form, vector) for vector in stencil]
        counts["proper_keys"] += len(keys)
        for vector, key in zip(stencil, keys):
            for shift in (-3, -1, 1, 4):
                complement = tuple(x + shift * y for x, y in zip(key["complement"], vector))
                assert _det(vector, complement) == 1
                assert _polar(form, vector, complement) % (2 * key["norm"]) == key["root"]
                counts["completion_shift_controls"] += 1
        for i, first in enumerate(keys):
            for j in range(i):
                second = keys[j]
                determinant = _det(stencil[i], stencil[j])
                alpha = _det(stencil[j], first["complement"])
                beta = _det(stencil[i], second["complement"])
                assert determinant * (first["middle"] + second["middle"]) == 2 * beta * second["norm"] - 2 * alpha * first["norm"]
                counts["polarization_identities"] += 1
                common = _remove_support(gcd(first["norm"], second["norm"]), abs(2 * discriminant))
                if common > 1 and gcd(common, determinant) == 1:
                    _pair_identity(form, stencil[i], stencil[j], first, second, common)
                    counts["nontrivial_unramified_pair_controls"] += 1
                    label = f"{first['sign']},{second['sign']}"
                    signs[label] = signs.get(label, 0) + 1
        for vectors, remove_determinants, counter in ((FOUR, False, "four_stencil_graph_identities"), (stencil, True, "larger_stencil_graph_identities")):
            record = _graph_control(form, vectors, remove_determinants=remove_determinants)
            counts[counter] += 1
            counts["independent_prime_incidence_controls"] += record["prime_occurrences"]
    assert len(signs) == 4
    power_records = []
    for prime in (5, 13, 17, 29, 37, 41, 53, 61):
        for exponent in (1, 2, 3, 4, 8, 16):
            for other_exponent in sorted({1, exponent, exponent + 1, 2 * exponent + 1}):
                record = _prime_power_pair(prime, exponent, other_exponent)
                counts["prime_power_pairs"] += 1
                counts["determinant_counterexamples"] += 1
                power_records.append({key: record[key] for key in ("prime", "exponent", "other_exponent", "first_norm_bits", "second_norm_bits", "gcd_bits")})
    large_power = _prime_power_pair(5, 128, 257)
    counts["prime_power_pairs"] += 1
    counts["determinant_counterexamples"] += 1
    explicit_graph = _graph_control((65, 3, 7), FOUR, remove_determinants=False)
    assert explicit_graph["residuals"] == [65, 7, 75, 69]
    assert explicit_graph["leftovers"] == [13, 7, 5, 23]
    assert explicit_graph["integer_content"] == 15
    assert hashes == {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in sources}
    return {"status": "passed", "counts": counts, "canonical_sign_patterns": signs,
            "four_stencil": FOUR, "larger_stencil": stencil,
            "prime_power_controls": power_records, "large_prime_power_control": large_power,
            "explicit_graph": explicit_graph, "source_sha256": hashes,
            "scope": "Exact proper-root correlations and rational-content cancellation; no relation-supply, rank-per-budget, or factoring-score bound."}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run_checks()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], **result["counts"], "largest_gcd_bits": result["large_prime_power_control"]["gcd_bits"], "largest_norm_bits": result["large_prime_power_control"]["second_norm_bits"]}))


if __name__ == "__main__":
    main()
