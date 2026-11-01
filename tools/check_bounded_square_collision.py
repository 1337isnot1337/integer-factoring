#!/usr/bin/env python3
"""Independent word and ambiguity replay for bounded square collisions."""

from itertools import combinations, product
from math import comb, gcd, prod
from pathlib import Path
import json

from bounded_square_collision_factoring import run
from check_square_relation_factoring import trial_prime
from quadratic_form_arithmetic import (
    compose_primitive, inverse_primitive, power_primitive, reduce_primitive,
)
from universal_smooth_square import representation_from_square_word


def verify(record, certify_primes=False):
    odd_n = record["odd_n"]
    primes = record["factor_base_primes"]
    roots = {int(prime): root for prime, root in record["roots"].items()}
    count, identity = len(primes), (1, 0, odd_n)
    assert count > 0 and len(set(primes)) == count
    assert all(trial_prime(prime) and gcd(prime, 4*odd_n) == 1 for prime in primes)
    assert all((roots[prime]**2+4*odd_n) % (4*prime) == 0 for prime in primes)
    assert not any(record[field] for field in ("class_number_supplied", "annihilator_supplied",
                                               "input_factors_supplied", "complete_generation_claimed",
                                               "uniform_success_guarantee_claimed"))
    stats = record["statistics"]
    assert 1 <= stats["stopped_depth"] <= record["depth_bound"]
    assert stats["maximum_frontier"] <= record["width"]
    assert stats["retained_keys"] <= 1+record["width"]*record["depth_bound"]
    assert stats["move_compositions"] == stats["squaring_compositions"]
    assert stats["comparison_compositions"] <= stats["move_compositions"]
    assert stats["move_compositions"] <= 2*count*record["width"]*record["depth_bound"]
    if record["status"] != "factored":
        assert record["status"] == "bounded_search_found_no_helpful_collision"
        assert not record["negative_answer_is_exact"]
        return {"odd_n": odd_n, "status": record["status"],
                "miss_is_inconclusive": True, "retained_keys": stats["retained_keys"]}

    generators = [reduce_primitive((prime, roots[prime],
                                    (roots[prime]**2+4*odd_n)//(4*prime)))
                  for prime in primes]
    collision = record["collision"]
    assert collision["kind"] in ("two_generated_roots", "self_inverse_square_key")
    key, oriented_roots = tuple(collision["key"]), []
    word_lengths = []
    for side in ("left", "right"):
        entry = collision[side]
        word, root = entry["word"], tuple(entry["root"])
        assert len(word) == count and entry["orientation"] in (-1, 1)
        length = sum(abs(exponent) for exponent in word)
        assert length <= entry["depth"] <= record["depth_bound"]
        value = identity
        for generator, exponent in zip(generators, word):
            value = compose_primitive(value, power_primitive(generator, exponent))
        assert value == root
        assert power_primitive(root, 2*entry["orientation"]) == key
        oriented_roots.append(power_primitive(root, entry["orientation"]))
        word_lengths.append(length)
    difference = compose_primitive(oriented_roots[0], inverse_primitive(oriented_roots[1]))
    assert difference == tuple(collision["ambiguous_difference"])
    assert power_primitive(difference, 2) == identity and difference != identity
    if collision["kind"] == "self_inverse_square_key":
        assert difference == key
        assert tuple(collision["right"]["root"]) == inverse_primitive(tuple(collision["left"]["root"]))
        assert list(collision["right"]["word"]) == [-value for value in collision["left"]["word"]]
    a, b, c = difference
    if b == 0:
        candidates = (a, c)
    elif b == a:
        assert a % 2 == 0
        candidates = (a//2, 2*c-a//2)
    else:
        assert a == c and b % 2 == 0
        candidates = (a-b//2, a+b//2)
    assert candidates[0]*candidates[1] == odd_n and gcd(*candidates) == 1
    factors = sorted(candidates)
    assert 1 < factors[0] <= factors[1] < odd_n
    assert factors == record["factorization"] == collision["factorization"]
    assert sorted(collision["factor_data"]["coprime_factors"]) == factors
    if certify_primes:
        # Finite post-recovery verification, outside the search algorithm.
        assert max(factors) < 1 << 41
        assert all(trial_prime(factor) for factor in factors)
    result = {"odd_n": odd_n, "input_bits": odd_n.bit_length(), "status": "factored",
            "word_lengths": word_lengths, "factorization": factors,
            "returned_factors_checked_prime": certify_primes,
            "retained_keys": stats["retained_keys"],
            "stopped_depth": stats["stopped_depth"],
            "move_compositions": stats["move_compositions"],
            "squaring_compositions": stats["squaring_compositions"],
            "comparison_compositions": stats["comparison_compositions"],
            "search_rerun": False}
    if certify_primes:
        # Convert the already verified helpful word to a purely integer
        # congruence-of-squares certificate.  The norm constructor receives
        # only the prime forms and that word, never the recovered factors.
        left, right = collision["left"], collision["right"]
        word = [left["orientation"]*a-right["orientation"]*b
                for a, b in zip(left["word"], right["word"])]
        cap = max(1, max(abs(value) for value in word))
        representation = representation_from_square_word(identity, primes, roots, word, cap)
        root_norm = prod(prime**abs(exponent) for prime, exponent in zip(primes, word))
        x, y = representation["primitive_vector"]
        assert representation["primitive_norm"] == root_norm*root_norm
        assert x*x+odd_n*y*y == root_norm*root_norm and gcd(x, y) == 1
        assert gcd(root_norm, odd_n) == 1
        recovered = sorted((gcd(root_norm-x, odd_n), gcd(root_norm+x, odd_n)))
        assert recovered == factors
        result["integer_congruence_certificate"] = {
            "root_norm": root_norm, "root_norm_bits": root_norm.bit_length(),
            "x": x, "y": y, "primitive": True,
            "factorization_of_root_norm": [[prime, abs(exponent)]
                                           for prime, exponent in zip(primes, word) if exponent],
            "gcd_minus_and_plus_sorted": recovered,
            "input_factors_supplied_to_constructor": False}
    return result


def distinct_prime_count(value):
    result, prime = 0, 2
    while prime*prime <= value:
        if value % prime == 0:
            result += 1
            while value % prime == 0:
                value //= prime
        prime += 1
    return result+(value > 1)


def finite_controls():
    successes = prime_power_misses = composite_misses = 0
    records = []
    for odd_n in range(3, 258, 2):
        result = run(odd_n, 4, 32, 16)
        checked = verify(result)
        distinct = distinct_prime_count(odd_n)
        if result["status"] == "factored":
            assert distinct > 1
            successes += 1
        elif distinct == 1:
            prime_power_misses += 1
        else:
            composite_misses += 1
        records.append(checked)
    subgroup = run(87, 1, 32, 16)
    assert subgroup["factor_base_primes"] == [7]
    assert subgroup["status"] != "factored"
    verify(subgroup)
    narrow = run(1621982935400141, 9, 1, 32)
    assert narrow["status"] != "factored"
    verify(narrow)
    tie = run(6557, 1, 1, 2)
    assert tie["status"] == "factored"
    assert tie["collision"]["kind"] == "self_inverse_square_key"
    assert tie["factorization"] == [79, 83]
    verify(tie)
    return {"odd_inputs_checked": len(records), "successful_coprime_splits": successes,
            "prime_power_misses": prime_power_misses,
            "other_composite_misses": composite_misses,
            "proper_subgroup_control": subgroup, "narrow_frontier_control": narrow,
            "self_inverse_key_control": tie,
            "records": records}


def injectivity_and_boundary_controls():
    primes, radius = [3, 5, 7, 11, 13, 17, 19, 23], 4
    support = prod(primes)
    threshold = max(primes)**(4*radius)
    multiplier = threshold//(8*support)+1
    odd_n = (8*multiplier+4)*support-1
    assert threshold < odd_n+1 and odd_n % 8 == 3
    identity = (1, 0, odd_n)
    generators = [reduce_primitive((prime, 2, (odd_n+1)//prime)) for prime in primes]
    keys = set()
    for positions in combinations(range(len(primes)), radius):
        # Fix the first sign to remove the unavoidable e <-> -e duplicate.
        for tail in product((-1, 1), repeat=radius-1):
            word = [0]*len(primes)
            for index, sign in zip(positions, (1,)+tail):
                word[index] = sign
            root = identity
            for generator, exponent in zip(generators, word):
                root = compose_primitive(root, power_primitive(generator, exponent))
            square = power_primitive(root, 2)
            key = min(square, inverse_primitive(square))
            assert key not in keys
            keys.add(key)
    assert len(keys) == (1 << (radius-1))*comb(len(primes), radius)
    boundaries = []
    for radius in range(1, 9):
        root_norm = 3**(2*radius)
        boundary_n = root_norm*root_norm-4
        record = run(boundary_n, 1, 1, radius)
        assert record["status"] == "factored"
        assert record["statistics"]["stopped_depth"] == radius
        assert record["collision"]["kind"] == "self_inverse_square_key"
        assert record["factorization"] == [root_norm-2, root_norm+2]
        assert record["statistics"]["pruned_key_occurrences"] == 0
        boundaries.append({"radius": radius, "odd_n": boundary_n,
                           "threshold": 3**(4*radius), **verify(record)})
    return {"injective_ball_subset": {"odd_n": odd_n, "factor_base_primes": primes,
                                       "radius": 4, "prime_bound_four_radius_power": threshold,
                                       "distinct_square_keys": len(keys),
                                       "exact_support_count": 560},
            "boundary_family": "N=3^(4r)-4; the helpful square key appears exactly at radius r",
            "boundary_controls": boundaries}


def main():
    root = Path(__file__).resolve().parent.parent
    certificates = []
    for bits in (51, 60, 80):
        path = root/f"evidence/bounded-square-collision-{bits}bit.json"
        record = json.loads(path.read_text())
        certificates.append({"certificate": path.name, **verify(record, certify_primes=True)})
    print(json.dumps({"scope": "Independent signed-word, square-key and ambiguous-factor replay",
                      "certificates": certificates, "finite_controls": finite_controls(),
                      "injectivity_and_boundary_controls": injectivity_and_boundary_controls(),
                      "uniform_success_guarantee_claimed": False}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
