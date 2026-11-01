#!/usr/bin/env python3
"""Concrete relation construction and factor extraction from square-box search.

The solver receives only N and the chosen finite search parameters.  No
class number, annihilator, complete class-group enumeration, or factors of
N are supplied.  Failure to fill the relation matrix is reported exactly;
the chosen small factor base has no unconditional generation certificate.
"""

from math import gcd
import argparse
import json

from check_positive_relation_kernel import determinant
from early_stop_class_factoring import (
    ambiguous_factorization, factor_from_annihilated_generators,
)
from factored_square_box_search import SquareBoxSearch
from quadratic_form_arithmetic import compose_primitive, inverse_primitive, power_primitive
from universal_smooth_square import dominant_square_relation


def first_split_primes(discriminant, count):
    primes, selected, roots = [], [], {}
    candidate = 2
    while len(selected) < count:
        if all(candidate % prime for prime in primes if prime*prime <= candidate):
            primes.append(candidate)
            if gcd(candidate, discriminant) == 1:
                root = next((value for value in range(candidate+1)
                             if (value*value-discriminant) % (4*candidate) == 0), None)
                if root is not None:
                    selected.append(candidate)
                    roots[candidate] = root
        candidate += 1
    return selected, roots


def run(odd_n, count, cap):
    if odd_n <= 1 or odd_n % 2 == 0 or count < 1 or cap < 1:
        raise ValueError("odd N>1, a positive factor-base size, and a positive cap are required")
    primes, roots = first_split_primes(-4*odd_n, count)
    search = SquareBoxSearch(-4*odd_n, primes, roots, [cap]*count)
    half_diagonal = cap*count+1
    rows, targets, missing = [], [], []
    for index, generator in enumerate(search.generators):
        known_root = power_primitive(generator, half_diagonal)
        target = compose_primitive(known_root, known_root)
        answer = search.solve(target)
        if not answer["represented"]:
            missing.append(index)
            targets.append({"index": index, "target": target, "represented": False,
                            "query_group_multiplications": answer["query_group_multiplications"]})
            continue
        word = answer["half_word"]
        row = dominant_square_relation(index, primes, word, cap)
        check, recovered_root = search.identity, search.identity
        for base, exponent, half_exponent in zip(search.generators, row, word):
            check = compose_primitive(check, power_primitive(base, exponent))
            recovered_root = compose_primitive(recovered_root, power_primitive(base, half_exponent))
        assert check == search.identity
        root_difference = compose_primitive(known_root, inverse_primitive(recovered_root))
        direct = ambiguous_factorization(root_difference, odd_n)
        rows.append(row)
        targets.append({"index": index, "target": target, **answer,
                        "relation_row": row, "root_difference": root_difference,
                        "direct_factor_from_root_difference": direct["proper_factor"]})
    output = {"odd_n": odd_n, "input_bits": odd_n.bit_length(), "factor_base_primes": primes,
              "roots": roots, "cap_fixed_before_targets": cap,
              "half_diagonal": half_diagonal, "index_statistics": search.statistics(),
              "targets": targets, "missing_relation_indices": missing,
              "class_number_supplied": False, "annihilator_supplied": False,
              "complete_generation_claimed": False}
    if missing:
        output["status"] = "bounded_search_did_not_fill_relation_matrix"
        return output
    annihilator = abs(determinant(rows))
    assert annihilator > 0
    for generator in search.generators:
        assert power_primitive(generator, annihilator) == search.identity
    extraction = factor_from_annihilated_generators(odd_n, search.generators, annihilator,
                                                   with_trace=True)
    output.update({"status": extraction["status"], "relation_matrix": rows,
                   "annihilator": annihilator, "annihilator_bits": annihilator.bit_length(),
                   "extraction": extraction})
    if extraction["factor"] is not None:
        factor = extraction["factor"]
        assert 1 < factor < odd_n and odd_n % factor == 0
        output["factorization"] = sorted((factor, odd_n//factor))
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--primes", type=int, default=9)
    parser.add_argument("--cap", type=int, default=4)
    args = parser.parse_args()
    print(json.dumps(run(args.n, args.primes, args.cap), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
