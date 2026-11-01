#!/usr/bin/env python3
"""Finite dominant-relation construction with a bounded-width word frontier."""

import argparse
import json

from bounded_class_word_search import BoundedClassWordSearch
from check_positive_relation_kernel import determinant
from early_stop_class_factoring import factor_from_annihilated_generators
from quadratic_form_arithmetic import compose_primitive, power_primitive, reduce_primitive
from run_square_relation_factoring import first_split_primes


def run(odd_n, count, width, depth, selection="stratified"):
    if odd_n <= 1 or odd_n % 2 == 0 or count < 1 or width < 1 or depth < 1:
        raise ValueError("odd N>1 and positive finite search parameters are required")
    discriminant = -4*odd_n
    primes, roots = first_split_primes(discriminant, count)
    generators = [reduce_primitive((prime, roots[prime],
                                    (roots[prime]**2-discriminant)//(4*prime)))
                  for prime in primes]
    search = BoundedClassWordSearch(discriminant, generators, width, depth//2,
                                    depth-depth//2, selection)
    diagonal = depth+1
    rows, entries, missing = [], [], []
    for index, generator in enumerate(generators):
        target = power_primitive(generator, diagonal)
        answer = search.solve(target)
        entry = {"index": index, "target": target, **answer}
        if not answer["found"]:
            missing.append(index)
        else:
            row = [-value for value in answer["word"]]
            row[index] += diagonal
            assert row[index] > sum(abs(value) for j, value in enumerate(row) if j != index)
            value = search.identity
            for base, exponent in zip(generators, row):
                value = compose_primitive(value, power_primitive(base, exponent))
            assert value == search.identity
            rows.append(row)
            entry["relation_row"] = row
        entries.append(entry)
    result = {"odd_n": odd_n, "input_bits": odd_n.bit_length(),
              "factor_base_primes": primes, "roots": roots,
              "diagonal_fixed_before_targets": diagonal,
              "search": search.statistics(), "targets": entries,
              "missing_relation_indices": missing,
              "class_number_supplied": False, "annihilator_supplied": False,
              "input_factors_supplied": False, "complete_generation_claimed": False,
              "uniform_success_guarantee_claimed": False}
    if missing:
        result["status"] = "bounded_frontier_did_not_fill_relation_matrix"
        return result
    annihilator = abs(determinant(rows))
    assert annihilator > 0
    assert all(power_primitive(generator, annihilator) == search.identity for generator in generators)
    extraction = factor_from_annihilated_generators(odd_n, generators, annihilator, with_trace=True)
    result.update({"relation_matrix": rows, "annihilator": annihilator,
                   "annihilator_bits": annihilator.bit_length(),
                   "extraction": extraction, "status": extraction["status"]})
    if extraction["factor"] is not None:
        factor = extraction["factor"]
        assert 1 < factor < odd_n and odd_n % factor == 0
        result["factorization"] = sorted((factor, odd_n//factor))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--primes", type=int, default=9)
    parser.add_argument("--width", type=int, default=1024)
    parser.add_argument("--depth", type=int, default=64)
    parser.add_argument("--selection", choices=["height", "stratified"], default="stratified")
    args = parser.parse_args()
    print(json.dumps(run(args.n, args.primes, args.width, args.depth, args.selection),
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
