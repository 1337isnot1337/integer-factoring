#!/usr/bin/env python3
"""Factor from a helpful square collision in an explicit bounded frontier.

No class number, relation determinant, annihilator, input factors, or
random samples are supplied.  Failure within the finite width/depth
budget is inconclusive.  The loop costs are recorded explicitly.
"""

from bisect import bisect_left
from heapq import merge
import argparse
import json

from bounded_class_word_search import choose_frontier
from early_stop_class_factoring import ambiguous_factorization
from quadratic_form_arithmetic import (
    compose_primitive, inverse_primitive, power_primitive, reduce_primitive,
)
from run_square_relation_factoring import first_split_primes


def run(odd_n, count, width, depth):
    if odd_n <= 1 or odd_n % 2 == 0 or min(count, width, depth) < 1:
        raise ValueError("odd N>1 and positive search parameters are required")
    discriminant = -4*odd_n
    primes, roots = first_split_primes(discriminant, count)
    generators = [reduce_primitive((prime, roots[prime],
                                    (roots[prime]**2-discriminant)//(4*prime)))
                  for prime in primes]
    identity = (1, 0, odd_n)
    exceptional = (2, 2, (odd_n+1)//2) if odd_n % 4 == 1 else identity
    moves = [(index, sign, generator if sign == 1 else inverse_primitive(generator))
             for index, generator in enumerate(generators) for sign in (-1, 1)]
    # Entry: canonical square key, root, signed generator word, depth,
    # and sign s satisfying root^(2s)=key.
    initial = (identity, identity, (0,)*count, 0, 1)
    seen, seen_keys, frontier = [initial], [identity], [initial]
    stats = {"move_compositions": 0, "squaring_compositions": 0,
             "comparison_compositions": 0, "unhelpful_collisions": 0,
             "pruned_key_occurrences": 0, "maximum_frontier": 1, "layers": []}

    def collision(left, right):
        assert left[0] == right[0]
        left_oriented = left[1] if left[4] == 1 else inverse_primitive(left[1])
        right_oriented = right[1] if right[4] == 1 else inverse_primitive(right[1])
        difference = compose_primitive(left_oriented, inverse_primitive(right_oriented))
        stats["comparison_compositions"] += 1
        data = ambiguous_factorization(difference, odd_n)
        factor = data["proper_factor"]
        if factor is None:
            assert difference in (identity, exceptional)
            stats["unhelpful_collisions"] += 1
            return None
        assert 1 < factor < odd_n and odd_n % factor == 0
        return {"key": left[0],
                "kind": "two_generated_roots",
                "left": {"root": left[1], "word": left[2], "depth": left[3],
                         "orientation": left[4]},
                "right": {"root": right[1], "word": right[2], "depth": right[3],
                          "orientation": right[4]},
                "ambiguous_difference": difference, "factor_data": data,
                "factorization": sorted((factor, odd_n//factor))}

    def finish(answer, stopped_depth):
        stats["retained_keys"] = len(seen)
        stats["stopped_depth"] = stopped_depth
        result = {"odd_n": odd_n, "input_bits": odd_n.bit_length(),
                  "factor_base_primes": primes, "roots": roots,
                  "width": width, "depth_bound": depth,
                  "selection": "stratified", "statistics": stats,
                  "class_number_supplied": False, "annihilator_supplied": False,
                  "input_factors_supplied": False, "complete_generation_claimed": False,
                  "uniform_success_guarantee_claimed": False,
                  "count_scope": "Move, squaring, and collision-comparison loops; setup and certificate replay are separate"}
        if answer is None:
            result.update({"status": "bounded_search_found_no_helpful_collision",
                           "negative_answer_is_exact": False})
        else:
            result.update({"status": "factored", "collision": answer,
                           "factorization": answer["factorization"]})
        return result

    for current_depth in range(1, depth+1):
        pending = []
        for _, form, word, _, _ in frontier:
            for index, sign, generator in moves:
                root = compose_primitive(form, generator)
                stats["move_compositions"] += 1
                new_word = list(word)
                new_word[index] += sign
                new_word = tuple(new_word)
                square = compose_primitive(root, root)
                stats["squaring_compositions"] += 1
                inverse_square = inverse_primitive(square)
                key, orientation = ((square, 1) if square <= inverse_square
                                    else (inverse_square, -1))
                entry = (key, root, new_word, current_depth, orientation)
                if square == inverse_square:
                    # An order-four root can expose a helpful ambiguous
                    # square before a second root has been retained.  In
                    # this tie, identifying h and h^-1 would otherwise
                    # discard the helpful quotient h^2.  The two explicit
                    # words e and -e give an independently replayable pair.
                    data = ambiguous_factorization(square, odd_n)
                    if data["proper_factor"] is not None:
                        factor = data["proper_factor"]
                        answer = {"key": square, "kind": "self_inverse_square_key",
                                  "left": {"root": root, "word": new_word,
                                           "depth": current_depth, "orientation": 1},
                                  "right": {"root": inverse_primitive(root),
                                            "word": tuple(-value for value in new_word),
                                            "depth": current_depth, "orientation": 1},
                                  "ambiguous_difference": square, "factor_data": data,
                                  "factorization": sorted((factor, odd_n//factor))}
                        return finish(answer, current_depth)
                position = bisect_left(seen_keys, key)
                if position < len(seen_keys) and seen_keys[position] == key:
                    answer = collision(seen[position], entry)
                    if answer is not None:
                        return finish(answer, current_depth)
                    continue
                pending.append(entry)
        pending.sort()
        distinct = []
        for entry in pending:
            if distinct and distinct[-1][0] == entry[0]:
                answer = collision(distinct[-1], entry)
                if answer is not None:
                    return finish(answer, current_depth)
            else:
                distinct.append(entry)
        ranked = sorted((root, word, entry_depth, key, orientation)
                        for key, root, word, entry_depth, orientation in distinct)
        selected = choose_frontier(ranked, width, "stratified")
        chosen = [(key, root, word, entry_depth, orientation)
                  for root, word, entry_depth, key, orientation in selected]
        pruned = len(distinct)-len(chosen)
        stats["pruned_key_occurrences"] += pruned
        stats["maximum_frontier"] = max(stats["maximum_frontier"], len(chosen))
        stats["layers"].append({"depth": current_depth, "input_frontier": len(frontier),
                                "distinct_new_keys": len(distinct),
                                "retained_frontier": len(chosen), "pruned_new_keys": pruned})
        if not chosen:
            return finish(None, current_depth)
        seen = list(merge(seen, sorted(chosen)))
        seen_keys = [entry[0] for entry in seen]
        frontier = chosen
    return finish(None, depth)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--primes", type=int, default=12)
    parser.add_argument("--width", type=int, default=1024)
    parser.add_argument("--depth", type=int, default=64)
    args = parser.parse_args()
    print(json.dumps(run(args.n, args.primes, args.width, args.depth), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
