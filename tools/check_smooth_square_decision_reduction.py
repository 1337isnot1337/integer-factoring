#!/usr/bin/env python3
"""Exact finite decision oracles for the smooth-square search reduction.

The oracle enumerates complete small subgroups and is not an efficient
implementation of the decision contract for large discriminants.
"""

from collections import Counter, deque
from functools import lru_cache
from pathlib import Path
import json

from binary_quadratic_forms import reduced_forms
from check_smooth_norm_relation import factor_base
from quadratic_form_arithmetic import compose_primitive, power_primitive
from smooth_square_decision_reduction import recover_uniform_square_word
from universal_smooth_square import word_from_square_representation


def check(record):
    discriminant = record["discriminant"]
    cap = record["exponent_cap_fixed_before_targets"]
    forms = reduced_forms(discriminant)
    identity = next(form for form in forms if form[0] == 1)
    base = factor_base(discriminant, record["smooth_bound"])
    primes = [prime for prime, _ in base]
    generators = [form for _, form in base]
    roots = {prime: next(value for value in range(prime+1)
                         if (value*value-discriminant) % (4*prime) == 0)
             for prime in primes}
    assert primes == record["factor_base_primes"]
    multiply = lru_cache(None)(compose_primitive)
    subgroup = {identity}
    queue = deque([identity])
    while queue:
        value = queue.popleft()
        for generator in generators:
            neighbor = multiply(value, generator)
            if neighbor not in subgroup:
                subgroup.add(neighbor)
                queue.append(neighbor)
    targets = {multiply(value, value) for value in subgroup}
    assert len(subgroup) == record["subgroup_order"]
    assert len(targets) == record["square_subgroup_order"]
    cyclic_powers = []
    for generator in generators:
        square = multiply(generator, generator)
        powers, value = [identity], square
        while value != identity:
            assert len(powers) <= len(subgroup)
            powers.append(value)
            value = multiply(value, square)
        cyclic_powers.append(powers)

    @lru_cache(None)
    def reachable(caps):
        # Saturating each cyclic interval is justified here by full finite
        # order enumeration; that is part of this finite oracle's cost.
        states = {identity}
        for powers, radius in zip(cyclic_powers, caps):
            order = len(powers)
            choices = (powers if 2*radius+1 >= order else
                       [powers[exponent % order] for exponent in range(-radius, radius+1)])
            states = {multiply(left, right) for left in states for right in choices}
        return states

    def decide(target, caps):
        assert all(0 <= radius <= cap for radius in caps)
        return target in reachable(caps)

    query_count = maximum_queries = negative_answers = 0
    kinds = Counter()
    example = None
    for target in sorted(targets):
        assert decide(target, (cap,)*len(primes))
        recovered = recover_uniform_square_word(target, primes, roots, cap, decide, with_trace=True)
        query_count += recovered["query_count"]
        maximum_queries = max(maximum_queries, recovered["query_count"])
        witness = recovered["representation"]
        normalized = word_from_square_representation(target, *witness["vector"], primes, cap)
        assert normalized["half_word"] == recovered["half_word"]
        for query in recovered["queries"]:
            assert decide(tuple(query["target"]), tuple(query["caps"])) == query["answer"]
            negative_answers += not query["answer"]
            kinds[query["kind"]] += 1
        if example is None or recovered["query_count"] > example["queries"]:
            example = {"target": target, "half_word": recovered["half_word"],
                       "queries": recovered["query_count"],
                       "primitive_norm_bits": witness["primitive_norm"].bit_length()}
    return {"discriminant": discriminant, "factor_base_primes": primes,
            "exponent_cap": cap, "subgroup_order": len(subgroup),
            "targets_checked": len(targets), "decision_calls": query_count,
            "negative_answers": negative_answers, "maximum_calls_per_target": maximum_queries,
            "query_bound_per_target": len(primes)*((cap-1).bit_length()+2),
            "distinct_cap_vectors_in_finite_oracle": reachable.cache_info().currsize,
            "query_kinds": dict(kinds), "example": example}


def main():
    source = Path(__file__).resolve().parent.parent / "evidence/universal-smooth-square.json"
    records = [check(record) for record in json.loads(source.read_text())["records"]]
    # An always-YES callback is not a valid replacement for the varying-cap
    # decision contract, even on a promised universal-square input.
    invalid_rejected = False
    try:
        recover_uniform_square_word((1, 1, 6), [2], {2: 1}, 1,
                                    lambda target, caps: True)
    except ValueError as error:
        assert "inconsistent" in str(error)
        invalid_rejected = True
    assert invalid_rejected
    assert sum(record["negative_answers"] for record in records) > 0
    print(json.dumps({
        "scope": "Search-to-decision reduction for variable-cap smooth square norms; every decider in this checker uses exhaustive finite group information",
        "discriminants_checked": len(records),
        "target_representations_recovered": sum(record["targets_checked"] for record in records),
        "decision_calls_checked": sum(record["decision_calls"] for record in records),
        "negative_decisions_used": sum(record["negative_answers"] for record in records),
        "invalid_always_yes_oracle_rejected": invalid_rejected,
        "records": records,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
