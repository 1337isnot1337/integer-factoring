#!/usr/bin/env python3
"""Exact early-stop controls in product groups and quadratic class groups."""

import json
from math import lcm, prod
import random

from binary_quadratic_forms import reduced_forms
from early_ambiguous_factor import factor_from_annihilator, first_useful_top


def product_group_control(moduli, generators, bad_top, annihilator=None):
    identity = (0,)*len(moduli)
    exponents = [(modulus & -modulus).bit_length()-1 for modulus in moduli]
    order_bound = max(exponents, default=0)
    if annihilator is None:
        annihilator = (lcm(*moduli) if moduli else 1)*2**17*45

    def multiply(left, right):
        return tuple((a+b) % modulus for a, b, modulus in zip(left, right, moduli))

    def inverse(value):
        return tuple(-a % modulus for a, modulus in zip(value, moduli))

    top, certificate, stats = first_useful_top(
        generators, annihilator, identity, multiply, inverse, order_bound,
        lambda value: None if value == bad_top else "outside exceptional subgroup",
    )
    rank = sum(exponent > 0 for exponent in exponents)
    useful_exists = rank > 1 or (rank == 1 and bad_top is None)
    assert (top is not None) == useful_exists
    if top is not None:
        assert top != identity and top != bad_top
        assert multiply(top, top) == identity
        assert certificate is not None
    assert stats["peak_retained_cyclic_directions"] <= 1
    bound = len(generators)*(2*annihilator.bit_length()+4*(order_bound+1)**2)
    assert stats["group_multiplications"] <= bound
    assert stats["group_inversions"] <= len(generators)*order_bound
    return {
        "moduli": moduli,
        "group_size": prod(moduli),
        "two_rank": rank,
        "bad_top": bad_top,
        "returned_top": top,
        "annihilator_bits": annihilator.bit_length(),
        "polynomial_operation_bound": bound,
        "stats": stats,
    }


def has_multiple_prime_factors(integer):
    remaining = integer
    distinct = 0
    divisor = 2
    while divisor*divisor <= remaining:
        if remaining % divisor == 0:
            distinct += 1
            while remaining % divisor == 0:
                remaining //= divisor
        divisor += 1
    distinct += remaining > 1
    return distinct >= 2


def class_group_control(integer):
    forms = reduced_forms(-4*integer)
    identity = (1, 0, integer)
    # All-class enumeration only supplies a fully certified small fixture.
    # The algorithm under test never enumerates classes or computes their size.
    generators = [form for form in forms if form != identity]
    tau = (2, 2, (integer+1)//2) if integer % 4 == 1 else None
    if tau in generators:
        generators.remove(tau)
        generators.insert(0, tau)
    annihilator = len(forms)*2**37*45
    result = factor_from_annihilator(integer, generators, annihilator)
    factor = result["factor"]
    assert (factor is not None) == has_multiple_prime_factors(integer)
    if factor is not None:
        assert 1 < factor < integer and integer % factor == 0
    stats = result["stats"]
    bound = len(generators)*(2*annihilator.bit_length()+4*((4*integer).bit_length()+3)**2)
    assert stats["group_multiplications"] <= bound
    assert stats["peak_retained_cyclic_directions"] <= 1
    return {
        "integer": integer,
        "class_number": len(forms),
        "generator_fixture": "all nonprincipal reduced classes, with exceptional top first",
        "annihilator_bits": annihilator.bit_length(),
        "returned_factor": factor,
        "ambiguous_form": result["ambiguous_form"],
        "stats": stats,
    }


def main():
    mixed = product_group_control([8, 2], [(1, 0), (1, 1)], (4, 0))
    assert mixed["stats"]["cancellations"] == 1
    assert mixed["returned_top"] == (0, 1)
    swap = product_group_control([8, 2], [(2, 1), (1, 0)], (4, 0))
    assert swap["stats"]["order_swaps"] == 1
    assert swap["stats"]["cancellations"] == 1
    deep = product_group_control([2**128], [(1,), (2**128-1,)], (2**127,))
    assert deep["returned_top"] is None
    assert deep["stats"]["cancellations"] == 128

    rng = random.Random(20261031)
    records = []
    for dimension in range(25):
        for trial in range(12):
            moduli = [2**rng.randrange(0, 13)*rng.choice((1, 3, 5, 9, 15))
                      for _ in range(dimension)]
            generators = [tuple(int(i == j) % modulus for i, modulus in enumerate(moduli))
                          for j in range(dimension)]
            if dimension > 1:
                for _ in range(4*dimension):
                    i, j = rng.sample(range(dimension), 2)
                    generators[i] = tuple((a+b) % modulus
                                          for a, b, modulus in zip(generators[i], generators[j], moduli))
            rng.shuffle(generators)
            generators += [(0,)*dimension]
            even_indices = [i for i, modulus in enumerate(moduli) if modulus % 2 == 0]
            bad_top = None
            if even_indices and trial % 3:
                selected = {i for i in even_indices if rng.randrange(2)}
                selected.add(rng.choice(even_indices))
                bad_top = tuple(modulus//2 if i in selected else 0
                                for i, modulus in enumerate(moduli))
            records.append(product_group_control(moduli, generators, bad_top))

    integers = (3, 5, 7, 9, 11, 13, 15, 21, 25, 27, 33, 35, 45, 49,
                55, 65, 75, 77, 81, 105, 121, 125, 165, 231, 385, 1155,
                1365, 10403, 10807, 16637, 15015, 19635)
    classes = [class_group_control(integer) for integer in integers]
    print(json.dumps({
        "scope": "Exact polynomial early-extraction controls; class-group generators and annihilators are supplied by exhaustive small fixtures",
        "generic_controls": len(records)+3,
        "maximum_generic_two_rank": max(row["two_rank"] for row in records),
        "maximum_generic_group_size": max(row["group_size"] for row in records),
        "maximum_retained_cyclic_directions": max(row["stats"]["peak_retained_cyclic_directions"] for row in records),
        "mixed_generator_control": mixed,
        "order_swap_control": swap,
        "deep_cancellation_control": deep,
        "class_group_controls": len(classes),
        "successful_factor_targets": sum(row["returned_factor"] is not None for row in classes),
        "class_group_records": classes,
    }, indent=2))


if __name__ == "__main__":
    main()
