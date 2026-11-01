#!/usr/bin/env python3
"""Exact controls for the obstruction to decreasing reduced leading coefficients.

The small case enumerates its complete class group to certify a minimum
escape weight.  This finite enumeration is not an efficient general
navigation algorithm.  The large case verifies the elementary product-
norm argument on a bounded exponent box without enumerating its group.
"""

from heapq import heappop, heappush
from itertools import product
from math import prod
import json

from binary_quadratic_forms import reduced_forms
from quadratic_form_arithmetic import (
    compose_primitive, inverse_primitive, power_primitive, reduce_primitive,
)
from run_square_relation_factoring import first_split_primes


def prime_form(odd_n, prime, root):
    return reduce_primitive((prime, root, (root*root+4*odd_n)//(4*prime)))


def finite_escape_certificate():
    odd_n = 1022117
    primes, roots = first_split_primes(-4*odd_n, 4)
    omitted_prime = primes.pop()
    generators = [prime_form(odd_n, prime, roots[prime]) for prime in primes]
    exceptional = (2, 2, (odd_n+1)//2)
    start = prime_form(odd_n, omitted_prime, roots[omitted_prime])
    identity = (1, 0, odd_n)
    assert primes == [3, 7, 11] and omitted_prime == 17
    assert start == (17, 10, 60126)
    assert power_primitive(exceptional, 2) == identity
    generators.append(exceptional)
    weights = primes+[2]
    moves = [(index, sign, weight, generator if sign == 1 else inverse_primitive(generator))
             for index, (weight, generator) in enumerate(zip(weights, generators))
             for sign in (-1, 1)]
    distances, predecessors = {start: 1}, {start: None}
    queue, settled = [(1, start)], set()
    while queue:
        distance, form = heappop(queue)
        if form in settled:
            continue
        settled.add(form)
        for index, sign, weight, generator in moves:
            neighbor = compose_primitive(form, generator)
            proposal = distance*weight
            if proposal < distances.get(neighbor, proposal+1):
                distances[neighbor] = proposal
                predecessors[neighbor] = (form, index, sign)
                heappush(queue, (proposal, neighbor))

    # Independent completeness and optimality certificate: the labels
    # satisfy every edge inequality and have strictly descending equality
    # predecessors.  These properties characterize exact path costs.
    forms = reduced_forms(-4*odd_n)
    assert settled == set(forms)
    for form in forms:
        assert distances[form] >= 1
        if form != start:
            parent, index, sign = predecessors[form]
            generator = power_primitive(generators[index], sign)
            assert compose_primitive(parent, generator) == form
            assert distances[form] == weights[index]*distances[parent]
            assert distances[parent] < distances[form]
        for index, sign, weight, generator in moves:
            neighbor = compose_primitive(form, generator)
            assert distances[neighbor] <= weight*distances[form]
    assert distances[start] == 1

    # Verify that the odd unramified prime forms already generate the whole
    # group, so the obstruction is not explained by a missing generator.
    reached, pending = {identity}, [identity]
    while pending:
        form = pending.pop()
        for generator in generators[:-1]:
            neighbor = compose_primitive(form, generator)
            if neighbor not in reached:
                reached.add(neighbor)
                pending.append(neighbor)
    assert reached == set(forms)

    escape = min((distances[form], form) for form in forms if form[0] < start[0])
    escape_weight, escape_form = escape
    assert escape == (33957, (14, 6, 73009))
    assert distances[identity] == 475398
    path, word, form = [], [0]*len(generators), escape_form
    while form != start:
        parent, index, sign = predecessors[form]
        word[index] += sign
        path.append({"prime_or_exception_weight": weights[index], "sign": sign,
                     "resulting_form": form})
        form = parent
    path.reverse()
    assert word == [-2, -3, -1, 0]
    neighbors = [{"prime_or_exception_weight": weight, "sign": sign,
                  "form": compose_primitive(start, generator)}
                 for index, sign, weight, generator in moves]
    assert all(entry["form"][0] > start[0] for entry in neighbors)
    assert all(distance >= escape_weight or form[0] >= start[0]
               for form, distance in distances.items())
    labels = []
    for form in forms:
        predecessor = predecessors[form]
        labels.append({"form": form, "minimum_path_weight": distances[form],
                       "predecessor": (None if predecessor is None else
                                       {"form": predecessor[0], "generator_index": predecessor[1],
                                        "sign": predecessor[2]})})
    return {"odd_n": odd_n, "discriminant": -4*odd_n,
            "complete_class_count": len(forms), "full_generation_checked": True,
            "factor_base_primes": primes, "roots": roots,
            "omitted_prime": omitted_prime, "start": start,
            "exceptional_generator": exceptional, "neighbors": neighbors,
            "minimum_descent_weight": escape_weight, "first_descent_form": escape_form,
            "descent_half_word_including_exception": word, "descent_path": path,
            "minimum_identity_weight": distances[identity],
            "classes_with_smaller_weight_than_first_descent":
                sum(distance < escape_weight for distance in distances.values()),
            "edge_inequalities_checked": len(forms)*len(moves),
            "equality_predecessors_checked": len(forms)-1,
            "optimality_labels": labels}


def large_box_control():
    primes, omitted_prime = [3, 5, 7, 11], 13
    support_product = prod(primes+[omitted_prime])
    odd_n = (8*(1 << 180)+4)*support_product-1
    assert odd_n % 8 == 3
    identity = (1, 0, odd_n)
    generators = [prime_form(odd_n, prime, 2) for prime in primes]
    start = prime_form(odd_n, omitted_prime, 2)
    exceptional = (4, 2, (odd_n+1)//4)
    exceptional_classes = [identity, exceptional, inverse_primitive(exceptional)]
    assert power_primitive(exceptional, 3) == identity
    assert len(set(exceptional_classes)) == 3
    cap, checked, maximum_oriented_norm = 3, 0, 0
    for word in product(range(-cap, cap+1), repeat=len(primes)):
        oriented_norm = prod(prime**abs(exponent) for prime, exponent in zip(primes, word))
        maximum_oriented_norm = max(maximum_oriented_norm, oriented_norm)
        assert (4*omitted_prime*oriented_norm)**2 < odd_n
        value = start
        for generator, exponent in zip(generators, word):
            value = compose_primitive(value, power_primitive(generator, exponent))
        for exceptional_class in exceptional_classes:
            result = compose_primitive(value, exceptional_class)
            expected = exceptional_class[0]*omitted_prime*oriented_norm
            assert result[0] == expected
            checked += 1
    return {"odd_n": odd_n, "input_bits": odd_n.bit_length(),
            "discriminant_bits": (4*odd_n).bit_length(),
            "factor_base_primes": primes, "all_splitting_roots": 2,
            "omitted_prime": omitted_prime, "start": start,
            "exceptional_classes": exceptional_classes, "cap": cap,
            "signed_words_checked": (2*cap+1)**len(primes),
            "exceptional_coset_products_checked": checked,
            "maximum_oriented_norm": maximum_oriented_norm,
            "full_generation_claimed": False, "class_group_enumerated": False,
            "factorization_of_input_used": False}


def main():
    print(json.dumps({
        "scope": "Exact arithmetic obstruction to strict reduced-leading-coefficient descent",
        "general_navigation_lower_bound_claimed": False,
        "improved_complete_factoring_bound_claimed": False,
        "finite_complete_group": finite_escape_certificate(),
        "large_product_norm_control": large_box_control(),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
