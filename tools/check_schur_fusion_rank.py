#!/usr/bin/env python3
"""Exact finite checks for the saturated-filtration Schur fusion rank proof.

The monomial expansion uses set-partition Mobius inversion. Products are
also checked directly in finite variables. Fusion uses the independent
pinned source power-sum recursion. No even denominator is reduced mod 2.
"""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, factorial
from pathlib import Path

from check_partition_channels import operator, partitions
from check_modular_crossing import rank_two


@lru_cache(None)
def set_partitions(size):
    if not size:
        return ((),)
    result = []
    for blocks in set_partitions(size-1):
        result.append(blocks+((size-1,),))
        for i in range(len(blocks)):
            result.append(blocks[:i]+(blocks[i]+(size-1,),)+blocks[i+1:])
    return tuple(result)


def multiplicity_factor(partition):
    result = 1
    for count in Counter(partition).values():
        result *= factorial(count)
    return result


@lru_cache(None)
def monomial_power_sums(partition):
    """m_partition in rational power-sum coordinates, by distinct indices."""
    result = defaultdict(Fraction)
    denominator = multiplicity_factor(partition)
    for blocks in set_partitions(len(partition)):
        coefficient = 1
        for block in blocks:
            coefficient *= (-1)**(len(block)-1)*factorial(len(block)-1)
        term = tuple(sorted(sum(partition[i] for i in block) for block in blocks))
        result[term] += Fraction(coefficient, denominator)
    return tuple((term, value) for term, value in result.items() if value)


def power_sums_to_monomials(polynomial, weight):
    remaining = defaultdict(Fraction, polynomial)
    result = {}
    for term in sorted(partitions(weight), key=lambda p: (-len(p), p)):
        value = remaining[term]*multiplicity_factor(term)
        if not value:
            continue
        result[term] = value
        for output, coefficient in monomial_power_sums(term):
            remaining[output] -= value*coefficient
    assert not any(remaining.values())
    return result


def power_sum_product(left, right):
    result = defaultdict(Fraction)
    for a, x in left:
        for b, y in right:
            result[tuple(sorted(a+b))] += x*y
    return {term: value for term, value in result.items() if value}


def finite_variable_monomial(partition, values):
    """Sum distinct exponent assignments, independently of power sums."""
    if len(partition) > len(values):
        return 0
    counts = Counter(partition+(0,)*(len(values)-len(partition)))
    exponents = tuple(counts)

    def recurse(index):
        if index == len(values):
            return 1
        total = 0
        for exponent in exponents:
            if counts[exponent]:
                counts[exponent] -= 1
                total += values[index]**exponent*recurse(index+1)
                counts[exponent] += 1
        return total

    return recurse(0)


def direct_product_coefficient(first, second, output):
    """Coefficient of one representative monomial x^output in m_A m_B."""
    if len(first) > len(output) or len(second) > len(output):
        return 0
    counts = Counter(first+(0,)*(len(output)-len(first)))
    exponents = tuple(counts)

    def recurse(index, remainder):
        if index == len(output):
            return int(tuple(sorted(k for k in remainder if k)) == second)
        total = 0
        for exponent in exponents:
            if counts[exponent] and exponent <= output[index]:
                counts[exponent] -= 1
                total += recurse(index+1, remainder+(output[index]-exponent,))
                counts[exponent] += 1
        return total

    return recurse(0, ())


def graded_coefficient(first, second):
    a, b = Counter(first), Counter(second)
    result = 1
    for part in a.keys() | b.keys():
        result *= comb(a[part]+b[part], a[part])
    return result


def carry_free(first, second):
    a, b = Counter(first), Counter(second)
    return all(not (a[part] & b[part]) for part in a.keys() | b.keys())


def fusion_monomials(first, second):
    assert sum(first) == sum(second)
    weight = sum(first)
    a, b = monomial_power_sums(first), monomial_power_sums(second)
    result = {}
    for output in partitions(2*weight):
        value = sum(x*y*operator(mu, output, nu)
                    for mu, x in a for nu, y in b)
        if value:
            result[output] = value
    return power_sums_to_monomials(result, 2*weight)


def stream_partitions(weight, least=1):
    if weight == 0:
        yield ()
        return
    for part in range(least, weight+1):
        for tail in stream_partitions(weight-part, part):
            yield (part,)+tail


def split_atoms(partition, target):
    """Subset-sum test on the tagged binary multiplicity atoms."""
    reachable = 1
    for part, count in Counter(partition).items():
        power = 1
        while count:
            if count & 1:
                reachable |= reachable << (part*power)
            count >>= 1
            power <<= 1
    return bool((reachable >> target) & 1)


def split_ordinary_parts(partition, target):
    reachable = 1
    for part in partition:
        reachable |= reachable << part
    return bool((reachable >> target) & 1)


def check_expansions():
    cases = leading = evaluations = 0
    values_list = ((1, -1, 2, -2, 0, 1, 3, -1),
                   (0, 1, 1, 2, 2, 3, -1, -2))
    for weight in range(9):
        for term in partitions(weight):
            polynomial = dict(monomial_power_sums(term))
            assert polynomial.get(term) == Fraction(1, multiplicity_factor(term))
            assert all(len(other) < len(term) for other in polynomial if other != term)
            leading += 1
            assert power_sums_to_monomials(polynomial, weight) == {term: Fraction(1)}
            for values in values_list:
                power_sums = {k: sum(x**k for x in values) for k in range(1, weight+1)}
                evaluated = Fraction(0)
                for mu, coefficient in polynomial.items():
                    value = coefficient
                    for k in mu:
                        value *= power_sums[k]
                    evaluated += value
                assert evaluated == finite_variable_monomial(term, values)
                evaluations += 1
            cases += len(polynomial)
    return {"rational_expansion_terms": cases,
            "unique_leading_terms": leading,
            "independent_finite_variable_evaluations": evaluations}


def check_products():
    labels = tuple(term for weight in range(5) for term in partitions(weight))
    entries = graded = 0
    for first in labels:
        for second in labels:
            weight = sum(first)+sum(second)
            product_ = power_sums_to_monomials(power_sum_product(
                monomial_power_sums(first), monomial_power_sums(second)), weight)
            assert all(c.denominator == 1 for c in product_.values())
            for output in partitions(weight):
                actual = direct_product_coefficient(first, second, output)
                assert product_.get(output, 0) == actual
                entries += 1
            union = tuple(sorted(first+second))
            top = graded_coefficient(first, second)
            assert product_.get(union, 0) == top
            assert all(len(p) < len(union) for p in product_ if p != union)
            assert (top % 2 == 1) == carry_free(first, second)
            graded += 1
    return {"independent_monomial_product_entries": entries,
            "carry_free_graded_products": graded}


def check_fusions():
    records = []
    top_entries = contraction_entries = 0
    for weight in range(1, 5):
        labels, output = partitions(weight), partitions(2*weight)
        rows = [0]*len(output)
        selected = {}
        column = 0
        for first in labels:
            for second in labels:
                fusion = fusion_monomials(first, second)
                product_ = power_sums_to_monomials(power_sum_product(
                    monomial_power_sums(first), monomial_power_sums(second)), 2*weight)
                assert all(c.denominator == 1 for c in fusion.values())
                length = len(first)+len(second)
                for j, term in enumerate(output):
                    value = fusion.get(term, 0)
                    difference = value-product_.get(term, 0)
                    assert not difference or len(term) <= length-2
                    contraction_entries += 1
                    if int(value) % 2:
                        rows[j] |= 1 << column
                union = tuple(sorted(first+second))
                top = graded_coefficient(first, second)
                assert fusion.get(union, 0) == top
                assert all(len(p) < length for p in fusion if p != union)
                top_entries += 1
                if carry_free(first, second):
                    # Keep one fusion vector for each distinct leading output.
                    vector = sum((int(fusion.get(p, 0)) % 2) << j
                                 for j, p in enumerate(output))
                    selected.setdefault(union, vector)
                column += 1
        atom_outputs = {p for p in output if split_atoms(p, weight)}
        assert set(selected) == atom_outputs
        assert rank_two(list(selected.values())) == len(selected)
        actual = rank_two(rows)
        assert actual >= len(selected)
        expected = {1: 1, 2: 3, 4: 16}
        if weight in expected:
            assert actual == expected[weight]
        records.append({"h": weight, "labels": len(labels), "ambient": len(output),
                        "carry_free_splittable_outputs": len(selected),
                        "independent_triangular_vectors": len(selected),
                        "actual_fusion_rank": actual})
    return {"fusion_top_identities": top_entries,
            "integral_two_length_drop_entries": contraction_entries,
            "small_fusion_ranks": records}


def count_splittable(max_half_weight):
    records = []
    for h in range(1, max_half_weight+1):
        total = ordinary = binary = 0
        for partition in stream_partitions(2*h):
            total += 1
            ordinary += split_ordinary_parts(partition, h)
            binary += split_atoms(partition, h)
        assert binary <= ordinary <= total
        records.append({"h": h, "partitions_2h": total,
                        "ordinary_splittable": ordinary,
                        "carry_free_splittable": binary})
    return records


def verify(max_half_weight=16):
    result = {}
    result.update(check_expansions())
    result.update(check_products())
    result.update(check_fusions())
    counts = count_splittable(max_half_weight)
    result["partition_counts"] = counts
    assert split_ordinary_parts((1, 1), 1)
    assert not split_atoms((1, 1), 1)
    assert graded_coefficient((1,), (1,)) == 2
    assert dict(monomial_power_sums((1, 1))) == {
        (1, 1): Fraction(1, 2), (2,): Fraction(-1, 2)}

    source_path = Path(__file__).resolve().parents[1]/"evidence/schur-fusion-gram.json"
    source = json.loads(source_path.read_text())
    recorded = next(item for item in source["weights"] if item["h"] == 8)
    assert (recorded["dimension"], recorded["ambient_dimension"],
            recorded["fusion_rank"], recorded["central_rank"]) == (22, 231, 229, 228)
    if max_half_weight >= 8:
        assert counts[7]["carry_free_splittable"] <= recorded["fusion_rank"]
    result["reused_weight_eight_checkpoint"] = {
        "source": "evidence/schur-fusion-gram.json",
        "source_sha256": sha256(source_path.read_bytes()).hexdigest(),
        "recomputed_by_this_checker": False,
        "record": recorded,
    }
    result["negative_controls"] = [
        "m_(1,1) has a nonintegral leading power-sum coefficient before saturation",
        "ordinary subset splitting does not imply carry-free splitting",
    ]
    result["scope"] = ("finite checks of the saturated integral filtration and "
                       "carry-free fusion lower bound; asymptotic counting is proved separately")
    result["new_factoring_exponent_claimed"] = False
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-half-weight", type=int, default=16)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify(args.max_half_weight)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": "pass",
                      "leading_terms": result["unique_leading_terms"],
                      "finite_variable_evaluations": result["independent_finite_variable_evaluations"],
                      "product_entries": result["independent_monomial_product_entries"],
                      "fusion_entries": result["integral_two_length_drop_entries"],
                      "max_counted_half_weight": args.max_half_weight}))


if __name__ == "__main__":
    main()
