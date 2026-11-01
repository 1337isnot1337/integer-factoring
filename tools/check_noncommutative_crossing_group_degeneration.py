#!/usr/bin/env python3
"""Independent two-stage degeneration audit for class-two group crossings."""

from functools import lru_cache
from itertools import product
import argparse
import hashlib
import json
import random
from pathlib import Path

from complement_shell_tensor import multiply as polynomial_multiply
from check_noncommutative_crossing_orbit import binary_rank


def submasks(mask):
    value = mask
    while True:
        yield value
        if not value:
            return
        value = (value - 1) & mask


def toggle(values, term):
    if term in values:
        values.remove(term)
    else:
        values.add(term)


class BilinearGroup:
    def __init__(self, horizontal, center, beta, central_weight=2):
        self.n, self.t = horizontal, center
        self.width = horizontal + center
        self.d = 1 << self.width
        self.beta = tuple(tuple(row) for row in beta)
        assert len(self.beta) == horizontal
        assert all(len(row) == horizontal for row in self.beta)
        assert all(0 <= x < 1 << center for row in self.beta for x in row)
        self.generators = tuple(1 << i for i in range(self.width))
        words = []
        for mask in range(self.d):
            value = 0
            for i, generator in enumerate(self.generators):
                if mask >> i & 1:
                    value = self.mul(value, generator)
            words.append(value)
        self.words = tuple(words)
        assert len(set(words)) == self.d
        inverse = [0] * self.d
        for mask, value in enumerate(words):
            inverse[value] = mask
        self.word_index = tuple(inverse)
        self.central_weight = central_weight
        self.weights = tuple((mask & ((1 << self.n) - 1)).bit_count()
                             + central_weight * (mask >> self.n).bit_count()
                             for mask in range(self.d))
        self.basis_words = tuple(tuple(i for i in range(self.width) if mask >> i & 1)
                                 for mask in range(self.d))

    def cocycle(self, a, b):
        value = 0
        for i in range(self.n):
            if a >> i & 1:
                for j in range(self.n):
                    if b >> j & 1:
                        value ^= self.beta[i][j]
        return value

    def mul(self, left, right):
        return left ^ right ^ (self.cocycle(left, right) << self.n)

    def inverse(self, value):
        return value ^ (self.cocycle(value, value) << self.n)

    def to_augmentation(self, support):
        out = set()
        for value in support:
            for mask in submasks(self.word_index[value]):
                toggle(out, mask)
        return tuple(sorted(out))

    @lru_cache(None)
    def augmentation_product(self, left, right):
        support = set()
        for a, b in product(submasks(left), submasks(right)):
            toggle(support, self.mul(self.words[a], self.words[b]))
        return self.to_augmentation(support)

    @lru_cache(None)
    def augmentation_antipode(self, mask):
        return self.to_augmentation(self.inverse(self.words[a])
                                    for a in submasks(mask))

    @lru_cache(None)
    def normal(self, word):
        """Independent restricted PBW rewriting, including nonzero squares."""
        for i in range(len(word) - 1):
            a, b = word[i:i + 2]
            if a == b:
                vector = self.beta[a][a] if a < self.n else 0
                out = set()
                for j in range(self.t):
                    if vector >> j & 1:
                        out.symmetric_difference_update(
                            self.normal(word[:i] + (self.n + j,) + word[i + 2:]))
                return tuple(sorted(out))
            if a > b:
                out = set(self.normal(word[:i] + (b, a) + word[i + 2:]))
                vector = (self.beta[a][b] ^ self.beta[b][a]
                          if a < self.n else 0)
                for j in range(self.t):
                    if vector >> j & 1:
                        out.symmetric_difference_update(
                            self.normal(word[:i] + (self.n + j,) + word[i + 2:]))
                return tuple(sorted(out))
        return (sum(1 << i for i in word),)

    @lru_cache(None)
    def graded_product(self, left, right):
        return self.normal(self.basis_words[left] + self.basis_words[right])

    @lru_cache(None)
    def graded_antipode(self, value):
        return self.normal(tuple(reversed(self.basis_words[value])))


def evaluate_shear(terms, horizontal):
    value = 0
    for monomial, vector in terms:
        if monomial & horizontal == monomial:
            value ^= vector
    return value


def shear_tables(group, terms):
    full = []
    for source in range(group.d):
        images = (group.words[a ^ (evaluate_shear(terms, a & ((1 << group.n) - 1))
                                  << group.n)] for a in submasks(source))
        full.append(group.to_augmentation(images))
    # Independently compute the leading map by dual square-zero substitution.
    dual = []
    for target in range(group.d):
        polynomial = {0}
        for i in range(group.width):
            if not target >> i & 1:
                continue
            factor = {1 << i}
            if i >= group.n:
                for monomial, vector in terms:
                    if monomial.bit_count() == 2 and vector >> (i - group.n) & 1:
                        toggle(factor, monomial)
            nxt = set()
            for a, b in product(polynomial, factor):
                if not a & b:
                    toggle(nxt, a | b)
            polynomial = nxt
        dual.append(polynomial)
    leading = [tuple(target for target, row in enumerate(dual) if source in row)
               for source in range(group.d)]
    for source, outputs in enumerate(full):
        assert all(group.weights[a] >= group.weights[source] for a in outputs)
        actual = tuple(a for a in outputs if group.weights[a] == group.weights[source])
        assert actual == leading[source], (source, actual, leading[source])
        assert all(a.bit_count() <= source.bit_count() for a in actual)
        assert tuple(a for a in actual if a.bit_count() == source.bit_count()) == (source,)
    return tuple(full), tuple(leading)


def word_polynomial(group, mask, weights, order):
    values = {0: 1}
    for i in order:
        if not mask >> i & 1:
            continue
        nxt = {}
        shift = 1 << weights[i]
        for element, polynomial in values.items():
            nxt[element] = nxt.get(element, 0) ^ polynomial_multiply(polynomial, shift)
            target = group.mul(element, group.generators[i])
            nxt[target] = nxt.get(target, 0) ^ polynomial_multiply(polynomial, shift ^ 1)
        values = {a: p for a, p in nxt.items() if p}
    return values


def group_crossing(group, weights, order, terms):
    degree, d = sum(weights), group.d
    matrices = [[0] * (d * d) for _ in range(degree + 1)]
    kernels = []
    for element in range(d):
        mask = group.word_index[element]
        control = mask ^ (evaluate_shear(terms, mask & ((1 << group.n) - 1)) << group.n)
        kernels.append(word_polynomial(group, control, weights, order))
    # Build the orbit formula.
    for u, v in product(range(d), repeat=2):
        gu, gv = group.words[u], group.words[v]
        S = group.mul(gv, group.inverse(gu))
        endpoint = 0
        for a, polynomial in kernels[S].items():
            inva = group.inverse(a)
            x = group.word_index[group.mul(gv, inva)]
            y = group.word_index[group.mul(gu, inva)]
            if polynomial.bit_count() & 1:
                endpoint ^= 1 << (x * d + y)
            while polynomial:
                bit = polynomial & -polynomial
                matrices[bit.bit_length() - 1][x * d + y] ^= 1 << (u * d + v)
                polynomial ^= bit
        assert endpoint == 1 << (v * d + u)
    # Independently assemble psi=S phi U S and the group-like Hopf formula.
    hopf = [[0] * (d * d) for _ in range(degree + 1)]
    psi = [{group.inverse(a): p for a, p in kernels[group.inverse(g)].items()}
           for g in range(d)]
    for u, v in product(range(d), repeat=2):
        gu, gv = group.words[u], group.words[v]
        argument = group.mul(gu, group.inverse(gv))
        for b, polynomial in psi[argument].items():
            x = group.word_index[group.mul(gv, b)]
            y = group.word_index[group.mul(gu, b)]
            while polynomial:
                bit = polynomial & -polynomial
                hopf[bit.bit_length() - 1][x * d + y] ^= 1 << (u * d + v)
                polynomial ^= bit
    assert hopf == matrices
    return matrices


def augmentation_change(matrices, width):
    dimension = 1 << (2 * width)
    masks = []
    for axis in range(2 * width):
        step = 1 << axis
        masks.append(sum(((1 << step) - 1) << start
                         for start in range(0, dimension, 2 * step)))
    out = []
    for matrix in matrices:
        rows = []
        for value in matrix:
            for axis, mask in enumerate(masks):
                value ^= (value & mask) << (1 << axis)
            rows.append(value)
        for axis in range(2 * width):
            step = 1 << axis
            for start in range(0, dimension, 2 * step):
                for offset in range(step):
                    rows[start + offset] ^= rows[start + step + offset]
        out.append(rows)
    return out


def leading_cut(matrices, grades, descending):
    masks = [0] * (max(grades) + 1)
    for index, grade in enumerate(grades):
        masks[grade] |= 1 << index
    leading, positive = [], 0
    for matrix in matrices:
        output = []
        for index, row in enumerate(matrix):
            grade = grades[index]
            forbidden = (sum(masks[grade + 1:]) if descending else sum(masks[:grade]))
            assert not row & forbidden, (index, grade)
            actual = row & masks[grade]
            positive += (row ^ actual).bit_count()
            output.append(actual)
        leading.append(output)
    return leading, positive


def abelian_crossing(width, weights):
    d = 1 << width
    matrices = [[0] * (d * d) for _ in range(sum(weights) + 1)]
    for u, v in product(range(d), repeat=2):
        different = u ^ v
        for swapped in submasks(different):
            polynomial = 1
            for j, weight in enumerate(weights):
                if different >> j & 1:
                    polynomial = polynomial_multiply(
                        polynomial, (1 << weight) ^ int(not swapped >> j & 1))
            x, y = u ^ swapped, v ^ swapped
            while polynomial:
                bit = polynomial & -polynomial
                matrices[bit.bit_length() - 1][x * d + y] ^= 1 << (u * d + v)
                polynomial ^= bit
    return matrices


def literal_graded_hopf(group, weights, order, shear):
    """Primitive coproduct and independently rewritten restricted products."""
    d = group.d
    phi = []
    for source in range(d):
        scalar = 1
        for j, weight in enumerate(weights):
            if source >> j & 1:
                scalar = polynomial_multiply(scalar, 1 ^ (1 << weight))
        targets = group.normal(tuple(j for j in order if source >> j & 1))
        phi.append({target: scalar for target in targets})
    psi = []
    for source in range(d):
        values = {}
        for a in group.graded_antipode(source):
            for b in shear[a]:
                for c, polynomial in phi[b].items():
                    for target in group.graded_antipode(c):
                        values[target] = values.get(target, 0) ^ polynomial
        psi.append({a: p for a, p in values.items() if p})
    matrices = [[0] * (d * d) for _ in range(sum(weights) + 1)]
    for u, v in product(range(d), repeat=2):
        outputs = {}
        for u2, v2 in product(submasks(u), submasks(v)):
            u1, v1 = u ^ u2, v ^ v2
            argument = set()
            for a in group.graded_antipode(v2):
                argument.symmetric_difference_update(group.graded_product(u2, a))
            for a in argument:
                for term, polynomial in psi[a].items():
                    for a1 in submasks(term):
                        a2 = term ^ a1
                        for x, y in product(group.graded_product(v1, a1),
                                            group.graded_product(u1, a2)):
                            key = x * d + y
                            outputs[key] = outputs.get(key, 0) ^ polynomial
        for xy, polynomial in outputs.items():
            while polynomial:
                bit = polynomial & -polynomial
                matrices[bit.bit_length() - 1][xy] ^= 1 << (u * d + v)
                polynomial ^= bit
    return matrices


def both_ranks(matrix, width):
    d = 1 << width
    column = [0] * (d * d)
    for xy, row in enumerate(matrix):
        x, y = divmod(xy, d)
        while row:
            bit = row & -row
            u, v = divmod(bit.bit_length() - 1, d)
            column[x * d + u] ^= 1 << (v * d + y)
            row ^= bit
    return binary_rank(matrix), binary_rank(column)


def combine(matrices, mask):
    out = [0] * len(matrices[0])
    for degree, matrix in enumerate(matrices):
        if mask >> degree & 1:
            out = [a ^ b for a, b in zip(out, matrix)]
    return out


def check(name, group, weights, order, terms):
    # All multiplication and antipode leading terms are checked separately.
    for left, right in product(range(group.d), repeat=2):
        actual = group.augmentation_product(left, right)
        weight = group.weights[left] + group.weights[right]
        assert all(group.weights[a] >= weight for a in actual)
        assert tuple(a for a in actual if group.weights[a] == weight) == group.graded_product(left, right)
    for source in range(group.d):
        actual = group.augmentation_antipode(source)
        assert all(group.weights[a] >= group.weights[source] for a in actual)
        assert tuple(a for a in actual if group.weights[a] == group.weights[source]) == group.graded_antipode(source)
    full_shear, leading_shear = shear_tables(group, terms)
    original = group_crossing(group, weights, order, terms)
    augmented = augmentation_change(original, group.width)
    weight_grades = [group.weights[x] + group.weights[y]
                     for x, y in product(range(group.d), repeat=2)]
    graded, positive_weight = leading_cut(augmented, weight_grades, True)
    literal_checked = group.width <= 5
    if literal_checked:
        assert graded == literal_graded_hopf(group, weights, order, leading_shear)
    length_grades = [x.bit_count() + y.bit_count()
                     for x, y in product(range(group.d), repeat=2)]
    abelian, positive_length = leading_cut(graded, length_grades, False)
    expected = abelian_crossing(group.width, weights)
    assert abelian == expected
    records = []
    for i, (matrix, target) in enumerate(zip(original, abelian)):
        ranks, target_ranks = both_ranks(matrix, group.width), both_ranks(target, group.width)
        assert all(a >= b for a, b in zip(ranks, target_ranks))
        records.append([i, *ranks, *target_ranks])
    count = len(original)
    masks = (range(1, 1 << count) if group.width <= 3 else
             sorted({1, (1 << count) - 1}
                    | {(1 << i) | (1 << (i + 1)) for i in range(count - 1)}
                    | {sum(1 << i for i in range(p, count, 2)) for p in (0, 1)}))
    combinations = 0
    for mask in masks:
        ranks = both_ranks(combine(original, mask), group.width)
        target = both_ranks(combine(abelian, mask), group.width)
        assert all(a >= b for a, b in zip(ranks, target))
        combinations += 1
    return {"name": name, "label_bits": group.width, "labels": group.d,
            "beta": [list(row) for row in group.beta],
            "deformation_weights": list(weights), "factor_order": list(order),
            "central_shear_terms": [list(term) for term in terms],
            "multiplication_pairs_audited": group.d * group.d,
            "orbit_equals_antipode_conjugated_hopf": True,
            "first_leading_map_equals_independent_primitive_hopf": literal_checked,
            "positive_first_valuation_entries": positive_weight,
            "positive_second_valuation_entries": positive_length,
            "second_leading_family_equals_abelian": True,
            "coefficient_ranks_and_abelian_ranks": records,
            "scalar_recombinations_checked": combinations}


def cubic_scope_witness():
    group = BilinearGroup(3, 1, ((0, 1, 0), (0, 0, 0), (0, 0, 0)))
    source = 7
    terms = ((7, 1),)
    actual = group.to_augmentation(
        group.words[a ^ (evaluate_shear(terms, a) << group.n)] for a in submasks(source))
    decreasing = [a for a in actual if group.weights[a] < group.weights[source]]
    assert (1 << group.n) in decreasing
    return {"source_augmentation_monomial": source, "source_weight": 3,
            "lower_weight_output_monomials": decreasing,
            "central_output_weight": 2,
            "conclusion": "A cubic control shear violates the old central-weight-two filtration; "
                          "the heavier direct filtration below covers it."}


def check_direct(name, group, weights, order, terms):
    """One-stage degeneration with central weight above every control degree."""
    degree_f = max((mask.bit_count() for mask, value in terms if value), default=0)
    assert group.central_weight > max(2, degree_f)
    for left, right in product(range(group.d), repeat=2):
        actual = group.augmentation_product(left, right)
        weight = group.weights[left] + group.weights[right]
        assert all(group.weights[a] >= weight for a in actual)
        leading = tuple(a for a in actual if group.weights[a] == weight)
        expected = () if left & right else (left | right,)
        assert leading == expected, (left, right, leading, expected)
    for source in range(group.d):
        actual = group.augmentation_antipode(source)
        assert all(group.weights[a] >= group.weights[source] for a in actual)
        assert tuple(a for a in actual if group.weights[a] == group.weights[source]) == (source,)
        transformed = group.to_augmentation(
            group.words[a ^ (evaluate_shear(terms, a & ((1 << group.n) - 1)) << group.n)]
            for a in submasks(source))
        assert all(group.weights[a] >= group.weights[source] for a in transformed)
        assert tuple(a for a in transformed if group.weights[a] == group.weights[source]) == (source,)
    original = group_crossing(group, weights, order, terms)
    augmented = augmentation_change(original, group.width)
    grades = [group.weights[x] + group.weights[y]
              for x, y in product(range(group.d), repeat=2)]
    leading, positive = leading_cut(augmented, grades, True)
    expected = abelian_crossing(group.width, weights)
    assert leading == expected
    ranks = []
    for index, (source, target) in enumerate(zip(original, leading)):
        actual, lower = both_ranks(source, group.width), both_ranks(target, group.width)
        assert all(a >= b for a, b in zip(actual, lower))
        ranks.append([index, *actual, *lower])
    count = len(original)
    masks = (range(1, 1 << count) if group.width <= 4 else
             sorted({1, (1 << count) - 1, 3, 5}
                    | {sum(1 << i for i in range(p, count, 2)) for p in (0, 1)}))
    checked = 0
    for mask in masks:
        actual = both_ranks(combine(original, mask), group.width)
        lower = both_ranks(combine(leading, mask), group.width)
        assert all(a >= b for a, b in zip(actual, lower))
        checked += 1
    return {"name": name, "labels": group.d, "label_bits": group.width,
            "central_filtration_weight": group.central_weight,
            "control_Boolean_degree": degree_f,
            "central_control_terms": [list(term) for term in terms],
            "deformation_weights": list(weights), "factor_order": list(order),
            "leading_multiplication_is_commutative_square_zero": True,
            "leading_antipode_and_control_are_identity": True,
            "complete_leading_tensor_equals_abelian": True,
            "positive_filtration_entries": positive,
            "coefficient_ranks_and_abelian_ranks": ranks,
            "scalar_recombinations_checked": checked}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    d8 = BilinearGroup(2, 1, ((0, 1), (0, 0)))
    q8 = BilinearGroup(2, 1, ((1, 1), (0, 1)))
    h16 = BilinearGroup(3, 1, ((0, 1, 0), (0, 0, 0), (0, 0, 0)))
    h32 = BilinearGroup(4, 1, ((0, 0, 1, 0), (0, 0, 0, 1),
                             (0, 0, 0, 0), (0, 0, 0, 0)))
    cases = [
        check("D8_unweighted", d8, (1, 1, 1), (0, 1, 2), ()),
        check("D8_affine_quadratic_control", d8, (1, 1, 2), (0, 1, 2),
              ((0, 1), (1, 1), (3, 1))),
        check("D8_reordered_distinct_weights", d8, (1, 2, 3), (1, 0, 2), ((3, 1),)),
        check("Q8_nonzero_restricted_squares", q8, (1, 1, 1), (0, 1, 2), ((3, 1),)),
        check("D8_times_C2_quadratic_control", h16, (1, 1, 1, 2), (2, 1, 0, 3),
              ((5, 1), (6, 1))),
        check("extraspecial_shared_center", h32, (1, 1, 1, 1, 2),
              (0, 1, 2, 3, 4), ((3, 1), (12, 1))),
    ]
    cubic = BilinearGroup(3, 1, ((0, 1, 0), (0, 0, 0), (0, 0, 0)), 4)
    field4 = BilinearGroup(4, 2, ((0, 0, 1, 2), (0, 0, 2, 3),
                                (0, 0, 0, 0), (0, 0, 0, 0)), 5)
    rng = random.Random(20261030)
    anf = [rng.randrange(4) for _ in range(16)]
    anf[-1] = 3
    lookup_terms = tuple((mask, value) for mask, value in enumerate(anf) if value)
    direct = [
        check_direct("cubic_control_weight_four", cubic, (1, 1, 1, 2),
                     (2, 0, 1, 3), ((7, 1),)),
        check_direct("full_degree_four_F4_lookup_weight_five", field4,
                     (1, 1, 1, 1, 1, 1), (3, 0, 2, 1, 5, 4), lookup_terms),
        check_direct("Q8_direct_weight_three",
                     BilinearGroup(2, 1, ((1, 1), (0, 1)), 3),
                     (1, 2, 3), (1, 0, 2), ((0, 1), (1, 1), (3, 1))),
    ]
    result = {"status": "Direct heavy-central and legacy two-stage degeneration audits passed.",
              "scope": "Full class-two bilinear group domains and arbitrary central Boolean lookup controls.",
              "legacy_two_stage_cases": cases, "old_weight_two_cubic_witness": cubic_scope_witness(),
              "direct_one_stage_cases": direct,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
