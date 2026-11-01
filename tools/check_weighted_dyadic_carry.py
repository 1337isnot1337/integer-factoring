#!/usr/bin/env python3
"""Exact padded-carry, weighted Fourier, and common-quotient controls."""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path
import random

from check_carry_fourier_elimination import construct_field_and_root
from check_frobenius_sparse_recovery import TensorAlgebra, construct_field, degree_family


def digits(value, base, length):
    result = []
    for _ in range(length):
        result.append(value % base)
        value //= base
    assert value == 0
    return result


def convolution(left, right, length):
    result = [0] * length
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def genuine_carries(left, right, base, target_digits):
    coefficients = convolution(left, right, len(target_digits))
    carry, internal = 0, []
    for j, (coefficient, target) in enumerate(zip(coefficients, target_digits)):
        assert (coefficient + carry) % base == target
        carry = (coefficient + carry) // base
        if j + 1 < len(target_digits):
            internal.append(carry)
    assert carry == 0
    return internal


def check_assignment_bijection(counts):
    for base, d, m, p_digits, q_digits in ((2, 2, 3, 1, 1),
                                         (4, 2, 5, 1, 1),
                                         (8, 2, 7, 1, 1),
                                         (4, 2, 5, 2, 1)):
        length = p_digits + q_digits
        bound, old_bound, modulus = d * base, d * (base - 1), (1 << m) + 1
        assert d * base * base < modulus
        for p, q in product(range(base ** p_digits), range(base ** q_digits)):
            left, right = digits(p, base, p_digits), digits(q, base, q_digits)
            coefficients = convolution(left, right, length)
            assert max(coefficients) <= d * (base - 1) ** 2
            for target in range(1, base ** length):
                n = digits(target, base, length)
                accepted = []
                for internal in product(range(bound), repeat=length - 1):
                    carries = (0, *internal, 0)
                    exact = all(coefficients[j] + carries[j] == n[j] + base * carries[j + 1]
                                for j in range(length))
                    modular = all((coefficients[j] + carries[j] - n[j] - base * carries[j + 1])
                                  % modulus == 0 for j in range(length))
                    assert exact == modular
                    if modular:
                        accepted.append(internal)
                    counts["padded_modular_assignments"] += 1
                assert bool(accepted) == (p * q == target)
                if accepted:
                    assert len(accepted) == 1
                    assert all(carry <= old_bound for carry in accepted[0])
                    assert list(accepted[0]) == genuine_carries(left, right, base, n)
                    counts["unique_accepted_factor_assignments"] += 1
                    counts["accepted_assignments_with_nonzero_carry"] += any(accepted[0])

    # A too-small modulus really admits a wrong integer product.
    base, modulus, p, q, target, carry = 4, 5, 3, 3, 4, 1
    n = digits(target, base, 2)
    assert (p * q - n[0] - base * carry) % modulus == 0
    assert (carry - n[1]) % modulus == 0
    assert p * q != target
    counts["insufficient_modulus_negative_controls"] += 1


def tensor_character(algebra, root, modulus, exponent):
    return algebra.root_scalar(algebra.root_field.power(root, exponent % modulus))


def check_carry_transform(counts):
    rng = random.Random(11022026)
    records = []
    for base, d, m, length, ell in ((2, 2, 3, 3, 3),
                                   (4, 2, 5, 3, 3),
                                   (8, 2, 7, 3, 3)):
        field, root, construction = construct_field_and_root(m)
        weights = construct_field(ell)
        algebra = TensorAlgebra(field, weights)
        modulus, bound, gamma = (1 << m) + 1, d * base, 2
        gamma_tensor = algebra.weight_scalar(gamma)
        n = [rng.randrange(base) for _ in range(length)]
        frequencies = [(0,) * length]
        frequencies += [tuple(rng.randrange(modulus) for _ in range(length)) for _ in range(20)]
        for frequency in frequencies:
            phase = sum(s * digit for s, digit in zip(frequency, n))
            claimed = tensor_character(algebra, root, modulus, phase)
            for j in range(1, length):
                c = base * frequency[j - 1] - frequency[j]
                y = algebra.mul(gamma_tensor, tensor_character(algebra, root, modulus, c))
                factor = algebra.geometric(y, bound)
                denominator_inverse = algebra.mul(
                    algebra.geometric(y, modulus),
                    algebra.weight_scalar(weights.inverse(1 ^ weights.power(gamma, modulus))))
                inverse = algebra.mul(algebra.add(algebra.one, y),
                                      algebra.power(denominator_inverse, bound))
                assert algebra.mul(factor, inverse) == algebra.one
                claimed = algebra.mul(claimed, factor)
                counts["carry_factor_unit_inverses"] += 1
            literal = algebra.zero
            for internal in product(range(bound), repeat=length - 1):
                carries = (0, *internal, 0)
                coefficients = [(n[j] + base * carries[j + 1] - carries[j]) % modulus
                                for j in range(length)]
                character = tensor_character(algebra, root, modulus,
                                             sum(s * a for s, a in zip(frequency, coefficients)))
                weight = algebra.weight_scalar(weights.power(gamma, sum(internal)))
                literal = algebra.add(literal, algebra.mul(weight, character))
                counts["literal_weighted_carry_terms"] += 1
            assert literal == claimed
            counts["whole_carry_fourier_products"] += 1
        records.append({**construction, "base": base, "D": d,
                        "T": length, "weight_degree": ell,
                        "width_is_multiple_of_weight_degree": (base.bit_length() - 1) % ell == 0})
    return records


def check_full_fourier(counts):
    records = []
    for base, d, m, ell in ((4, 2, 5, 3), (8, 2, 7, 3)):
        field, root, _ = construct_field_and_root(m)
        weights = construct_field(ell)
        algebra = TensorAlgebra(field, weights)
        modulus, bound = (1 << m) + 1, d * base
        alpha, beta, gamma = 2, 3, 2
        characters = [1]
        for _ in range(1, modulus):
            characters.append(field.mul(characters[-1], root))
        assert len(set(characters)) == modulus
        carry_values = []
        for c in range(modulus):
            y = algebra.mul(algebra.weight_scalar(gamma), algebra.root_scalar(characters[c]))
            carry_values.append(algebra.geometric(y, bound))
        input_values = []
        for s0 in range(modulus):
            value = algebra.zero
            for p in range(base):
                y = algebra.mul(algebra.weight_scalar(beta),
                                algebra.root_scalar(characters[(-p * s0) % modulus]))
                box = algebra.geometric(y, base)
                denominator_inverse = algebra.mul(
                    algebra.geometric(y, modulus),
                    algebra.weight_scalar(weights.inverse(1 ^ weights.power(beta, modulus))))
                box_inverse = algebra.mul(algebra.add(algebra.one, y),
                                          algebra.power(denominator_inverse, base))
                assert algebra.mul(box, box_inverse) == algebra.one
                value = algebra.add(value, algebra.mul(algebra.weight_scalar(weights.power(alpha, p)), box))
                counts["coupled_input_box_units"] += 1
            input_values.append(value)
        transformed = [algebra.mul(carry_values[(base * s0 - s1) % modulus], input_values[s0])
                       for s0, s1 in product(range(modulus), repeat=2)]
        target_results = []
        for target in sorted({1, base + 1, (base - 1) ** 2, base ** 2 - 1}):
            n0, n1 = digits(target, base, 2)
            fourier = algebra.zero
            for value, (s0, s1) in zip(transformed, product(range(modulus), repeat=2)):
                phase = characters[(s0 * n0 + s1 * n1) % modulus]
                fourier = algebra.add(fourier, tuple(field.mul(phase, coefficient) for coefficient in value))
                counts["full_fourier_frequency_terms"] += 1
            direct = 0
            factor_pairs = []
            for p, q in product(range(base), repeat=2):
                if p * q == target:
                    carry = p * q // base
                    direct ^= weights.product(weights.power(alpha, p), weights.power(beta, q),
                                              weights.power(gamma, carry))
                    factor_pairs.append([p, q])
            assert fourier == algebra.weight_scalar(direct)
            counts["full_coupled_fourier_contractions"] += 1
            target_results.append({"N": target, "factor_pairs_in_digit_boxes": factor_pairs,
                                   "weighted_field_value": direct})
        records.append({"base": base, "D": d, "M": modulus, "weight_degree": ell,
                        "width_is_multiple_of_weight_degree": (base.bit_length() - 1) % ell == 0,
                        "targets": target_results})
    return records


class CyclicAlgebra:
    def __init__(self, field, modulus, dimensions):
        self.field, self.modulus, self.dimensions = field, modulus, dimensions
        self.zero_exponent = (0,) * dimensions
        self.one = {self.zero_exponent: 1}

    def add(self, left, right):
        result = dict(left)
        for exponent, coefficient in right.items():
            result[exponent] = result.get(exponent, 0) ^ coefficient
            if not result[exponent]:
                del result[exponent]
        return result

    def mul(self, left, right):
        result = {}
        for u, a in left.items():
            for v, b in right.items():
                exponent = tuple((x + y) % self.modulus for x, y in zip(u, v))
                result[exponent] = result.get(exponent, 0) ^ self.field.mul(a, b)
        return {exponent: coefficient for exponent, coefficient in result.items() if coefficient}

    def scalar(self, polynomial, value):
        return {exponent: self.field.mul(coefficient, value) for exponent, coefficient in polynomial.items()
                if self.field.mul(coefficient, value)}

    def binomial(self, scalar, exponent):
        exponent = tuple(x % self.modulus for x in exponent)
        return self.add(self.one, {exponent: scalar})

    def geometric(self, scalar, exponent, length):
        result, coefficient = {}, 1
        for r in range(length):
            key = tuple(r * x % self.modulus for x in exponent)
            result[key] = result.get(key, 0) ^ coefficient
            coefficient = self.field.mul(coefficient, scalar)
        return {key: value for key, value in result.items() if value}

    def binomial_inverse(self, scalar, exponent):
        denominator = 1 ^ self.field.power(scalar, self.modulus)
        assert denominator != 0
        inverse = self.scalar(self.geometric(scalar, exponent, self.modulus), self.field.inverse(denominator))
        assert self.mul(self.binomial(scalar, exponent), inverse) == self.one
        return inverse

    def frobenius(self, polynomial, base):
        return {tuple(base * x % self.modulus for x in exponent): self.field.power(coefficient, base)
                for exponent, coefficient in polynomial.items()}


def carry_vectors(base, length):
    result = []
    for j in range(1, length):
        exponent = [0] * length
        exponent[j - 1], exponent[j] = base, -1
        result.append(tuple(exponent))
    return result


def check_common_quotient(counts):
    records = []
    for base, d, modulus, length, p in ((2, 2, 9, 3, (1, 1)),
                                       (4, 2, 33, 2, (0,)),
                                       (4, 2, 33, 2, (3,)),
                                       (8, 2, 129, 2, (1,)),
                                       (8, 2, 129, 2, (3,)),
                                       (8, 2, 129, 2, (7,))):
        field = construct_field(3)
        algebra = CyclicAlgebra(field, modulus, length)
        beta, gamma = 3, 2
        carry, low, denominator, inverse = (dict(algebra.one) for _ in range(4))
        for vector in carry_vectors(base, length):
            carry = algebra.mul(carry, algebra.geometric(gamma, vector, d * base))
            low = algebra.mul(low, algebra.geometric(gamma, vector, d))
            high_vector = tuple(d * x for x in vector)
            high_scalar = field.power(gamma, d)
            denominator = algebra.mul(denominator, algebra.binomial(high_scalar, high_vector))
            inverse = algebra.mul(inverse, algebra.binomial_inverse(high_scalar, high_vector))
        input_vector = tuple(-x for x in p) + (0,) * (length - len(p))
        qbox = algebra.geometric(beta, input_vector, base)
        denominator = algebra.mul(denominator, algebra.binomial(beta, input_vector))
        inverse = algebra.mul(inverse, algebra.binomial_inverse(beta, input_vector))
        assert algebra.mul(denominator, inverse) == algebra.one
        left = algebra.mul(carry, qbox)
        right = algebra.mul(low, algebra.mul(algebra.frobenius(denominator, base), inverse))
        assert left == right
        assert len(low) <= d ** (length - 1)
        if (base.bit_length() - 1) % 3 == 0:
            substituted = {tuple(base * x % modulus for x in exponent): coefficient
                           for exponent, coefficient in denominator.items()}
            assert substituted == algebra.frobenius(denominator, base)
            counts["recovery_field_fixed_frobenius_quotients"] += 1
        counts["exact_cyclic_common_quotients"] += 1
        records.append({"base": base, "D": d, "M": modulus, "T": length,
                        "first_operand_digits": list(p), "denominator_terms": len(denominator),
                        "inverse_terms": len(inverse), "low_carry_terms": len(low),
                        "joint_box_terms": len(left)})

    # Recover every target coefficient from a full first-operand sum.
    for base, modulus in ((4, 33), (8, 129)):
        field, length, d = construct_field(3), 2, 2
        algebra = CyclicAlgebra(field, modulus, length)
        alpha, beta, gamma = 2, 3, 2
        carry = algebra.geometric(gamma, carry_vectors(base, length)[0], d * base)
        whole = {}
        for p in range(base):
            qbox = algebra.geometric(beta, (-p, 0), base)
            whole = algebra.add(whole, algebra.scalar(algebra.mul(carry, qbox), field.power(alpha, p)))
        for target in range(1, base ** 2):
            target_digits = digits(target, base, length)
            coefficient = whole.get(tuple(-x % modulus for x in target_digits), 0)
            expected = 0
            for p, q in product(range(base), repeat=2):
                if p * q == target:
                    expected ^= field.product(field.power(alpha, p), field.power(beta, q),
                                              field.power(gamma, p * q // base))
            assert coefficient == expected
            counts["cyclic_weighted_target_coefficients"] += 1
    return records


def check_weighted_prefix_recovery(counts):
    fields = {degree: construct_field(degree) for degree in (3, 5, 7)}
    records = []
    for target in (9, 25, 49, 121):
        n, d = target.bit_length(), 2
        divisors = [value for value in range(1, target + 1) if target % value == 0]
        family = degree_family(n, len(divisors))
        assert set(family) <= set(fields)
        terms = {}
        carry_weights = []
        for ell in family:
            width0 = (n + d - 1) // d
            width = ell * ((width0 + ell - 1) // ell)
            base, field, gamma = 1 << width, fields[ell], 2
            degree_terms = []
            degree_record = []
            for divisor in divisors:
                carries = genuine_carries(digits(divisor, base, d),
                                           digits(target // divisor, base, d), base,
                                           digits(target, base, 2 * d))
                rho = sum(carries)
                coefficient = field.power(gamma, rho)
                assert coefficient != 0
                if divisor == 1:
                    assert rho == 0 and coefficient == 1
                else:
                    degree_terms.append((coefficient, divisor, target // divisor))
                degree_record.append({"divisor": divisor, "carry_sum": rho,
                                      "nonzero_coefficient": coefficient})
            terms[ell] = degree_terms
            carry_weights.append({"degree": ell, "digit_width": width, "weights": degree_record})

        def nonzero_prefix(k, prefix):
            for ell in family:
                field = fields[ell]
                selected = [(c, u, v) for c, u, v in terms[ell] if u >> k == prefix]
                for alpha, beta in product(range(2, field.size), repeat=2):
                    value = 0
                    for coefficient, u, v in selected:
                        value ^= field.product(coefficient, field.power(alpha, u), field.power(beta, v))
                    counts["weighted_prefix_field_evaluations"] += 1
                    if value:
                        return True
            return False

        prefix = 0
        for k in range(n - 1, -1, -1):
            prefix *= 2
            if not nonzero_prefix(k, prefix):
                prefix += 1
        assert prefix == min(value for value in divisors if value > 1)
        counts["degree_dependent_carry_weight_recoveries"] += 1
        records.append({"N": target, "recovered_smallest_prime": prefix,
                        "degree_dependent_carry_weights": carry_weights})
    return records


def check_aligned_common_denominator(counts):
    records = []
    for base, d, modulus, p in ((8, 2, 129, 3), (4, 8, 129, 1)):
        field = construct_field(3)
        algebra = CyclicAlgebra(field, modulus, 2)
        beta = 3
        t = d.bit_length() - 1
        inverse_frobenius_exponent = 1 << ((-t) % field.degree)
        gamma = field.power(beta, inverse_frobenius_exponent)
        assert gamma != 1 and field.power(gamma, d) == beta
        v = carry_vectors(base, 2)[0]
        columns = [tuple(d * x for x in v), (-p, 0)]
        denominator, full_cycle_product = dict(algebra.one), dict(algebra.one)
        for column in columns:
            denominator = algebra.mul(denominator, algebra.binomial(beta, column))
            full_cycle_product = algebra.mul(full_cycle_product,
                                             algebra.geometric(beta, column, modulus))
        scalar = field.power(1 ^ field.power(beta, modulus), len(columns))
        assert scalar != 0
        carry = algebra.geometric(gamma, v, d * base)
        qbox = algebra.geometric(beta, (-p, 0), base)
        left = algebra.scalar(algebra.mul(carry, qbox), scalar)
        low = algebra.geometric(gamma, v, d)
        right = algebra.mul(low, algebra.mul(algebra.frobenius(denominator, base), full_cycle_product))
        assert left == right
        counts["aligned_common_scalar_denominators"] += 1
        records.append({"base": base, "D": d, "M": modulus, "beta": beta, "gamma": gamma,
                        "inverse_frobenius_exponent": inverse_frobenius_exponent,
                        "clearing_scalar": scalar, "cleared_joint_terms": len(left),
                        "full_cycle_product_terms": len(full_cycle_product)})
    return records


def check_aligned_univariate_recovery(counts):
    fields = {degree: construct_field(degree) for degree in (3, 5, 7)}
    records = []
    for target in (9, 25, 49, 121):
        n, d, t = target.bit_length(), 2, 1
        divisors = [value for value in range(1, target + 1) if target % value == 0]
        family = degree_family(n, len(divisors))
        assert set(family) <= set(fields)
        terms = {}
        coefficient_records = []
        for ell in family:
            width0 = (n + d - 1) // d
            width = ell * ((width0 + ell - 1) // ell)
            m = d * ((2 * width + t + d - 1) // d)
            base, modulus, field, beta = 1 << width, (1 << m) + 1, fields[ell], 3
            gamma = field.power(beta, 1 << ((-t) % ell))
            assert field.power(gamma, d) == beta and gamma != 1
            degree_terms = []
            coefficient_record = []
            for divisor in divisors:
                rho = sum(genuine_carries(digits(divisor, base, d),
                                           digits(target // divisor, base, d), base,
                                           digits(target, base, 2 * d)))
                coefficient = field.mul(field.power(beta, target // divisor), field.power(gamma, rho))
                assert coefficient != 0
                q_digit_count = (n - divisor.bit_length() + 1 + width - 1) // width
                nu = 2 * d - 1 + q_digit_count
                clearing_scalar = field.power(1 ^ field.power(beta, modulus), nu)
                assert clearing_scalar != 0
                cleared = field.mul(coefficient, clearing_scalar)
                assert cleared != 0
                if divisor == 1:
                    assert rho == 0 and coefficient == field.power(beta, target)
                else:
                    degree_terms.append((divisor, coefficient, cleared))
                coefficient_record.append({"divisor": divisor, "carry_sum": rho,
                                           "coefficient": coefficient, "clearing_scalar": clearing_scalar,
                                           "cleared_coefficient": cleared})
            terms[ell] = degree_terms
            coefficient_records.append({"degree": ell, "beta": beta, "gamma": gamma,
                                        "digit_width": width, "coefficients": coefficient_record})

        recovered = []
        for coefficient_index in (1, 2):
            def nonzero_prefix(k, prefix):
                counts["aligned_univariate_prefix_queries"] += 1
                for ell in family:
                    field = fields[ell]
                    selected = [term for term in terms[ell] if term[0] >> k == prefix]
                    for alpha in range(2, field.size):
                        value = 0
                        for term in selected:
                            value ^= field.mul(term[coefficient_index], field.power(alpha, term[0]))
                        counts["aligned_univariate_field_evaluations"] += 1
                        if value:
                            return True
                return False

            prefix = 0
            for k in range(n - 1, -1, -1):
                prefix *= 2
                if not nonzero_prefix(k, prefix):
                    prefix += 1
            assert prefix == min(value for value in divisors if value > 1)
            recovered.append(prefix)
            counts["aligned_univariate_recoveries"] += 1
            if coefficient_index == 2:
                counts["instance_rescaled_univariate_recoveries"] += 1
        records.append({"N": target, "recovered_smallest_prime": recovered[0],
                        "recovered_after_instance_rescaling": recovered[1],
                        "degree_dependent_coefficients": coefficient_records})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    counts = Counter()
    check_assignment_bijection(counts)
    carries = check_carry_transform(counts)
    fourier = check_full_fourier(counts)
    quotients = check_common_quotient(counts)
    recovery = check_weighted_prefix_recovery(counts)
    common_scalars = check_aligned_common_denominator(counts)
    aligned_recovery = check_aligned_univariate_recovery(counts)
    result = {"status": "all exact controls passed", "score_improvement_claimed": False,
              "summary": dict(sorted(counts.items())), "carry_transform_models": carries,
              "coupled_fourier_models": fourier, "cyclic_quotient_models": quotients,
              "weighted_recovery_models": recovery,
              "aligned_scalar_denominator_models": common_scalars,
              "aligned_univariate_recovery_models": aligned_recovery}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "summary": result["summary"]}, indent=2))


if __name__ == "__main__":
    main()
