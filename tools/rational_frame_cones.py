#!/usr/bin/env python3
"""Fixed rational functions expanded in a common positive Laurent cone.

Polynomials are dictionaries from exponent tuples to prime-field scalars.
The distinguished first variable is x. Every denominator must have a
nonzero x=0 part. Families and their supports are fixed before word length
is varied; this module makes no claim for arbitrary formal-series frames.
"""


def clean(polynomial, prime):
    return {tuple(exponent): value % prime for exponent, value in polynomial.items()
            if value % prime}


def add(left, right, prime):
    result = left.copy()
    for exponent, value in right.items():
        result[exponent] = (result.get(exponent, 0) + value) % prime
    return {exponent: value for exponent, value in result.items() if value}


def multiply(left, right, prime, weights=None, upper_weight=None):
    result = {}
    for a, x in left.items():
        for b, y in right.items():
            exponent = tuple(u + v for u, v in zip(a, b))
            if upper_weight is not None and sum(w * e for w, e in zip(weights, exponent)) > upper_weight:
                continue
            result[exponent] = (result.get(exponent, 0) + x * y) % prime
    return {exponent: value for exponent, value in result.items() if value}


class FixedRationalCone:
    def __init__(self, prime, variables, functions):
        if variables < 1 or prime < 2:
            raise ValueError('Require positive dimension and a prime coefficient modulus')
        self.prime = prime
        self.variables = variables
        self.zero = (0,) * variables
        self.functions = [(clean(p, prime), clean(q, prime)) for p, q in functions]
        for numerator, denominator in self.functions:
            if not denominator:
                raise ValueError('Zero denominator')
            for exponent in tuple(numerator) + tuple(denominator):
                if len(exponent) != variables or min(exponent) < 0:
                    raise ValueError('Input numerator and denominator must be ordinary polynomials')
            if not any(exponent[0] == 0 for exponent in denominator):
                raise ValueError('Denominators must be x-adic units after cancellation')
        bound = max((max(exponent[1:], default=0)
                     for _, denominator in self.functions for exponent in denominator), default=0)
        base = max(2, bound + 1)
        other_weights = tuple(base ** i for i in range(variables - 1))
        x_weight = 1 + max((sum(w * e for w, e in zip(other_weights, exponent[1:]))
                            for _, denominator in self.functions for exponent in denominator), default=0)
        self.weights = (x_weight,) + other_weights
        self.normalized = []
        generators = set()
        all_shifts = []
        for numerator, denominator in self.functions:
            leading = min(denominator, key=self.weight)
            assert leading[0] == 0
            assert sum(self.weight(exponent) == self.weight(leading) for exponent in denominator) == 1
            inverse = pow(denominator[leading], -1, prime)
            shifted = {tuple(a - b for a, b in zip(exponent, leading)): value * inverse % prime
                       for exponent, value in numerator.items()}
            correction = {}
            for exponent, value in denominator.items():
                if exponent == leading:
                    continue
                difference = tuple(a - b for a, b in zip(exponent, leading))
                assert difference[0] >= 0 and self.weight(difference) > 0
                correction[difference] = value * inverse % prime
                generators.add(difference)
            self.normalized.append((shifted, correction))
            all_shifts.extend(shifted)
        self.generators = tuple(sorted(generators))
        self.shifts = tuple(sorted(set(all_shifts)))
        self._inverse_cache = {}
        self._product_cache = {}

    def weight(self, exponent):
        return sum(w * e for w, e in zip(self.weights, exponent))

    def lower_weight(self, indices):
        total = 0
        for index in indices:
            shifted, _ = self.normalized[index]
            if not shifted:
                return None
            total += min(map(self.weight, shifted))
        return total

    def inverse_unit_series(self, index, budget):
        key = (index, budget)
        if key in self._inverse_cache:
            return self._inverse_cache[key]
        if budget < 0:
            return {}
        _, correction = self.normalized[index]
        result = {self.zero: 1}
        term = result
        if correction:
            negative = {exponent: -value % self.prime for exponent, value in correction.items()}
            minimum = min(map(self.weight, negative))
            for _ in range(budget // minimum):
                term = multiply(term, negative, self.prime, self.weights, budget)
                if not term:
                    break
                result = add(result, term, self.prime)
        self._inverse_cache[key] = result
        return result

    def product_series(self, indices, upper_weight):
        """All coefficients of a product with weight at most upper_weight."""
        indices = tuple(sorted(indices))
        key = (indices, upper_weight)
        if key in self._product_cache:
            return self._product_cache[key]
        lower = self.lower_weight(indices)
        if lower is None or lower > upper_weight:
            return {}
        budget = upper_weight - lower
        numerator = {self.zero: 1}
        unit = {self.zero: 1}
        for index in indices:
            shifted, _ = self.normalized[index]
            numerator = multiply(numerator, shifted, self.prime)
            unit = multiply(unit, self.inverse_unit_series(index, budget),
                            self.prime, self.weights, budget)
        result = multiply(numerator, unit, self.prime, self.weights, upper_weight)
        self._product_cache[key] = result
        return result

    def coefficient(self, indices, exponent):
        return self.product_series(indices, self.weight(exponent)).get(tuple(exponent), 0)
