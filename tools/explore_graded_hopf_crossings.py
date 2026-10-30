#!/usr/bin/env python3
"""Exploratory exact crossings from graded restricted Heisenberg Hopf algebras.

PBW generators x_i,y_i,z_j square to zero, z_j are central, and
[y_i,x_i]=sum_j L[j,i] z_j. All other pairs commute. Degrees are 1,1,2.
This is an exploratory construction, not an improved factoring theorem.
"""

from functools import lru_cache
import argparse
import json
from pathlib import Path

from check_noncommutative_crossing_orbit import binary_rank


def submasks(value):
    sub = value
    while True:
        yield sub
        if not sub:
            return
        sub = (sub - 1) & value


def add_support(target: set[int], values):
    target.symmetric_difference_update(values)


class RestrictedHeisenberg:
    def __init__(self, horizontal: int, center: int, columns: tuple[int, ...]):
        self.m, self.t = horizontal, center
        if len(columns) != horizontal or any(x <= 0 or x >= 1 << center for x in columns):
            raise ValueError("one nonzero center vector is required for each pair")
        self.columns = columns
        self.width = 2 * horizontal + center
        self.d = 1 << self.width
        self.hmask = (1 << horizontal) - 1
        self.cshift = 2 * horizontal

    def weight(self, value):
        return (value & ((1 << self.cshift) - 1)).bit_count() + 2 * (value >> self.cshift).bit_count()

    @lru_cache(None)
    def central_product(self, selected: int):
        result = {0}
        for i in range(self.m):
            if selected >> i & 1:
                nxt = set()
                for monomial in result:
                    for j in range(self.t):
                        bit = 1 << j
                        if self.columns[i] & bit and not monomial & bit:
                            value = monomial | bit
                            if value in nxt:
                                nxt.remove(value)
                            else:
                                nxt.add(value)
                result = nxt
        return tuple(sorted(result))

    @lru_cache(None)
    def multiply(self, left: int, right: int):
        a, b, c = left & self.hmask, (left >> self.m) & self.hmask, left >> self.cshift
        A, B, C = right & self.hmask, (right >> self.m) & self.hmask, right >> self.cshift
        if c & C:
            return ()
        out = set()
        for selected in submasks(b & A):
            aa, bb = A ^ selected, b ^ selected
            if a & aa or bb & B:
                continue
            xy = (a | aa) | ((bb | B) << self.m)
            for z in self.central_product(selected):
                if not z & (c | C):
                    term = xy | ((c | C | z) << self.cshift)
                    if term in out:
                        out.remove(term)
                    else:
                        out.add(term)
        return tuple(sorted(out))

    @lru_cache(None)
    def antipode(self, value: int):
        a, b, c = value & self.hmask, (value >> self.m) & self.hmask, value >> self.cshift
        return self.multiply((b << self.m) | (c << self.cshift), a)

    def product_support(self, left, right):
        result = set()
        for a in left:
            for b in right:
                add_support(result, self.multiply(a, b))
        return result

    def check_algebra(self):
        # Full associativity for dimension <= 32, plus antipode/coproduct
        # identities. This independently checks the normal-order construction.
        for a in range(self.d):
            assert self.multiply(a, 0) == self.multiply(0, a) == (a,)
            assert self.product_support(self.antipode(a), (0,)) == set(self.antipode(a))
            convolution = set()
            for left in submasks(a):
                add_support(convolution, self.product_support((left,), self.antipode(a ^ left)))
            assert convolution == ({0} if a == 0 else set())
            assert self.product_support(self.antipode(a), (0,)) == set(self.antipode(a))
        if self.d <= 32:
            for a in range(self.d):
                for b in range(self.d):
                    for c in range(self.d):
                        assert self.product_support(self.multiply(a, b), (c,)) == self.product_support((a,), self.multiply(b, c))

    def coefficient_matrices(self):
        d = self.d
        maximum_degree = 2 * (2 * self.m + 2 * self.t)
        coefficients = [[0] * (d * d) for _ in range(maximum_degree + 1)]
        # R_t(u,v)=sum (v1 a1) tensor (u1 a2), where
        # Delta(a)=Delta(phi_t(u2 S(v2))) and Delta(u),Delta(v) are used.
        # At t=0 this is exactly the swap. Substitute t=1+z.
        for u in range(d):
            for v in range(d):
                outputs = {}
                for u2 in submasks(u):
                    u1 = u ^ u2
                    for v2 in submasks(v):
                        v1 = v ^ v2
                        w = self.weight(u2) + self.weight(v2)
                        binomial = sum(1 << j for j in submasks(w))
                        a = self.product_support((u2,), self.antipode(v2))
                        for term in a:
                            assert self.weight(term) == w
                            for a1 in submasks(term):
                                a2 = term ^ a1
                                for x in self.multiply(v1, a1):
                                    for y in self.multiply(u1, a2):
                                        key = x * d + y
                                        outputs[key] = outputs.get(key, 0) ^ binomial
                for xy, polynomial in outputs.items():
                    assert (polynomial.bit_count() & 1) == int(xy == v * d + u)
                    x, y = divmod(xy, d)
                    assert self.weight(x) + self.weight(y) == self.weight(u) + self.weight(v)
                    while polynomial:
                        bit = polynomial & -polynomial
                        coefficients[bit.bit_length() - 1][xy] ^= 1 << (u * d + v)
                        polynomial ^= bit
                assert (outputs.get(v * d + u, 0).bit_count() & 1) == 1
        return coefficients

    def generator_images(self):
        """Dual coalgebra-map generator substitution after t=1+z."""
        images = []
        for side in range(2):
            for variable in range(self.width):
                if variable < 2 * self.m:
                    images.append({1 << (side * self.width + variable): 3,
                                   1 << ((1 - side) * self.width + variable): 2})
                    continue
                center = variable - 2 * self.m
                result = {1 << (side * self.width + variable): 5,
                          1 << ((1 - side) * self.width + variable): 4}
                for i, column in enumerate(self.columns):
                    if not (column >> center) & 1:
                        continue
                    xL, yL = 1 << i, 1 << (self.m + i)
                    xR, yR = xL << self.width, yL << self.width
                    terms = ((xL | yR, 3), (yL | xR, 5), (xR | yR, 6)) if side == 0 else (
                            (xL | yL, 3), (yL | xR, 6), (xR | yR, 5))
                    for monomial, polynomial in terms:
                        result[monomial] = result.get(monomial, 0) ^ polynomial
                images.append({m: p for m, p in result.items() if p})
        return images

    def substitution_matrices(self, images=None):
        from check_noncommutative_crossing_orbit import polynomial_mul
        if images is None:
            images = self.generator_images()
        d = self.d
        degree = sum(max(p.bit_length() - 1 for p in image.values()) for image in images)
        coefficients = [[0] * (d * d) for _ in range(degree + 1)]

        def visit(variable, output, polynomial):
            if variable == 2 * self.width:
                # Our basis bit convention has the first tensor factor in
                # the low bits, whereas cut rows use first*d+second.
                x, y = output & (d - 1), output >> self.width
                row = x * d + y
                for monomial, p in polynomial.items():
                    u, v = monomial & (d - 1), monomial >> self.width
                    while p:
                        bit = p & -p
                        coefficients[bit.bit_length() - 1][row] ^= 1 << (u * d + v)
                        p ^= bit
                return
            visit(variable + 1, output, polynomial)
            product = {}
            for a, p in polynomial.items():
                for b, q in images[variable].items():
                    if not a & b:
                        monomial = a | b
                        product[monomial] = product.get(monomial, 0) ^ polynomial_mul(p, q)
            visit(variable + 1, output | (1 << variable), {m: p for m, p in product.items() if p})

        visit(0, 0, {0: 1})
        return coefficients

    def rank_profile(self, coefficients):
        d = self.d
        ranks = []
        for degree, rows in enumerate(coefficients):
            if not any(rows):
                continue
            cols = [0] * (d * d)
            for xy, row in enumerate(rows):
                x, y = divmod(xy, d)
                while row:
                    bit = row & -row
                    u, v = divmod(bit.bit_length() - 1, d)
                    cols[x * d + u] ^= 1 << (y * d + v)
                    row ^= bit
            r, c = binary_rank(rows), binary_rank(cols)
            ranks.append([degree, r, c, r * c])
        return {"horizontal_pairs": self.m, "center_bits": self.t,
                "columns": list(self.columns), "labels": self.d,
                "label_bits": self.width, "coefficient_profiles": ranks,
                "maximum_rank_product": max(row[3] for row in ranks)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--horizontal", type=int, default=1)
    parser.add_argument("--center", type=int, default=1)
    parser.add_argument("--columns", default="1")
    parser.add_argument("--method", choices=("literal", "substitution"), default="substitution")
    parser.add_argument("--cross-check", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    algebra = RestrictedHeisenberg(args.horizontal, args.center, tuple(map(int, args.columns.split(','))))
    algebra.check_algebra()
    matrices = algebra.substitution_matrices() if args.method == "substitution" else algebra.coefficient_matrices()
    if args.cross_check:
        expected = algebra.coefficient_matrices() if args.method == "substitution" else algebra.substitution_matrices()
        assert matrices == expected
    result = algebra.rank_profile(matrices)
    result["literal_substitution_cross_check"] = args.cross_check
    encoded = json.dumps(result, sort_keys=True, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")
