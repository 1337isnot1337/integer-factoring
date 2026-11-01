"""Polynomial cyclic-carry sums over unramified 2-adic coefficient rings.

The carry construction is adapted from joint_frobenius_carries.py,
blob 1a2392a35185746f8f659496ae8763bfa7de69b7 at base b283e97.
New content: actual Teichmuller moment weights and full ring arithmetic.
The field-only inverse expansion is deliberately not part of this interface.
"""
from itertools import product
from math import comb


class UnramifiedRing:
    """Z/(2^precision)[z]/f for a given monic basic irreducible f.

    Elements pack degree-many coefficient residues in precision-bit blocks.
    Irreducibility modulo two is an input assumption, as in the field model.
    """
    def __init__(self, precision, polynomial):
        if precision < 1 or len(polynomial) < 2 or polynomial[-1] != 1:
            raise ValueError('Positive precision and monic polynomial required')
        self.precision = precision
        self.polynomial = tuple(polynomial)
        self.degree = len(polynomial)-1
        self.modulus = 1 << precision
        self.field_size = 1 << self.degree
        self.generator = self.from_polynomial([0, 1])
        root = self.mod_two(self.mul(self.generator, self.generator))
        derivative = [i*c for i, c in enumerate(polynomial)][1:]
        bits = 1
        while bits < precision:
            error = self.evaluate(polynomial, root)
            slope = self.evaluate(derivative, root)
            root = self.sub(root, self.mul(error, self.unit_inverse(slope)))
            bits *= 2
        assert self.evaluate(polynomial, root) == 0
        self.frobenius_generator = root

    def coefficients(self, value):
        return [(value >> (self.precision*i)) & (self.modulus-1)
                for i in range(self.degree)]

    def from_polynomial(self, values):
        work = list(values)+[0]*max(0, self.degree-len(values))
        for k in range(len(work)-1, self.degree-1, -1):
            c = work[k] % self.modulus
            for j in range(self.degree):
                work[k-self.degree+j] -= c*self.polynomial[j]
        return sum((work[i] % self.modulus) << (self.precision*i)
                   for i in range(self.degree))

    def integer(self, value):
        return value % self.modulus

    def add(self, a, b):
        return self.from_polynomial([x+y for x, y in zip(self.coefficients(a), self.coefficients(b))])

    def neg(self, a):
        return self.from_polynomial([-x for x in self.coefficients(a)])

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        if not a or not b:
            return 0
        left, right = self.coefficients(a), self.coefficients(b)
        values = [0]*(2*self.degree-1)
        for i, x in enumerate(left):
            for j, y in enumerate(right):
                values[i+j] += x*y
        return self.from_polynomial(values)

    def power(self, a, exponent):
        if exponent < 0:
            raise ValueError('Nonnegative exponent required')
        result = 1
        while exponent:
            if exponent & 1:
                result = self.mul(result, a)
            a = self.mul(a, a)
            exponent >>= 1
        return result

    def mod_two(self, a):
        return self.from_polynomial([x & 1 for x in self.coefficients(a)])

    def unit_inverse(self, a):
        residue = self.mod_two(a)
        if not residue:
            raise ValueError('Nonunit has no ring inverse')
        inverse = self.power(residue, self.field_size-2)
        bits = 1
        while bits < self.precision:
            inverse = self.mul(inverse, self.sub(self.integer(2), self.mul(a, inverse)))
            bits *= 2
        assert self.mul(a, inverse) == 1
        return inverse

    def evaluate(self, integer_coefficients, a):
        result = 0
        for c in reversed(integer_coefficients):
            result = self.add(self.mul(result, a), self.integer(c))
        return result

    def frobenius(self, a, power=1):
        for _ in range(power % self.degree):
            result = 0
            for c in reversed(self.coefficients(a)):
                result = self.add(self.mul(result, self.frobenius_generator), self.integer(c))
            a = result
        return a

    def teichmuller(self, residue):
        a = self.from_polynomial([(residue >> i) & 1 for i in range(self.degree)])
        # If a=tau(residue) modulo 2, Rabinoff Lemma 1.4 gives equality
        # modulo 2^precision after this many squarings; q-powers fix tau.
        return self.power(a, 1 << (self.degree*self.precision))


def polynomial_sum(ring, width, row_modulus, variables, factors, nonzero=False):
    """Sum explicit synchronized polynomial factors on Teichmuller inputs.

    A factor is (row, [(coefficient, counts), ...]); counts use variable,
    then relative column, exactly as in the credited field interface.
    No inverse factors, root-indicator compilation, or field enumeration.
    """
    t = row_modulus
    if width < 1 or ring.degree % width or t < 1 or width % t or variables < 1:
        raise ValueError('Incompatible ring, domain, or row dimensions')
    columns = width//t
    dimension = variables*columns
    q = 1 << width
    punctured = [nonzero]*variables if isinstance(nonzero, bool) else list(nonzero)
    if len(punctured) != variables:
        raise ValueError('Bad domain mask')
    zero = (0,)*dimension
    local = []
    for rho in range(t):
        values = {zero: 1}
        for row, terms in factors:
            if not 0 <= row < t:
                raise ValueError('Unsynchronized row')
            if row != rho:
                continue
            nxt = {}
            for h, a in values.items():
                for c, k in terms:
                    if len(k) != dimension or min(k) < 0:
                        raise ValueError('Bad monomial')
                    out = tuple(x+y for x, y in zip(h, k))
                    nxt[out] = ring.add(nxt.get(out, 0), ring.mul(a, c))
            values = {k: c for k, c in nxt.items() if c}
        local.append(values)
    bounds = []
    for v in range(variables):
        for j in range(columns):
            rotated = sum(max([0]+[
                sum(k[v*columns+a] << ((row+t*a-t*j) % width) for a in range(columns))
                for c, k in terms if c]) for row, terms in factors)
            bounds.append(rotated//(q-1))
    answer, cache = 0, {}
    for initial in product(*(range(b+1) for b in bounds)):
        weight = 1
        for v in range(variables):
            moment = q if not punctured[v] and initial[v*columns] == 0 else q-1
            weight = ring.mul(weight, ring.integer(moment))
        if not weight:
            continue
        final = tuple(initial[v*columns+(j+1) % columns]
                      for v in range(variables) for j in range(columns))
        states = {initial: weight}
        for rho in range(t):
            nxt = {}
            for incoming, a in states.items():
                key = (rho, incoming)
                if key not in cache:
                    transitions = []
                    for k, c in local[rho].items():
                        total = tuple(x+y for x, y in zip(incoming, k))
                        if all(x % 2 == 0 for x in total):
                            transitions.append((tuple(x//2 for x in total), c))
                    cache[key] = transitions
                for outgoing, c in cache[key]:
                    nxt[outgoing] = ring.add(nxt.get(outgoing, 0), ring.mul(a, c))
            states = {k: c for k, c in nxt.items() if c}
        answer = ring.add(answer, states.get(final, 0))
    return answer


def parity_lift_coefficients(precision):
    """Integral H with H=0 mod T^h and H=1 mod (1-T)^h, degree 2h-1."""
    if precision < 1:
        raise ValueError('Positive precision required')
    h = precision
    return [0]*h+[(-1)**j*comb(2*h-1, h+j)*comb(h+j-1, j) for j in range(h)]


def residue_root_count(ring, variables, terms, nonzero=False):
    """Count zeros of the reduction modulo 2, modulo 2^ring.precision.

    Terms use ordinary exponents on variables (one relative column).
    Each acquired norm moment has power at most 2*precision-1.
    This counts residue-field roots, not ring congruence roots.
    """
    norm_factors = [(r, [(ring.frobenius(c, r), h) for c, h in terms])
                    for r in range(ring.degree)]
    answer = polynomial_sum(ring, ring.degree, ring.degree, variables, [], nonzero)
    for exponent, coefficient in enumerate(parity_lift_coefficients(ring.precision)):
        if coefficient:
            moment = polynomial_sum(ring, ring.degree, ring.degree, variables,
                                    norm_factors*exponent, nonzero)
            answer = ring.sub(answer, ring.mul(ring.integer(coefficient), moment))
    assert not any(ring.coefficients(answer)[1:])
    return answer


def system_root_count(precision, polynomial, variables, equations):
    """Count a residue-field polynomial system modulo 2^precision.

    Input coefficients are binary-packed residue-field values. All variables
    have the full residue field as domain. The added linear-combination
    variables and the nonunit precision shift are explicit, not free.
    """
    degree = len(polynomial)-1
    q, modulus = 1 << degree, 1 << precision
    count = len(equations)
    if not count:
        return pow(q, variables, modulus)
    shift = degree*(count-1)
    ring = UnramifiedRing(precision+shift, polynomial)
    combined = []
    for i, equation in enumerate(equations):
        for coefficient, exponents in equation:
            h = list(exponents)+[0]*count
            h[variables+i] = 1
            combined.append((ring.teichmuller(coefficient), tuple(h)))
    total = residue_root_count(ring, variables+count, combined)
    baseline = ring.integer(1 << (degree*(variables+count-1)))
    numerator = ring.sub(total, baseline)
    assert numerator % (1 << shift) == 0
    return ((numerator >> shift)*pow(q-1, -1, modulus)) % modulus
