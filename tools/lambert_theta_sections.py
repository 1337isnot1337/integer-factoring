"""Exact binary-index signed divisor coefficients from a theta quotient.

This is an exploratory coefficient engine, without divisor enumeration. Its
explicit worst-case cost is exponential in the input bit length; it supplies
an exact section representation, not a stronger factoring bound.
"""

from math import isqrt

from cyclotomic_character_decoder import modular_product


class AdditivePhaseAlgebra:
    """F_2[Z]/(1+Z+...+Z**(q-1)), for an odd prime q."""

    def __init__(self, prime):
        if (not isinstance(prime, int) or prime < 3 or prime % 2 == 0
                or any(prime % d == 0 for d in range(3, isqrt(prime)+1, 2))):
            raise ValueError("an odd prime additive phase order is required")
        self.prime = prime
        self.modulus = (1 << prime)-1
        self.phases = tuple((1 << k) if k < prime-1 else (1 << (prime-1))-1
                            for k in range(prime))

    def phase(self, exponent):
        return self.phases[exponent % self.prime]

    def product(self, left, right):
        return modular_product(left, right, self.modulus)

    def inverse_theta_constant(self, phase):
        # (1+Z) * (Z+Z^3+...+Z^(q-2)) = 1 modulo Phi_q.
        inverse = 0
        for k in range(1, self.prime-1, 2):
            inverse ^= self.phase(phase*k)
        assert self.product(1 ^ self.phase(phase), inverse) == 1
        return inverse

    def signed_residues(self, value):
        # Two degree-<q lifts differ by Phi_q. Exactly one has augmentation
        # zero because q is odd. Genuine signed outputs also have c_0=0.
        lifted = value ^ (self.modulus if value.bit_count() & 1 else 0)
        if lifted & 1 or any((lifted >> a & 1) != (lifted >> (self.prime-a) & 1)
                             for a in range(1, self.prime)):
            raise ValueError("the phase value is not a signed residue vector")
        return tuple(lifted >> a & 1 for a in range(self.prime))


def theta_series(algebra, limit, phase=1, derivative=False):
    """Theta or its formal-argument Euler derivative, at Z**phase.

    Differentiate in an independent argument W, then set W=Z**phase.
    The derivative is taken BEFORE imposing Z**q=1; differentiating a
    reduced phase representative would give a different operation.
    """
    values = {}
    bound = (isqrt(8*limit+1)-1)//2
    for s in range(bound+1):
        exponent = s*(s+1)//2
        if derivative:
            value = algebra.phase(phase*(s+1 if s % 2 == 0 else -s))
        else:
            value = algebra.phase(phase*(s+1)) ^ algebra.phase(-phase*s)
        if value:
            values[exponent] = value
    return values


def coefficient(integer, prime, retain_numerators=False):
    """Return the exact signed divisor residue value and charged section sizes."""
    if not isinstance(integer, int) or integer < 1:
        raise ValueError("a positive coefficient index is required")
    algebra = AdditivePhaseAlgebra(prime)
    numerator = theta_series(algebra, integer, derivative=True)
    target, phase = integer, 1
    records = []
    while target:
        denominator = theta_series(algebra, target, phase)
        parity, child_limit = target & 1, target//2
        child = {}
        pairs = products = 0
        for i, a in numerator.items():
            for j, b in denominator.items():
                if i+j > target:
                    break
                pairs += 1
                if (i+j) & 1 != parity:
                    continue
                products += 1
                position = (i+j)//2
                value = child.get(position, 0) ^ algebra.product(a, b)
                if value:
                    child[position] = value
                else:
                    child.pop(position, None)
        record = {"target": target, "bit": parity, "phase": phase,
                  "numerator_terms": len(numerator), "denominator_terms": len(denominator),
                  "candidate_pairs": pairs, "ring_products": products,
                  "child_terms": len(child), "child_limit": child_limit}
        if retain_numerators:
            record["numerator"] = [[k, v] for k, v in sorted(numerator.items())]
        records.append(record)
        numerator, target, phase = child, child_limit, 2*phase % prime
    result = algebra.product(numerator.get(0, 0), algebra.inverse_theta_constant(phase))
    return {"value": result, "signed_residues": algebra.signed_residues(result),
            "sections": records,
            "peak_numerator_terms": max(record["numerator_terms"] for record in records),
            "ring_products": sum(record["ring_products"] for record in records)}


def q3_coefficient(integer):
    """Closed high-index control: the 3-free part is a square iff the value is 1."""
    if not isinstance(integer, int) or integer < 1:
        raise ValueError("a positive coefficient index is required")
    while integer % 3 == 0:
        integer //= 3
    return int(isqrt(integer)**2 == integer)


def square_section_rank_witness(split_bits):
    """An identity submatrix in [a_(2**k*u+v)] at additive phase order 3.

    For k>=6, its size is at least floor(2**k/36). This is a lower
    bound for a common linear section space, not for arbitrary algorithms.
    """
    if not isinstance(split_bits, int) or split_bits < 3:
        raise ValueError("at least three split bits are required")
    modulus = 1 << split_bits
    classes = {a: [] for a in range(3)}
    for root in range(modulus//2+1, 3*modulus//4, 2):
        if root % 3:
            classes[(root*root//modulus) % 3].append(root)
    selected = max(classes, key=lambda a: len(classes[a]))
    roots = classes[selected]
    return {"split_bits": split_bits, "modulus": modulus,
            "roots": roots, "rows": [root*root//modulus for root in roots],
            "columns": [root*root % modulus for root in roots],
            "row_class_mod_3": selected, "dimension": len(roots),
            "uniform_lower_bound": modulus//36}
