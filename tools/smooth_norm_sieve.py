"""Deterministic, complete-at-a-supplied-cutoff quadratic-form norm sieve.

This is a concrete search baseline for the short-norm contract, not a
subexponential factoring theorem. Its running time includes the coordinate
region actually traversed. All arithmetic and square-root bounds are exact.
"""

from dataclasses import dataclass
from math import gcd, isqrt

from prime_form_factor_base import modular_square_root


@dataclass(frozen=True)
class NormSieveResult:
    vector: tuple | None
    norm: int | None
    multiplicity: int | None
    maximum_norm: int
    rows_visited: int
    coordinate_points_visited: int
    primitive_coprime_candidates: int
    modular_strikes: int
    prime_divisions: int
    exhausted: bool


def primes_through(bound):
    if bound < 2:
        return ()
    sieve = bytearray(b"\x01") * (bound + 1)
    sieve[0:2] = b"\x00\x00"
    for prime in range(2, isqrt(bound) + 1):
        if sieve[prime]:
            start = prime * prime
            sieve[start:bound + 1:prime] = b"\x00" * (((bound - start) // prime) + 1)
    return tuple(value for value in range(2, bound + 1) if sieve[value])


class SmoothNormSieve:
    def __init__(self, discriminant, smooth_bound):
        if discriminant >= 0 or discriminant % 4 not in (0, 1):
            raise ValueError("expected a negative discriminant")
        if smooth_bound < 2:
            raise ValueError("smooth bound must be at least two")
        self.discriminant = discriminant
        self.smooth_bound = smooth_bound
        prime_roots = []
        for prime in primes_through(smooth_bound):
            if discriminant % prime == 0:
                continue
            if prime == 2:
                if discriminant % 8 == 1:
                    prime_roots.append((2, 1))
                continue
            # Prime enumeration supplies the root helper's known-prime
            # contract. Canonicalize the orientation for reproducibility.
            root = modular_square_root(discriminant, prime)
            if root is not None:
                root = min(root, prime - root)
                prime_roots.append((prime, root))
        self.prime_roots = tuple(prime_roots)
        self._scheduled_form = None
        self._root_schedules = None

    @staticmethod
    def _row_roots(form, y, prime, square_root):
        """None means every x is a root; otherwise give all root residues."""
        a, b, c = form
        if prime == 2:
            return tuple(x for x in (0, 1) if (a * x * x + b * x * y + c * y * y) % 2 == 0)
        if a % prime:
            inverse = pow(2 * a, -1, prime)
            return tuple(sorted({((-b + sign * square_root) * y * inverse) % prime
                                 for sign in (-1, 1)}))
        if y % prime == 0:
            return None
        # p does not divide Delta, so B is nonzero when p divides A.
        return ((-c * y * pow(b, -1, prime)) % prime,)

    def search(self, form, maximum_norm, multiplicity_budget):
        """Find a primitive coprime smooth norm, or exhaust the exact ellipse.

        Vectors (x,y) and (-x,-y) are identified. The scan takes y>=0;
        at y=0 only x=1 can be primitive. A successful return need not be
        the smallest represented norm. An exhausted return covers every
        primitive vector of norm at most maximum_norm.
        """
        a, b, c = form
        if (a <= 0 or b * b - 4 * a * c != self.discriminant
                or gcd(gcd(a, b), c) != 1):
            raise ValueError("form does not match this primitive discriminant")
        if maximum_norm < 1 or multiplicity_budget < 0:
            raise ValueError("invalid norm or multiplicity budget")
        if self._scheduled_form != tuple(form):
            root_schedules = []
            for prime, square_root in self.prime_roots:
                if prime == 2:
                    root_schedules.append((prime, (
                        self._row_roots(form, 0, prime, square_root),
                        self._row_roots(form, 1, prime, square_root),
                    ), False))
                else:
                    # Every nonexceptional row root is y times a fixed
                    # slope. Retain one schedule across this target's caps.
                    slopes = self._row_roots(form, 1, prime, square_root)
                    assert slopes is not None
                    root_schedules.append((prime, slopes, a % prime == 0))
            self._scheduled_form = tuple(form)
            self._root_schedules = tuple(root_schedules)
        root_schedules = self._root_schedules
        rows = points = candidates = strikes = divisions = 0
        y_limit = isqrt((4 * a * maximum_norm) // (-self.discriminant))
        for y in range(y_limit + 1):
            if y == 0:
                if a > maximum_norm:
                    continue
                lower = upper = 1
            else:
                radicand = 4 * a * maximum_norm + self.discriminant * y * y
                radius = isqrt(radicand)
                lower = -((radius + b * y) // (2 * a))
                upper = (radius - b * y) // (2 * a)
                if lower > upper:
                    continue
            rows += 1
            width = upper - lower + 1
            points += width
            remainders = [0] * width
            multiplicities = [0] * width
            norm = a * lower * lower + b * lower * y + c * y * y
            step = a * (2 * lower + 1) + b * y
            for index in range(width):
                x = lower + index
                assert 0 < norm <= maximum_norm
                if gcd(x, y) == 1 and gcd(norm, self.discriminant) == 1:
                    remainders[index] = norm
                    candidates += 1
                norm += step
                step += 2 * a
            for prime, slopes, leading_divisible in root_schedules:
                if prime == 2:
                    roots = slopes[y & 1]
                elif leading_divisible and y % prime == 0:
                    roots = None
                else:
                    roots = tuple(sorted({slope * y % prime for slope in slopes}))
                if roots is None:
                    progressions = (range(width),)
                else:
                    progressions = tuple(range((root - lower) % prime, width, prime) for root in roots)
                for progression in progressions:
                    for index in progression:
                        remainder = remainders[index]
                        if not remainder or remainder == 1:
                            continue
                        strikes += 1
                        while remainder % prime == 0:
                            remainder //= prime
                            multiplicities[index] += 1
                            divisions += 1
                            if multiplicities[index] > multiplicity_budget:
                                remainder = 0
                                break
                        remainders[index] = remainder
            for index, remainder in enumerate(remainders):
                if remainder == 1:
                    x = lower + index
                    norm = a * x * x + b * x * y + c * y * y
                    return NormSieveResult(
                        (x, y), norm, multiplicities[index], maximum_norm,
                        rows, points, candidates, strikes, divisions, False,
                    )
        return NormSieveResult(
            None, None, None, maximum_norm,
            rows, points, candidates, strikes, divisions, True,
        )
