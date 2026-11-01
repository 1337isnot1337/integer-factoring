"""Deterministic finite-control factoring with recursive Lucas certificates.

The fixed Pollard schedule is only an accelerator. Trial division is the
complete fallback. No subexponential worst-case bound is claimed for this
helper; a recursive asymptotic factoring routine may replace it.
"""

from math import gcd, isqrt, prod


class CertifiedFactorizer:
    def __init__(self):
        self.certificates = {2: {'factors': [], 'witnesses': []}}
        self.cache = {1: ()}
        self.stats = {'rho_steps': 0, 'primality_certificates': 1,
                      'trial_division_fallbacks': 0}

    @staticmethod
    def _probable(n):
        if n < 2:
            return False
        for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
            if n % p == 0:
                return n == p
        d, s = n-1, 0
        while d % 2 == 0:
            d //= 2
            s += 1
        for a in (2, 3, 5, 7, 11, 13, 17):
            if a >= n:
                continue
            value = pow(a, d, n)
            if value in (1, n-1):
                continue
            for _ in range(s-1):
                value = value*value % n
                if value == n-1:
                    break
            else:
                return False
        return True

    def _certificate(self, n, exhaustive=False):
        if n in self.certificates:
            return True
        predecessor = self.factor(n-1)
        witnesses = []
        for q, _ in predecessor:
            for a in range(2, n if exhaustive else min(n, 258)):
                if (pow(a, n-1, n) == 1
                        and gcd(pow(a, (n-1)//q, n)-1, n) == 1):
                    witnesses.append((q, a))
                    break
            else:
                return False
        self.certificates[n] = {'factors': predecessor, 'witnesses': witnesses}
        self.stats['primality_certificates'] += 1
        return True

    def _divisor(self, n):
        for c in range(1, 33):
            y, r, divisor, steps = 2+c, 1, 1, 0
            saved, x = y, y
            while divisor == 1 and steps < 1 << 20:
                x = y
                for _ in range(r):
                    y = (y*y+c) % n
                    steps += 1
                offset = 0
                while offset < r and divisor == 1:
                    saved, batch = y, 1
                    for _ in range(min(64, r-offset)):
                        y = (y*y+c) % n
                        batch = batch*abs(x-y) % n
                        steps += 1
                    divisor = gcd(batch, n)
                    offset += 64
                r *= 2
            if divisor == n:
                divisor = 1
                while divisor == 1 and steps < 1 << 21:
                    saved = (saved*saved+c) % n
                    divisor = gcd(abs(x-saved), n)
                    steps += 1
            self.stats['rho_steps'] += steps
            if 1 < divisor < n:
                return divisor
        self.stats['trial_division_fallbacks'] += 1
        for divisor in range(3, isqrt(n)+1, 2):
            if n % divisor == 0:
                return divisor
        return n

    def factor(self, integer):
        if not isinstance(integer, int) or integer < 1:
            raise ValueError('A positive integer is required.')
        if integer in self.cache:
            return self.cache[integer]
        remaining, factors = integer, {}

        def add(p, exponent):
            factors[p] = factors.get(p, 0)+exponent

        for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
            exponent = 0
            while remaining % p == 0:
                remaining //= p
                exponent += 1
            if exponent:
                assert self._certificate(p)
                add(p, exponent)
        if remaining > 1:
            if self._probable(remaining) and self._certificate(remaining):
                add(remaining, 1)
            else:
                divisor = self._divisor(remaining)
                if divisor == remaining:
                    # Exhaustive division has established primality. The full
                    # witness search now terminates and supplies a certificate.
                    assert self._certificate(remaining, exhaustive=True)
                    add(remaining, 1)
                else:
                    for p, exponent in self.factor(divisor)+self.factor(remaining//divisor):
                        add(p, exponent)
        assert prod(p**e for p, e in factors.items()) == integer
        answer = tuple(sorted(factors.items()))
        self.cache[integer] = answer
        return answer


def verify_prime_certificates(certificates):
    """Verify every supplied Lucas certificate, without a probable-prime test."""
    records = {int(p): record for p, record in certificates.items()}
    verified = set()
    for n in sorted(records):
        record = records[n]
        factors = [(int(q), int(e)) for q, e in record['factors']]
        witnesses = [(int(q), int(a)) for q, a in record['witnesses']]
        if n == 2:
            if factors or witnesses:
                raise ValueError('The base certificate for two must be empty.')
            verified.add(2)
            continue
        if (n < 3 or any(q not in verified or e < 1 for q, e in factors)
                or len({q for q, _ in factors}) != len(factors)
                or prod(q**e for q, e in factors) != n-1
                or {q for q, _ in witnesses} != {q for q, _ in factors}
                or len(witnesses) != len(factors)):
            raise ValueError('Invalid factored predecessor in a prime certificate.')
        for q, a in witnesses:
            if (not 1 < a < n or pow(a, n-1, n) != 1
                    or gcd(pow(a, (n-1)//q, n)-1, n) != 1):
                raise ValueError('Invalid Lucas witness.')
        verified.add(n)
    return verified
