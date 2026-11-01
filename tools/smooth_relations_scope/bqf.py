"""Exact positive-definite binary quadratic form arithmetic (self-contained).

Forms are integer triples (a, b, c) with b^2 - 4ac = D < 0, a > 0.
All arithmetic is exact Python integers; no floating point is used.
"""
from math import gcd, isqrt


def xgcd(a, b):
    """Return (g, x, y) with x*a + y*b = g = gcd(a, b) >= 0."""
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        q, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    if a < 0:
        a, x0, y0 = -a, -x0, -y0
    return a, x0, y0


def disc(f):
    a, b, c = f
    return b * b - 4 * a * c


def reduce_form(f):
    a, b, c = f
    assert a > 0 and c > 0
    while True:
        # normalize -a < b <= a
        if not (-a < b <= a):
            k = (a - b) // (2 * a)  # b + 2ak in (-a, a]
            c = c + b * k + a * k * k
            b = b + 2 * a * k
        if a > c:
            a, b, c = c, -b, a
            continue
        if a == c and b < 0:
            b = -b
        return (a, b, c)


def is_reduced(f):
    a, b, c = f
    if not (-a < b <= a and a <= c):
        return False
    if a == c and b < 0:
        return False
    return True


def compose(f1, f2):
    """Gaussian composition (Cohen, A Course in Computational Algebraic
    Number Theory, Alg. 5.4.7), followed by reduction."""
    D = disc(f1)
    assert disc(f2) == D
    a1, b1, c1 = f1
    a2, b2, c2 = f2
    if a1 > a2:
        a1, b1, c1, a2, b2, c2 = a2, b2, c2, a1, b1, c1
    s = (b1 + b2) // 2
    n = b2 - s
    if a2 % a1 == 0:
        y1 = 0
        d = a1
    else:
        d, u, v = xgcd(a2, a1)
        y1 = u
    if s % d == 0:
        y2 = -1
        x2 = 0
        d1 = d
    else:
        d1, x2, y2 = xgcd(s, d)
        y2 = -y2
    v1 = a1 // d1
    v2 = a2 // d1
    r = (y1 * y2 * n - x2 * c2) % v1
    b3 = b2 + 2 * v2 * r
    a3 = v1 * v2
    num = b3 * b3 - D
    assert num % (4 * a3) == 0
    c3 = num // (4 * a3)
    return reduce_form((a3, b3, c3))


def identity(D):
    if D % 4 == 0:
        return (1, 0, -D // 4)
    return (1, 1, (1 - D) // 4)


def inverse(f):
    a, b, c = f
    return reduce_form((a, -b, c))


def power(f, k):
    D = disc(f)
    r = identity(D)
    base = f
    while k:
        if k & 1:
            r = compose(r, base)
        base = compose(base, base)
        k >>= 1
    return r


def reduced_forms(D):
    """All primitive reduced forms of discriminant D < 0."""
    out = []
    amax = isqrt(-D // 3)
    for a in range(1, amax + 1):
        for b in range(-a + 1, a + 1):
            if (b - D) % 2:
                continue
            num = b * b - D
            if num % (4 * a):
                continue
            c = num // (4 * a)
            if c < a:
                continue
            if (a == c or b == a) and b < 0:
                continue
            if gcd(gcd(a, b), c) != 1:
                continue
            out.append((a, b, c))
    return out


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    t = pow(a, (p - 1) // 2, p)
    return 1 if t == 1 else -1


def sqrt_mod_prime(n, p):
    """Deterministic exhaustive square root for small odd prime p (tests only)."""
    n %= p
    for x in range(p):
        if x * x % p == n:
            return x
    return None


def prime_form(D, p):
    """Reduced prime form of norm p, or None if p not split (p odd prime, or p=2)."""
    if p == 2:
        if D % 8 == 1:
            return reduce_form((2, 1, (1 - D) // 8))
        return None
    if D % p == 0 or legendre(D, p) != 1:
        return None
    b = sqrt_mod_prime(D, p)
    if (b - D) % 2:
        b = p - b
    num = b * b - D
    assert num % (4 * p) == 0
    return reduce_form((p, b, num // (4 * p)))


def is_prime(n):
    if n < 2:
        return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:  # deterministic for n < 3.3e24
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, isqrt(n) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(sieve[i * i :: i]))
    return [i for i in range(n + 1) if sieve[i]]


def genus_char(f, q):
    """chi_q(f) = (m/q) for any m represented by f with q not dividing m.
    q must be an odd prime dividing disc(f), f primitive."""
    a, b, c = f
    for m in (a, c, a + b + c, a - b + c):
        if m % q:
            return legendre(m, q)
    raise ValueError("no represented value coprime to q found")


def ambiguous_leading(f):
    """For a reduced ambiguous form return the positive represented divisor m
    of |D| used for splitting: a if b in {0, a}; 2a-b = f(1,-1) if a == c."""
    a, b, c = f
    if b == 0 or b == a:
        return a
    if a == c:
        return 2 * a - b
    raise ValueError("not ambiguous")


def is_ambiguous_reduced(f):
    a, b, c = f
    return b == 0 or b == a or a == c
