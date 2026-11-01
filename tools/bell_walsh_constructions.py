"""Exact Bell paths from arbitrary quadratic-parity Walsh energies.

The energies need not be weights of a binary code.  The constructor uses
only a positive energy on each nonzero frequency and a prescribed parity.
Search is optional and uses NumPy; the algebraic constructor is standard
library only.  No factoring improvement is asserted by this module.
"""

from collections import Counter
from fractions import Fraction
from math import comb, exp, log, log2

from bell_crossing import transform


def quadratic(value, bits):
    mask = (1 << bits) - 1
    return ((value & mask) & (value >> bits)).bit_count() & 1


def parity(value):
    return value.bit_count() & 1


def validate_energies(energies, bits, shift=0):
    q = 1 << (2 * bits)
    if bits < 1 or len(energies) != q or not 0 <= shift < q:
        raise ValueError("wrong alphabet or shift")
    if energies[0] != 0 or any(e < 1 for e in energies[1:]):
        raise ValueError("only the zero frequency may have energy zero")
    for value, energy in enumerate(energies):
        if not isinstance(energy, int) or energy & 1 != (
                quadratic(value, bits) ^ parity(value & shift)):
            raise ValueError("energies must have the specified quadratic parity")


def walsh(values, prime):
    result = list(values)
    step = 1
    while step < len(result):
        for block in range(0, len(result), 2 * step):
            for j in range(block, block + step):
                a, b = result[j], result[j + step]
                result[j] = (a + b) % prime
                result[j + step] = (a - b) % prime
        step *= 2
    return tuple(result)


def energy_polynomials(degree, prime):
    """Rows are (1+z)^(degree-e) (1-z)^e, for 0 <= e <= degree."""
    rows = [tuple(comb(degree, j) % prime for j in range(degree + 1))]
    for _ in range(degree):
        old, new = rows[-1], []
        for j in range(degree + 1):
            new.append((old[j] - (old[j - 1] if j else 0)
                        - (new[j - 1] if j else 0)) % prime)
        rows.append(tuple(new))
    return tuple(rows)


def coefficients(energies, bits, prime=3, shift=0):
    validate_energies(energies, bits, shift)
    if prime == 2:
        raise ValueError("requires odd characteristic")
    degree = max(energies)
    polynomials = energy_polynomials(degree, prime)
    scale = pow(pow(2, degree, prime), -1, prime)
    return tuple(tuple(scale * x % prime for x in walsh(
        tuple(polynomials[e][j] for e in energies), prime))
        for j in range(degree + 1))


def shifted_dual(values, bits, prime, shift):
    other = transform(tuple(values), bits, prime)
    return tuple(other[j ^ shift] for j in range(len(values)))


def valuation_histogram(rows):
    degree = len(rows) - 1
    counts = Counter(next(i for i in range(degree + 1) if rows[i][j])
                     for j in range(len(rows[0])))
    return tuple(counts[i] for i in range(degree + 1))


def pressure(histogram, scale):
    scale = Fraction(scale)
    if scale < 1:
        raise ValueError("scale must be at least one")
    degree = len(histogram) - 1
    total = sum(Fraction(n) * scale ** (-i) for i, n in enumerate(histogram))
    return scale ** degree * total * total


def optimum(histogram):
    """Exploratory floating optimum; use pressure() for a certificate."""
    degree = len(histogram) - 1
    active = [(i, n) for i, n in enumerate(histogram) if n]

    def evaluate(t):
        terms = [(log(n) - i * t, i) for i, n in active]
        maximum = max(a for a, _ in terms)
        weights = [(exp(a - maximum), i) for a, i in terms]
        total = sum(a for a, _ in weights)
        derivative = degree / 2 - sum(a * i for a, i in weights) / total
        return degree * t / 2 + maximum + log(total), derivative

    if evaluate(0)[1] >= 0:
        return evaluate(0)[0], 1.0
    left, right = 0.0, 1.0
    while evaluate(right)[1] < 0:
        right *= 2
    for _ in range(42):
        middle = (left + right) / 2
        if evaluate(middle)[1] < 0:
            left = middle
        else:
            right = middle
    t = (left + right) / 2
    return evaluate(t)[0], exp(t)


def receipt(energies, bits, prime=3, shift=0):
    energies = tuple(int(x) for x in energies)
    rows = coefficients(energies, bits, prime, shift)
    degree = len(rows) - 1
    assert all(sum(row[j] for row in rows) % prime == 1
               for j in range(len(energies)))
    assert tuple(shifted_dual(row, bits, prime, shift) for row in rows) == rows[::-1]
    histogram = valuation_histogram(rows)
    _, approximate_scale = optimum(histogram)
    scale = Fraction(approximate_scale).limit_denominator(65536)
    bound = pressure(histogram, scale)
    return {
        "physical_bits": bits,
        "prime": prime,
        "quadratic_shift": shift,
        "energies": list(energies),
        "degree": degree,
        "valuation_histogram": list(histogram),
        "scale": [scale.numerator, scale.denominator],
        "pressure": [bound.numerator, bound.denominator],
        "rho_pressure": (log2(bound.numerator) - log2(bound.denominator)) / (4 * bits),
        "strict_sufficient_target": bound ** 577 < 2 ** (4 * bits * 408),
        "coefficients_by_degree": [list(row) for row in rows],
    }


def folded_weights(bits):
    """Put the even-(2b+1)-cube polar form in symplectic coordinates."""
    dimension = 2 * bits

    def polar(a, b):
        return parity(a & b) ^ (parity(a) & parity(b))

    remaining = [1 << j for j in range(dimension)]
    pairs = []
    while remaining:
        a = remaining.pop(0)
        j = next(j for j, b in enumerate(remaining) if polar(a, b))
        b = remaining.pop(j)
        pairs.append((a, b))
        remaining = [v ^ (a if polar(v, b) else 0)
                     ^ (b if polar(v, a) else 0) for v in remaining]
    basis = [a for a, _ in pairs] + [b for _, b in pairs]
    weights = []
    for value in range(1 << dimension):
        image = 0
        for j, vector in enumerate(basis):
            if value >> j & 1:
                image ^= vector
        weights.append(image.bit_count())
    return tuple(weights)


def folded_energies(bits):
    """The canonical folded-cube energies, with their linear parity shift."""
    dimension = 2 * bits
    energies = [(weight + 1) // 2 for weight in folded_weights(bits)]
    shift = sum((energies[1 << j] & 1) << j for j in range(dimension))
    validate_energies(energies, bits, shift)
    return tuple(energies), shift


def radial_energies(levels, bits):
    if len(levels) != 2 * bits + 1 or levels[0] != 0:
        raise ValueError("one energy is required for each Hamming weight")
    if any(e < 1 or e % 2 != ((weight + 1) // 2) % 2
           for weight, e in enumerate(levels[1:], 1)):
        raise ValueError("radial energies violate positivity or quadratic parity")
    _, shift = folded_energies(bits)
    return tuple(levels[weight] for weight in folded_weights(bits)), shift


def krawtchouk_matrix(dimension, prime):
    return tuple(tuple(sum((-1) ** j * comb(r, j) * comb(dimension - r, w - j)
                           for j in range(max(0, w - dimension + r), min(r, w) + 1)) % prime
                       for w in range(dimension + 1))
                 for r in range(dimension + 1))


def radial_histogram(levels, bits, prime=3):
    degree, dimension = max(levels), 2 * bits
    polys = energy_polynomials(degree, prime)
    matrix = krawtchouk_matrix(dimension, prime)
    counts = Counter()
    for r, row in enumerate(matrix):
        order = next(j for j in range(degree + 1)
                     if sum(a * polys[e][j] for a, e in zip(row, levels)) % prime)
        counts[order] += comb(dimension, r)
    return tuple(counts[j] for j in range(degree + 1))


def quadratic_two_level(bits):
    """An explicit degree-two path beyond positive-code-weight quotients."""
    return tuple(0 if value == 0 else 1 if quadratic(value, bits) else 2
                 for value in range(1 << (2 * bits)))


def symplectic_relabel_energies(source, bits):
    """Convert any nondegenerate quadratic parity to standard Bell indices."""
    dimension = 2 * bits
    if len(source) != 1 << dimension or source[0] != 0:
        raise ValueError("wrong source energy table")

    def polar(a, b):
        return (source[a] ^ source[b] ^ source[a ^ b]) & 1

    remaining = [1 << j for j in range(dimension)]
    pairs = []
    while remaining:
        a = remaining.pop(0)
        j = next((j for j, b in enumerate(remaining) if polar(a, b)), None)
        if j is None:
            raise ValueError("degenerate quadratic parity")
        b = remaining.pop(j)
        pairs.append((a, b))
        remaining = [v ^ (a if polar(v, b) else 0)
                     ^ (b if polar(v, a) else 0) for v in remaining]
    basis = [a for a, _ in pairs] + [b for _, b in pairs]
    energies = []
    for value in range(1 << dimension):
        image = 0
        for j, vector in enumerate(basis):
            if value >> j & 1:
                image ^= vector
        energies.append(source[image])
    shift = sum((energies[1 << j] & 1) << j for j in range(dimension))
    validate_energies(energies, bits, shift)
    return tuple(energies), shift


def norm_energies(bits, polynomial, output_basis=None):
    """Energies from an anisotropic quadratic norm over F_(2^bits).

    The supplied binary modulus must be irreducible; the checker verifies
    each modulus it uses.  output_basis is an optional invertible binary
    coordinate change on the norm value.
    """
    size = 1 << bits

    def multiply(a, b):
        value = 0
        while b:
            if b & 1:
                value ^= a
            b >>= 1
            a <<= 1
            if a & size:
                a ^= polynomial
        return value

    def trace(a):
        total = 0
        for _ in range(bits):
            total ^= a
            a = multiply(a, a)
        return total

    c = next(a for a in range(1, size) if trace(a) == 1)
    squares = [multiply(a, a) for a in range(size)]
    source = []
    for y in range(size):
        for x in range(size):
            value = squares[x] ^ multiply(x, y) ^ multiply(c, squares[y])
            if output_basis is not None:
                image = 0
                for j, vector in enumerate(output_basis):
                    if value >> j & 1:
                        image ^= vector
                value = image
            source.append(value.bit_count())
    return symplectic_relabel_energies(source, bits)


def reciprocal_retime(rows, degrees):
    """Move coefficient indices while preserving a reciprocal endpoint."""
    old_degree = len(rows) - 1
    if len(degrees) != len(rows) or degrees[0] != 0:
        raise ValueError("wrong degree schedule")
    if any(a >= b for a, b in zip(degrees, degrees[1:])):
        raise ValueError("degree schedule must be strictly increasing")
    degree = degrees[-1]
    if any(degrees[i] + degrees[old_degree - i] != degree
           for i in range(old_degree + 1)):
        raise ValueError("degree schedule must preserve reciprocal pairs")
    result = [(0,) * len(rows[0]) for _ in range(degree + 1)]
    for index, row in zip(degrees, rows):
        result[index] = tuple(row)
    return tuple(result)


def retiming_blocks(low, high):
    """Exact pooled blocks for optimal monotone coefficient retiming.

    The infimum pressure is (sum sqrt(A_B B_B))**2. The returned integer
    blocks are also a first-order optimality certificate; no numerical
    optimization is used. Zero bins are allowed, but the path must have
    a nonzero constant coefficient and a nonzero highest opposite degree.
    """
    low, high = tuple(low), tuple(high)
    if len(low) != len(high) or not low or low[0] <= 0 or high[-1] <= 0:
        raise ValueError("normalize the common monomial before retiming")
    if any(not isinstance(n, int) or n < 0 for n in low + high):
        raise ValueError("rank multiplicities must be nonnegative integers")
    blocks = []
    for index, (a, b) in enumerate(zip(low, high)):
        if a + b == 0:
            continue
        blocks.append([index, index, a, b])
        while len(blocks) > 1 and blocks[-2][2] * blocks[-1][3] > blocks[-1][2] * blocks[-2][3]:
            right = blocks.pop()
            left = blocks[-1]
            left[1] = right[1]
            left[2] += right[2]
            left[3] += right[3]
    for i, block in enumerate(blocks):
        block[1] = blocks[i + 1][0] - 1 if i + 1 < len(blocks) else len(low) - 1
        assert block[2] > 0 and block[3] > 0
        # The gradient prefix has the sign required by monotonicity.
        a = b = 0
        for index in range(block[0], block[1] + 1):
            a += low[index]
            b += high[index]
            assert a * block[3] >= block[2] * b
    return tuple(tuple(block) for block in blocks)


def _nearest_half_log(ratio, scale):
    """Nearest integer to log(ratio)/(2 log(scale)), using exact powers."""
    ratio, scale = Fraction(ratio), Fraction(scale)
    if ratio < 1 or scale <= 1:
        raise ValueError("requires ratio >= 1 and scale > 1")
    low, high = -1, 1
    while scale ** (2 * high + 1) < ratio:
        high *= 2
    while high - low > 1:
        middle = (low + high) // 2
        if scale ** (2 * middle + 1) < ratio:
            low = middle
        else:
            high = middle
    return high


def monotone_retiming(low, high, accuracy=64):
    """An exact rational schedule within x**(degree+1) of the infimum.

    Here x=1+1/accuracy. All coefficient indices are strictly increasing;
    reciprocity is not required. The integer coefficient vectors themselves
    do not change. Higher accuracy can produce large exact rational receipts.
    """
    if not isinstance(accuracy, int) or accuracy < 1:
        raise ValueError("accuracy must be a positive integer")
    blocks = retiming_blocks(low, high)
    scale = Fraction(accuracy + 1, accuracy)
    first_ratio = Fraction(blocks[0][2], blocks[0][3])
    degrees = [None] * len(low)
    for start, stop, a, b in blocks:
        rounded = _nearest_half_log(Fraction(a, b) / first_ratio, scale)
        for index in range(start, stop + 1):
            degrees[index] = rounded + index
    assert degrees[0] == 0
    assert all(a < b for a, b in zip(degrees, degrees[1:]))
    a = sum(Fraction(n) * scale ** (-i) for n, i in zip(low, degrees))
    b = sum(Fraction(n) * scale ** i for n, i in zip(high, degrees))
    bound = a * b
    factor = scale ** len(low)
    return {
        "blocks": [list(block) for block in blocks],
        "retimed_degrees": degrees,
        "scale": [scale.numerator, scale.denominator],
        "pressure": [bound.numerator, bound.denominator],
        "approximation_factor": [factor.numerator, factor.denominator],
    }


def folded_retiming(bits):
    """Exact rational approximation to all independently optimal exponents.

    The schedule is chosen entirely by rational power comparisons. It has
    degree O(bits^3), and uses scale 1+1/bits^2.
    """
    counts = tuple(comb(2 * bits + 1, i) for i in range(bits + 1))
    scale = Fraction(bits * bits + 1, bits * bits)
    positive = []
    for i in range((bits + 1) // 2):
        ratio = Fraction(counts[bits - i], counts[i])
        positive.append(_nearest_half_log(ratio, scale))
    signed = positive + ([0] if bits % 2 == 0 else []) + [-x for x in reversed(positive)]
    assert all(a > b for a, b in zip(signed, signed[1:]))
    degrees = tuple(signed[0] - x for x in signed)
    root_pressure = sum(Fraction(n) * scale ** exponent
                        for n, exponent in zip(counts, signed))
    return {
        "physical_bits": bits,
        "original_histogram": list(counts),
        "retimed_degrees": list(degrees),
        "signed_exponents": list(signed),
        "scale": [scale.numerator, scale.denominator],
        "pressure": [root_pressure.numerator ** 2, root_pressure.denominator ** 2],
        "rho_pressure": (log2(root_pressure.numerator) - log2(root_pressure.denominator)) / (2 * bits),
    }


def _numpy_histogram(energies, prime, polynomial_cache):
    import numpy as np
    degree = int(max(energies))
    key = (degree, prime)
    if key not in polynomial_cache:
        polynomial_cache[key] = np.array(energy_polynomials(degree, prime), dtype=np.int32)
    work = polynomial_cache[key][energies].copy()
    q = len(energies)
    step = 1
    while step < q:
        blocks = work.reshape(-1, 2 * step, degree + 1)
        left = blocks[:, :step, :].copy()
        right = blocks[:, step:, :].copy()
        blocks[:, :step, :] = (left + right) % prime
        blocks[:, step:, :] = (left - right) % prime
        step *= 2
    orders = np.argmax(work != 0, axis=1)
    return tuple(int(n) for n in np.bincount(orders, minlength=degree + 1))


def search(bits, prime=3, steps=10000, seed=0, dilation=1, max_degree=None):
    """Completed finite energy walk, with exact endpoint-preserving moves.

    Moves add/subtract even multiples of affine-subspace indicators or
    differences of linear indicators.  Negative intermediate energies are
    rejected, as are new zero frequencies.  This is not an exhaustive search.
    """
    import numpy as np
    import random
    rng = random.Random(seed)
    base, shift = folded_energies(bits)
    energies = np.array([e * dilation for e in base], dtype=np.int32)
    if dilation % 2 == 0:
        raise ValueError("dilation must be odd")
    q = len(energies)
    if max_degree is None:
        max_degree = max(energies) + 2 * bits
    masks = np.array([[parity(v & a) for v in range(q)]
                      for a in range(q)], dtype=np.int32)
    cache = {}
    histogram = _numpy_histogram(energies, prime, cache)
    current = optimum(histogram)[0]
    best, best_energies = current, energies.copy()
    distinct = {histogram}
    accepted = evaluated = 0
    for step in range(steps):
        choice = rng.randrange(5)
        if choice < 2:
            mask = masks[rng.randrange(1, q)].copy()
        elif choice == 2:
            mask = masks[rng.randrange(1, q)] - masks[rng.randrange(1, q)]
        else:
            mask = masks[rng.randrange(1, q)].copy()
            for _ in range(rng.randrange(1, min(4, 2 * bits))):
                other = masks[rng.randrange(1, q)]
                if rng.randrange(2):
                    other = 1 - other
                mask *= other
        amount = 2 * (prime ** rng.randrange(2)) * (1 if rng.randrange(2) else -1)
        candidate = energies + amount * mask
        if candidate[0] != 0 or np.any(candidate[1:] <= 0) or max(candidate) > max_degree:
            continue
        target = _numpy_histogram(candidate, prime, cache)
        value = optimum(target)[0]
        evaluated += 1
        distinct.add(target)
        temperature = (0.015 + 0.12 * (1 - step / max(1, steps))) * log(q)
        if value <= current or rng.random() < exp((current - value) / temperature):
            energies, current = candidate, value
            accepted += 1
        if value < best - 1e-12:
            best, best_energies = value, candidate.copy()
        if (step + 1) % 2048 == 0:
            energies, current = best_energies.copy(), best
    result = receipt(tuple(int(x) for x in best_energies), bits, prime, shift)
    result["search"] = {
        "seed": seed, "steps_completed": steps, "evaluated": evaluated,
        "accepted": accepted, "distinct_valuation_histograms": len(distinct),
        "dilation": dilation, "maximum_degree": int(max_degree),
        "scope": "finite heuristic walk; no global optimality assertion",
    }
    return result


def main():
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bits", type=int, default=3)
    parser.add_argument("--prime", type=int, default=3)
    parser.add_argument("--steps", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261101)
    parser.add_argument("--dilation", type=int, default=1)
    parser.add_argument("--max-degree", type=int)
    parser.add_argument("--output")
    args = parser.parse_args()
    result = search(args.bits, args.prime, args.steps, args.seed,
                    args.dilation, args.max_degree)
    if args.output:
        Path(args.output).write_text(json.dumps(result, indent=2) + "\n")
    summary = {key: value for key, value in result.items()
               if key not in ("energies", "coefficients_by_degree")}
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
