#!/usr/bin/env python3
"""Construct polynomial Bell paths and minimize a uniform rank-pressure bound.

Search is heuristic. For every retained path the rational pressure bound
is exact, and applies to every tensor length by the proof in the note.
"""

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import json
from math import ceil, exp, log, log2
import random
import time

from bell_crossing import sign, transform


def profiles(coefficients: tuple[tuple[int, ...], ...]
             ) -> tuple[tuple[int, ...], tuple[int, ...]]:
    degree, q = len(coefficients) - 1, len(coefficients[0])
    low, high = [0] * (degree + 1), [0] * (degree + 1)
    for j in range(q):
        low[next(i for i in range(degree + 1) if coefficients[i][j])] += 1
        high[next(i for i in range(degree, -1, -1) if coefficients[i][j])] += 1
    return tuple(low), tuple(high)


@lru_cache(None)
def optimum(low: tuple[int, ...], high: tuple[int, ...]) -> tuple[float, float]:
    """Floating search only; retained bounds are recomputed over Q."""
    def value(y: float) -> float:
        a = sum(n * exp(-i * y) for i, n in enumerate(low))
        b = sum(n * exp(i * y) for i, n in enumerate(high))
        return log(a) + log(b)
    left, right = 0.0, 20.0
    for _ in range(60):
        a = (2 * left + right) / 3
        b = (left + 2 * right) / 3
        if value(a) <= value(b):
            right = b
        else:
            left = a
    middle = (left + right) / 2
    if value(0) <= value(middle):
        return value(0), 1.0
    return value(middle), exp(middle)


def pressure(low: tuple[int, ...], high: tuple[int, ...], x: Fraction) -> Fraction:
    return (sum(Fraction(n) / x ** i for i, n in enumerate(low))
            * sum(Fraction(n) * x ** i for i, n in enumerate(high)))


def tensor_start(bits: int, degree: int) -> tuple[tuple[int, ...], ...]:
    if degree < bits:
        raise ValueError("tensor initialization needs degree >= bits")
    d, q = 1 << bits, 1 << (2 * bits)
    result = [[0] * q for _ in range(degree + 1)]
    for a in range(d):
        for b in range(d):
            polynomial = [1]
            for j in range(bits):
                local = (2, 2) if (a >> j & 1) and (b >> j & 1) else (0, 1)
                new = [0] * (len(polynomial) + 1)
                for i, value in enumerate(polynomial):
                    new[i] = (new[i] + value * local[0]) % 3
                    new[i + 1] = (new[i + 1] + value * local[1]) % 3
                polynomial = new
            for i, value in enumerate(polynomial):
                result[i][a * d + b] = value
    return tuple(map(tuple, result))


def record(coefficients: tuple[tuple[int, ...], ...], bits: int) -> dict:
    other = tuple(transform(row, bits, 3) for row in coefficients)
    assert all(sum(col) % 3 == 1 for col in zip(*coefficients))
    assert all(sum(col) % 3 == 1 for col in zip(*other))
    row_low, row_high = profiles(coefficients)
    col_low, col_high = profiles(other)
    options = []
    for orientation, low, high in (("row-low/column-high", row_low, col_high),
                                   ("column-low/row-high", col_low, row_high)):
        _, approximate_x = optimum(low, high)
        x = Fraction(approximate_x).limit_denominator(256)
        bound = pressure(low, high, x)
        rho = (log2(bound.numerator) - log2(bound.denominator)) / (4 * bits)
        upper_numerator = ceil(1000 * rho)
        while bound ** 1000 > 2 ** (4 * bits * upper_numerator):
            upper_numerator += 1
        options.append({"orientation": orientation, "low_counts": low,
                        "high_counts": high, "scale": [x.numerator, x.denominator],
                        "pressure": [bound.numerator, bound.denominator],
                        "rho_approx": rho,
                        "rho_upper_exact": [upper_numerator, 1000],
                        "strict_below_critical_certified":
                            bound ** 577 < 2 ** (4 * bits * 408)})
    return {"coefficients_by_degree": coefficients,
            "transformed_coefficients_by_degree": other,
            "row_low": row_low, "row_high": row_high,
            "column_low": col_low, "column_high": col_high,
            "bounds": options}


def run(bits: int, degree: int, steps: int, seed: int) -> dict:
    start = time.monotonic()
    q, d = 1 << (2 * bits), 1 << bits
    rng = random.Random(seed)
    signs = [[pow(d, -1, 3) * sign((a ^ c) & (b ^ e)) % 3
              for a in range(d) for b in range(d)]
             for c in range(d) for e in range(d)]
    initial = tensor_start(bits, degree)
    coefficients = initial
    other = tuple(transform(row, bits, 3) for row in coefficients)
    evaluations = 0

    def objective(source, dual):
        nonlocal evaluations
        evaluations += 1
        rl, rh = profiles(source)
        cl, ch = profiles(dual)
        return min(optimum(rl, ch)[0], optimum(cl, rh)[0])

    current = best = objective(coefficients, other)
    witness = coefficients
    improvements = []
    for step in range(steps):
        if step % 1024 == 0:
            coefficients = initial if step % 4096 == 0 else witness
            other = tuple(transform(row, bits, 3) for row in coefficients)
            current = objective(coefficients, other)
        source = [list(row) for row in coefficients]
        dual = [list(row) for row in other]
        if rng.random() < 0.55:
            a, b = rng.sample(range(q), 2)
            for i in range(degree + 1):
                delta = source[i][b] - source[i][a]
                source[i][a], source[i][b] = source[i][b], source[i][a]
                if delta:
                    dual[i] = [(x + delta * (row[a] - row[b])) % 3
                               for x, row in zip(dual[i], signs)]
        else:
            i, j = rng.randrange(degree), rng.randrange(q)
            delta = rng.choice((1, 2))
            source[i][j] = (source[i][j] + delta) % 3
            source[degree][j] = (source[degree][j] - delta) % 3
            dual[i] = [(x + delta * row[j]) % 3 for x, row in zip(dual[i], signs)]
            dual[degree] = [(x - delta * row[j]) % 3
                            for x, row in zip(dual[degree], signs)]
        source, dual = tuple(map(tuple, source)), tuple(map(tuple, dual))
        value = objective(source, dual)
        temperature = 0.03 * (1 - (step % 1024) / 1024) + 0.0001
        if value < current or rng.random() < exp(min(0, (current - value) / temperature)):
            coefficients, other, current = source, dual, value
        if value < best - 1e-12:
            best, witness = value, source
            improvements.append({"step": step, "rho_search_approx": best / (4 * bits * log(2))})
    assert tuple(transform(row, bits, 3) for row in witness) == tuple(
        tuple(sum(r[j] * values[j] for j in range(q)) % 3 for r in signs)
        for values in witness)
    return {"status": "heuristic search completed", "bits": bits,
            "local_dimension": d, "degree": degree, "coefficient_prime": 3,
            "steps": steps, "seed": seed, "candidate_evaluations": evaluations,
            "profile_cache_entries": optimum.cache_info().currsize,
            "initial": record(initial, bits), "best": record(witness, bits),
            "improvements": improvements, "elapsed_seconds": time.monotonic() - start,
            "scope": "Heuristic path search; rational bounds are uniform in tensor length"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bits", type=int, default=3)
    parser.add_argument("--degree", type=int, default=3)
    parser.add_argument("--steps", type=int, default=32768)
    parser.add_argument("--seed", type=int, default=20261031)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = run(args.bits, args.degree, args.steps, args.seed)
    with open(args.output, "w") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in (
        "status", "local_dimension", "degree", "candidate_evaluations", "elapsed_seconds")}))
    print(json.dumps({"bounds": result["best"]["bounds"], "improvements": result["improvements"]}))


if __name__ == "__main__":
    main()
