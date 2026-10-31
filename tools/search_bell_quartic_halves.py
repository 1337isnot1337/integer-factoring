#!/usr/bin/env python3
"""Search coordinate-sparse halves of a self-reciprocal quartic Bell path.

The free vector b determines the midpoint 1-b-Sb. This is the matrix
decomposition in research/bell-quartic-rigidity.md, with exact endpoint.
"""

import argparse
import json
from math import exp
import os
from pathlib import Path
import random
import time

from bell_crossing import sign, transform
from folded_bell import coefficient_vectors, coordinates
from search_bell_profiles import optimum, record


def balanced_pool(bits, shift, rng):
    """Affine Fourier moves whose two support sizes have product d^2."""
    d, q = 1 << bits, 1 << (2 * bits)
    pool = []
    for _ in range(min(q, 256)):
        row, frequency = rng.randrange(d), rng.randrange(d)
        dimension = rng.randrange(bits + 1)
        mask = (1 << dimension) - 1
        source = [(row * d + b, sign(frequency & b) % 3) for b in range(1 << dimension)]
        scale = pow(2, dimension, 3) * pow(d, -1, 3) % 3
        dual = [(c * d + e, scale * sign((row ^ c) & e) % 3)
                for c in range(d) if not ((row ^ c ^ frequency) & mask)
                for e in range(d)]
        assert len(source) * len(dual) == q
        for _ in range(4):
            t = rng.randrange(q)
            while sign((t // d) & (t % d)) != -1:
                t = rng.randrange(q)
            def move(index):
                polar = ((index // d & t % d).bit_count()
                         + (index % d & t // d).bit_count()) % 2
                return index ^ t if polar else index
            source = [(move(j), c) for j, c in source]
            dual = [(move(j), c) for j, c in dual]
        dual = [(j ^ shift, c) for j, c in dual]
        pool.append((tuple(source), tuple(dual)))
    return tuple(pool)


def run(bits, steps, seed, shift=None, outer=None, checkpoint=None, balanced=False):
    start = time.monotonic()
    d, q = 1 << bits, 1 << (2 * bits)
    _, folded_shift = coordinates(bits)
    shift = folded_shift if shift is None else shift
    assert 0 <= shift < q
    folded = coefficient_vectors(bits, 3)
    base = [tuple(values[j ^ folded_shift] for j in range(q)) for values in folded]
    if bits == 3:
        # Multiply the known cubic by (1+z)/2 to obtain a quartic start.
        initial = tuple((base[0][j] + 2 * base[1][j]) % 3 for j in range(q))
        default_outer = 2 * base[0][0] % 3
    else:
        initial = tuple((base[0][j] + base[1][j]) % 3 for j in range(q))
        default_outer = base[0][0]
    outer = default_outer if outer is None else outer
    assert outer in (1, 2)
    a0 = (outer,) + (0,) * (q - 1)
    def dual(values):
        return transform(tuple(values[j ^ shift] for j in range(q)), bits, 3)
    a4 = dual(a0)
    columns = [tuple(pow(d, -1, 3) * sign(((j ^ shift) // d ^ c) & ((j ^ shift) % d ^ e)) % 3
                     for c in range(d) for e in range(d)) for j in range(q)]
    rng = random.Random(seed)
    moves = balanced_pool(bits, shift, rng) if balanced else ()
    for source, target in moves[:8]:
        vector, expected = [0] * q, [0] * q
        for j, value in source:
            vector[j] = value
        for j, value in target:
            expected[j] = value
        assert dual(vector) == tuple(expected)
    evaluations = 0
    seen = set()
    improvements = []

    def evaluate(b, sb):
        nonlocal evaluations
        counts = [1, 0, 0, 0, 0]
        for j in range(1, q):
            if b[j]:
                counts[1] += 1
            elif sb[j] != 1:
                counts[2] += 1
            elif a4[j] != 1:
                counts[3] += 1
            else:
                counts[4] += 1
        low = tuple(counts)
        seen.add(low)
        evaluations += 1
        value = optimum(low, tuple(reversed(low)))[0]
        mean_order = sum(j * n for j, n in enumerate(low)) / q
        # Convex pressure is flat at rho=1 when mean order <= degree/2.
        # Guide that plateau by its exact necessary mean-order deficit.
        return value, low, value + max(0., 2 - mean_order)

    def receipt(b, completed, status):
        sb = dual(b)
        rows = (a0, tuple((b[j] - a0[j]) % 3 for j in range(q)),
                tuple((1 - b[j] - sb[j]) % 3 for j in range(q)),
                tuple((sb[j] - a4[j]) % 3 for j in range(q)), a4)
        assert tuple(dual(values) for values in rows) == tuple(reversed(rows))
        physical = tuple(tuple(values[j ^ shift] for j in range(q)) for values in rows)
        return {"status": status, "bits": bits, "degree": 4, "seed": seed,
                "steps": steps, "steps_completed": completed, "shift": shift,
                "outer_value": outer, "candidate_evaluations": evaluations,
                "balanced_moves": balanced,
                "move_model": "multiscale affine Fourier moves and coordinate moves; mean-deficit-guided annealing",
                "distinct_valuation_profiles": len(seen), "sparse_half": b,
                "improvements": improvements, "best": record(physical, bits),
                "elapsed_seconds": time.monotonic() - start,
                "scope": "Specified heuristic sparse-half schedule; every retained path has exact endpoint and rational uniform pressure"}

    b, sb = initial, dual(initial)
    value, _, current = evaluate(b, sb)
    best, witness = value, b
    for step in range(steps):
        if step % 8192 == 0:
            b = initial if step % 32768 == 0 else witness
            sb = dual(b)
            _, _, current = evaluate(b, sb)
        delta = rng.choice((1, 2))
        candidate = list(b)
        if balanced and rng.random() < .65:
            source, target = rng.choice(moves)
            other = list(sb)
            for j, coefficient in source:
                candidate[j] = (candidate[j] + delta * coefficient) % 3
            for j, coefficient in target:
                other[j] = (other[j] + delta * coefficient) % 3
            other = tuple(other)
        else:
            index = rng.randrange(q)
            candidate[index] = (candidate[index] + delta) % 3
            other = tuple((x + delta * y) % 3 for x, y in zip(sb, columns[index]))
        candidate = tuple(candidate)
        value, _, energy = evaluate(candidate, other)
        temperature = .3 * (1 - (step % 8192) / 8192) ** 2 + .0001
        if energy <= current or rng.random() < exp(min(0, (current - energy) / temperature)):
            b, sb, current = candidate, other, energy
        if value < best - 1e-12:
            best, witness = value, candidate
            improvements.append({"step": step, "log_pressure_approx": value})
        if checkpoint and (step == 0 or (step + 1) % 1024 == 0):
            checkpoint(receipt(witness, step + 1, "sparse-half search in progress"))
    return receipt(witness, steps, "sparse-half schedule completed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bits", type=int, default=4)
    parser.add_argument("--steps", type=int, default=65536)
    parser.add_argument("--seed", type=int, default=20261101)
    parser.add_argument("--shift", type=int)
    parser.add_argument("--outer", type=int, choices=(1, 2))
    parser.add_argument("--output", required=True)
    parser.add_argument("--fourier-moves", "--balanced-moves", dest="balanced_moves", action="store_true")
    args = parser.parse_args()
    output = Path(args.output)
    def save(value):
        temporary = output.with_name(output.name + ".tmp")
        temporary.write_text(json.dumps(value, indent=2) + "\n")
        os.replace(temporary, output)
    result = run(args.bits, args.steps, args.seed, args.shift, args.outer, save,
                 args.balanced_moves)
    save(result)
    print(json.dumps({"status": result["status"], "completed": result["steps_completed"],
                      "profiles": result["distinct_valuation_profiles"],
                      "bounds": result["best"]["bounds"], "improvements": result["improvements"]}))


if __name__ == "__main__":
    main()
