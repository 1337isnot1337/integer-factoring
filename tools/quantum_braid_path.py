#!/usr/bin/env python3
"""Height-three braid path quotient and a strict F4 crossing polynomial."""

from __future__ import annotations

from bounded_spin_crossing import binary_path, paths
from check_fourier_multiplication_gates import BinaryField
from quantum_block_braid import braid_word

FIELD = BinaryField(0b111)
MULTIPLY = tuple(tuple(FIELD.multiply(a, b) for b in range(4)) for a in range(4))
# alpha=2, alpha^2=alpha+1. In the ambient F16 model alpha=q+q^-1=7.
EMBED = (0, 1, 7, 6)


def local_action(path, site):
    if path[site] == path[site + 1]:
        return ()
    h = sum(path[:site])
    scalar = {(0, 1): 2, (1, 1): 1, (1, -1): 3, (2, 1): 3, (2, -1): 1, (3, -1): 2}[
        h, path[site]
    ]
    result = [(path, scalar)]
    if (path[site] == 1 and h > 0) or (path[site] == -1 and h < 3):
        result.append(
            (path[:site] + (path[site + 1], path[site]) + path[site + 2 :], scalar)
        )
    return tuple(result)


def crossing_weights(length):
    phase = length * length // 2 % 5
    # Choose the F4-linear projection pi(1)=1, pi(q)=tau with pi(q^phase)=1.
    tau = (0, 1, 0, 2, 3)[phase]
    constants, linear = (1, 0, 1, 2, 2), (0, 1, 2, 2, 1)
    values = [a ^ MULTIPLY[tau][b] for a, b in zip(constants, linear)]
    result = tuple(values[(phase - j) % 5] for j in range(length * length + 1))
    assert result[0] == 1
    return result


def coefficient_entry(four_paths):
    """Coefficients of the middle-zero block of product(I+s e_i), before descent."""
    n = len(four_paths[0])
    assert n > 0 and n % 2 == 0 and all(len(p) == n and sum(p) == 0 for p in four_paths)
    initial, final = four_paths[0] + four_paths[1], four_paths[2] + four_paths[3]
    states, degree = {initial: [1]}, 0
    for site in braid_word(n, n):
        degree += 1
        following = {p: values + [0] for p, values in states.items()}
        for path, values in states.items():
            for target, scalar in local_action(path, site):
                out = following.setdefault(target, [0] * (degree + 1))
                for j, value in enumerate(values):
                    out[j + 1] ^= MULTIPLY[scalar][value]
        states = following
    return states.get(final, [0] * (n * n + 1))


def polynomial_entry(four_paths):
    values = coefficient_entry(four_paths)
    return [
        MULTIPLY[w][v] for w, v in zip(crossing_weights(len(four_paths[0])), values)
    ]


def binary_entry(width, labels):
    return polynomial_entry(tuple(binary_path(width, label) for label in labels))


def coefficient_tensor(length):
    """Optional NumPy batch of every coefficient, axes degree,P,Q,R,S."""
    import numpy as np

    ps, qs = paths(length), paths(2 * length)
    indices = {p: i for i, p in enumerate(qs)}
    d = len(ps)
    boundary = np.array([indices[p + q] for p in ps for q in ps])
    multiply = np.array(MULTIPLY, dtype=np.uint8)
    states = np.zeros((len(qs), d * d, 1), dtype=np.uint8)
    states[boundary, np.arange(d * d), 0] = 1
    for site in braid_word(length, length):
        diagonal = np.zeros(len(qs), dtype=np.uint8)
        source, target = [], []
        for i, p in enumerate(qs):
            for q, c in local_action(p, site):
                if p == q:
                    diagonal[i] = c
                else:
                    source.append(i)
                    target.append(indices[q])
        source, target = np.array(source, dtype=np.int64), np.array(
            target, dtype=np.int64
        )
        following = np.zeros((len(qs), d * d, states.shape[2] + 1), dtype=np.uint8)
        following[:, :, :-1] = states
        following[:, :, 1:] ^= multiply[diagonal[:, None, None], states]
        # Each swapped target has one source, so indexed XOR has no collisions.
        following[target, :, 1:] ^= multiply[
            diagonal[source, None, None], states[source]
        ]
        states = following
    return states[boundary].transpose(2, 1, 0).reshape(length * length + 1, d, d, d, d)


def crossing_tensor(length):
    import numpy as np

    coefficients = coefficient_tensor(length)
    weights = np.array(crossing_weights(length), dtype=np.uint8)
    return np.array(MULTIPLY, dtype=np.uint8)[
        weights[:, None, None, None, None], coefficients
    ]
