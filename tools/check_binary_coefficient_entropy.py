#!/usr/bin/env python3
"""Exact inclusion ranks and binary crossing coefficient entropy checks."""
from functools import lru_cache
from itertools import combinations, product
from math import comb, log2
import json


def binary_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            index = row.bit_length() - 1
            if index not in pivots:
                pivots[index] = row
                break
            row ^= pivots[index]
    return len(pivots)


@lru_cache(None)
def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def inclusion_formula(n, low, high):
    gap = high - low
    bound = min(low, n - high)
    return sum(choose(n, index) - choose(n, index - 1)
               for index in range(bound + 1)
               if not (gap & (bound - index)))


def inclusion_direct(n, low, high):
    inputs = [sum(1 << j for j in subset)
              for subset in combinations(range(n), low)]
    outputs = [sum(1 << j for j in subset)
               for subset in combinations(range(n), high)]
    return binary_rank([sum(1 << j for j, target in enumerate(outputs)
                            if source & target == source)
                        for source in inputs])


@lru_cache(None)
def elementary_rank(n, degree):
    if not 0 <= degree <= n:
        return 0
    return sum(inclusion_formula(n, low, low + degree)
               for low in range(n - degree + 1))


def coefficient_rank(t, index):
    return sum(choose(t, k) * (1 << (t-k)) * elementary_rank(k, t-index)
               for k in range(t-index, t+1))


def direct_crossing(t):
    size = 1 << t
    rows = [[0] * (size*size) for _ in range(t+1)]
    cols = [[0] * (size*size) for _ in range(t+1)]
    for x, y, u, v in product(range(size), repeat=4):
        polynomial = 1
        for bit in range(t):
            a, b, c, d = ((value >> bit) & 1 for value in (x, y, u, v))
            rr = int(a == b and c == d)
            ee = int(a == d and b == c)
            local = rr | ((rr ^ ee) << 1)
            polynomial = ((polynomial if local & 1 else 0)
                          ^ (polynomial << 1 if local & 2 else 0))
        for index in range(t+1):
            if (polynomial >> index) & 1:
                rows[index][size*x+y] |= 1 << (size*v+u)
                cols[index][size*x+u] |= 1 << (size*v+y)
    return [(binary_rank(row), binary_rank(col)) for row, col in zip(rows, cols)]


def verify():
    inclusion_cases = 0
    for n in range(11):
        for low in range(n+1):
            for high in range(low, n+1):
                assert inclusion_formula(n, low, high) == inclusion_direct(n, low, high)
                inclusion_cases += 1
    complete_cases = 0
    for t in range(1, 5):
        actual = direct_crossing(t)
        expected = [(coefficient_rank(t, index), coefficient_rank(t, t-index))
                    for index in range(t+1)]
        assert actual == expected, (t, actual, expected)
        complete_cases += len(actual)

    integer_bounds = 0
    profiles = []
    for t in range(1, 129):
        ranks = [coefficient_rank(t, index) for index in range(t+1)]
        tail = 0
        for index, value in enumerate(ranks):
            tail += choose(2*t, index)
            assert tail <= (t+1)*value <= 2*(t+1)*tail
            assert 1 <= value <= 1 << (2*t)
            integer_bounds += 1
        if t & (t-1) == 0:
            maximum = max(a*b for a, b in zip(ranks, reversed(ranks)))
            profiles.append({'bits': t, 'maximum_rank_product': maximum,
                             'finite_rho': log2(maximum)/(4*t),
                             'central_ranks': [ranks[t//2], ranks[t-t//2]]})

    # Ignoring characteristic-two cancellations gives the wrong rank.
    assert inclusion_direct(4, 1, 2) == 3
    assert inclusion_direct(4, 1, 2) < choose(4, 1)
    entropy = -0.25*log2(0.25)-0.75*log2(0.75)
    return {'status': 'pass',
            'scope': 'exact ranks for the literal full-cube coefficient family only',
            'direct_inclusion_rank_cases': inclusion_cases,
            'complete_crossing_profiles': complete_cases,
            'exact_tail_bound_cases': integer_bounds,
            'limiting_rho': entropy,
            'selected_profiles': profiles,
            'negative_control': 'characteristic-zero full inclusion rank fails over F2',
            'new_factoring_exponent_claimed': False}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2, sort_keys=True))
