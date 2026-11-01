#!/usr/bin/env python3
"""Exact finite checks for streamed type and polynomial seed factors."""

import argparse
from functools import lru_cache
from itertools import product
import json
from math import factorial
from pathlib import Path


def pmul(a, b):
    """Multiply binary polynomials represented by coefficient bit sets."""
    answer = 0
    while b:
        if b & 1:
            answer ^= a
        a <<= 1
        b >>= 1
    return answer


def pproduct(values):
    answer = 1
    for value in values:
        answer = pmul(answer, value)
    return answer


def mmul(a, b):
    return [[xor(pmul(a[i][k], b[k][j]) for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def xor(values):
    answer = 0
    for value in values:
        answer ^= value
    return answer


def inverse(a):
    n = len(a)
    rows = [sum(a[i][j] << j for j in range(n)) | (1 << (n+i))
            for i in range(n)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if (rows[i] >> j) & 1)
        rows[j], rows[pivot] = rows[pivot], rows[j]
        for i in range(n):
            if i != j and (rows[i] >> j) & 1:
                rows[i] ^= rows[j]
    return [[(rows[i] >> (n+j)) & 1 for j in range(n)] for i in range(n)]


def diagonal(values):
    return [[value if i == j else 0 for j in range(len(values))]
            for i, value in enumerate(values)]


def basis_matrix(columns):
    return [[(column >> i) & 1 for column in columns]
            for i in range(len(columns))]


def flatten_column(matrix, d):
    result = [[0]*(d*d) for _ in range(d*d)]
    for x, y, u, v in product(range(d), repeat=4):
        result[x*d+u][y*d+v] = matrix[x*d+y][u*d+v]
    return result


def seed_polynomial():
    result = [[0]*4 for _ in range(4)]
    for x, y, u, v in product(range(2), repeat=4):
        result[2*x+y][2*u+v] = (2 if x == v and y == u else 0) ^ (
            3 if x == y and u == v else 0)
    return result


def polynomial_certificates():
    valuations = (0, 2, 1, 1)
    right_block = [[2, 1, 0, 0], [1, 0, 0, 0],
                   [0, 0, 1, 0], [0, 0, 0, 1]]
    row_left = [[3, 1, 0, 0], [2, 1, 0, 0],
                [0, 0, 1, 0], [0, 0, 0, 1]]
    col_left = [[1, 0, 0, 0], [2, 1, 0, 0],
                [0, 0, 1, 0], [0, 0, 0, 1]]
    row_basis = basis_matrix((9, 1, 2, 4))
    col_basis = basis_matrix((6, 2, 1, 8))
    swap = [[int(i//2 == j % 2 and i % 2 == j//2)
             for j in range(4)] for i in range(4)]
    row_l = mmul(row_basis, row_left)
    row_r = mmul(mmul(right_block, inverse(row_basis)), swap)
    col_l = mmul(col_basis, col_left)
    col_r = mmul(right_block, inverse(col_basis))
    row = seed_polynomial()
    col = flatten_column(row, 2)
    reverse_col = [[((v & 1) << 1) | ((v >> 1) & 1) for v in r]
                   for r in col]
    diag = diagonal([1 << v for v in valuations])
    assert mmul(mmul(row_l, diag), row_r) == row
    assert mmul(mmul(col_l, diag), col_r) == reverse_col
    return row, (row_l, row_r), (col_l, col_r), valuations


def words(label, t, base):
    result = []
    for _ in range(t):
        label, digit = divmod(label, base)
        result.append(digit)
    assert label == 0
    return tuple(reversed(result))


@lru_cache(None)
def budget_count(m, budget, valuations):
    if budget < 0:
        return 0
    if m == 0:
        return budget+1
    return sum(budget_count(m-1, budget-v, valuations) for v in valuations)


def budget_unrank(t, budget, valuations, index):
    assert 0 <= index < budget_count(t, budget, valuations)
    state = []
    for m in range(t-1, -1, -1):
        for i, value in enumerate(valuations):
            count = budget_count(m, budget-value, valuations)
            if index < count:
                state.append(i)
                budget -= value
                break
            index -= count
        else:
            raise AssertionError("channel unranking failed")
    assert 0 <= index <= budget
    return tuple(state), index


def polynomial_factor_tables(t, budget, left, right, valuations):
    """Query each small example channel independently; no factor elimination."""
    d = 2**t
    channels = budget_count(t, budget, valuations)
    ltab = [[0]*channels for _ in range(d*d)]
    rtab = [[0]*(d*d) for _ in range(channels)]
    direct_channels = [(state, k)
                       for state in product(range(len(valuations)), repeat=t)
                       for k in range(max(0, budget-sum(valuations[i]
                                                       for i in state)+1))]
    for channel in range(channels):
        state, k = budget_unrank(t, budget, valuations, channel)
        assert (state, k) == direct_channels[channel]
        residual = budget-sum(valuations[i] for i in state)-k
        for a, b in product(range(d), repeat=2):
            da, db = words(a, t, 2), words(b, t, 2)
            ltab[a*d+b][channel] = (pproduct(
                left[2*da[j]+db[j]][state[j]] for j in range(t)) >> k) & 1
            rtab[channel][a*d+b] = (pproduct(
                right[state[j]][2*da[j]+db[j]] for j in range(t)) >> residual) & 1
    return ltab, rtab


def check_polynomial():
    seed, row_factors, col_factors, vals = polynomial_certificates()
    checks = 0
    profiles = []
    for t in range(1, 4):
        d = 2**t
        profile = []
        for j in range(t+1):
            rl, rr = polynomial_factor_tables(t, j, *row_factors, vals)
            cl, cr = polynomial_factor_tables(t, t-j, *col_factors, vals)
            row_matrix, col_matrix = mmul(rl, rr), mmul(cl, cr)
            for x, y, u, v in product(range(d), repeat=4):
                dx, dy, du, dv = (words(a, t, 2) for a in (x, y, u, v))
                expected = (pproduct(seed[2*dx[l]+dy[l]][2*du[l]+dv[l]]
                                     for l in range(t)) >> j) & 1
                assert row_matrix[x*d+y][u*d+v] == expected
                assert col_matrix[x*d+u][y*d+v] == expected
                checks += 2
            r, c = len(rr), len(cr)
            # B=256/27 at x=1/3. Clear denominators exactly.
            assert r*c*27**t <= (t+1)**2*256**t
            profile.append([r, c])
        profiles.append({"copies": t, "channel_counts": profile})
    # The worked example fails the strict rho<7/10 sufficient test.
    assert 256**5 > 2**14*27**5
    return {"polynomial_identities": 2, "factor_entry_checks": checks,
            "tilt": [1, 3], "bound_B": [256, 27],
            "strict_seven_tenths_test": False, "profiles": profiles}


def rank_factor(matrix):
    """Exact row-basis factorization over F2, independent of seed formulas."""
    pivots, basis = {}, []
    for row in matrix:
        value = sum(v << j for j, v in enumerate(row))
        reduced, coordinate = value, 1 << len(basis)
        while reduced:
            p = reduced.bit_length()-1
            if p not in pivots:
                pivots[p] = reduced, coordinate
                basis.append(value)
                break
            other, mask = pivots[p]
            reduced ^= other
            coordinate ^= mask
    left = []
    for row in matrix:
        value, coordinate = sum(v << j for j, v in enumerate(row)), 0
        while value:
            other, mask = pivots[value.bit_length()-1]
            value ^= other
            coordinate ^= mask
        left.append([(coordinate >> j) & 1 for j in range(len(basis))])
    right = [[(row >> j) & 1 for j in range(len(matrix[0]))] for row in basis]
    assert mmul(left, right) == matrix
    return left, right


def multinomial(counts):
    answer = factorial(sum(counts))
    for n in counts:
        answer //= factorial(n)
    return answer


def type_unrank(counts, ranks, index):
    local_count = 1
    for n, rank in zip(counts, ranks):
        local_count *= rank**n
    word_index, local_index = divmod(index, local_count)
    remaining, word = list(counts), []
    for _ in range(sum(counts)):
        for i, n in enumerate(remaining):
            if n == 0:
                continue
            remaining[i] -= 1
            count = multinomial(remaining)
            if word_index < count:
                word.append(i)
                break
            word_index -= count
            remaining[i] += 1
        else:
            raise AssertionError("type unranking failed")
    local = []
    for i in reversed(word):
        local_index, digit = divmod(local_index, ranks[i])
        local.append(digit)
    assert local_index == word_index == 0
    return tuple(word), tuple(reversed(local))


def type_entry(seed, counts, coords):
    # Independent coefficient DP in two variables, using one count coordinate.
    t = sum(counts)
    digits = [words(a, t, 3) for a in coords]
    coefficients = [1]+[0]*t
    for l in range(t):
        x, y, u, v = (a[l] for a in digits)
        values = [a[3*x+y][3*u+v] for a in seed]
        nxt = [0]*(t+1)
        for j in range(l+1):
            nxt[j] ^= coefficients[j]*values[0]
            nxt[j+1] ^= coefficients[j]*values[1]
        coefficients = nxt
    return coefficients[counts[1]]


def check_types():
    seed = [[[0]*9 for _ in range(9)] for _ in range(2)]
    for x, y, u, v in product(range(3), repeat=4):
        r = int(x == y and u == v)
        seed[0][3*x+y][3*u+v] = r
        seed[1][3*x+y][3*u+v] = r ^ int(x == v and y == u)
    row_factors = [rank_factor(a) for a in seed]
    col_factors = [rank_factor(flatten_column(a, 3)) for a in seed]
    checks, capped = 0, 0
    profiles = []
    d = 8  # Binary prefix of the full nine-label two-copy space.
    for counts in ((2, 0), (1, 1), (0, 2)):
        matrices = []
        profile = []
        for cut, factors in enumerate((row_factors, col_factors)):
            ranks = [len(r) for _, r in factors]
            channel_count = multinomial(counts)
            for n, r in zip(counts, ranks):
                channel_count *= r**n
            profile.append(channel_count)
            if channel_count > d*d:
                # The identity-factor cap is checked entry by entry below.
                left = diagonal([1]*(d*d))
                right = [[0]*(d*d) for _ in range(d*d)]
                for x, y, u, v in product(range(d), repeat=4):
                    a, b = (x*d+y, u*d+v) if cut == 0 else (x*d+u, y*d+v)
                    right[a][b] = type_entry(seed, counts, (x, y, u, v))
                capped += 1
            else:
                left = [[0]*channel_count for _ in range(d*d)]
                right = [[0]*(d*d) for _ in range(channel_count)]
                seen = set()
                for channel in range(channel_count):
                    word, local = type_unrank(counts, ranks, channel)
                    assert (word, local) not in seen
                    seen.add((word, local))
                    assert tuple(word.count(i) for i in range(2)) == counts
                    for a, b in product(range(d), repeat=2):
                        da, db = words(a, 2, 3), words(b, 2, 3)
                        lv = rv = 1
                        for l in range(2):
                            fl, fr = factors[word[l]]
                            lv *= fl[3*da[l]+db[l]][local[l]]
                            rv *= fr[local[l]][3*da[l]+db[l]]
                        left[a*d+b][channel] = lv
                        right[channel][a*d+b] = rv
            matrices.append(mmul(left, right))
        for x, y, u, v in product(range(d), repeat=4):
            expected = type_entry(seed, counts, (x, y, u, v))
            assert matrices[0][x*d+y][u*d+v] == expected
            assert matrices[1][x*d+u][y*d+v] == expected
            checks += 2
        profiles.append({"type": counts, "uncapped_channels": profile})
    for x, y, u, v in product(range(d), repeat=4):
        assert xor(type_entry(seed, n, (x, y, u, v))
                   for n in ((2, 0), (1, 1), (0, 2))) == int(x == v and y == u)
    return {"seed_labels": 3, "copies": 2, "restricted_labels": d,
            "seed_rank_pairs": [[len(row_factors[i][1]), len(col_factors[i][1])]
                                for i in range(2)],
            "factor_entry_checks": checks, "identity_caps": capped,
            "crossing_entries": d**4, "profiles": profiles}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert 2**3544 > 17**867
    assert 2*577**2 > 816**2
    report = {"status": "PASS", "polynomial_certificate": check_polynomial(),
              "type_certificate": check_types(),
              "four_label_smith_exclusion": {
                  "rational_separation": [577, 816],
                  "power_inequality": "2^3544 > 17^867",
                  "scope": "usual cumulative Smith envelope only"},
              "scope": "Finite exact controls; no improved time exponent claimed."}
    encoded = json.dumps(report, indent=2)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
