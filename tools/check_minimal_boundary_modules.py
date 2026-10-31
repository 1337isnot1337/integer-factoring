#!/usr/bin/env python3
"""Exact reachable/observable reduction before module rate certification.

All arithmetic is over F_2(tau).  Physical choices remain grouped by channel:
the partition function counts a channel word once, using the minimum order
over its physical choices.  The checker exercises the representation change,
not a floating-point estimate of an asymptotic rate.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

import check_cancellation_module_automata as module


ZERO, ONE = module.ZERO, module.ONE


def scalar_product(left, right):
    answer = ZERO
    for a, b in zip(left, right):
        answer = module.radd(answer, module.rmul(a, b))
    return answer


def transpose(matrix):
    return tuple(tuple(row[j] for row in matrix) for j in range(len(matrix[0])))


def matrix_product(left, right):
    columns = transpose(right)
    return tuple(tuple(scalar_product(row, column) for column in columns) for row in left)


def row_product(row, matrix):
    return tuple(scalar_product(row, column) for column in transpose(matrix))


def column_product(matrix, column):
    return tuple(scalar_product(row, column) for row in matrix)


def row_basis(rows, width):
    """Canonical reduced row basis over F_2(tau), with leftmost pivots."""
    pivots = {}
    for source in rows:
        row = list(source)
        assert len(row) == width
        for column in sorted(pivots):
            if row[column] != ZERO:
                multiple = row[column]
                row = [module.radd(value, module.rmul(multiple, old))
                       for value, old in zip(row, pivots[column])]
        if all(value == ZERO for value in row):
            continue
        column = next(j for j, value in enumerate(row) if value != ZERO)
        divisor = row[column]
        row = tuple(module.rdiv(value, divisor) for value in row)
        for old_column in list(pivots):
            multiple = pivots[old_column][column]
            if multiple != ZERO:
                pivots[old_column] = tuple(
                    module.radd(value, module.rmul(multiple, new))
                    for value, new in zip(pivots[old_column], row))
        pivots[column] = row
    columns = tuple(sorted(pivots))
    basis = tuple(pivots[column] for column in columns)
    assert all(basis[i][column] == (ONE if i == j else ZERO)
               for i in range(len(basis)) for j, column in enumerate(columns))
    return basis, columns


def closure(seed, matrices, column=False):
    """Reachable row span, or observable column span represented as rows."""
    width = len(seed)
    basis, pivots = row_basis((seed,), width)
    rounds = 0
    while basis:
        images = [column_product(matrix, vector) if column else row_product(vector, matrix)
                  for vector in basis for matrix in matrices]
        following, new_pivots = row_basis((*basis, *images), width)
        rounds += 1
        assert rounds <= width
        if len(following) == len(basis):
            assert following == basis
            return basis, pivots
        basis, pivots = following, new_pivots
    return basis, pivots


def linear_combination(coefficients, rows):
    if not rows:
        return ()
    return tuple(scalar_product(coefficients, column) for column in transpose(rows))


def minimize(channels, initial, terminal):
    """Preserve every scalar word, then remove unobservable row directions.

    Two exact intertwining identities are asserted for every physical matrix.
    They certify equality of all lengths, independently of finite word tests.
    """
    dimension = len(initial)
    assert dimension == len(terminal)
    matrices = tuple(matrix for choices in channels for matrix in choices)
    reachable, pivots = closure(initial, matrices)
    if not reachable:
        return {"dimension": 0, "reachable_dimension": 0, "original_dimension": dimension,
                "channels": (), "initial": (), "terminal": (), "identities": 0}
    reduced_initial = tuple(initial[j] for j in pivots)
    reduced_terminal = column_product(reachable, terminal)
    assert linear_combination(reduced_initial, reachable) == tuple(initial)
    reachable_channels = []
    identities = 1
    for choices in channels:
        changed = []
        for matrix in choices:
            images = matrix_product(reachable, matrix)
            following = tuple(tuple(row[j] for j in pivots) for row in images)
            assert matrix_product(following, reachable) == images
            identities += 1
            changed.append(following)
        reachable_channels.append(tuple(changed))
    reachable_channels = tuple(reachable_channels)
    observable, observable_pivots = closure(
        reduced_terminal, tuple(matrix for choices in reachable_channels for matrix in choices),
        column=True)
    if not observable:
        return {"dimension": 0, "reachable_dimension": len(reachable), "original_dimension": dimension,
                "channels": (), "initial": (), "terminal": (), "identities": identities}
    quotient = transpose(observable)
    final_initial = row_product(reduced_initial, quotient)
    final_terminal = tuple(reduced_terminal[j] for j in observable_pivots)
    assert column_product(quotient, final_terminal) == reduced_terminal
    identities += 1
    final_channels = []
    for choices in reachable_channels:
        changed = []
        for matrix in choices:
            images = matrix_product(matrix, quotient)
            following = tuple(images[j] for j in observable_pivots)
            assert matrix_product(quotient, following) == images
            identities += 1
            changed.append(following)
        final_channels.append(tuple(changed))
    final_channels = tuple(final_channels)
    final_matrices = tuple(matrix for choices in final_channels for matrix in choices)
    assert len(closure(final_initial, final_matrices)[0]) == len(observable)
    assert len(closure(final_terminal, final_matrices, column=True)[0]) == len(observable)
    return {"dimension": len(observable), "reachable_dimension": len(reachable),
            "original_dimension": dimension, "channels": final_channels,
            "initial": final_initial, "terminal": final_terminal, "identities": identities}


def regularize(representation, exponents):
    """Remove only matrix poles and normalize the two boundary vectors.

    A channel with no nonzero reduced physical matrix contributes zero to
    every positive-length word containing it, and can be deleted.
    """
    if not representation["dimension"]:
        return {**representation, "exponents": (), "kept_channels": (),
                "matrix_shifts": (), "boundary_shift": 0}
    channels, changed_exponents, kept, shifts = [], [], [], []
    for symbol, (choices, exponent) in enumerate(zip(representation["channels"], exponents)):
        orders = [module.rorder(value) for matrix in choices for row in matrix for value in row
                  if value != ZERO]
        if not orders:
            continue
        shift = min(0, min(orders))
        channels.append(tuple(tuple(tuple(module.rshift(value, -shift) for value in row)
                                    for row in matrix) for matrix in choices))
        changed_exponents.append(exponent - shift)
        shifts.append(shift)
        kept.append(symbol)
    initial_order = min(module.rorder(value) for value in representation["initial"] if value != ZERO)
    terminal_order = min(module.rorder(value) for value in representation["terminal"] if value != ZERO)
    return {**representation, "channels": tuple(channels), "exponents": tuple(changed_exponents),
            "initial": tuple(module.rshift(value, -initial_order) for value in representation["initial"]),
            "terminal": tuple(module.rshift(value, -terminal_order) for value in representation["terminal"]),
            "kept_channels": tuple(kept), "matrix_shifts": tuple(shifts),
            "boundary_shift": initial_order + terminal_order}


def literal_output(channels, initial, terminal, word, physical):
    vector = tuple(initial)
    for symbol, choice in zip(word, physical):
        vector = row_product(vector, channels[symbol][choice])
    return scalar_product(vector, terminal)


def word_order(channels, initial, terminal, word, free=False):
    dimension = len(initial)
    orders = []
    for physical in product(*(range(len(channels[symbol])) for symbol in word)):
        if free:
            matrix = tuple(tuple(ONE if a == b else ZERO for b in range(dimension))
                           for a in range(dimension))
            for symbol, choice in zip(word, physical):
                matrix = matrix_product(matrix, channels[symbol][choice])
            orders.extend(module.rorder(value) for row in matrix for value in row if value != ZERO)
        else:
            value = literal_output(channels, initial, terminal, word, physical)
            if value != ZERO:
                orders.append(module.rorder(value))
    return min(orders) if orders else None


def partition(channels, initial, terminal, exponents, point, length, free=False):
    answer = Fraction(0)
    for word in product(range(len(channels)), repeat=length):
        order = word_order(channels, initial, terminal, word, free)
        if order is not None:
            answer += point ** (sum(exponents[symbol] for symbol in word) - order)
    return answer


def raw_word_basis(channels, seed, column=False):
    """A basis of literal prefix rows or suffix columns, with their words."""
    dimension = len(seed)
    selected, pending = [], [(tuple(seed), (), ())]
    cursor = 0
    while cursor < len(pending) and len(selected) < dimension:
        vector, word, physical = pending[cursor]
        cursor += 1
        old = [record[0] for record in selected]
        if len(row_basis((*old, vector), dimension)[0]) == len(old):
            continue
        selected.append((vector, word, physical))
        for symbol, choices in enumerate(channels):
            for choice, matrix in enumerate(choices):
                if column:
                    pending.append((column_product(matrix, vector), (symbol,) + word,
                                    (choice,) + physical))
                else:
                    pending.append((row_product(vector, matrix), word + (symbol,),
                                    physical + (choice,)))
    assert len(selected) == dimension
    assert max(map(lambda record: len(record[1]), selected)) <= dimension - 1
    return selected


def inverse(matrix):
    dimension = len(matrix)
    augmented = tuple(tuple(row) + tuple(ONE if i == j else ZERO for j in range(dimension))
                      for i, row in enumerate(matrix))
    basis, pivots = row_basis(augmented, 2 * dimension)
    assert pivots == tuple(range(dimension))
    answer = tuple(row[dimension:] for row in basis)
    identity = tuple(tuple(ONE if i == j else ZERO for j in range(dimension)) for i in range(dimension))
    assert matrix_product(matrix, answer) == matrix_product(answer, matrix) == identity
    return answer


def window_certificate(representation, exponents, point):
    """Exact finite-window constants in Q_free(t) <= sum_k H_k Q_b(t+k)."""
    channels, initial, terminal = (representation[key] for key in ("channels", "initial", "terminal"))
    prefixes = raw_word_basis(channels, initial)
    suffixes = raw_word_basis(channels, terminal, column=True)
    row_inverse = inverse(tuple(record[0] for record in prefixes))
    column_inverse = inverse(transpose(tuple(record[0] for record in suffixes)))
    offset = (min(module.rorder(value) for row in row_inverse for value in row if value != ZERO)
              + min(module.rorder(value) for row in column_inverse for value in row if value != ZERO))
    window = {}
    for _, prefix, _ in prefixes:
        for _, suffix, _ in suffixes:
            length = len(prefix) + len(suffix)
            weight = point ** (-offset - sum(exponents[symbol] for symbol in prefix + suffix))
            window[length] = window.get(length, Fraction(0)) + weight
    boundary_offset = (min(module.rorder(value) for value in initial if value != ZERO)
                       + min(module.rorder(value) for value in terminal if value != ZERO))
    return window, point ** (-boundary_offset)


def verify_words(channels, initial, terminal, result, maximum_length):
    comparisons = 0
    for length in range(maximum_length + 1):
        for word in product(range(len(channels)), repeat=length):
            for physical in product(*(range(len(channels[symbol])) for symbol in word)):
                before = literal_output(channels, initial, terminal, word, physical)
                after = (literal_output(result["channels"], result["initial"], result["terminal"], word, physical)
                         if result["dimension"] else ZERO)
                assert before == after
                comparisons += 1
    return comparisons


def hidden_component_control():
    tau = module.rat(2)
    channels = ((((ONE, ZERO), (ZERO, tau)),),)
    initial, terminal = (ONE, ONE), (ZERO, ONE)
    reduced = minimize(channels, initial, terminal)
    assert reduced["dimension"] == 1
    assert reduced["channels"] == ((((tau,),),),)
    rows = []
    for precision in (1, 2, 3):
        original = module.ModuleAutomaton(channels, precision, initial=initial, terminal=terminal)
        minimized = module.ModuleAutomaton(reduced["channels"], precision, free_boundaries=True)
        for length in range(1, 9):
            word = (0,) * length
            assert original.word_bound(word) == min(length, precision)
            assert minimized.word_bound(word) == length
        transfer, final = minimized.weighted_matrix((0,), Fraction(2))
        assert transfer == [{0: Fraction(1, 2)}]
        assert final == (Fraction(1),)
        rows.append({"precision": precision, "unreduced_rate": "1", "minimal_module_rate": "1/2"})
    return {"original_dimension": 2, "minimal_dimension": 1, "records": rows,
            "word_comparisons": verify_words(channels, initial, terminal, reduced, 8)}


def parity_control():
    channels = ((((ZERO, ONE), (ONE, ZERO)),),)
    initial = terminal = (ONE, ZERO)
    reduced = minimize(channels, initial, terminal)
    assert reduced["dimension"] == 2
    values = [partition(channels, initial, terminal, (0,), Fraction(2), length)
              for length in range(9)]
    assert values == [Fraction(int(length % 2 == 0)) for length in range(9)]
    assert all(partition(channels, initial, terminal, (0,), Fraction(2), length, free=True) == 1
               for length in range(9))
    return {"dimension": 2, "boundary_partition": list(map(str, values)),
            "limsup_rate": "1", "liminf_rate": "0"}


def grouped_physical_control():
    # Two identical physical choices still give one channel-word contribution.
    channels = ((((ONE,),), ((ONE,),)),)
    for length in range(8):
        assert partition(channels, (ONE,), (ONE,), (0,), Fraction(2), length) == 1
    automaton = module.ModuleAutomaton(channels, 1, free_boundaries=True)
    transfer, final = automaton.weighted_matrix((0,), Fraction(2))
    assert transfer == [{0: Fraction(1)}]
    return {"channels": 1, "physical_choices": 2, "rate": "1"}


def change_gauge(channels, initial, terminal, rng):
    dimension = len(initial)
    for _ in range(3):
        a, b = rng.sample(range(dimension), 2)
        coefficient = module.rat(rng.randrange(1, 8), rng.randrange(1, 8))
        gauge = [[ONE if i == j else ZERO for j in range(dimension)] for i in range(dimension)]
        gauge[a][b] = coefficient
        gauge = tuple(map(tuple, gauge))
        # In characteristic two this elementary transvection is its own inverse.
        channels = tuple(tuple(matrix_product(matrix_product(gauge, matrix), gauge)
                               for matrix in choices) for choices in channels)
        initial, terminal = row_product(initial, gauge), column_product(gauge, terminal)
    heights = [rng.randrange(-2, 3) for _ in range(dimension)]
    channels = tuple(tuple(tuple(tuple(module.rshift(value, heights[b] - heights[a])
                                      for b, value in enumerate(row))
                                for a, row in enumerate(matrix)) for matrix in choices)
                     for choices in channels)
    initial = tuple(module.rshift(value, height) for value, height in zip(initial, heights))
    terminal = tuple(module.rshift(value, -height) for value, height in zip(terminal, heights))
    return channels, initial, terminal


def embedded_controls(count):
    rng = Random(0xB0A0D)
    main = (
        (((module.rat(2), ONE), (ZERO, module.rat(2))),
         ((module.rat(3), ONE), (ZERO, module.rat(3)))),
        (((ZERO, ZERO), (ONE, ZERO)),
         ((ZERO, ZERO), (module.rat(2), ZERO))),
    )
    records, comparisons, identities, window_checks = [], 0, 0, 0
    for case in range(count):
        channels = []
        for choices in main:
            grown = []
            for matrix in choices:
                whole = [[ZERO] * 4 for _ in range(4)]
                for a, b in product(range(2), repeat=2):
                    whole[a][b] = matrix[a][b]
                # The third coordinate is reachable and unobservable; the
                # fourth coordinate can be observed but is unreachable.
                for a in range(2):
                    whole[a][2] = module.rat(rng.randrange(1, 8))
                whole[2][2] = module.rat(rng.randrange(1, 8))
                whole[3][2] = module.rat(rng.randrange(1, 8))
                whole[3][3] = module.rat(rng.randrange(1, 8))
                grown.append(tuple(map(tuple, whole)))
            channels.append(tuple(grown))
        initial, terminal = (ONE, ZERO, ONE, ZERO), (ONE, ZERO, ZERO, ONE)
        channels, initial, terminal = change_gauge(tuple(channels), initial, terminal, rng)
        reduced = minimize(channels, initial, terminal)
        assert reduced["reachable_dimension"] == 3 and reduced["dimension"] == 2
        comparisons += verify_words(channels, initial, terminal, reduced, 3)
        identities += reduced["identities"]
        regular = regularize(reduced, (1, -1))
        assert regular["kept_channels"] == (0, 1)
        for length in range(4):
            before = partition(channels, initial, terminal, (1, -1), Fraction(3, 2), length)
            after = partition(regular["channels"], regular["initial"], regular["terminal"],
                              regular["exponents"], Fraction(3, 2), length)
            assert after == Fraction(3, 2) ** regular["boundary_shift"] * before
        point = Fraction(3, 2)
        window, boundary_constant = window_certificate(regular, regular["exponents"], point)
        counts = {length: partition(regular["channels"], regular["initial"], regular["terminal"],
                                    regular["exponents"], point, length)
                  for length in range(3 + max(window))}
        for length in range(3):
            free = partition(regular["channels"], regular["initial"], regular["terminal"],
                             regular["exponents"], point, length, free=True)
            assert counts[length] <= boundary_constant * free
            assert free <= sum(weight * counts[length + shift] for shift, weight in window.items())
            window_checks += 2
        records.append({"case": case, "reachable_dimension": 3, "minimal_dimension": 2,
                        "matrix_shifts": list(regular["matrix_shifts"]),
                        "boundary_shift": regular["boundary_shift"],
                        "window": {str(shift): str(weight) for shift, weight in sorted(window.items())}})
    return {"cases": records, "literal_scalar_word_comparisons": comparisons,
            "exact_intertwining_identities": identities, "finite_window_inequalities": window_checks}


def crossing_controls():
    records, identities = [], 0
    for name, constructor in (("auxiliary elementary", lambda: module.elementary(2)),
                              ("ordered pairing", module.pairing),
                              ("mixed nilpotent loop", module.mixed_loop)):
        local, _ = constructor()
        for cut, side, infinity in product(range(2), range(2), (False, True)):
            factors = module.local_factors(local, cut, side, (0, 0), infinity)
            reduced = minimize(factors["channels"], (ONE, ZERO), (ONE, ZERO))
            identities += reduced["identities"]
            assert reduced["dimension"] >= 1
            if name == "auxiliary elementary":
                assert reduced["dimension"] == 1
            regular = regularize(reduced, factors["exponents"])
            records.append({"name": name, "cut": cut, "side": side, "infinity": infinity,
                            "minimal_dimension": reduced["dimension"],
                            "retained_channels": len(regular["channels"]),
                            "original_channels": len(factors["channels"])})
    return {"cases": records, "exact_intertwining_identities": identities}


def check(random_cases=8):
    zero = minimize(((((ZERO, ZERO), (ZERO, ZERO)),),), (ONE, ZERO), (ZERO, ONE))
    assert zero["dimension"] == 0
    return {"status": "pass", "field": "F_2(tau)",
            "scope": "Constructive boundary-rate completeness for a fixed grouped-channel factorization; no new factoring exponent.",
            "hidden_component": hidden_component_control(), "parity_gap": parity_control(),
            "grouped_physical_choices": grouped_physical_control(),
            "rational_gauge_embeddings": embedded_controls(random_cases),
            "crossing_reductions": crossing_controls(), "zero_series_control": True,
            "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "module_checker_sha256": sha256(Path(module.__file__).read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--random-cases", type=int, default=8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    evidence = check(args.random_cases)
    text = json.dumps(evidence, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
