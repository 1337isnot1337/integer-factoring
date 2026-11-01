#!/usr/bin/env python3
"""Exact adjugate controls and an all-degree seven-letter entropy witness."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path

import check_block_smith_crossings as smith
from check_cancellation_module_automata import pmul


ROOT = Path(__file__).resolve().parent.parent


def polynomial_power(value, exponent):
    answer = 1
    while exponent:
        if exponent & 1:
            answer = pmul(answer, value)
        value = pmul(value, value)
        exponent //= 2
    return answer


def determinant(matrix):
    """Independent small determinant by permutations in characteristic two."""
    size = len(matrix)
    if not size:
        return 1
    answer = 0
    for selected in permutations(range(size)):
        value = 1
        for row, column in enumerate(selected):
            value = pmul(value, matrix[row][column])
        answer ^= value
    return answer


def adjugate(matrix):
    size = len(matrix)
    return [[determinant([[matrix[i][j] for j in range(size) if j != row]
                          for i in range(size) if i != column])
             for column in range(size)] for row in range(size)]


def matrix_product(left, right):
    output = []
    for row in left:
        output_row = []
        for column in zip(*right):
            value = 0
            for a, b in zip(row, column):
                value ^= pmul(a, b)
            output_row.append(value)
        output.append(output_row)
    return output


def reshuffle(operator, q):
    return [[operator[q*x+u][q*v+y] for u, v in product(range(q), repeat=2)]
            for x, y in product(range(q), repeat=2)]


def transpose_channel_controls():
    records = []
    for q in (2, 3):
        left = [[int(i == j) for j in range(q)] for i in range(q)]
        right = [[int(i == j) for j in range(q)] for i in range(q)]
        left[0][0] = 2  # singular at zero, but generically invertible
        for i in range(q-1):
            left[i][i+1] = 3
            right[i+1][i] = 2
        adj_left, adj_right = adjugate(left), adjugate(right)
        det_pair = pmul(determinant(left), determinant(right))
        factor = polynomial_power(det_pair, q-1)
        operator = [[pmul(left[i][l], right[k][j])
                     for k, l in product(range(q), repeat=2)]
                    for i, j in product(range(q), repeat=2)]
        claimed = [[pmul(factor, pmul(adj_right[l][i], adj_left[j][k]))
                    for k, l in product(range(q), repeat=2)]
                   for i, j in product(range(q), repeat=2)]
        det_operator = polynomial_power(det_pair, q)
        expected = [[det_operator*int(i == j) for j in range(q*q)] for i in range(q*q)]
        assert matrix_product(operator, claimed) == expected
        assert matrix_product(claimed, operator) == expected
        if q == 2:
            assert claimed == adjugate(operator)
        flat = reshuffle(claimed, q)
        minors = 0
        for i in range(q*q):
            for j in range(i+1, q*q):
                for k in range(q*q):
                    for l in range(k+1, q*q):
                        assert pmul(flat[i][k], flat[j][l]) == pmul(flat[i][l], flat[j][k])
                        minors += 1
        records.append({"alphabet": q, "operator_dimension": q*q,
                        "singular_leading_channel_factors": True,
                        "operator_adjugate_product_entries": 2*q**4,
                        "reshuffled_two_minors_zero": minors})
    return records


def elementary_fixture(function, vacuum_defect=0):
    local = {}
    for x, y, u, v in product(range(2), repeat=4):
        straight, crossing = int(x == u and y == v), int(x == v and y == u)
        value = straight ^ (function if straight ^ crossing else 0)
        if x == y == u == v == 0:
            value ^= vacuum_defect
        assert value.bit_count() % 2 == crossing
        local[x, y, u, v] = value
    return local


def crossing_matrices(local, q):
    row = [[local[x, y, u, v] for u, v in product(range(q), repeat=2)]
           for x, y in product(range(q), repeat=2)]
    column = [[local[x, y, u, v] for y, v in product(range(q), repeat=2)]
              for x, u in product(range(q), repeat=2)]
    return row, column


def orientation_control(local, name):
    q, n = 2, 4
    degree = max(value.bit_length()-1 for value in local.values())
    row, column = crossing_matrices(local, q)
    reverse_row = [[smith.reverse_polynomial(value, degree) for value in line] for line in row]
    mu, _ = smith.smith_profile(reverse_row, degree)
    nu, _ = smith.smith_profile(column, degree)
    operator = [[smith.reverse_polynomial(column[q*x+u][q*y+v], degree)
                 for v, y in product(range(q), repeat=2)]
                for x, u in product(range(q), repeat=2)]
    assert all(value.bit_count() % 2 == int(i == j)
               for i, line in enumerate(operator) for j, value in enumerate(line))
    adj = adjugate(operator)
    det = determinant(operator)
    expected = [[det*int(i == j) for j in range(n)] for i in range(n)]
    assert matrix_product(operator, adj) == matrix_product(adj, operator) == expected

    adj_degree = max(value.bit_length()-1 for line in adj for value in line)
    assert adj_degree == (n-1)*degree-sum(nu[:-1])
    adj_row = reshuffle(adj, q)
    eta, _ = smith.smith_profile(adj_row, adj_degree)
    reversed_adj = [[smith.reverse_polynomial(value, adj_degree) for value in line] for line in adj]
    theta, _ = smith.smith_profile(reversed_adj, adj_degree)
    assert theta == tuple(nu[-1]-value for value in reversed(nu))
    assert all(value >= mu[1] for value in eta[1:])
    capped = tuple(min(value, 2*mu[1]) for value in mu)
    assert all(actual >= lower for actual, lower in zip(eta, capped))
    assert (n-1)*mu[1] <= n*adj_degree
    assert sum(capped) <= n*adj_degree
    for selected in range(2, q+1):
        assert sum(eta[:selected])+sum(theta[:selected]) <= selected*adj_degree
        assert (selected-1)*mu[1]+selected*sum(nu)-sum(nu[-selected:]) \
            <= selected*(n-1)*degree
        assert sum(capped[:selected])+selected*sum(nu)-sum(nu[-selected:]) \
            <= selected*(n-1)*degree
    leading_rank = smith.binary_rank([
        sum((value & 1) << j for j, value in enumerate(line)) for line in operator])
    return {
        "name": name, "degree": degree, "row_infinity": mu, "column_zero": nu,
        "adjugate_degree": adj_degree,
        "adjugate_row_zero": eta, "adjugate_column_infinity": theta,
        "universal_two_jet_lower_profile": capped,
        "leading_column_operator_rank": leading_rank,
        "determinant_is_monomial": det.bit_count() == 1,
        "positive_second_row_order": mu[1] > 0,
        "literal_adjugate_product_entries": 2*n*n,
    }


def profile_controls():
    # A strict cap control: the original last order is seven, whereas
    # the adjugate last order is six.  A theorem claiming eta_i>=mu_i
    # without truncation would therefore be false.
    cap_values = (8, 0, 0, 12, 3, 3, 2, 0, 3, 2, 3, 5, 9, 0, 0, 13)
    cap_fixture = dict(zip(product(range(2), repeat=4), cap_values))
    fixtures = (
        ("elementary", elementary_fixture(2)),
        ("ramified_cubic", elementary_fixture(8)),
        ("nonmonomial_cubic", elementary_fixture(14)),
        ("singular_leading_vacuum_defect", elementary_fixture(2, 10)),
        ("strict_two_jet_cap", cap_fixture),
    )
    records = []
    for name, local in fixtures:
        records.append(orientation_control(local, name))
        degree = max(value.bit_length()-1 for value in local.values())
        exchanged = {(x, y, u, v): smith.reverse_polynomial(local[x, u, y, v], degree)
                     for x, y, u, v in local}
        reverse_record = orientation_control(exchanged, name+"_exchanged")
        assert reverse_record["row_infinity"] == records[-1]["column_zero"]
        assert reverse_record["column_zero"] == records[-1]["row_infinity"]
        records.append(reverse_record)
    assert any(record["leading_column_operator_rank"] < 4
               and not record["determinant_is_monomial"]
               and record["positive_second_row_order"] for record in records)
    cap_record = next(record for record in records if record["name"] == "strict_two_jet_cap")
    assert cap_record["row_infinity"] == (0, 1, 2, 7)
    assert cap_record["adjugate_row_zero"] == (0, 1, 2, 6)
    return records


def entropy_certificate():
    q, n = 7, 49
    prefix = [Fraction(0)]+[Fraction(int(i < q)) for i in range(1, n)]
    adj_det = [Fraction(0)]+[Fraction(n)]*(n-2)+[Fraction(0)]
    adj_det[1] += n-1
    adj_shared = [Fraction(0)]+[Fraction(q-int(i >= n-q)) for i in range(1, n)]
    adj_shared[1] += q-1
    multipliers = (Fraction(34, 875), Fraction(1, 28000), Fraction(1, 1200))
    right_sides = (Fraction(7, 2), Fraction(2352), Fraction(336))
    assert sum(a*b for a, b in zip(multipliers, right_sides)) == Fraction(1, 2)
    probabilities = [sum(multiplier*row[i] for multiplier, row in zip(
        multipliers, (prefix, adj_det, adj_shared))) for i in range(n)]
    probabilities[0] = 1-sum(probabilities[1:])
    counts = [33913, 4465]+[3901]*5+[637]*35+[567]*6+[420]
    denominator = 84000
    assert probabilities == [Fraction(value, denominator) for value in counts]
    assert sum(counts) == denominator and all(value > 0 for value in counts)
    assert Fraction(119112, denominator) == Fraction(709, 500)
    left = denominator**denominator
    right = 7**119112
    for value, multiplicity in Counter(counts).items():
        right *= value**(value*multiplicity)
    assert left > right
    assert 2*709**2 > 1000**2
    digest = lambda value: sha256(value.to_bytes((value.bit_length()+7)//8, "big")).hexdigest()
    return {
        "multipliers": [str(value) for value in multipliers],
        "weighted_rhs": "1/2",
        "probability_denominator": denominator,
        "probability_counts": sorted(Counter(counts).items(), reverse=True),
        "exact_comparison_seven_exponent": 119112,
        "left_integer_bits": left.bit_length(), "right_integer_bits": right.bit_length(),
        "left_sha256": digest(left), "right_sha256": digest(right),
        "strict_integer_entropy_inequality": True,
        "rho_lower_bound": "709/1000",
        "strictly_above_one_over_sqrt_two": 2*709**2-1000**2,
    }


def eight_letter_relaxation_witness():
    """A feasible formal profile pair; this is not a tensor construction."""
    q, n, degree = 8, 64, 140
    mu = nu = (0,)+(80,)*7+(147,)*55+(224,)
    assert len(mu) == n and mu[0] == 0 and tuple(sorted(mu)) == mu
    assert sum(mu) <= n*degree
    capped = tuple(min(value, 2*mu[1]) for value in mu)
    old_det = (n-1)*mu[1]+n*sum(nu[:-1])
    new_det = sum(capped)+n*sum(nu[:-1])
    assert old_det <= n*(n-1)*degree
    assert new_det <= n*(n-1)*degree
    prefixes = []
    for selected in range(2, q+1):
        original = sum(mu[:selected])+sum(nu[:selected])
        old_adj = (selected-1)*mu[1]+selected*sum(nu)-sum(nu[-selected:])
        new_adj = sum(capped[:selected])+selected*sum(nu)-sum(nu[-selected:])
        assert original <= selected*degree
        assert old_adj <= selected*(n-1)*degree
        assert new_adj <= selected*(n-1)*degree
        prefixes.append({
            "size": selected, "original_lhs": original,
            "adjugate_lhs": old_adj, "two_jet_lhs": new_adj,
            "adjugate_rhs": selected*(n-1)*degree,
        })
    polynomial = 36**224+7*35**80*36**144+55*35**147*36**77+35**224
    left = polynomial**500
    right = 8**707*36**77000*35**35000
    assert left < right
    assert 2*707**2 < 1000**2
    digest = lambda value: sha256(value.to_bytes((value.bit_length()+7)//8, "big")).hexdigest()
    return {
        "alphabet": q, "degree": degree,
        "both_profiles": sorted(Counter(mu).items()),
        "profile_total": sum(mu), "profile_total_budget": n*degree,
        "capped_profile_total": sum(capped),
        "adjugate_determinant_lhs": old_det,
        "two_jet_determinant_lhs": new_det,
        "adjugate_determinant_rhs": n*(n-1)*degree,
        "prefix_checks": prefixes,
        "scale": "36/35", "strict_rho_upper_bound": "707/1000",
        "left_integer_bits": left.bit_length(), "right_integer_bits": right.bit_length(),
        "left_sha256": digest(left), "right_sha256": digest(right),
        "strict_integer_pressure_inequality": True,
        "strictly_below_one_over_sqrt_two": 1000**2-2*707**2,
        "crossing_realization_claimed": False,
    }


def run():
    return {
        "scope": "all-degree bare seven-letter Smith pressure, all fields and fixed auxiliary dimension",
        "scalar_profile_constraints_extend_by_constant_contraction_only_at_pressure_level": True,
        "transpose_channel_controls": transpose_channel_controls(),
        "small_polynomial_adjugate_controls": profile_controls(),
        "entropy_certificate": entropy_certificate(),
        "eight_letter_formal_relaxation_witness": eight_letter_relaxation_witness(),
        "floating_arithmetic_used": False,
        "factoring_exponent_improved": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "evidence" / "adjugate-seven-letter-floor.json")
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
