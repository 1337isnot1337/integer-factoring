#!/usr/bin/env python3
"""Independent full-state controls and exact multilevel gain certificates."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import random

from charged_gaussian_floor import Dyadic
from symmetric_binary_filter_synthesis import synthesize, transform


F = Fraction


def xor(values):
    result = 0
    for value in values:
        result ^= value
    return result


def product_polynomial(left, right):
    # This independent tensor control uses prime-field coefficients.
    result = {}
    for a, c in left.items():
        for b, d in right.items():
            result[a + b] = result.get(a + b, 0) ^ (c * d)
    return {degree: coefficient for degree, coefficient in result.items() if coefficient}


def full_transform(values, modes, m):
    """Expand each tensor coordinate, without symmetry or Lucas compression."""
    result = []
    for state in range(1 << modes):
        value = {}
        for source in range(1 << modes):
            if source & state != source:
                continue
            power = state.bit_count() - source.bit_count()
            for shift in range(power + 1):
                if comb(power, shift) % 2:
                    for degree, coefficient in values[source].items():
                        at = degree + m * shift
                        value[at] = value.get(at, 0) ^ coefficient
        result.append({degree: coefficient for degree, coefficient in value.items() if coefficient})
    return tuple(result)


def uncompressed_feasible(modes, m, degrees, orders):
    variables = [(weight, degree) for weight, bound in enumerate(degrees)
                 for degree in range(bound + 1)]
    equations = {}
    endpoint_masks = [0] * (modes + 1)
    for column, (weight, degree) in enumerate(variables):
        bit = 1 << column
        endpoint_masks[weight] |= bit
        for state in range(1 << modes):
            for source in range(1 << modes):
                if source.bit_count() != weight or source & state != source:
                    continue
                power = state.bit_count() - weight
                for shift in range(power + 1):
                    at = degree + m * shift
                    if comb(power, shift) % 2 and at < orders[state.bit_count()]:
                        equations[state, at] = equations.get((state, at), 0) ^ bit
    pivots = {}
    for value in equations.values():
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    # A coordinate is forced to vanish precisely when its endpoint linear
    # functional belongs to the coefficient-equation row space.
    for mask in endpoint_masks:
        for pivot in sorted(pivots, reverse=True):
            if (mask >> pivot) & 1:
                mask ^= pivots[pivot]
        if not mask:
            return False
    return True


def exact_ratio(n, energies, orders, degrees, base, charge):
    zero = infinity = norm_zero = norm_infinity = F(0)
    for state in range(1 << len(energies)):
        weight = state.bit_count()
        energy = sum(value for bit, value in enumerate(energies) if (state >> bit) & 1)
        wz = base ** (energy + charge * weight)
        wi = base ** energy
        norm_zero += wz
        norm_infinity += wi
        zero += wz * base ** (2 * (orders[state] - n * weight))
        infinity += wi * base ** (-2 * degrees[state])
    return zero * infinity / (norm_zero * norm_infinity)


def tensor_controls():
    n = m = 4
    left = ({0: 1, 4: 1, 5: 1}, {0: 1}, {0: 1})
    right = ({0: 1}, {0: 1, 8: 1, 11: 1}, {0: 1, 4: 1, 6: 1})
    full_left = tuple(left[state.bit_count()] for state in range(4))
    full_right = tuple(right[state.bit_count()] for state in range(4))
    beta_left = full_transform(full_left, 2, m)
    beta_right = full_transform(full_right, 2, m)
    full = tuple(product_polynomial(full_left[state & 3], full_right[state >> 2])
                 for state in range(16))
    beta = full_transform(full, 4, m)
    for state in range(16):
        assert beta[state] == product_polynomial(beta_left[state & 3], beta_right[state >> 2])
        assert max(full[state]) == max(full_left[state & 3]) + max(full_right[state >> 2])
        assert min(beta[state]) == min(beta_left[state & 3]) + min(beta_right[state >> 2])
        assert xor(full[state].values()) == 1
    ratios = 0
    for base in (F(49, 50), F(4, 5), F(2, 3)):
        for charge in (-12, -8, 0, 8, 12):
            first = exact_ratio(n, (7, 1), tuple(map(min, beta_left)),
                                tuple(map(max, full_left)), base, charge)
            second = exact_ratio(n, (5, 3), tuple(map(min, beta_right)),
                                 tuple(map(max, full_right)), base, charge)
            combined = exact_ratio(n, (7, 1, 5, 3), tuple(map(min, beta)),
                                   tuple(map(max, full)), base, charge)
            assert combined == first * second
            ratios += 1
    return {'different_block_polynomial_states': 16,
            'exact_two_block_cost_factorizations': ratios,
            'block_mode_energies': [[7, 1], [5, 3]],
            'different_block_polynomials': [left, right]}


def gaussian_at_rho(arithmetic, n, scale, charge, rho):
    result = arithmetic.zero
    for i in range(n):
        odd = arithmetic.rational(2 * i + 1)
        for energy in (arithmetic.subtract(odd, charge), arithmetic.add(odd, charge)):
            result = arithmetic.add(result, arithmetic.softplus(arithmetic.multiply(energy, scale)))
    for i in range(2 * n):
        result = arithmetic.add(result, arithmetic.softplus(arithmetic.multiply(
            arithmetic.rational(2 * i + 1 - 2 * n), scale)))
    ratio = arithmetic.rational(2 * rho.numerator, rho.denominator)
    physical_scale = arithmetic.divide(scale, ratio)
    physical = arithmetic.zero
    for i in range(n):
        center = arithmetic.rational(4 * i - 2 * n + 2)
        for energy in (arithmetic.subtract(center, charge), arithmetic.add(center, charge)):
            physical = arithmetic.add(physical, arithmetic.softplus(
                arithmetic.multiply(energy, physical_scale)))
    return arithmetic.subtract(result, arithmetic.multiply(ratio, physical))


def multilevel_relative_gain_certificate():
    n = m = 4
    values = ({0: 1, 8: 1, 16: 1}, {0: 1}, {0: 1}, {0: 1}, {0: 1})
    beta = transform(values, m)
    expected_orders, expected_degrees = (0, 4, 24, 8, 8), (16, 0, 0, 0, 0)
    assert tuple(map(min, beta)) == expected_orders
    assert tuple(map(max, values)) == expected_degrees
    full = tuple(values[state.bit_count()] for state in range(16))
    independently_transformed = full_transform(full, 4, m)
    assert independently_transformed == tuple(beta[state.bit_count()] for state in range(16))
    assert all(xor(polynomial.values()) == 1 for polynomial in values + beta)

    arithmetic = Dyadic(80)
    # exp(-t)=49/50, so every occupation expectation is rational.
    scale = arithmetic.log_one_plus(arithmetic.rational(1, 49))
    cases = []
    for charge in (0, 8):
        ratio = exact_ratio(n, (7, 5, 3, 1),
                            tuple(expected_orders[state.bit_count()] for state in range(16)),
                            tuple(expected_degrees[state.bit_count()] for state in range(16)),
                            F(49, 50), charge)
        assert 0 < ratio < 1
        excess = 1 / ratio - 1
        log_inverse = arithmetic.log_one_plus(arithmetic.rational(excess.numerator, excess.denominator))
        log_ratio = -log_inverse[1], -log_inverse[0]
        full_pressure = arithmetic.add(arithmetic.finite_pressure(
            n, scale, arithmetic.rational(charge)), log_ratio)
        assert full_pressure[0] > 0
        root_low, root_high = F(966, 1000), F(967, 1000)
        low_pressure = arithmetic.add(gaussian_at_rho(
            arithmetic, n, scale, arithmetic.rational(charge), root_low), log_ratio)
        high_pressure = arithmetic.add(gaussian_at_rho(
            arithmetic, n, scale, arithmetic.rational(charge), root_high), log_ratio)
        assert low_pressure[0] > 0 > high_pressure[1]
        cases.append({
            'charge': charge, 'exp_minus_t': '49/50',
            'exact_relative_cost': str(ratio), 'relative_cost_is_below_one': True,
            'log_relative_cost_interval_numerators': [str(x) for x in log_ratio],
            'full_critical_pressure_interval_numerators': [str(x) for x in full_pressure],
            'genuine_rho_root_bracket': [str(root_low), str(root_high)],
            'lower_rho_pressure_interval_numerators': [str(x) for x in low_pressure],
            'upper_rho_pressure_interval_numerators': [str(x) for x in high_pressure],
        })
    return {
        'n': n, 'ramification': 1, 'modes': 4, 'f_by_weight': values,
        'beta_by_weight': beta, 'full_degrees': expected_degrees,
        'exact_beta_orders': expected_orders, 'all_endpoint_values': 1,
        'common_denominator': str(arithmetic.scale), 'cases': cases,
        'scope': 'Exact relative gains at natural and interior charges; complete pressure stays positive and rho remains above the factoring threshold.',
    }


def main():
    directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=directory.parent / 'evidence' /
                        'symmetric-binary-filter-synthesis.json')
    args = parser.parse_args()
    generator = random.Random(720491)
    comparisons = witnesses = full_states = 0
    for modes in (1, 2, 3, 4):
        for m in (1, 2, 3):
            for _ in range(12):
                degrees = tuple(generator.randrange(6) for _ in range(modes + 1))
                orders = tuple(generator.randrange(2 * m + 6) for _ in range(modes + 1))
                expected = uncompressed_feasible(modes, m, degrees, orders)
                witness = synthesize(modes, m, degrees, orders)
                assert bool(witness) == expected
                comparisons += 1
                if witness:
                    full = tuple(witness['f'][state.bit_count()] for state in range(1 << modes))
                    beta = full_transform(full, modes, m)
                    assert beta == tuple(witness['beta'][state.bit_count()] for state in range(1 << modes))
                    assert all(xor(polynomial.values()) for polynomial in full + beta)
                    witnesses += 1
                    full_states += 1 << modes
    for modes in (1, 2, 4, 8, 16, 32, 64):
        witness = synthesize(modes, modes, (0,) * (modes + 1),
                             tuple(modes * weight for weight in range(modes + 1)))
        assert witness and witness['actual_degrees'] == (0,) * (modes + 1)
        assert witness['actual_orders'] == tuple(modes * weight for weight in range(modes + 1))
    sources = ('symmetric_binary_filter_synthesis.py', Path(__file__).name,
               'charged_gaussian_floor.py', 'polynomial_state_filter_synthesis.py',
               'weyl_lattice_filter_synthesis.py')
    record = {
        'status': 'Exact full-state feasibility, polynomial tensorization, and relative-gain certificates passed.',
        'scope': 'Constructive symmetric multilevel filter interface; no improved factoring exponent.',
        'uncompressed_feasibility_comparisons': comparisons,
        'feasible_literal_witnesses': witnesses, 'full_state_transform_controls': full_states,
        'baseline_lengths_through_64': [1, 2, 4, 8, 16, 32, 64],
        'tensor_product_controls': tensor_controls(),
        'natural_and_interior_relative_gain': multilevel_relative_gain_certificate(),
        'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in sources},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
