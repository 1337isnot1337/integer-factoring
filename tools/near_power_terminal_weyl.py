#!/usr/bin/env python3
"""Exact dummy-mode Weyl reduction for n=2^k-1 terminal operator cores.

One independent YU Gaussian completes the rectangular core. The square
Weyl transform is computed with complete polynomials, then the known
dummy Gaussian factor is divided out of the charge-tagged Smith profile.
The dummy is local profile arithmetic, not a new physical crossing mode.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache

from compiled_f4_identity_operators import _compile_matrices, _orders
from coupled_identity_mode_filters import add as add2, multiply as multiply2, shift, algebra_multiply
from coupled_identity_operator_filters import states
from coupled_identity_weyl_profiles import determinant, adjugate
from f4_identity_operator_filters import normalized_kernel
from f4_repaired_identity_filters import MUL, add, multiply, evaluate_one, laurent_orders
from near_power_terminal_operators import geometry_data


@lru_cache(maxsize=16)
def weyl_data(n, modes, initial_modes=0):
    data = geometry_data(n, modes, initial_modes)
    selected = data['augmented']
    fx, fv = data['augmented_x'], data['augmented_v']
    size = modes+1
    dx, dv = determinant(fx, 2), determinant(fv, 2)
    delta = multiply2({0: 1, 1: 1}, multiply2(dx, dv, 2), 2)
    assert delta.get(0) == 1
    ax, av = adjugate(fx, 2), adjugate(fv, 2)
    numerator = []
    for i in range(size):
        row = []
        for j in range(size):
            value = dict(delta) if i == j else {}
            for k in range(size):
                term = multiply2(multiply2(ax[i][k], av[j][k], 2), data['unit'], 2)
                value = add2(value, shift(term, selected[i]+data['retained'][k]+1, 2), 2)
            assert value.get(0, 0) == int(i == j)
            if value:
                assert min(value) >= max(0, selected[i]-selected[j])
            row.append(value)
        numerator.append(tuple(row))
    return {'selected': selected, 'retained': data['retained'],
            'denominator_unit': delta, 'inverse_parameter_numerator': tuple(numerator)}


@lru_cache(maxsize=16)
def reduced_components(n, modes, initial_modes=0):
    """Complete polynomial sources Delta^(r+1) Psi(F tensor I)."""
    data = weyl_data(n, modes, initial_modes)
    size = modes+1
    alphabet = states(2, size)
    lookup = {state: i for i, state in enumerate(alphabet)}
    delta, numerator = data['denominator_unit'], data['inverse_parameter_numerator']
    delta_powers = [{0: 1}]
    for exponent in range(size):
        delta_powers.append(multiply2(delta_powers[-1], delta, 2))
    forms = []
    for j in range(size):
        form = {}
        for i in range(size):
            powers = [0]*size
            powers[i] = 1
            if numerator[j][i]:
                form[tuple(powers)] = numerator[j][i]
        forms.append(form)
    cache = {(0,)*size: {(0,)*size: {0: 1}}}

    def monomial(powers):
        if powers not in cache:
            j = next(i for i, digit in enumerate(powers) if digit)
            previous = list(powers)
            previous[j] = 0
            cache[powers] = algebra_multiply(forms[j], monomial(tuple(previous)), 2)
        return cache[powers]

    result = {}
    for b in states(2, modes):
        for a in states(2, modes):
            if sum(a) != sum(b):
                continue
            matrix = {}
            for ghost in (0, 1):
                aa, bb = a+(ghost,), b+(ghost,)
                exponents = tuple(1-digit for digit in bb)
                polynomial = {powers: multiply2(value, delta_powers[sum(bb)], 2)
                              for powers, value in monomial(exponents).items()}
                derivatives = tuple(1-digit for digit in aa)
                for column, source in enumerate(alphabet):
                    for powers, value in polynomial.items():
                        total = tuple(x+y for x, y in zip(source, powers))
                        if any(digit >= 2 for digit in total):
                            continue
                        if any(digit < derivative for digit, derivative in zip(total, derivatives)):
                            continue
                        target = tuple(digit-derivative for digit, derivative in zip(total, derivatives))
                        entry = (lookup[target], column)
                        matrix[entry] = add2(matrix.get(entry, {}), value, 2)
            result[a, b] = {entry: value for entry, value in matrix.items() if value}
    return result


def remove_dummy(augmented, dummy_order, modes):
    """Divide exact charge/order multiplicities by 1+u*z^dummy_order."""
    previous = Counter()
    result = []
    for charge in range(-modes-1, modes+2):
        current = Counter(augmented.get(charge, ()))
        for order, count in previous.items():
            current[order+dummy_order] -= count
        if any(count < 0 for count in current.values()):
            raise ArithmeticError('The augmented profile does not contain its promised dummy Gaussian.')
        current = +current
        if charge <= modes:
            result.append({'charge': charge,
                           'orders': tuple(sorted(current.elements()))})
        elif current:
            raise ArithmeticError('A residual profile remained after exact dummy division.')
        previous = current
    assert sum(len(block['orders']) for block in result) == 2**(2*modes+1)
    return tuple(result)


def _source_matrices(n, ramification, modes, kernel, initial_modes=0):
    alphabet = states(2, modes+1)
    selected = weyl_data(n, modes, initial_modes)['selected']
    components = reduced_components(n, modes, initial_modes)
    matrix = [[{} for column in alphabet] for row in alphabet]
    for index, scalar in kernel.items():
        if index not in components:
            raise ValueError('The kernel must preserve total binary occupation.')
        for (row, column), polynomial in components[index].items():
            weight = ramification*sum(j*digit for j, digit in zip(selected, alphabet[row]))
            value = {ramification*degree+weight: coefficient for degree, coefficient in polynomial.items()}
            matrix[row][column] = add(matrix[row][column], multiply(value, scalar))
    return matrix


def zero_profiles(n, ramification, modes, kernel, initial_modes=0):
    selected = geometry_data(n, modes, initial_modes)['selected']
    kernel, inverse = normalized_kernel(modes, kernel, selected)
    alphabet = states(2, modes+1)
    augmented_selected = selected+(n,)
    matrix = _source_matrices(n, ramification, modes, kernel, initial_modes)
    result = defaultdict(list)
    for occupation in range(modes+2):
        indices = [i for i, state in enumerate(alphabet) if sum(state) == occupation]
        block = [[matrix[i][j] for j in indices] for i in indices]
        orders = tuple(2*value for value in laurent_orders(block)[0])
        for auxiliary in alphabet:
            shift_order = -2*ramification*sum(j*digit for j, digit in zip(augmented_selected, auxiliary))
            result[occupation-sum(auxiliary)].extend(order+shift_order for order in orders)
    return remove_dummy(result, 2*ramification*n, modes)


def compile_family(n, ramification, modes, constant, perturbations, initial_modes=0):
    """Compile complete F4 sources, including exact endpoint validation."""
    selected = geometry_data(n, modes, initial_modes)['selected']
    augmented_selected = selected+(n,)
    alphabet = states(2, modes)
    augmented_alphabet = states(2, modes+1)
    zero_sources, infinity_sources, endpoint_sources = [], [], []
    for kernel in (constant,)+tuple(perturbations):
        zero_sources.append(_source_matrices(n, ramification, modes, kernel, initial_modes))
        infinity = []
        for row in alphabet:
            entries = []
            for column in alphabet:
                weight = ramification*sum((2*j+1-2*n)*digit for j, digit in zip(selected, column))
                entries.append({weight-2*degree: value for degree, value in kernel.get((row, column), {}).items()})
            infinity.append(entries)
        infinity_sources.append(infinity)
        endpoint_sources.append(tuple(evaluate_one(kernel.get((row, column), {}))
                                      for row in alphabet for column in alphabet))
    zero_compiled, infinity_compiled = [], []
    for occupation in range(modes+2):
        indices = [i for i, state in enumerate(augmented_alphabet) if sum(state) == occupation]
        zero_compiled.append((occupation, _compile_matrices(
            [[[matrix[i][j] for j in indices] for i in indices] for matrix in zero_sources])))
    for occupation in range(modes+1):
        indices = [i for i, state in enumerate(alphabet) if sum(state) == occupation]
        infinity_compiled.append((occupation, _compile_matrices(
            [[[matrix[i][j] for j in indices] for i in indices] for matrix in infinity_sources])))
    auxiliary = [(sum(state), -2*ramification*sum(j*digit for j, digit in zip(augmented_selected, state)))
                 for state in augmented_alphabet]

    @lru_cache(maxsize=4096)
    def validate_endpoint(values):
        endpoint = {(row, column): {0: values[i*len(alphabet)+j]}
                    for i, row in enumerate(alphabet) for j, column in enumerate(alphabet)
                    if values[i*len(alphabet)+j]}
        normalized_kernel(modes, endpoint, selected)

    def evaluate(colors):
        if len(colors) != len(perturbations) or any(color not in (0, 1, 2, 3) for color in colors):
            raise ValueError('Use one F4 color per affine perturbation.')
        endpoint = list(endpoint_sources[0])
        for source, color in zip(endpoint_sources[1:], colors):
            if color:
                for i, value in enumerate(source):
                    endpoint[i] ^= MUL[color][value]
        validate_endpoint(tuple(endpoint))
        zero = defaultdict(list)
        for occupation, compiled in zero_compiled:
            orders = tuple(2*value for value in _orders(compiled, colors))
            for auxiliary_occupation, shift_order in auxiliary:
                zero[occupation-auxiliary_occupation].extend(order+shift_order for order in orders)
        infinity = tuple(sorted(value for occupation, compiled in infinity_compiled
                                for value in _orders(compiled, colors)))
        return remove_dummy(zero, 2*ramification*n, modes), infinity

    return evaluate
