"""Search consecutive prime-form powers for base-smooth leading norms.

Candidates are tested by a modular primorial power and then decoded by
division over a fixed known-prime base. No
integer-factorization oracle, coordinate search, or new prime is used.
The window length is a finite budget, not a proved uniform success bound.
"""

from math import gcd, isqrt, prod

from arithmetic_two_group import power
from binary_quadratic_forms import ambiguous_factor, inverse_form, reduce_form
from dynamic_leading_norm_relations import leading_norm_word
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top
from weighted_norm_relations import least_power_above


def compose_with_prime_form(form, prime, root, reduced_prime_form):
    """Coprime CRT composition, with exact HNF fallback at a shared prime."""
    leading, middle, trailing = form
    if leading % prime == 0:
        return compose_forms_hnf(form, reduced_prime_form)
    delta = middle * middle - 4 * leading * trailing
    correction = ((root - middle) // 2) * pow(leading, -1, prime) % prime
    new_middle = middle + 2 * leading * correction
    new_leading = leading * prime
    return reduce_form((new_leading, new_middle,
                        (new_middle * new_middle - delta) // (4 * new_leading)))


def factor_with_smooth_leading_powers(integer, cutoff, window, *, use_primorial=True):
    if not isinstance(integer, int) or integer < 3 or integer % 2 == 0:
        raise ValueError('An odd integer of at least three is required.')
    if (not isinstance(cutoff, int) or cutoff < 2
            or not isinstance(window, int) or window < 1):
        raise ValueError('The cutoff must be at least two and the window positive.')
    delta = -integer if integer % 4 == 3 else -4 * integer
    parity = delta % 2
    identity = (1, parity, (parity - delta) // 4)
    height = isqrt((-delta) // 3)
    base = build_prime_form_factor_base(delta, cutoff)
    records = {r.prime: (r.root, r.positive) for r in base['records']}
    if integer % 4 == 1:
        records[2] = (2, reduce_form((2, 2, (integer + 1) // 2)))
    primes = sorted(records)
    prime_product = prod(primes)
    smoothness_exponent = height.bit_length()
    indices = {p: i for i, p in enumerate(primes)}
    rows, witnesses, attempts, parity_pivots = [], [], [], {}
    queries = trial_divisions = maximum_leading_bits = maximum_exponent = 0
    modular_smoothness_tests = 0
    last_extraction = None
    last_extraction_rows = -1

    def extract():
        nonlocal last_extraction, last_extraction_rows
        _, factor, last_extraction = first_useful_relation_top(
            [records[p][1] for p in primes], rows, identity,
            compose_forms_hnf, inverse_form,
            lambda form: ambiguous_factor(form, integer), verify_relations=True)
        last_extraction_rows = len(rows)
        return factor

    def finish(factor, method, terminal=None):
        return {
            'integer': integer, 'input_bits': integer.bit_length(),
            'discriminant': delta, 'height_bound': height,
            'cutoff': cutoff, 'window': window, 'factor': factor, 'method': method,
            'base_size': len(primes), 'relation_rank': len(rows),
            'norm_queries': queries, 'trial_divisions': trial_divisions,
            'modular_smoothness_tests': modular_smoothness_tests,
            'smoothness_method': 'primorial power' if use_primorial else 'division',
            'prime_product_bits': prime_product.bit_length(),
            'maximum_leading_bits': maximum_leading_bits,
            'maximum_exponent_tested': maximum_exponent,
            'factoring_oracle_calls': 0, 'new_primes_introduced': 0,
            'prime_records': records, 'active_primes': primes,
            'attempts': attempts, 'witnesses': witnesses, 'relation_matrix': rows,
            'last_extraction': last_extraction, 'terminal_witness': terminal,
        }

    for prime in base['excluded_primes']:
        divisor = gcd(prime, integer)
        if 1 < divisor < integer:
            return finish(divisor, 'factor-base gcd', {'prime': prime})

    for index, prime in enumerate(primes):
        first = least_power_above(prime, height)
        target = power(records[prime][1], first, identity, compose_forms_hnf)
        selected = None
        for offset in range(window):
            exponent = first + offset
            queries += 1
            maximum_exponent = max(maximum_exponent, exponent)
            leading = target[0]
            maximum_leading_bits = max(maximum_leading_bits, leading.bit_length())
            divisor = gcd(leading, integer)
            if 1 < divisor < integer:
                attempts.append({'prime': prime, 'queries': offset + 1,
                                 'selected_exponent': None, 'terminal_gcd': True})
                return finish(divisor, 'leading-norm gcd', {
                    'prime': prime, 'exponent': exponent, 'target': target})
            if divisor != 1 or not 1 <= leading <= height:
                raise ValueError('The reduced leading-norm invariant failed.')
            remaining, parts = leading, []
            if use_primorial:
                modular_smoothness_tests += 1
            smooth_candidate = (not use_primorial or
                                pow(prime_product, smoothness_exponent, leading) == 0)
            if smooth_candidate:
                for divisor_prime in primes:
                    count = 0
                    while True:
                        trial_divisions += 1
                        if remaining % divisor_prime:
                            break
                        remaining //= divisor_prime
                        count += 1
                    if count:
                        parts.append((divisor_prime, count))
                    if remaining == 1:
                        break
            if remaining == 1:
                terms, new = leading_norm_word(integer, target, parts, records)
                if new:
                    raise ValueError('A fixed-base smooth norm introduced a prime.')
                row = [0] * len(primes)
                row[index] = exponent
                for divisor_prime, count in terms.items():
                    row[indices[divisor_prime]] -= count
                rows.append(row)
                selected = exponent
                witnesses.append({'prime': prime, 'exponent': exponent,
                                  'target': target, 'factors': parts, 'word': terms})
                attempts.append({'prime': prime, 'queries': offset + 1,
                                 'selected_exponent': exponent, 'terminal_gcd': False})
                # This is only a cheap trigger for trying extraction early.
                # A final extraction is also made, so identity half-relations
                # cannot hide a useful top behind an unchanged raw kernel.
                bits = sum((entry & 1) << j for j, entry in enumerate(row))
                while bits:
                    column = bits.bit_length() - 1
                    if column not in parity_pivots:
                        parity_pivots[column] = bits
                        break
                    bits ^= parity_pivots[column]
                if bits == 0:
                    factor = extract()
                    if factor is not None:
                        return finish(factor, 'smooth-power relation halving')
                break
            if offset + 1 < window:
                target = compose_with_prime_form(target, prime, *records[prime])
        if selected is None:
            attempts.append({'prime': prime, 'queries': window,
                             'selected_exponent': None, 'terminal_gcd': False})
    if rows and last_extraction_rows != len(rows):
        factor = extract()
        if factor is not None:
            return finish(factor, 'smooth-power relation halving')
    return finish(None, 'complete fixed-base presentation' if len(rows) == len(primes)
                  else 'power window exhausted')
