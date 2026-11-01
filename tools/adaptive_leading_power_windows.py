"""Extend only unfinished prime-form power windows, retaining exact rows.

Every earlier candidate remains charged. Accepted targets are not searched
again. This changes finite scheduling, not the unresolved uniform window bound.
"""

from copy import deepcopy
from math import gcd, isqrt, prod

from arithmetic_two_group import power
from binary_quadratic_forms import ambiguous_factor, inverse_form, reduce_form
from dynamic_leading_norm_relations import leading_norm_word
from leading_power_smooth_relations import (compose_with_prime_form,
                                           factor_with_smooth_leading_powers)
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top
from weighted_norm_relations import least_power_above


def _validated_state(record):
    """Validate a completed bounded phase before using any retained relation."""
    integer, cutoff = record['integer'], record['cutoff']
    if (not isinstance(integer, int) or integer < 3 or integer % 2 == 0
            or not isinstance(cutoff, int) or cutoff < 2):
        raise ValueError('Invalid integer or prime cutoff.')
    delta = -integer if integer % 4 == 3 else -4 * integer
    identity = (1, delta % 2, (delta % 2 - delta) // 4)
    height = isqrt((-delta) // 3)
    base = build_prime_form_factor_base(delta, cutoff)
    records = {r.prime: (r.root, tuple(r.positive)) for r in base['records']}
    if integer % 4 == 1:
        records[2] = (2, reduce_form((2, 2, (integer + 1) // 2)))
    supplied = {int(p): (r[0], tuple(r[1])) for p, r in record['prime_records'].items()}
    primes = sorted(records)
    if (supplied != records or record['active_primes'] != primes
            or record['discriminant'] != delta or record['height_bound'] != height
            or record['base_size'] != len(primes)
            or record['smoothness_method'] != 'primorial power'):
        raise ValueError('The retained factor base or height is inconsistent.')
    old_window = record['window']
    if not isinstance(old_window, int) or old_window < 1:
        raise ValueError('Invalid previous window.')
    attempts = {a['prime']: a for a in record['attempts']}
    if len(attempts) != len(record['attempts']) or set(attempts) != set(primes):
        raise ValueError('A resumable phase must have attempted every target once.')
    indices = {p: i for i, p in enumerate(primes)}
    witnesses = record['witnesses']
    rows = record['relation_matrix']
    if len(witnesses) != len(rows) or record['relation_rank'] != len(rows):
        raise ValueError('The retained row count is inconsistent.')
    seen = set()
    expected_divisions = 0
    for witness, row in zip(witnesses, rows):
        prime, exponent = witness['prime'], witness['exponent']
        if prime not in records or prime in seen:
            raise ValueError('Repeated or unknown relation target.')
        seen.add(prime)
        first = least_power_above(prime, height)
        attempt = attempts[prime]
        if (exponent != attempt['selected_exponent']
                or attempt['queries'] != exponent - first + 1
                or not first <= exponent < first + old_window):
            raise ValueError('The retained exponent and query count disagree.')
        target = power(records[prime][1], exponent, identity, compose_forms_hnf)
        parts = tuple((int(p), int(e)) for p, e in witness['factors'])
        if (tuple(witness['target']) != target or not 1 <= target[0] <= height
                or gcd(target[0], integer) != 1
                or any(p not in records or e < 1 or e > height.bit_length()
                       for p, e in parts)
                or len({p for p, _ in parts}) != len(parts)
                or prod(p ** e for p, e in parts) != target[0]):
            raise ValueError('A retained smooth norm is invalid.')
        terms, new = leading_norm_word(integer, target, parts, records)
        expected = [0] * len(primes)
        expected[indices[prime]] = exponent
        for p, e in terms.items():
            expected[indices[p]] -= e
        if (new or {int(p): e for p, e in witness['word'].items()} != terms
                or row != expected):
            raise ValueError('A retained relation does not match its norm.')
        expected_divisions += (indices[max(p for p, _ in parts)] + 1
                               + sum(e for _, e in parts)) if parts else 1
    for prime, attempt in attempts.items():
        if attempt.get('terminal_gcd', False):
            raise ValueError('A terminal factor record does not need resuming.')
        if prime not in seen and (attempt['selected_exponent'] is not None
                                 or attempt['queries'] != old_window):
            raise ValueError('An unfinished target has an incomplete prefix.')
    if (record['norm_queries'] != sum(a['queries'] for a in attempts.values())
            or record['modular_smoothness_tests'] != record['norm_queries']
            or record['trial_divisions'] != expected_divisions
            or record['maximum_exponent_tested'] != max(
                (least_power_above(p, height) + a['queries'] - 1
                 for p, a in attempts.items()), default=0)
            or record['factoring_oracle_calls'] != 0
            or record['new_primes_introduced'] != 0):
        raise ValueError('The retained work counters are inconsistent.')
    return records, primes, indices, attempts, identity


def resume_leading_power_windows(previous, new_window):
    """Extend an exhausted phase to a larger absolute window.

The stored identities are verified. The caller remains responsible for the
claimed history of rejected candidates; full replay regenerates that history.
"""
    if previous.get('factor') is not None:
        raise ValueError('The previous phase already returned a factor.')
    if not isinstance(new_window, int) or new_window <= previous['window']:
        raise ValueError('The new window must be strictly larger.')
    result = deepcopy(previous)
    records, primes, indices, attempts, identity = _validated_state(result)
    integer, height = result['integer'], result['height_bound']
    old_window = result['window']
    rows = result['relation_matrix']
    prime_product, smoothness_exponent = prod(primes), height.bit_length()
    completed = {w['prime'] for w in result['witnesses']}
    parity_pivots = {}

    def add_parity(row):
        bits = sum((e & 1) << j for j, e in enumerate(row))
        while bits:
            pivot = bits.bit_length() - 1
            if pivot not in parity_pivots:
                parity_pivots[pivot] = bits
                return False
            bits ^= parity_pivots[pivot]
        return True

    for row in rows:
        add_parity(row)
    starting_queries, starting_rows = result['norm_queries'], len(rows)
    result['window'] = new_window
    result.setdefault('resume_phases', [])

    def extract():
        _, factor, detail = first_useful_relation_top(
            [records[p][1] for p in primes], rows, identity,
            compose_forms_hnf, inverse_form,
            lambda form: ambiguous_factor(form, integer), verify_relations=True)
        result['last_extraction'] = detail
        return factor

    def finish(factor, method, terminal=None):
        result['factor'], result['method'] = factor, method
        result['terminal_witness'] = terminal
        result['relation_rank'] = len(rows)
        result['resume_phases'].append({
            'previous_window': old_window, 'new_window': new_window,
            'targets_already_completed': len(completed),
            'additional_norm_queries': result['norm_queries'] - starting_queries,
            'additional_relations': len(rows) - starting_rows,
            'total_norm_queries': result['norm_queries'],
        })
        return result

    for prime in primes:
        if prime in completed:
            continue
        first = least_power_above(prime, height)
        target = power(records[prime][1], first + old_window, identity, compose_forms_hnf)
        attempt = attempts[prime]
        for offset in range(old_window, new_window):
            exponent = first + offset
            result['norm_queries'] += 1
            attempt['queries'] += 1
            result['maximum_exponent_tested'] = max(result['maximum_exponent_tested'], exponent)
            leading = target[0]
            result['maximum_leading_bits'] = max(result['maximum_leading_bits'], leading.bit_length())
            divisor = gcd(leading, integer)
            if 1 < divisor < integer:
                attempt['terminal_gcd'] = True
                return finish(divisor, 'resumed leading-norm gcd',
                              {'prime': prime, 'exponent': exponent, 'target': target})
            if divisor != 1 or not 1 <= leading <= height:
                raise ValueError('The reduced leading-norm invariant failed.')
            result['modular_smoothness_tests'] += 1
            if pow(prime_product, smoothness_exponent, leading) == 0:
                remaining, parts = leading, []
                for p in primes:
                    count = 0
                    while True:
                        result['trial_divisions'] += 1
                        if remaining % p:
                            break
                        remaining //= p
                        count += 1
                    if count:
                        parts.append((p, count))
                    if remaining == 1:
                        break
                if remaining != 1:
                    raise ValueError('The primorial test and division disagree.')
                terms, new = leading_norm_word(integer, target, parts, records)
                if new:
                    raise ValueError('A fixed-base norm introduced a new prime.')
                row = [0] * len(primes)
                row[indices[prime]] = exponent
                for p, e in terms.items():
                    row[indices[p]] -= e
                rows.append(row)
                result['witnesses'].append({'prime': prime, 'exponent': exponent,
                                           'target': target, 'factors': parts, 'word': terms})
                attempt['selected_exponent'] = exponent
                if add_parity(row):
                    factor = extract()
                    if factor is not None:
                        return finish(factor, 'resumed smooth-power relation halving')
                break
            if offset + 1 < new_window:
                target = compose_with_prime_form(target, prime, *records[prime])
    if len(rows) > starting_rows:
        factor = extract()
        if factor is not None:
            return finish(factor, 'resumed smooth-power relation halving')
    return finish(None, 'complete fixed-base presentation' if len(rows) == len(primes)
                  else 'power window exhausted')


def factor_with_expanding_power_windows(integer, cutoff, windows):
    windows = tuple(windows)
    if (not windows or any(not isinstance(w, int) or w < 1 for w in windows)
            or any(a >= b for a, b in zip(windows, windows[1:]))):
        raise ValueError('A nonempty strictly increasing window list is required.')
    result = factor_with_smooth_leading_powers(integer, cutoff, windows[0])
    for window in windows[1:]:
        if result['factor'] is not None or result['relation_rank'] == result['base_size']:
            break
        result = resume_leading_power_windows(result, window)
    return result
