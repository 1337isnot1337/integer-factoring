"""Construct exact prime-form relations by cancelling residual ideals.

The residual integer need not be prime or factored.  A bounded run is a
finite search, not a uniform short-window theorem.  All norms, rejected
candidates, and retained residual keys are charged in the returned counters.
"""

from fractions import Fraction
from math import gcd, isqrt, prod

from arithmetic_two_group import power
from binary_quadratic_forms import ambiguous_factor, inverse_form, reduce_form
from leading_power_smooth_relations import compose_with_prime_form
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top
from weighted_norm_relations import least_power_above


class RationalRowBasis:
    """Select independent original rows by exact rational elimination."""

    def __init__(self, dimension):
        self.dimension = dimension
        self.pivots = {}
        self.rows = []

    def add(self, row):
        if len(row) != self.dimension or any(not isinstance(x, int) for x in row):
            raise ValueError('One integer entry per generator is required.')
        residual = {j: Fraction(x) for j, x in enumerate(row) if x}
        for pivot in sorted(self.pivots):
            coefficient = residual.get(pivot, 0)
            if not coefficient:
                continue
            for j, value in self.pivots[pivot].items():
                value = residual.get(j, 0) - coefficient * value
                if value:
                    residual[j] = value
                else:
                    residual.pop(j, None)
        if not residual:
            return False
        pivot = min(residual)
        denominator = residual[pivot]
        self.pivots[pivot] = {j: value / denominator for j, value in residual.items()}
        self.rows.append(list(row))
        return True


def partial_leading_word(integer, form, primes, records, prime_product=None, *,
                         known_smooth_part=None):
    """Separate all known-prime factors from a primitive leading norm.

    Return a signed base word and a canonically oriented residual ideal.
    No trial division by an unknown integer and no residual primality test
    occurs.  The caller handles a proper gcd with the factoring input first.
    """
    delta = -integer if integer % 4 == 3 else -4 * integer
    leading, middle, trailing = form
    if (leading < 1 or middle * middle - 4 * leading * trailing != delta
            or gcd(gcd(leading, middle), trailing) != 1
            or gcd(leading, integer) != 1 or tuple(form) != reduce_form(form)):
        raise ValueError('A reduced primitive form with coprime leading norm is required.')
    if prime_product is None:
        prime_product = prod(primes)
    smooth = (gcd(leading, pow(prime_product, leading.bit_length(), leading))
              if known_smooth_part is None else known_smooth_part)
    if not isinstance(smooth, int) or smooth < 1 or leading % smooth:
        raise ValueError('The supplied smooth part must divide the leading norm.')
    residual = leading // smooth
    remaining, parts, word = smooth, [], {}
    divisions = 0
    for prime in primes:
        exponent = 0
        while True:
            divisions += 1
            if remaining % prime:
                break
            remaining //= prime
            exponent += 1
        if exponent:
            parts.append((prime, exponent))
            if prime == 2 and integer % 4 == 1:
                if exponent != 1 or middle % 4 != 2:
                    raise ValueError('The ramified-two orientation invariant failed.')
                word[prime] = 1
            else:
                root = records[prime][0]
                residue = middle % (2 * prime)
                if residue not in (root, (-root) % (2 * prime)):
                    raise ValueError('The known prime orientation is inconsistent.')
                word[prime] = exponent if residue == root else -exponent
        if remaining == 1:
            break
    if remaining != 1 or gcd(smooth, residual) != 1:
        raise ArithmeticError('The modular smooth-part split failed.')
    if residual == 1:
        key, orientation, residual_form = None, 1, None
    else:
        if gcd(residual, 2 * delta) != 1:
            raise ValueError('The residual ideal must have odd unramified norm.')
        residue = middle % (2 * residual)
        canonical = min(residue, (-residue) % (2 * residual))
        orientation = 1 if residue == canonical else -1
        key = (residual, canonical)
        if (canonical * canonical - delta) % (4 * residual):
            raise ArithmeticError('The residual root does not define an ideal.')
        residual_form = reduce_form((residual, canonical,
                                    (canonical * canonical - delta) // (4 * residual)))
    return {'smooth_part': smooth, 'residual': residual, 'factors': parts,
            'word': word, 'key': key, 'orientation': orientation,
            'residual_form': residual_form, 'trial_divisions': divisions}


def word_value(word, records, identity):
    result = identity
    for prime, exponent in word.items():
        form = records[prime][1]
        if exponent < 0:
            form, exponent = inverse_form(form), -exponent
        result = compose_forms_hnf(result, power(form, exponent, identity, compose_forms_hnf))
    return result


def verify_partial_identity(integer, form, split, records):
    delta = -integer if integer % 4 == 3 else -4 * integer
    identity = (1, delta % 2, (delta % 2 - delta) // 4)
    value = word_value(split['word'], records, identity)
    if split['residual_form'] is not None:
        residual = split['residual_form']
        if split['orientation'] < 0:
            residual = inverse_form(residual)
        value = compose_forms_hnf(value, residual)
    return value == tuple(form)


def factor_with_cofactor_collisions(integer, cutoff, window, *,
                                   max_cofactor=None, max_queries=None,
                                   verify_candidates=False, progress=None):
    """Bounded round-robin prime-power scan with exact cofactor collisions.

    ``max_cofactor=1`` is the matched smooth-only traversal.  Set a larger
    bound to retain partially smooth norms.  A missing factor certifies
    only this finite search's outcome.  Initial factor-base generation is
    not claimed to satisfy a uniform analytic cutoff in finite controls.
    """
    if not isinstance(integer, int) or integer < 3 or integer % 2 == 0:
        raise ValueError('An odd integer of at least three is required.')
    if not isinstance(cutoff, int) or cutoff < 2 or not isinstance(window, int) or window < 1:
        raise ValueError('Require cutoff >= 2 and a positive integer window.')
    if max_queries is not None and (not isinstance(max_queries, int) or max_queries < 0):
        raise ValueError('The query cap must be a nonnegative integer.')
    delta = -integer if integer % 4 == 3 else -4 * integer
    identity = (1, delta % 2, (delta % 2 - delta) // 4)
    height = isqrt((-delta) // 3)
    if max_cofactor is None:
        max_cofactor = height
    if not isinstance(max_cofactor, int) or max_cofactor < 1:
        raise ValueError('The cofactor cap must be a positive integer.')
    base = build_prime_form_factor_base(delta, cutoff)
    records = {r.prime: (r.root, r.positive) for r in base['records']}
    if integer % 4 == 1:
        records[2] = (2, reduce_form((2, 2, (integer + 1) // 2)))
    primes = sorted(records)
    indices = {p: i for i, p in enumerate(primes)}
    prime_product = prod(primes)
    basis, residuals, witnesses, parity_pivots = RationalRowBasis(len(primes)), {}, [], {}
    counters = dict(norm_queries=0, modular_smooth_part_tests=0, trial_divisions=0,
                    cofactor_rejections=0, retained_cofactor_keys=0, cofactor_collisions=0,
                    smooth_candidates=0, candidate_relations=0, independent_collisions=0,
                    independent_smooth_rows=0, dependent_relations=0, extractions=0,
                    maximum_exponent_tested=0, maximum_leading_bits=0,
                    maximum_retained_cofactor_bits=0)
    last_extraction, last_extraction_rows = None, -1

    def extract():
        nonlocal last_extraction, last_extraction_rows
        counters['extractions'] += 1
        _, factor, last_extraction = first_useful_relation_top(
            [records[p][1] for p in primes], basis.rows, identity,
            compose_forms_hnf, inverse_form,
            lambda form: ambiguous_factor(form, integer), verify_relations=True)
        last_extraction_rows = len(basis.rows)
        return factor

    def finish(factor, method, terminal=None):
        return {'integer': integer, 'input_bits': integer.bit_length(),
                'discriminant': delta, 'height_bound': height, 'cutoff': cutoff,
                'window': window, 'max_cofactor': max_cofactor, 'max_queries': max_queries,
                'factor': factor, 'method': method, 'base_size': len(primes),
                'relation_rank': len(basis.rows), 'active_primes': primes,
                'prime_records': records, 'relation_matrix': basis.rows,
                'witnesses': witnesses, 'last_extraction': last_extraction,
                'terminal_witness': terminal, 'factoring_oracle_calls': 0,
                'residual_primality_tests': 0, 'new_primes_introduced': 0, **counters}

    for prime in base['excluded_primes']:
        divisor = gcd(prime, integer)
        if 1 < divisor < integer:
            return finish(divisor, 'factor-base gcd', {'prime': prime})
    starts = {p: least_power_above(p, height) for p in primes}
    targets = {p: power(records[p][1], starts[p], identity, compose_forms_hnf) for p in primes}
    for offset in range(window):
        for prime in primes:
            if max_queries is not None and counters['norm_queries'] >= max_queries:
                if basis.rows and last_extraction_rows != len(basis.rows):
                    factor = extract()
                    if factor is not None:
                        return finish(factor, 'cofactor relation halving')
                return finish(None, 'query cap exhausted')
            target, exponent = targets[prime], starts[prime] + offset
            counters['norm_queries'] += 1
            counters['maximum_exponent_tested'] = max(counters['maximum_exponent_tested'], exponent)
            counters['maximum_leading_bits'] = max(counters['maximum_leading_bits'], target[0].bit_length())
            terminal = {'prime': prime, 'exponent': exponent, 'target': target}
            if target != identity and target == inverse_form(target):
                divisor = ambiguous_factor(target, integer)
                if divisor is not None:
                    return finish(divisor, 'target involution', terminal)
            divisor = gcd(target[0], integer)
            if 1 < divisor < integer:
                return finish(divisor, 'leading-norm gcd', terminal)
            if divisor != 1 or not 1 <= target[0] <= height:
                raise ArithmeticError('The leading-norm invariant failed.')
            counters['modular_smooth_part_tests'] += 1
            smooth = gcd(target[0], pow(prime_product, target[0].bit_length(), target[0]))
            residual = target[0] // smooth
            if residual > max_cofactor:
                counters['cofactor_rejections'] += 1
            else:
                split = partial_leading_word(integer, target, primes, records, prime_product,
                                             known_smooth_part=smooth)
                counters['trial_divisions'] += split['trial_divisions']
                if verify_candidates and not verify_partial_identity(integer, target, split, records):
                    raise ArithmeticError('A partial ideal identity failed.')
                row = {indices[prime]: exponent}
                for p, count in split['word'].items():
                    j = indices[p]
                    row[j] = row.get(j, 0) - count
                    if not row[j]:
                        del row[j]
                witness = {**terminal, **split, 'vector': sorted(row.items())}
                previous = None
                if split['key'] is None:
                    counters['smooth_candidates'] += 1
                    candidate = row
                else:
                    counters['maximum_retained_cofactor_bits'] = max(
                        counters['maximum_retained_cofactor_bits'], residual.bit_length())
                    signed = {j: split['orientation'] * value for j, value in row.items()}
                    key = split['key']
                    if key not in residuals:
                        residuals[key] = (signed, witness)
                        counters['retained_cofactor_keys'] += 1
                        candidate = None
                    else:
                        counters['cofactor_collisions'] += 1
                        old, previous = residuals[key]
                        candidate = dict(signed)
                        for j, value in old.items():
                            candidate[j] = candidate.get(j, 0) - value
                if candidate is not None:
                    counters['candidate_relations'] += 1
                    dense = [candidate.get(j, 0) for j in range(len(primes))]
                    if basis.add(dense):
                        counters['independent_collisions' if previous is not None
                                 else 'independent_smooth_rows'] += 1
                        witnesses.append({'current': witness, 'previous': previous})
                        bits = sum((entry & 1) << j for j, entry in enumerate(dense))
                        while bits:
                            pivot = bits.bit_length() - 1
                            if pivot not in parity_pivots:
                                parity_pivots[pivot] = bits
                                break
                            bits ^= parity_pivots[pivot]
                        if (not bits or len(basis.rows) % 8 == 0
                                or len(basis.rows) == len(primes)):
                            factor = extract()
                            if factor is not None:
                                return finish(factor, 'cofactor relation halving')
                        if len(basis.rows) == len(primes):
                            return finish(None, 'complete fixed-base presentation')
                    else:
                        counters['dependent_relations'] += 1
            if offset + 1 < window:
                targets[prime] = compose_with_prime_form(target, prime, *records[prime])
        if progress is not None:
            progress({'round': offset + 1, 'relation_rank': len(basis.rows), **counters})
    if basis.rows and last_extraction_rows != len(basis.rows):
        factor = extract()
        if factor is not None:
            return finish(factor, 'cofactor relation halving')
    return finish(None, 'complete fixed-base presentation' if len(basis.rows) == len(primes)
                  else 'power window exhausted')
