"""Exact factor extraction from selected smooth integers with a height bound.

Each target exponent is chosen from its prime's weight. No multiplicity
promise is imposed on a candidate norm. Selecting the integers remains
an external search cost.
"""

from math import gcd

from arithmetic_two_group import power
from factored_norm_representation import represent_factored_norm
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from weighted_norm_relations import least_power_above, factor_from_weighted_relations


def factor_from_height_selected_norms(integer, cutoff, height_bound,
                                      selected_norms, *, extract_each=True,
                                      keep_certificates=False):
    """Consume (p,n) for targets g_p**K_p, where p**K_p > height_bound.

    K_p is the least such positive exponent. Each accepted n is positive,
    at most height_bound, coprime to -4N, and primitively represented by
    its target. Its prime support must be in the constructed base.
    The routine computes factorizations, orientations, and coordinates;
    none of those data are supplied by the caller.

    One weighted dominant relation is kept per successful target.
    Failures and retries are permitted. Every returned factor is exact.
    Failure describes only the accepted partial presentation's torsion
    image. Full generation plus a valid answer for every target gives
    the full guarantee for an integer with two distinct prime factors.
    """
    if not isinstance(integer, int) or integer <= 1 or integer % 2 == 0:
        raise ValueError('The target must be an odd integer greater than one.')
    if not isinstance(cutoff, int) or cutoff < 2:
        raise ValueError('The prime cutoff must be an integer at least two.')
    if not isinstance(height_bound, int) or height_bound < 1:
        raise ValueError('The norm height bound must be a positive integer.')
    discriminant = -4*integer
    base = build_prime_form_factor_base(discriminant, cutoff)
    records = base['records']
    positions = {record.prime: index for index, record in enumerate(records)}
    generators = tuple(record.positive for record in records)
    primes = tuple(record.prime for record in records)
    identity = (1, 0, integer)
    rows, pivots, completed = [], [], set()
    target_exponents, targets = {}, {}
    certificates = []
    last_extraction = None
    stats = {
        'factor_base': base['stats'], 'selected_integers_examined': 0,
        'accepted_relations': 0, 'duplicate_completed_targets': 0,
        'over_height_candidates': 0, 'noncoprime_candidates': 0,
        'nonsmooth_candidates': 0, 'unrepresented_candidates': 0,
        'division_tests': 0, 'prime_divisions': 0,
        'target_group_multiplications': 0,
        'decoder_group_multiplications': 0,
        'extraction_calls': 0, 'extraction_group_multiplications': 0,
        'maximum_norm_bits': 0, 'maximum_multiplicity': 0,
        'maximum_distinct_prime_support': 0, 'maximum_target_exponent': 0,
        'multiplicity_budget_imposed': False,
        'coordinates_supplied': False, 'orientations_supplied': False,
        'factorizations_supplied': False, 'determinant_computed': False,
        'complete_class_group_enumerated': False,
    }

    def finish(factor, method):
        if factor is not None:
            assert 1 < factor < integer and integer % factor == 0
        result = {
            'integer': integer, 'factor': factor, 'method': method,
            'prime_cutoff': cutoff, 'height_bound': height_bound,
            'factor_base_primes': list(primes),
            'accepted_target_primes': [primes[index] for index in pivots],
            'stats': stats, 'last_extraction': last_extraction,
        }
        if keep_certificates:
            result['certificates'] = certificates
        return result

    for prime in base['excluded_primes']:
        if 1 < prime < integer and integer % prime == 0:
            return finish(prime, 'factor-base trial division')

    def multiply(left, right):
        stats['target_group_multiplications'] += 1
        return compose_forms_hnf(left, right)

    def extract():
        nonlocal last_extraction
        last_extraction = factor_from_weighted_relations(
            integer, generators, rows, primes, pivots,
            height_bound=height_bound, verify_relations=False,
        )
        # The orientation decoder has already certified each word equals
        # its explicitly computed target, so all rows are genuine.
        stats['extraction_calls'] += 1
        stats['extraction_group_multiplications'] += (
            last_extraction['stats']['group_multiplications'])
        return last_extraction['factor']

    for prime, norm in selected_norms:
        if prime not in positions:
            raise ValueError('The selected target prime is absent from the factor base.')
        if not isinstance(norm, int) or norm < 1:
            raise ValueError('A selected norm must be a positive integer.')
        stats['selected_integers_examined'] += 1
        index = positions[prime]
        if index in completed:
            stats['duplicate_completed_targets'] += 1
            continue
        if norm > height_bound:
            stats['over_height_candidates'] += 1
            continue
        common = gcd(norm, integer)
        if 1 < common < integer:
            return finish(common, 'selected-norm gcd')
        if common != 1:
            stats['noncoprime_candidates'] += 1
            continue
        remaining, multiplicity, factors = norm, 0, []
        for record in records:
            if remaining == 1:
                break
            exponent = 0
            while True:
                stats['division_tests'] += 1
                quotient, remainder = divmod(remaining, record.prime)
                if remainder:
                    break
                remaining = quotient
                exponent += 1
                multiplicity += 1
                stats['prime_divisions'] += 1
            if exponent:
                factors.append((record.prime, exponent, record.root))
        if remaining != 1:
            stats['nonsmooth_candidates'] += 1
            continue
        if index not in targets:
            exponent = least_power_above(prime, height_bound)
            target_exponents[index] = exponent
            targets[index] = power(generators[index], exponent, identity, multiply)
            stats['maximum_target_exponent'] = max(
                stats['maximum_target_exponent'], exponent)
        decoded = represent_factored_norm(targets[index], factors)
        stats['decoder_group_multiplications'] += decoded['stats']['group_multiplications']
        represented = decoded['representation']
        if represented is None:
            stats['unrepresented_candidates'] += 1
            continue
        row = [0]*len(records)
        row[index] = target_exponents[index]
        for word_prime, exponent, _ in decoded['prime_terms']:
            row[positions[word_prime]] -= exponent
        assert prime**target_exponents[index] > height_bound >= norm
        rows.append(tuple(row))
        pivots.append(index)
        completed.add(index)
        stats['accepted_relations'] += 1
        stats['maximum_norm_bits'] = max(stats['maximum_norm_bits'], norm.bit_length())
        stats['maximum_multiplicity'] = max(stats['maximum_multiplicity'], multiplicity)
        stats['maximum_distinct_prime_support'] = max(
            stats['maximum_distinct_prime_support'], len(factors))
        if keep_certificates:
            certificates.append({
                'generator_prime': prime,
                'target_exponent': target_exponents[index],
                'target_form': list(targets[index]), 'norm': norm,
                'primitive_vector': list(represented.vector),
                'prime_terms': [list(term) for term in decoded['prime_terms']],
                'relation_row': list(row), 'pivot_column': index,
            })
        if extract_each:
            factor = extract()
            if factor is not None:
                return finish(factor, 'height-only selected norms and partial relation halving')
    if not extract_each or last_extraction is None:
        factor = extract()
        if factor is not None:
            return finish(factor, 'height-only selected norms and relation halving')
    return finish(None, 'no useful top in the accepted weighted relation torsion image')
