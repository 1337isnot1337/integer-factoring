"""Height-only factoring consumer allowing arbitrary represented smooth norms.

Unsigned candidates may have imprimitive representations and odd inert
prime factors. Complete fixed-norm decoding supplies a primitive norm
divisor to the existing weighted consumer. Norm selection is external.
"""

from math import gcd

from arbitrary_norm_representation import represent_arbitrary_factored_norm
from arithmetic_two_group import power
from height_selected_norm_factoring import factor_from_height_selected_norms
from prime_form_factor_base import build_prime_form_factor_base, primes_through
from quadratic_ideal_hnf import compose_forms_hnf
from weighted_norm_relations import least_power_above


def factor_from_arbitrary_selected_norms(integer, cutoff, height_bound,
                                         selected_norms, *, extract_each=True,
                                         keep_certificates=False):
    """Consume only (target_prime,n), with any represented odd smooth n<=H.

    The target is g_p**K_p for the least K_p with p**K_p>H. No
    factorization, representation, primitivity, or split-support promise
    is supplied by the caller. Odd norms coprime to N are decoded; a
    proper gcd is returned immediately through the existing consumer.
    Even candidates are outside this wrapper's shortcut and are rejected.

    The fixed-norm decoder strips integer content and sends a primitive
    norm divisor to the weighted consumer. Re-decoding that divisor adds
    at most another instance of the same worst-case processing bound.
    """
    if not isinstance(integer, int) or integer <= 1 or integer % 2 == 0:
        raise ValueError('The target must be an odd integer greater than one.')
    if not isinstance(cutoff, int) or cutoff < 2:
        raise ValueError('The prime cutoff must be an integer at least two.')
    if not isinstance(height_bound, int) or height_bound < 1:
        raise ValueError('The norm height bound must be a positive integer.')
    discriminant = -4*integer
    base = build_prime_form_factor_base(discriminant, cutoff)
    records = {record.prime: record for record in base['records']}
    all_primes = tuple(primes_through(cutoff))
    identity = (1, 0, integer)
    targets, exponents, completed = {}, {}, set()
    translations = []
    stats = {
        'factor_base': base['stats'], 'selected_integers_examined': 0,
        'candidate_bits_examined': 0, 'over_height_candidates': 0,
        'noncoprime_candidates': 0, 'even_candidates': 0,
        'nonsmooth_candidates': 0, 'unrepresented_candidates': 0,
        'odd_inert_exponent_candidates': 0, 'duplicate_completed_targets': 0,
        'primitive_norms_supplied': 0, 'norms_with_removed_content': 0,
        'norms_with_forced_inert_content': 0,
        'proper_candidate_gcds': 0,
        'division_tests': 0, 'prime_divisions': 0,
        'target_group_multiplications': 0,
        'decoder_group_multiplications': 0,
        'maximum_choice_product': 0, 'maximum_table_capacity': 0,
        'maximum_content_bits': 0,
        'coordinates_supplied': False, 'factorizations_supplied': False,
        'orientations_supplied': False, 'primitivity_promised': False,
        'split_support_promised': False,
    }

    def multiply(left, right):
        stats['target_group_multiplications'] += 1
        return compose_forms_hnf(left, right)

    def primitive_candidates():
        for prime, norm in selected_norms:
            if prime not in records:
                raise ValueError('The selected target prime is absent from the factor base.')
            if not isinstance(norm, int) or norm < 1:
                raise ValueError('A selected norm must be a positive integer.')
            stats['selected_integers_examined'] += 1
            stats['candidate_bits_examined'] += norm.bit_length()
            if prime in completed:
                stats['duplicate_completed_targets'] += 1
                continue
            if norm > height_bound:
                stats['over_height_candidates'] += 1
                continue
            common = gcd(norm, integer)
            if 1 < common < integer:
                stats['proper_candidate_gcds'] += 1
                yield prime, norm
                return
            if common != 1:
                stats['noncoprime_candidates'] += 1
                continue
            if norm % 2 == 0:
                stats['even_candidates'] += 1
                continue
            remaining, factors = norm, []
            for candidate_prime in all_primes:
                if remaining == 1:
                    break
                exponent = 0
                while True:
                    stats['division_tests'] += 1
                    quotient, remainder = divmod(remaining, candidate_prime)
                    if remainder:
                        break
                    remaining = quotient
                    exponent += 1
                    stats['prime_divisions'] += 1
                if exponent:
                    record = records.get(candidate_prime)
                    factors.append((candidate_prime, exponent, None if record is None else record.root))
            if remaining != 1:
                stats['nonsmooth_candidates'] += 1
                continue
            if prime not in targets:
                exponents[prime] = least_power_above(prime, height_bound)
                targets[prime] = power(records[prime].positive, exponents[prime], identity, multiply)
            decoded = represent_arbitrary_factored_norm(targets[prime], factors)
            work = decoded['stats']
            stats['decoder_group_multiplications'] += work['group_multiplications']
            stats['maximum_choice_product'] = max(stats['maximum_choice_product'], work['choice_product'])
            stats['maximum_table_capacity'] = max(
                stats['maximum_table_capacity'], work['left_capacity'], work['right_capacity'])
            if decoded['representation'] is None:
                stats['unrepresented_candidates'] += 1
                stats['odd_inert_exponent_candidates'] += work['odd_inert_exponent']
                continue
            primitive = decoded['primitive_representation']
            content = decoded['integer_content']
            assert content*content*primitive.norm == norm
            assert 1 <= primitive.norm <= norm <= height_bound
            stats['primitive_norms_supplied'] += 1
            stats['norms_with_removed_content'] += content > 1
            stats['norms_with_forced_inert_content'] += work['inert_primes'] > 0
            stats['maximum_content_bits'] = max(stats['maximum_content_bits'], content.bit_length())
            completed.add(prime)
            if keep_certificates:
                translations.append({
                    'generator_prime': prime, 'target_exponent': exponents[prime],
                    'target_form': list(targets[prime]), 'selected_norm': norm,
                    'primitive_norm': primitive.norm, 'integer_content': content,
                    'selected_norm_vector': list(decoded['representation'].vector),
                    'primitive_vector': list(primitive.vector),
                    'prime_terms': [list(term) for term in decoded['prime_terms']],
                })
            yield prime, primitive.norm

    result = factor_from_height_selected_norms(
        integer, cutoff, height_bound, primitive_candidates(),
        extract_each=extract_each, keep_certificates=keep_certificates,
    )
    result['arbitrary_norm_preprocessing'] = stats
    if keep_certificates:
        result['norm_translations'] = translations
    if (result['factor'] is not None and stats['primitive_norms_supplied']
            and not stats['proper_candidate_gcds']):
        result['method'] = 'arbitrary selected norms, primitive divisors, and weighted relation halving'
    return result
