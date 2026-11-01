"""Consume selected smooth integer norms and attempt exact factor extraction.

The caller supplies only (generator_prime, candidate_integer) pairs. Roots,
orientations, primitive coordinates, relations, and the final factor are
computed here. Selecting useful candidate integers is a separate cost.
"""

from math import gcd

from arithmetic_two_group import power
from factored_norm_representation import represent_factored_norm
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import factor_from_partial_relations


def factor_from_selected_norms(integer, cutoff, multiplicity_budget,
                               selected_norms, *, extract_each=True,
                               keep_certificates=False):
    """Try selected integers for targets g_p**(multiplicity_budget+1).

    ``integer`` is odd and greater than one. The factor base contains one
    chosen orientation of every split prime p<=cutoff coprime to -4N.
    ``selected_norms`` is an iterable of (p,n) pairs. A failed candidate
    may be followed by another for the same target. Successful targets
    occur only once among the independent relation rows.

    Every returned factor is exact. A missing factor concerns only the
    torsion image of the accepted partial relations. Full generation and
    eventual valid candidates for all targets give the theorem's full
    guarantee; this function does not find those candidates by itself.
    ``extract_each`` tries every accepted prefix. False tries only the
    final set, which can reduce repeated finite linear algebra.
    """
    if not isinstance(integer, int) or integer <= 1 or integer % 2 == 0:
        raise ValueError("The target must be an odd integer greater than one.")
    if not isinstance(cutoff, int) or cutoff < 2:
        raise ValueError("The prime cutoff must be an integer at least two.")
    if (not isinstance(multiplicity_budget, int)
            or multiplicity_budget < 0):
        raise ValueError("The multiplicity budget must be a nonnegative integer.")
    discriminant = -4 * integer
    base = build_prime_form_factor_base(discriminant, cutoff)
    records = base['records']
    positions = {record.prime: index for index, record in enumerate(records)}
    generators = tuple(record.positive for record in records)
    identity = (1, 0, integer)
    diagonal = multiplicity_budget + 1
    rows, pivots, accepted, targets = [], [], set(), {}
    stats = {
        'factor_base': base['stats'],
        'selected_integers_examined': 0,
        'duplicate_completed_targets': 0,
        'accepted_relations': 0,
        'noncoprime_candidates': 0,
        'over_budget_candidates': 0,
        'nonsmooth_candidates': 0,
        'unrepresented_candidates': 0,
        'division_tests': 0,
        'prime_divisions': 0,
        'target_group_multiplications': 0,
        'decoder_group_multiplications': 0,
        'extraction_calls': 0,
        'extraction_group_multiplications': 0,
        'maximum_norm_bits': 0,
        'maximum_multiplicity': 0,
        'maximum_distinct_prime_support': 0,
        'complete_class_group_enumerated': False,
        'coordinates_supplied': False,
        'orientations_supplied': False,
        'factorizations_supplied': False,
        'determinant_computed': False,
        'annihilator_supplied': False,
    }
    certificates = []
    last_extraction = None

    def finish(factor, method):
        if factor is not None:
            assert 1 < factor < integer and integer % factor == 0
        result = {
            'integer': integer, 'factor': factor, 'method': method,
            'prime_cutoff': cutoff, 'multiplicity_budget': multiplicity_budget,
            'diagonal_parameter': diagonal,
            'factor_base_primes': [record.prime for record in records],
            'accepted_target_primes': [records[index].prime for index in pivots],
            'stats': stats, 'last_extraction': last_extraction,
        }
        if keep_certificates:
            result['certificates'] = certificates
        return result

    # A discriminant prime encountered while constructing the base can
    # already be a proper factor. Excluding it from the base loses none.
    for prime in base['excluded_primes']:
        if 1 < prime < integer and integer % prime == 0:
            return finish(prime, 'factor-base trial division')

    def multiply(left, right):
        stats['target_group_multiplications'] += 1
        return compose_forms_hnf(left, right)

    def extract():
        nonlocal last_extraction
        last_extraction = factor_from_partial_relations(
            integer, generators, rows, pivots, verify_relations=False,
        )
        # Each identity is already certified by the orientation decoder:
        # its signed prime word equals the explicitly computed target.
        stats['extraction_calls'] += 1
        stats['extraction_group_multiplications'] += (
            last_extraction['stats']['group_multiplications'])
        return last_extraction['factor']

    for prime, norm in selected_norms:
        if prime not in positions:
            raise ValueError("The selected target prime is absent from the factor base.")
        if not isinstance(norm, int) or norm < 1:
            raise ValueError("A selected norm must be a positive integer.")
        stats['selected_integers_examined'] += 1
        index = positions[prime]
        if index in accepted:
            stats['duplicate_completed_targets'] += 1
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
                if multiplicity > multiplicity_budget:
                    break
            if exponent:
                factors.append((record.prime, exponent, record.root))
            if multiplicity > multiplicity_budget:
                break
        if multiplicity > multiplicity_budget:
            stats['over_budget_candidates'] += 1
            continue
        if remaining != 1:
            stats['nonsmooth_candidates'] += 1
            continue
        if index not in targets:
            targets[index] = power(generators[index], diagonal, identity, multiply)
        decoded = represent_factored_norm(targets[index], factors)
        stats['decoder_group_multiplications'] += (
            decoded['stats']['group_multiplications'])
        represented = decoded['representation']
        if represented is None:
            stats['unrepresented_candidates'] += 1
            continue
        row = [0] * len(records)
        row[index] = diagonal
        for word_prime, exponent, _ in decoded['prime_terms']:
            row[positions[word_prime]] -= exponent
        # A strict minor remains after any ordering or omission of targets.
        assert abs(row[index]) > sum(abs(value) for j, value in enumerate(row)
                                     if j != index)
        rows.append(tuple(row))
        pivots.append(index)
        accepted.add(index)
        stats['accepted_relations'] += 1
        stats['maximum_norm_bits'] = max(stats['maximum_norm_bits'], norm.bit_length())
        stats['maximum_multiplicity'] = max(stats['maximum_multiplicity'], multiplicity)
        stats['maximum_distinct_prime_support'] = max(
            stats['maximum_distinct_prime_support'], len(factors))
        if keep_certificates:
            certificates.append({
                'generator_prime': prime, 'target_form': list(targets[index]),
                'norm': norm, 'primitive_vector': list(represented.vector),
                'prime_terms': [list(term) for term in decoded['prime_terms']],
                'relation_row': list(row), 'pivot_column': index,
            })
        if extract_each:
            factor = extract()
            if factor is not None:
                return finish(factor, 'selected norms and partial relation halving')

    if not extract_each or last_extraction is None:
        factor = extract()
        if factor is not None:
            return finish(factor, 'selected norms and relation halving')
    return finish(None, 'no useful top in the accepted relation torsion image')
