"""Compare exact-key and oriented-atom rows on the same power-scan prefix."""

from fractions import Fraction
from math import gcd, isqrt, lcm, prod

from arithmetic_two_group import power
from binary_quadratic_forms import ambiguous_factor, inverse_form, reduce_form
from cofactor_collision_relations import RationalRowBasis, partial_leading_word, verify_partial_identity
from leading_power_smooth_relations import compose_with_prime_form
from oriented_coprime_basis import oriented_coprime_basis, transport_relation_rows
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top
from weighted_norm_relations import least_power_above


def eliminate_atom_columns(rows, base_size):
    """Emit integral base rows with integer-combination certificates.

    Clear denominators in the combination coefficients, not only the output
    row: the latter could silently turn a genuine relation into a half-word.
    """
    if not rows:
        return {'rows': [], 'combinations': [], 'atom_rank': 0}
    width = len(rows[0])
    if any(len(row) != width or any(not isinstance(x, int) for x in row) for row in rows):
        raise ValueError('An integer rectangular matrix is required.')
    if not isinstance(base_size, int) or not 0 <= base_size <= width:
        raise ValueError('The base size must lie within the matrix width.')
    pivots, output, combinations = {}, RationalRowBasis(base_size), []
    for index, original in enumerate(rows):
        row, combination = list(map(Fraction, original)), {index: Fraction(1)}
        for pivot in sorted(pivots):
            coefficient = row[pivot]
            if not coefficient:
                continue
            previous, old_combination = pivots[pivot]
            row = [x - coefficient * y for x, y in zip(row, previous)]
            for i, value in old_combination.items():
                new = combination.get(i, 0) - coefficient * value
                if new:
                    combination[i] = new
                else:
                    combination.pop(i, None)
        pivot = next((j for j in range(base_size, width) if row[j]), None)
        if pivot is not None:
            coefficient = row[pivot]
            pivots[pivot] = ([x / coefficient for x in row],
                             {i: x / coefficient for i, x in combination.items()})
        elif any(row[:base_size]):
            denominator = lcm(*(value.denominator for value in combination.values()))
            integral = [int(denominator * x) for x in row[:base_size]]
            multipliers = {i: int(denominator * value) for i, value in combination.items()}
            literal = [sum(c * rows[i][j] for i, c in multipliers.items()) for j in range(width)]
            assert literal == integral + [0] * (width - base_size)
            if output.add(integral):
                combinations.append(sorted(multipliers.items()))
    return {'rows': output.rows, 'combinations': combinations, 'atom_rank': len(pivots)}


def compare_power_prefix(integer, cutoff, max_queries, *, batch_size=8):
    """Retain every query from the existing round-robin prime-power stream."""
    delta = -integer if integer % 4 == 3 else -4 * integer
    identity = (1, delta % 2, (delta % 2 - delta) // 4)
    height = isqrt((-delta) // 3)
    base = build_prime_form_factor_base(delta, cutoff)
    if any(1 < gcd(p, integer) < integer for p in base['excluded_primes']):
        raise ValueError('The fixture is already split by its initial factor base.')
    records = {r.prime: (r.root, r.positive) for r in base['records']}
    if integer % 4 == 1:
        records[2] = (2, reduce_form((2, 2, (integer + 1) // 2)))
    primes = sorted(records)
    if not primes or max_queries < 1 or batch_size < 1:
        raise ValueError('Require a nonempty base and positive query and batch limits.')
    initial_size, prime_product = len(primes), prod(primes)
    first = {p: least_power_above(p, height) for p in primes}
    targets = {p: power(records[p][1], first[p], identity, compose_forms_hnf) for p in primes}
    keys, key_indices, rows, witnesses, checkpoints = [], {}, [], [], []
    total_gcd_decisions = trial_divisions = 0
    inspect = lambda form: ambiguous_factor(form, integer)

    def finish(extra=None):
        return {'integer': integer, 'cutoff': cutoff, 'norm_queries': query + 1,
                'batch_size': batch_size, 'max_queries': max_queries,
                'active_primes': primes, 'prime_records': records, 'checkpoints': checkpoints,
                'witnesses': witnesses, 'keys': keys, 'old_rows': rows,
                'trial_divisions': trial_divisions,
                'total_refinement_gcd_decisions': total_gcd_decisions, **(extra or {})}

    for query in range(max_queries):
        offset, column = divmod(query, initial_size)
        prime = primes[column]
        target, exponent = targets[prime], first[prime] + offset
        direct_factor = inspect(target) if target != identity and target == inverse_form(target) else None
        leading_gcd = gcd(target[0], integer)
        if 1 < leading_gcd < integer:
            direct_factor = leading_gcd
        if direct_factor is not None:
            return finish({'direct_factor': direct_factor,
                           'terminal': {'prime': prime, 'exponent': exponent, 'target': target}})
        split = partial_leading_word(integer, target, primes, records, prime_product)
        trial_divisions += split['trial_divisions']
        assert verify_partial_identity(integer, target, split, records)
        key = split['key']
        if key is not None and key not in key_indices:
            key_indices[key] = len(keys)
            keys.append(key)
            for row in rows:
                row.append(0)
        row = [0] * (initial_size + len(keys))
        row[column] = exponent
        for p, value in split['word'].items():
            row[primes.index(p)] -= value
        if key is not None:
            row[initial_size + key_indices[key]] -= split['orientation']
        rows.append(row)
        witnesses.append({'prime': prime, 'exponent': exponent, 'target': target, 'split': split})
        targets[prime] = compose_with_prime_form(target, prime, *records[prime])
        if (query + 1) % batch_size and query + 1 != max_queries:
            continue
        refined = oriented_coprime_basis(delta, keys)
        total_gcd_decisions += refined['decision_gcd_calls']
        transported = transport_relation_rows(rows, initial_size, refined['words'])
        old = eliminate_atom_columns(rows, initial_size)
        new = eliminate_atom_columns(transported, initial_size)
        assert len(new['rows']) >= len(old['rows'])
        checkpoint = {'queries': query + 1, 'distinct_residual_keys': len(keys),
                      'atom_count': len(refined['atoms']),
                      'old_base_rank': len(old['rows']), 'new_base_rank': len(new['rows'])}
        for name, emitted in (('old', old), ('new', new)):
            top, factor, stats = first_useful_relation_top(
                [records[p][1] for p in primes], emitted['rows'], identity,
                compose_forms_hnf, inverse_form, inspect)
            checkpoint[name + '_factor'] = factor
            checkpoint[name + '_stats'] = stats
        checkpoints.append(checkpoint)
        if checkpoint['new_factor'] is not None or query + 1 == max_queries:
            return finish({'refinement': refined, 'transported_rows': transported,
                           'old_emitted': old, 'new_emitted': new})
    return finish()
