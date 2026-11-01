"""Build an adaptive class-group presentation from reduced leading norms.

Every queried norm has at most half the input bits, up to a constant.
The finite-control factorizer is exact, but its running time and the
number of dynamically introduced primes have no subexponential bound.
"""

from collections import deque
from math import gcd, isqrt, prod

from arithmetic_two_group import power
from binary_quadratic_forms import ambiguous_factor, inverse_form, reduce_form
from certified_small_factor import CertifiedFactorizer
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top
from weighted_norm_relations import least_power_above


def _add(destination, source, scale):
    for prime, exponent in source.items():
        destination[prime] = destination.get(prime, 0) + scale * exponent
        if destination[prime] == 0:
            del destination[prime]


def leading_norm_word(integer, form, factors, records):
    """Decode a factored reduced leading coefficient into prime-form classes.

    ``records`` maps known primes to ``(root, reduced_form)``. The input
    factors must be distinct known primes, with positive exponents.
    The caller establishes primality. The product and form are checked.
    A proper gcd of the leading coefficient and integer must be handled
    before this call. The returned new records do not mutate the input.
    """
    delta = -integer if integer % 4 == 3 else -4 * integer
    leading, middle, trailing = form
    if (leading < 1 or gcd(gcd(leading, middle), trailing) != 1
            or middle * middle - 4 * leading * trailing != delta
            or form != reduce_form(form) or gcd(leading, integer) != 1):
        raise ValueError('A reduced primitive form with coprime leading norm is required.')
    factors = tuple(factors)
    if (any(not isinstance(p, int) or not isinstance(e, int) or p < 2 or e < 1
            for p, e in factors)
            or len({p for p, _ in factors}) != len(factors)
            or prod(p ** e for p, e in factors) != leading):
        raise ValueError('The supplied factorization does not equal the leading norm.')
    terms, new = {}, {}
    for prime, exponent in factors:
        if prime == 2 and integer % 4 == 1:
            if (exponent != 1 or leading % 4 != 2 or middle % 4 != 2
                    or trailing % 2 != 1 or 2 not in records):
                raise ValueError('The ramified-two leading-norm invariant failed.')
            terms[2] = 1
            continue
        if gcd(prime, delta) != 1:
            raise ValueError('A ramified odd prime should have supplied a proper gcd.')
        residue = middle % (2 * prime)
        if prime in records:
            root = records[prime][0]
        else:
            root = min(residue, 2 * prime - residue)
            if (root * root - delta) % (4 * prime):
                raise ValueError('The leading norm does not supply a prime-form root.')
            new[prime] = (root, reduce_form(
                (prime, root, (root * root - delta) // (4 * prime))))
        if residue not in (root, (-root) % (2 * prime)):
            raise ValueError('The stored prime root has an incompatible orientation.')
        terms[prime] = exponent if residue == root else -exponent
    return terms, new


def factor_with_dynamic_leading_norms(integer, cutoff, max_steps, *,
                                     lookahead=8, factor_norm=None):
    """Run a bounded deterministic adaptive presentation construction.

    A supplied ``factor_norm`` must return complete prime factorizations.
    The default returns recursive Lucas certificates for every prime used.
    All norm queries are recorded so a later replay can use certificates
    alone. A missing factor reports only this bounded attempt's outcome.
    It is not a primality certificate or a fast closure theorem.
    """
    if not isinstance(integer, int) or integer < 3 or integer % 2 == 0:
        raise ValueError('An odd input integer of at least three is required.')
    if (not isinstance(cutoff, int) or cutoff < 2
            or not isinstance(max_steps, int) or max_steps < 0
            or not isinstance(lookahead, int) or lookahead < 1):
        raise ValueError('Require cutoff >= 2, max_steps >= 0, and lookahead >= 1.')
    delta = -integer if integer % 4 == 3 else -4 * integer
    parity = delta % 2
    identity = (1, parity, (parity - delta) // 4)
    height = isqrt((-delta) // 3)
    certifier = CertifiedFactorizer() if factor_norm is None else None
    factor_norm = certifier.factor if certifier is not None else factor_norm
    base = build_prime_form_factor_base(delta, cutoff)
    records = {r.prime: (r.root, r.positive) for r in base['records']}
    if integer % 4 == 1:
        records[2] = (2, reduce_form((2, 2, (integer + 1) // 2)))
    active = sorted(records)
    pending = deque(active)
    expressions = {p: {p: 1} for p in active}
    rows, history, queries = [], [], {}
    calls = units = maximum_norm_bits = 0
    last_extraction = None
    initial_count = len(active)

    def word_value(word):
        value = identity
        for prime, exponent in word.items():
            generator = records[prime][1]
            if exponent < 0:
                generator, exponent = inverse_form(generator), -exponent
            value = compose_forms_hnf(value, power(
                generator, exponent, identity, compose_forms_hnf))
        return value

    def finish(factor, method, terminal=None):
        result = {
            'integer': integer, 'input_bits': integer.bit_length(),
            'discriminant': delta, 'height_bound': height,
            'cutoff': cutoff, 'max_steps': max_steps, 'lookahead': lookahead,
            'factor': factor, 'method': method, 'steps': len(history),
            'initial_generators': initial_count, 'known_prime_forms': len(records),
            'active_generators': len(active), 'relation_rank': len(rows),
            'maximum_active_prime': max(active, default=0),
            'new_active_prime_bound': isqrt(height),
            'maximum_new_active_prime': max(active[initial_count:], default=0),
            'deficiency': len(active) - len(rows), 'pending_primes': list(pending),
            'unit_eliminations': units, 'factor_calls': calls,
            'maximum_norm_bits': maximum_norm_bits,
            'history': history, 'norm_factorizations': queries,
            'prime_records': records, 'expressions': expressions,
            'active_primes': active, 'relation_rows': rows,
            'last_extraction': last_extraction, 'terminal_witness': terminal,
        }
        if certifier is not None:
            result['prime_certificates'] = certifier.certificates
            result['factoring_stats'] = certifier.stats
        return result

    for prime in base['excluded_primes']:
        divisor = gcd(prime, integer)
        if 1 < divisor < integer:
            return finish(divisor, 'factor-base gcd', {'prime': prime})

    for _ in range(max_steps):
        if not pending:
            break
        pivot = pending.popleft()
        start = least_power_above(pivot, height)
        candidates = []
        for exponent in range(start, start + lookahead):
            target = power(records[pivot][1], exponent, identity, compose_forms_hnf)
            leading = target[0]
            terminal = {'pivot_prime': pivot, 'exponent': exponent, 'target': target}
            if target == inverse_form(target) and target != identity:
                divisor = ambiguous_factor(target, integer)
                if divisor is not None:
                    return finish(divisor, 'target involution', terminal)
            divisor = gcd(leading, integer)
            if 1 < divisor < integer:
                return finish(divisor, 'leading-norm gcd', terminal)
            if divisor != 1:
                raise ValueError('The reduced leading coefficient unexpectedly contains the input.')
            if not 1 <= leading <= height:
                raise ValueError('The reduced leading-norm height bound failed.')
            parts = tuple((int(p), int(e)) for p, e in factor_norm(leading))
            terms, new = leading_norm_word(integer, target, parts, records)
            if leading in queries and queries[leading] != parts:
                raise ValueError('The factoring callback changed a previous factorization.')
            queries[leading] = parts
            calls += 1
            maximum_norm_bits = max(maximum_norm_bits, leading.bit_length())
            unit_candidates = [p for p in new if abs(terms[p]) == 1]
            score = (len(new) - bool(unit_candidates), len(new),
                     sum(p.bit_length() for p in new), leading, exponent)
            candidates.append((score, exponent, target, parts, terms, new, unit_candidates))
        _, exponent, target, parts, terms, new, unit_candidates = min(
            candidates, key=lambda item: item[0])
        eliminated = max(unit_candidates) if unit_candidates else None
        for prime in sorted(new):
            records[prime] = new[prime]
            pending.append(prime)
            if prime != eliminated:
                if prime * prime > target[0]:
                    raise ValueError('The fresh-unit active-prime bound failed.')
                active.append(prime)
                expressions[prime] = {prime: 1}
        if word_value(terms) != target:
            raise ValueError('A decoded leading norm did not equal its target class.')
        original = {p: -e for p, e in terms.items()}
        original[pivot] = original.get(pivot, 0) + exponent
        if (word_value(original) != identity
                or prod(p ** abs(e) for p, e in terms.items()) != target[0]
                or pivot ** exponent <= height):
            raise ValueError('A raw leading-norm relation failed exact verification.')
        row = {}
        _add(row, expressions[pivot], exponent)
        for prime, value in terms.items():
            if prime != eliminated:
                _add(row, expressions[prime], -value)
        stored = eliminated is None
        if eliminated is not None:
            expressions[eliminated] = {
                p: value * terms[eliminated] for p, value in row.items()}
            units += 1
        else:
            # Original rows have distinct weighted-dominant pivot columns.
            # Eliminating fresh unit columns preserves their independence.
            if not row:
                raise ValueError('An independent substituted relation became zero.')
            rows.append(row)
        history.append({
            'step': len(history) + 1, 'pivot_prime': pivot, 'exponent': exponent,
            'target': target, 'norm': target[0], 'factors': parts, 'word': terms,
            'new_primes': sorted(new), 'eliminated_prime': eliminated,
            'substituted_row': row, 'stored_relation': stored,
            'known_prime_forms': len(records), 'active_generators': len(active),
            'relation_rank': len(rows),
        })
        if stored:
            matrix = [[row.get(prime, 0) for prime in active] for row in rows]
            _, divisor, last_extraction = first_useful_relation_top(
                [records[p][1] for p in active], matrix, identity,
                compose_forms_hnf, inverse_form,
                lambda form: ambiguous_factor(form, integer), verify_relations=True)
            if divisor is not None:
                return finish(divisor, 'unit-eliminated relation halving')
    return finish(None, 'dynamic budget exhausted' if pending else 'closed relation set')
