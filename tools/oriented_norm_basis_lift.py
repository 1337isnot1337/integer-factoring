"""Attach exact quadratic orientations to a supplied coprime norm basis.

The norm-basis algorithm is a separate input-producing component. The
associated theorem may use Bernstein's classical fast construction; this
module does not label a supplied factorization as a free factoring oracle.
"""
from math import gcd

from binary_quadratic_forms import reduce_form
from oriented_coprime_basis import _form, _ideal_power, _key
from quadratic_ideal_hnf import IdealBasis, form_to_ideal, multiply_ideals


def _product_tree(values):
    levels = [list(values)]
    while len(levels[-1]) > 1:
        old = levels[-1]
        levels.append([old[i] * old[i+1] if i+1 < len(old) else old[i]
                       for i in range(0, len(old), 2)])
    return levels


def _remainders(value, levels):
    if not levels[0]:
        return []
    current = [value % levels[-1][0]]
    for moduli in reversed(levels[:-1]):
        current = [current[i // 2] % modulus for i, modulus in enumerate(moduli)]
    return current


def verify_coprime_norms(norms):
    """Batch-check all norm pairs without repeatedly scanning a huge product."""
    norms = list(norms)
    if any(not isinstance(a, int) or a <= 1 for a in norms):
        raise ValueError('Every coprime basis norm must be an integer greater than one.')
    if not norms:
        return {'norms': 0, 'product_bits': 0, 'batch_gcds': 0}
    product = _product_tree(norms)[-1][0]
    residues = _remainders(product, _product_tree([a*a for a in norms]))
    for a, remainder in zip(norms, residues):
        quotient, rest = divmod(remainder, a)
        assert rest == 0
        if gcd(a, quotient) != 1:
            raise ValueError('The supplied norm basis is not pairwise coprime.')
    return {'norms': len(norms), 'product_bits': product.bit_length(),
            'batch_gcds': len(norms)}


def _validated_factorization(keys, norms, words):
    if len(words) != len(keys):
        raise ValueError('Supply one sparse norm word for every input key.')
    incidence = [[] for _ in norms]
    for index, ((norm, _), word) in enumerate(zip(keys, words)):
        product, previous = 1, -1
        for entry in word:
            if not isinstance(entry, (tuple, list)) or len(entry) != 2:
                raise ValueError('A sparse norm entry must be a pair.')
            column, exponent = entry
            if (not isinstance(column, int) or not previous < column < len(norms)
                    or not isinstance(exponent, int) or not 0 < exponent <= norm.bit_length()):
                raise ValueError('Require ordered columns and positive bounded exponents.')
            previous = column
            factor = norms[column]
            if exponent * (factor.bit_length()-1) > norm.bit_length():
                raise ValueError('A claimed norm power exceeds its input norm.')
            product *= factor**exponent
            if product > norm:
                raise ValueError('A claimed norm product exceeds its input norm.')
            incidence[column].append((index, exponent))
        if product != norm:
            raise ValueError('A supplied norm word does not reconstruct its input norm.')
    if any(not rows for rows in incidence):
        raise ValueError('Discard unused norm-basis entries before orientation refinement.')
    return incidence


def verify_lifted_oriented_basis(record):
    """Check sparse words by independent literal integral-ideal arithmetic."""
    delta = record['discriminant']
    if not isinstance(delta, int) or delta >= 0 or delta % 4 not in (0, 1):
        raise ValueError('Require a negative quadratic-order discriminant.')
    keys = [_key(delta, *key) for key in record['keys']]
    atoms, words = record['atoms'], record['sparse_words']
    assert len(words) == len(keys)
    coprime = verify_coprime_norms([atom['norm'] for atom in atoms])
    for atom in atoms:
        norm, root = _key(delta, atom['norm'], atom['root'])
        assert norm > 1 and root == min(root, (-root) % (2*norm))
        assert tuple(atom['form']) == reduce_form(_form(delta, norm, root))
    entries = 0
    for (norm, root), word in zip(keys, words):
        product, literal, previous = 1, IdealBasis(1, 0, 1), -1
        for column, exponent in word:
            assert isinstance(column, int) and previous < column < len(atoms)
            assert isinstance(exponent, int) and 0 < abs(exponent) <= norm.bit_length()
            previous = column
            atom = atoms[column]
            oriented = atom['root'] if exponent > 0 else -atom['root']
            assert (root-oriented) % (2*atom['norm']) == 0
            product *= atom['norm']**abs(exponent)
            assert product <= norm
            ideal = form_to_ideal(_form(delta, atom['norm'], oriented))
            literal = multiply_ideals(literal, _ideal_power(ideal, abs(exponent), delta), delta)
            entries += 1
        assert product == norm
        assert literal == form_to_ideal(_form(delta, norm, root))
    return {'literal_ideal_identities': len(keys), 'pairwise_coprime_atoms': len(atoms),
            'sparse_usage_entries': entries, 'coprime_check': coprime}


def lift_coprime_norm_basis(discriminant, keys, norms, norm_words, *, verify=True):
    """Refine orientations of a validated, completely supplied norm factorization."""
    if (not isinstance(discriminant, int) or discriminant >= 0
            or discriminant % 4 not in (0, 1)):
        raise ValueError('Require a negative quadratic-order discriminant.')
    keys = [_key(discriminant, *key) for key in keys]
    norms = list(norms)
    coprime = verify_coprime_norms(norms)
    incidence = _validated_factorization(keys, norms, norm_words)
    pending, gcd_calls, splits, piece_visits = [], 0, 0, 0
    for norm, rows in zip(norms, incidence):
        reference = keys[rows[0][0]][1] % (2*norm)
        pieces = [norm]
        for index, _ in rows:
            difference = (keys[index][1]-reference) // 2
            following = []
            for piece in pieces:
                piece_visits += 1
                plus = gcd(piece, difference)
                minus = piece // plus
                gcd_calls += 1
                if plus == 1 or minus == 1:
                    following.append(piece)
                else:
                    # Odd unramified roots imply that these are coprime.
                    following.extend((plus, minus))
                    splits += 1
            pieces = following
        for piece in pieces:
            root = reference % (2*piece)
            root = min(root, (-root) % (2*piece))
            signed = []
            for index, exponent in rows:
                residue = keys[index][1] % (2*piece)
                if residue == root:
                    signed.append((index, exponent))
                else:
                    assert residue == (-root) % (2*piece)
                    signed.append((index, -exponent))
            pending.append((piece, root, signed))
    pending.sort(key=lambda item: item[0])
    atoms, words = [], [[] for _ in keys]
    for column, (norm, root, signed) in enumerate(pending):
        atoms.append({'norm': norm, 'root': root,
                      'form': reduce_form(_form(discriminant, norm, root))})
        for index, exponent in signed:
            words[index].append((column, exponent))
    input_bits = sum(norm.bit_length() for norm, _ in keys)
    usage = sum(len(word) for word in words)
    assert piece_visits == gcd_calls <= input_bits
    assert usage <= input_bits
    result = {'discriminant': discriminant, 'keys': keys, 'atoms': atoms,
              'sparse_words': words, 'input_bit_bound': input_bits,
              'orientation_splits': splits,
              'statistics': {'norm_basis_size': len(norms),
                             'norm_basis_incidences': sum(map(len, incidence)),
                             'orientation_gcd_calls': gcd_calls,
                             'sparse_usage_entries': usage,
                             'input_coprime_check': coprime},
              'scope': 'Orient a supplied validated norm factorization; norm construction is separate.',
              'factoring_oracle_calls': 0, 'residual_primality_tests': 0}
    if verify:
        result['verification'] = verify_lifted_oriented_basis(result)
    return result
