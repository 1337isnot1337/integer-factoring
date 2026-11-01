"""Refine supplied odd unramified quadratic ideals without factoring norms.

The output is a pairwise coprime integer basis with consistent root
orientations. Signed atom words describe ideal classes. Literal integral
ideal identities use conjugate atoms for negative signs and are checked
by independent two-dimensional HNF multiplication.
"""

from math import gcd

from binary_quadratic_forms import reduce_form
from quadratic_ideal_hnf import IdealBasis, form_to_ideal, multiply_ideals


def _key(discriminant, norm, root):
    if (not isinstance(norm, int) or norm < 1 or norm % 2 == 0
            or gcd(norm, discriminant) != 1 or not isinstance(root, int)
            or (root * root - discriminant) % (4 * norm)):
        raise ValueError('Require an odd unramified norm and a valid quadratic root.')
    return norm, root % (2 * norm)


def _form(discriminant, norm, root):
    return norm, root, (root * root - discriminant) // (4 * norm)


def _ideal_power(ideal, exponent, discriminant):
    result = IdealBasis(1, 0, 1)
    while exponent:
        if exponent & 1:
            result = multiply_ideals(result, ideal, discriminant)
        exponent >>= 1
        if exponent:
            ideal = multiply_ideals(ideal, ideal, discriminant)
    return result


def verify_oriented_basis(discriminant, keys, record):
    """Replay norm powers, orientations, and literal integral ideal identities."""
    keys = [_key(discriminant, *key) for key in keys]
    assert record['discriminant'] == discriminant
    assert [_key(discriminant, *key) for key in record['keys']] == keys
    atoms, words = record['atoms'], record['words']
    assert len(words) == len(keys)
    for index, atom in enumerate(atoms):
        norm, root = _key(discriminant, atom['norm'], atom['root'])
        assert norm > 1 and tuple(atom['form']) == reduce_form(_form(discriminant, norm, root))
        assert all(gcd(norm, other['norm']) == 1 for other in atoms[:index])
    for (norm, root), word in zip(keys, words):
        assert len(word) == len(atoms)
        reconstructed_norm = 1
        ideal = IdealBasis(1, 0, 1)
        for atom, exponent in zip(atoms, word):
            if not exponent:
                continue
            d, c = atom['norm'], atom['root']
            assert isinstance(exponent, int) and abs(exponent) <= norm.bit_length()
            oriented_root = c if exponent > 0 else -c
            assert (root - oriented_root) % (2 * d) == 0
            reconstructed_norm *= d ** abs(exponent)
            base = form_to_ideal(_form(discriminant, d, oriented_root))
            ideal = multiply_ideals(ideal, _ideal_power(base, abs(exponent), discriminant),
                                    discriminant)
        assert reconstructed_norm == norm
        assert ideal == form_to_ideal(_form(discriminant, norm, root))
    return {'literal_ideal_identities': len(keys), 'pairwise_coprime_atoms': len(atoms)}


def oriented_coprime_basis(discriminant, keys, *, verify=True):
    """Return gcd-free atoms and one exact signed atom word per supplied key."""
    if (not isinstance(discriminant, int) or discriminant >= 0
            or discriminant % 4 not in (0, 1)):
        raise ValueError('Require a negative quadratic-order discriminant.')
    keys = [_key(discriminant, *key) for key in keys]
    size = len(keys)
    pieces = [(norm, [int(j == i) for j in range(size)])
              for i, (norm, _) in enumerate(keys) if norm > 1]
    gcd_calls = norm_splits = orientation_splits = 0
    input_bit_bound = sum(norm.bit_length() for norm, _ in keys)
    while True:
        selected = None
        for i, (a, u) in enumerate(pieces):
            for j in range(i + 1, len(pieces)):
                b, v = pieces[j]
                common = gcd(a, b)
                gcd_calls += 1
                if common > 1:
                    selected = i, j, common
                    break
            if selected is not None:
                break
        if selected is None:
            break
        i, j, common = selected
        a, u = pieces[i]
        b, v = pieces[j]
        pieces = [piece for k, piece in enumerate(pieces) if k not in (i, j)]
        for norm, vector in ((common, [x + y for x, y in zip(u, v)]),
                             (a // common, u), (b // common, v)):
            if norm > 1:
                pieces.append((norm, vector))
        pieces.sort(key=lambda piece: (piece[0], piece[1]))
        norm_splits += 1
        assert norm_splits <= input_bit_bound

    while True:
        selected = None
        for index, (norm, vector) in enumerate(pieces):
            first = next(i for i, value in enumerate(vector) if value)
            c = keys[first][1] % (2 * norm)
            c = min(c, (-c) % (2 * norm))
            for i, multiplicity in enumerate(vector):
                if not multiplicity:
                    continue
                plus = gcd(norm, (keys[i][1] - c) // 2)
                minus = gcd(norm, (keys[i][1] + c) // 2)
                gcd_calls += 2
                assert plus * minus == norm and gcd(plus, minus) == 1
                if 1 < plus < norm:
                    selected = index, plus, minus
                    break
            if selected is not None:
                break
        if selected is None:
            break
        index, plus, minus = selected
        _, vector = pieces.pop(index)
        pieces.extend(((plus, vector), (minus, vector)))
        pieces.sort(key=lambda piece: (piece[0], piece[1]))
        orientation_splits += 1
        assert orientation_splits <= input_bit_bound

    atoms, words = [], [[] for _ in keys]
    for norm, vector in sorted(pieces):
        first = next(i for i, value in enumerate(vector) if value)
        root = keys[first][1] % (2 * norm)
        root = min(root, (-root) % (2 * norm))
        atoms.append({'norm': norm, 'root': root,
                      'form': reduce_form(_form(discriminant, norm, root))})
        for i, multiplicity in enumerate(vector):
            sign = 1 if keys[i][1] % (2 * norm) == root else -1
            if multiplicity:
                assert (keys[i][1] - sign * root) % (2 * norm) == 0
            words[i].append(sign * multiplicity)
    result = {'discriminant': discriminant, 'keys': keys,
              'atoms': atoms, 'words': words, 'decision_gcd_calls': gcd_calls,
              'norm_splits': norm_splits, 'orientation_splits': orientation_splits,
              'input_bit_bound': input_bit_bound,
              'factoring_oracle_calls': 0, 'residual_primality_tests': 0}
    if verify:
        result['verification'] = verify_oriented_basis(discriminant, keys, result)
    return result


def transport_relation_rows(rows, initial_size, words):
    """Substitute signed residual words, retaining all transported rows."""
    atom_count = len(words[0]) if words else 0
    if (not isinstance(initial_size, int) or initial_size < 0
            or any(len(word) != atom_count or any(not isinstance(x, int) for x in word)
                   for word in words)):
        raise ValueError('Require a nonnegative initial size and equal word dimensions.')
    output = []
    for row in rows:
        if len(row) != initial_size + len(words) or any(not isinstance(x, int) for x in row):
            raise ValueError('Require one integer entry per original generator.')
        new = list(row[:initial_size]) + [0] * atom_count
        for coefficient, word in zip(row[initial_size:], words):
            for index, exponent in enumerate(word):
                new[initial_size + index] += coefficient * exponent
        output.append(new)
    return output
