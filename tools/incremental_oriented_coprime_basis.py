"""Incremental oriented coprime atoms with sparse usage and a product tree.

This improves exact refinement of supplied quadratic ideals. It does not
claim Bernstein's essentially-linear batch cost. Only overlapping atoms
enter each local gcd-free refinement.
"""

from heapq import heapify, heappop, heappush
from math import gcd

from binary_quadratic_forms import reduce_form
from oriented_coprime_basis import _form, _ideal_power, _key
from quadratic_ideal_hnf import IdealBasis, form_to_ideal, multiply_ideals


class _ProductTree:
    """Stable leaf slots, unit-valued vacancies, and charged point updates."""

    def __init__(self):
        self.capacity = 1
        self.leaves = [1]
        self.products = [1, 1]
        self.free = [0]
        self.gcd_calls = self.gcd_operand_bits = 0
        self.multiplications = self.multiplication_operand_bits = 0
        self.leaf_updates = self.capacity_doublings = 0

    def _multiply(self, left, right):
        self.multiplications += 1
        self.multiplication_operand_bits += left.bit_length() + right.bit_length()
        return left * right

    def _grow(self):
        old = self.capacity
        self.capacity *= 2
        self.leaves.extend([1] * old)
        self.free.extend(range(old, self.capacity))
        heapify(self.free)
        self.products = [1] * (2 * self.capacity)
        self.products[self.capacity:] = self.leaves
        for node in range(self.capacity - 1, 0, -1):
            self.products[node] = self._multiply(
                self.products[2 * node], self.products[2 * node + 1])
        self.capacity_doublings += 1

    def _set(self, slot, norm):
        self.leaves[slot] = norm
        node = self.capacity + slot
        self.products[node] = norm
        self.leaf_updates += 1
        while node > 1:
            node //= 2
            self.products[node] = self._multiply(
                self.products[2 * node], self.products[2 * node + 1])

    def insert(self, norm):
        assert norm > 1
        if not self.free:
            self._grow()
        slot = heappop(self.free)
        self._set(slot, norm)
        return slot

    def remove(self, slot):
        assert self.leaves[slot] > 1
        self._set(slot, 1)
        heappush(self.free, slot)

    def overlaps(self, norm):
        """Return precisely the leaves with nontrivial gcd with norm."""
        found, stack = [], [(1, norm)]
        while stack:
            node, value = stack.pop()
            product = self.products[node]
            if product == 1:
                continue
            self.gcd_calls += 1
            self.gcd_operand_bits += value.bit_length() + product.bit_length()
            common = gcd(value, product)
            if common == 1:
                continue
            if node >= self.capacity:
                found.append(node - self.capacity)
            else:
                stack.append((2 * node + 1, common))
                stack.append((2 * node, common))
        return found


def _refine_overlaps(discriminant, keys):
    """The established gcd splits with sparse local exponent dictionaries."""
    pieces = [(norm, {index: 1}) for index, (norm, _) in enumerate(keys) if norm > 1]
    gcd_calls = norm_splits = orientation_splits = 0
    bit_bound = sum(norm.bit_length() for norm, _ in keys)
    sort_key = lambda piece: (piece[0], tuple(sorted(piece[1].items())))
    while True:
        selected = None
        for i, (a, _) in enumerate(pieces):
            for j in range(i + 1, len(pieces)):
                common = gcd(a, pieces[j][0])
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
        merged = dict(u)
        for index, exponent in v.items():
            merged[index] = merged.get(index, 0) + exponent
        pieces = [piece for k, piece in enumerate(pieces) if k not in (i, j)]
        pieces.extend((norm, word) for norm, word in
                      ((common, merged), (a // common, u), (b // common, v)) if norm > 1)
        pieces.sort(key=sort_key)
        norm_splits += 1
        assert norm_splits <= bit_bound
    while True:
        selected = None
        for column, (norm, word) in enumerate(pieces):
            root = keys[min(word)][1] % (2 * norm)
            root = min(root, (-root) % (2 * norm))
            for index in sorted(word):
                plus = gcd(norm, (keys[index][1] - root) // 2)
                minus = gcd(norm, (keys[index][1] + root) // 2)
                gcd_calls += 2
                assert plus * minus == norm and gcd(plus, minus) == 1
                if 1 < plus < norm:
                    selected = column, plus, minus
                    break
            if selected is not None:
                break
        if selected is None:
            break
        column, plus, minus = selected
        _, word = pieces.pop(column)
        pieces.extend(((plus, word), (minus, word)))
        pieces.sort(key=sort_key)
        orientation_splits += 1
        assert orientation_splits <= bit_bound
    atoms = []
    for norm, word in sorted(pieces, key=sort_key):
        root = keys[min(word)][1] % (2 * norm)
        root = min(root, (-root) % (2 * norm))
        signed = {index: exponent if keys[index][1] % (2 * norm) == root else -exponent
                  for index, exponent in word.items()}
        atoms.append({'norm': norm, 'root': root,
                      'form': reduce_form(_form(discriminant, norm, root)),
                      'local_word': signed})
    return {'atoms': atoms, 'decision_gcd_calls': gcd_calls,
            'norm_splits': norm_splits, 'orientation_splits': orientation_splits}


class IncrementalOrientedBasis:
    """Maintain exact signed atom words without a dense key-by-atom matrix."""

    def __init__(self, discriminant):
        if (not isinstance(discriminant, int) or discriminant >= 0
                or discriminant % 4 not in (0, 1)):
            raise ValueError('Require a negative quadratic-order discriminant.')
        self.discriminant = discriminant
        self.keys = []
        self.atoms = {}
        self.tree = _ProductTree()
        self.input_bit_bound = self.local_gcd_calls = 0
        self.norm_splits = self.orientation_splits = self.local_refinements = 0
        self.maximum_local_inputs = self.maximum_local_norm_bits = 0
        self.maximum_atoms = self.usage_entries = self.usage_writes = 0
        self.maximum_usage_entries = 0

    def add(self, norm, root):
        """Add one supplied ideal key; return its immutable input index."""
        norm, root = _key(self.discriminant, norm, root)
        index = len(self.keys)
        slots = sorted(self.tree.overlaps(norm),
                       key=lambda slot: self.atoms[slot]['norm']) if norm > 1 else []
        self.keys.append((norm, root))
        self.input_bit_bound += norm.bit_length()
        if norm == 1:
            return index
        if not slots:
            canonical = min(root, (-root) % (2 * norm))
            atom = {'norm': norm, 'root': canonical,
                    'form': reduce_form(_form(self.discriminant, norm, canonical)),
                    'usage': {index: 1 if canonical == root else -1}}
            self.atoms[self.tree.insert(norm)] = atom
            self.usage_entries += 1
            self.usage_writes += 1
        else:
            old_atoms = [self.atoms[slot] for slot in slots]
            local_keys = [(a['norm'], a['root']) for a in old_atoms] + [(norm, root)]
            refined = _refine_overlaps(self.discriminant, local_keys)
            self.local_refinements += 1
            self.maximum_local_inputs = max(self.maximum_local_inputs, len(local_keys))
            self.maximum_local_norm_bits = max(
                self.maximum_local_norm_bits, sum(a.bit_length() for a, _ in local_keys))
            self.local_gcd_calls += refined['decision_gcd_calls']
            self.norm_splits += refined['norm_splits']
            self.orientation_splits += refined['orientation_splits']
            old_by_key = {(a['norm'], a['root']): (i, slot)
                          for i, (slot, a) in enumerate(zip(slots, old_atoms))}
            retained, replacements = set(), []
            for new in refined['atoms']:
                old_coefficients = sorted((i, exponent) for i, exponent in new['local_word'].items()
                                          if i < len(old_atoms))
                reuse = old_by_key.get((new['norm'], new['root']))
                coefficient = new['local_word'].get(len(old_atoms), 0)
                if reuse is not None and old_coefficients == [(reuse[0], 1)]:
                    # A repeated use of an unchanged atom only appends one entry.
                    slot = reuse[1]
                    retained.add(slot)
                    assert coefficient and index not in self.atoms[slot]['usage']
                    self.atoms[slot]['usage'][index] = coefficient
                    self.usage_entries += 1
                    self.usage_writes += 1
                    continue
                usage = {}
                for old_index, multiplier in old_coefficients:
                    for key_index, exponent in old_atoms[old_index]['usage'].items():
                        assert key_index not in usage
                        usage[key_index] = multiplier * exponent
                        self.usage_writes += 1
                if coefficient:
                    usage[index] = coefficient
                    self.usage_writes += 1
                assert usage and all(usage.values())
                replacements.append({key: new[key] for key in ('norm', 'root', 'form')}
                                    | {'usage': usage})
            for slot in slots:
                if slot not in retained:
                    self.usage_entries -= len(self.atoms[slot]['usage'])
                    self.tree.remove(slot)
                    del self.atoms[slot]
            for atom in replacements:
                self.atoms[self.tree.insert(atom['norm'])] = atom
                self.usage_entries += len(atom['usage'])
        self.maximum_atoms = max(self.maximum_atoms, len(self.atoms))
        self.maximum_usage_entries = max(self.maximum_usage_entries, self.usage_entries)
        assert self.usage_entries <= self.input_bit_bound
        return index

    def result(self, *, dense=False, verify=True):
        """Export sparse words; request dense words only for legacy interfaces."""
        ordered = sorted(self.atoms.values(), key=lambda atom: atom['norm'])
        atoms = [{key: atom[key] for key in ('norm', 'root', 'form')} for atom in ordered]
        words = [[] for _ in self.keys]
        for column, atom in enumerate(ordered):
            for index, exponent in atom['usage'].items():
                words[index].append((column, exponent))
        tree = self.tree
        statistics = {
            'tree_gcd_calls': tree.gcd_calls,
            'tree_gcd_operand_bits': tree.gcd_operand_bits,
            'local_gcd_calls': self.local_gcd_calls,
            'local_refinements': self.local_refinements,
            'maximum_local_inputs': self.maximum_local_inputs,
            'maximum_local_norm_bits': self.maximum_local_norm_bits,
            'tree_multiplications': tree.multiplications,
            'tree_multiplication_operand_bits': tree.multiplication_operand_bits,
            'tree_leaf_updates': tree.leaf_updates,
            'tree_capacity_doublings': tree.capacity_doublings,
            'tree_capacity': tree.capacity,
            'tree_product_bits': tree.products[1].bit_length(),
            'maximum_atoms': self.maximum_atoms,
            'usage_entries': self.usage_entries,
            'usage_writes': self.usage_writes,
            'maximum_usage_entries': self.maximum_usage_entries,
        }
        record = {'discriminant': self.discriminant, 'keys': list(self.keys),
                  'atoms': atoms, 'sparse_words': words,
                  'decision_gcd_calls': tree.gcd_calls + self.local_gcd_calls,
                  'norm_splits': self.norm_splits, 'orientation_splits': self.orientation_splits,
                  'input_bit_bound': self.input_bit_bound, 'statistics': statistics,
                  'factoring_oracle_calls': 0, 'residual_primality_tests': 0}
        if dense:
            record['words'] = []
            for word in words:
                row = [0] * len(atoms)
                for column, exponent in word:
                    row[column] = exponent
                record['words'].append(row)
        if verify:
            record['verification'] = verify_incremental_oriented_basis(record)
        return record


def incremental_oriented_coprime_basis(discriminant, keys, *, dense=False, verify=True):
    """Process keys in the supplied order and return the final exact basis."""
    state = IncrementalOrientedBasis(discriminant)
    for norm, root in keys:
        state.add(norm, root)
    return state.result(dense=dense, verify=verify)


def verify_incremental_oriented_basis(record):
    """Replay sparse norm powers, roots, and literal integral HNF products."""
    discriminant = record['discriminant']
    if (not isinstance(discriminant, int) or discriminant >= 0
            or discriminant % 4 not in (0, 1)):
        raise ValueError('Require a negative quadratic-order discriminant.')
    keys = [_key(discriminant, *key) for key in record['keys']]
    atoms, words = record['atoms'], record['sparse_words']
    assert len(words) == len(keys)
    product = 1
    for atom in atoms:
        norm, root = _key(discriminant, atom['norm'], atom['root'])
        assert norm > 1 and gcd(norm, product) == 1
        assert tuple(atom['form']) == reduce_form(_form(discriminant, norm, root))
        product *= norm
    entries = 0
    for (norm, root), word in zip(keys, words):
        reconstructed, literal = 1, IdealBasis(1, 0, 1)
        previous = -1
        for column, exponent in word:
            assert isinstance(column, int) and previous < column < len(atoms)
            assert isinstance(exponent, int) and 0 < abs(exponent) <= norm.bit_length()
            previous = column
            atom = atoms[column]
            oriented_root = atom['root'] if exponent > 0 else -atom['root']
            assert (root - oriented_root) % (2 * atom['norm']) == 0
            reconstructed *= atom['norm'] ** abs(exponent)
            ideal = form_to_ideal(_form(discriminant, atom['norm'], oriented_root))
            literal = multiply_ideals(literal, _ideal_power(ideal, abs(exponent), discriminant),
                                     discriminant)
            entries += 1
        assert reconstructed == norm
        assert literal == form_to_ideal(_form(discriminant, norm, root))
    assert entries <= sum(norm.bit_length() for norm, _ in keys)
    if 'words' in record:
        assert len(record['words']) == len(words)
        for sparse, dense in zip(words, record['words']):
            assert len(dense) == len(atoms)
            assert list(map(tuple, sparse)) == [(i, value) for i, value in enumerate(dense) if value]
    return {'literal_ideal_identities': len(keys), 'pairwise_coprime_atoms': len(atoms),
            'sparse_usage_entries': entries}
