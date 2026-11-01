"""Prepared exact overlap queries with frozen bucket remainders and deltas."""
from coprime_bucket_forest import BucketOrientedBasis, CoprimeBucketForest
from oriented_coprime_basis import _key


class PreparedCoprimeBucketForest(CoprimeBucketForest):
    def __init__(self, bucket_size=512):
        super().__init__(bucket_size)
        self.prepared = None
        self.batch_statistics = {
            'prepared_batches': 0, 'prepared_candidates': 0,
            'prepared_product_multiplications': 0, 'prepared_product_operand_bits': 0,
            'prepared_remainders': 0, 'prepared_remainder_operand_bits': 0,
            'delta_multiplications': 0, 'delta_multiplication_operand_bits': 0,
            'delta_divisions': 0, 'delta_division_operand_bits': 0,
            'maximum_prepared_remainders': 0, 'maximum_batch_size': 0,
        }

    def _multiply_delta(self, left, right):
        self.batch_statistics['delta_multiplications'] += 1
        self.batch_statistics['delta_multiplication_operand_bits'] += left.bit_length() + right.bit_length()
        return left * right

    def _divide_delta(self, left, right):
        self.batch_statistics['delta_divisions'] += 1
        self.batch_statistics['delta_division_operand_bits'] += left.bit_length() + right.bit_length()
        quotient, remainder = divmod(left, right)
        assert remainder == 0
        return quotient

    def prepare_batch(self, norms):
        assert self.prepared is None
        norms = tuple(norm for norm in norms if norm > 1)
        stats = self.batch_statistics
        stats['prepared_batches'] += 1
        stats['prepared_candidates'] += len(norms)
        stats['maximum_batch_size'] = max(stats['maximum_batch_size'], len(norms))
        levels = [list(norms)] if norms else []
        while levels and len(levels[-1]) > 1:
            layer, following = levels[-1], []
            for i in range(0, len(layer), 2):
                if i + 1 == len(layer):
                    following.append(layer[i])
                else:
                    stats['prepared_product_multiplications'] += 1
                    stats['prepared_product_operand_bits'] += layer[i].bit_length() + layer[i + 1].bit_length()
                    following.append(layer[i] * layer[i + 1])
            levels.append(following)

        def remainder(value, modulus):
            if value < modulus:
                return value
            stats['prepared_remainders'] += 1
            stats['prepared_remainder_operand_bits'] += value.bit_length() + modulus.bit_length()
            return value % modulus

        residues, stored = [], 0
        for tree in self.buckets:
            value = tree.products[1]
            if value == 1 or not norms:
                residues.append(None)
                continue
            layer = [remainder(value, levels[-1][0])]
            for moduli in reversed(levels[:-1]):
                layer = [remainder(layer[i // 2], modulus) for i, modulus in enumerate(moduli)]
            residues.append(layer)
            stored += len(layer)
        stats['maximum_prepared_remainders'] = max(stats['maximum_prepared_remainders'], stored)
        old_slots = {i: norm for i, norm in enumerate(self.leaves) if norm > 1}
        self.prepared = {'norms': norms, 'index': 0, 'residues': residues,
                         'old_slots': old_slots, 'removed': [1] * len(self.buckets),
                         'added': [1] * len(self.buckets)}

    def finish_batch(self):
        assert self.prepared is not None
        assert self.prepared['index'] == len(self.prepared['norms'])
        self.prepared = None

    def insert(self, norm):
        slot = super().insert(norm)
        if self.prepared is not None:
            state = self.prepared
            while len(state['added']) < len(self.buckets):
                state['added'].append(1)
                state['removed'].append(1)
                state['residues'].append(None)
            bucket = slot // self.bucket_size
            state['added'][bucket] = self._multiply_delta(state['added'][bucket], norm)
        return slot

    def remove(self, slot):
        if self.prepared is not None:
            state = self.prepared
            bucket, local = divmod(slot, self.bucket_size)
            norm = self.buckets[bucket].leaves[local]
            if slot in state['old_slots']:
                assert state['old_slots'].pop(slot) == norm
                state['removed'][bucket] = self._multiply_delta(state['removed'][bucket], norm)
            else:
                state['added'][bucket] = self._divide_delta(state['added'][bucket], norm)
        super().remove(slot)

    def overlaps(self, norm):
        if self.prepared is None:
            return super().overlaps(norm)
        state = self.prepared
        offset = state['index']
        assert offset < len(state['norms']) and state['norms'][offset] == norm
        state['index'] += 1
        found, remaining = [], norm
        for index, tree in enumerate(self.buckets):
            if remaining == 1:
                break
            if tree.products[1] == 1:
                continue
            values = state['residues'][index]
            old_common = self._gcd(remaining, values[offset]) if values is not None else 1
            removed, added = state['removed'][index], state['added'][index]
            if old_common > 1 and removed > 1:
                old_common = self._divide_delta(old_common, self._gcd(old_common, removed))
            new_common = self._gcd(remaining, added) if added > 1 else 1
            common = self._multiply_delta(old_common, new_common) if old_common > 1 and new_common > 1 else max(old_common, new_common)
            if common > 1:
                found.extend(self._route(index, common))
                remaining //= common
                self.routing_divisions += 1
        return found


class BufferedBucketOrientedBasis:
    """Buffer supplied keys without changing their order or public indices."""
    def __init__(self, discriminant, *, batch_size=512, bucket_size=512):
        if not isinstance(batch_size, int) or batch_size < 1:
            raise ValueError('The batch size must be positive.')
        self.discriminant, self.batch_size = discriminant, batch_size
        self.engine = BucketOrientedBasis(discriminant, bucket_size=bucket_size)
        self.engine.tree = PreparedCoprimeBucketForest(bucket_size)
        self.keys, self.pending = [], []

    def add(self, norm, root):
        key = _key(self.discriminant, norm, root)
        index = len(self.keys)
        self.keys.append(key)
        self.pending.append(key)
        if len(self.pending) >= self.batch_size:
            self.flush()
        return index

    def flush(self):
        if not self.pending:
            return
        self.engine.tree.prepare_batch([norm for norm, _ in self.pending])
        start = len(self.engine.keys)
        for offset, (norm, root) in enumerate(self.pending):
            assert self.engine.add(norm, root) == start + offset
        self.engine.tree.finish_batch()
        self.pending.clear()

    def result(self, *, dense=False, verify=True):
        self.flush()
        record = self.engine.result(dense=dense, verify=verify)
        assert record['keys'] == self.keys
        record['tree_update_method'] = 'prepared coprime bucket remainders'
        record['statistics'].update(self.engine.tree.batch_statistics)
        record['statistics']['configured_batch_size'] = self.batch_size
        return record


def batched_bucket_oriented_coprime_basis(discriminant, keys, *, batch_size=512, bucket_size=512, verify=True):
    state = BufferedBucketOrientedBasis(discriminant, batch_size=batch_size, bucket_size=bucket_size)
    for norm, root in keys:
        state.add(norm, root)
    return state.result(verify=verify)
