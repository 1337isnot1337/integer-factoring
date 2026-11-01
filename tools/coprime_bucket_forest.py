"""Local product trees for pairwise-coprime oriented atoms."""
from heapq import heapify, heappop, heappush
from math import gcd
from incremental_oriented_coprime_basis import IncrementalOrientedBasis
from quotient_update_oriented_basis import _QuotientUpdateTree


class CoprimeBucketForest:
    def __init__(self, bucket_size=128):
        if not isinstance(bucket_size, int) or bucket_size < 2 or bucket_size & (bucket_size - 1):
            raise ValueError('The bucket size must be a power of two at least two.')
        self.bucket_size = bucket_size
        self.capacity = bucket_size
        self.buckets = [self._new_bucket()]
        self.free = list(range(bucket_size))
        self.capacity_doublings = 0
        self.gcd_calls = self.gcd_operand_bits = self.routing_divisions = 0
        self._export_multiplications = self._export_operand_bits = 0
        self._cached_product = 1

    def _new_bucket(self):
        tree = _QuotientUpdateTree()
        tree.capacity = self.bucket_size
        tree.leaves = [1] * self.bucket_size
        tree.products = [1] * (2 * self.bucket_size)
        tree.free = []
        return tree

    def _grow(self):
        old = self.capacity
        self.buckets.extend(self._new_bucket() for _ in range(len(self.buckets)))
        self.capacity *= 2
        self.free.extend(range(old, self.capacity))
        heapify(self.free)
        self.capacity_doublings += 1

    def insert(self, norm):
        assert norm > 1
        if not self.free:
            self._grow()
        slot = heappop(self.free)
        bucket, local = divmod(slot, self.bucket_size)
        self.buckets[bucket]._set(local, norm)
        self._cached_product = None
        return slot

    def remove(self, slot):
        bucket, local = divmod(slot, self.bucket_size)
        assert self.buckets[bucket].leaves[local] > 1
        self.buckets[bucket]._set(local, 1)
        heappush(self.free, slot)
        self._cached_product = None

    def _gcd(self, left, right):
        self.gcd_calls += 1
        self.gcd_operand_bits += left.bit_length() + right.bit_length()
        return gcd(left, right)

    def _route(self, bucket, common):
        tree = self.buckets[bucket]
        found, stack = [], [(1, common)]
        while stack:
            node, value = stack.pop()
            if value == 1:
                continue
            if node >= self.bucket_size:
                found.append(bucket * self.bucket_size + node - self.bucket_size)
                continue
            left, right = tree.products[2 * node], tree.products[2 * node + 1]
            if left == 1:
                stack.append((2 * node + 1, value))
            elif right == 1:
                stack.append((2 * node, value))
            else:
                if left.bit_length() <= right.bit_length():
                    first = self._gcd(value, left)
                    second = value // first
                else:
                    second = self._gcd(value, right)
                    first = value // second
                self.routing_divisions += 1
                if second > 1:
                    stack.append((2 * node + 1, second))
                if first > 1:
                    stack.append((2 * node, first))
        return found

    def overlaps(self, norm):
        found, remaining = [], norm
        for index, tree in enumerate(self.buckets):
            if remaining == 1:
                break
            product = tree.products[1]
            if product == 1:
                continue
            common = self._gcd(remaining, product)
            if common > 1:
                found.extend(self._route(index, common))
                remaining //= common
                self.routing_divisions += 1
        return found

    def prepare_export(self):
        if self._cached_product is None:
            layer = [tree.products[1] for tree in self.buckets if tree.products[1] != 1]
            while len(layer) > 1:
                following = []
                for i in range(0, len(layer), 2):
                    if i + 1 == len(layer):
                        following.append(layer[i])
                    else:
                        self._export_multiplications += 1
                        self._export_operand_bits += layer[i].bit_length() + layer[i + 1].bit_length()
                        following.append(layer[i] * layer[i + 1])
                layer = following
            self._cached_product = layer[0] if layer else 1

    @property
    def products(self):
        self.prepare_export()
        return {1: self._cached_product}

    @property
    def leaves(self):
        return [norm for tree in self.buckets for norm in tree.leaves]

    @property
    def leaf_updates(self):
        return sum(tree.leaf_updates for tree in self.buckets)

    @property
    def multiplications(self):
        return self._export_multiplications + sum(tree.multiplications for tree in self.buckets)

    @property
    def multiplication_operand_bits(self):
        return self._export_operand_bits + sum(tree.multiplication_operand_bits for tree in self.buckets)


class BucketOrientedBasis(IncrementalOrientedBasis):
    def __init__(self, discriminant, *, bucket_size=128):
        super().__init__(discriminant)
        self.tree = CoprimeBucketForest(bucket_size)

    def result(self, *, dense=False, verify=True):
        self.tree.prepare_export()
        record = super().result(dense=dense, verify=verify)
        record['tree_update_method'] = 'coprime local product forest'
        record['statistics'].update(
            bucket_size=self.tree.bucket_size, bucket_count=len(self.tree.buckets),
            tree_exact_divisions=sum(tree.exact_divisions for tree in self.tree.buckets),
            tree_division_operand_bits=sum(tree.division_operand_bits for tree in self.tree.buckets),
            tree_routing_divisions=self.tree.routing_divisions,
            tree_export_multiplications=self.tree._export_multiplications,
            tree_export_multiplication_operand_bits=self.tree._export_operand_bits)
        return record


def bucket_oriented_coprime_basis(discriminant, keys, *, bucket_size=128, dense=False, verify=True):
    state = BucketOrientedBasis(discriminant, bucket_size=bucket_size)
    for norm, root in keys:
        state.add(norm, root)
    return state.result(dense=dense, verify=verify)
