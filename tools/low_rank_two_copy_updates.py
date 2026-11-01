#!/usr/bin/env python3
"""Certified multiple pure corrections to the complete two-copy crossing."""
from __future__ import annotations

from entangled_two_copy_frames import EntangledTwoCopyFrames
from f4_recentered_lattices import pack, unpack
from native_certified_low_rank_smith_f4 import certified_low_rank_profile
from two_copy_crossing_couplings import reverse
from weighted_smith_oracle_f4 import scalar_multiply


class LowRankTwoCopyFrames(EntangledTwoCopyFrames):
    def profiles(self, updates, work_precision=512):
        if not 1 <= len(updates) <= 32:
            raise ValueError('Use between one and32 complete pure corrections.')
        profiles, certificates = [], []
        for side in (0, 1):
            vectors = []
            for update in updates:
                physical = update['vectors']; polynomial = update['polynomial']; coefficient = update.get('coefficient', 1)
                if coefficient not in (1, 2, 3) or polynomial < 0 or polynomial.bit_count() % 2:
                    raise ValueError('Use nonzero field scalars and endpoint-zero polynomials.')
                if polynomial.bit_length() > 2*self.compiler.degree+1:
                    raise ValueError('The complete degree bound is exceeded.')
                left, right = self.transformed_vectors(physical, side)
                local = reverse(polynomial, 2*self.compiler.degree) if side == 0 else polynomial
                local = scalar_multiply(pack((local, 0)), coefficient)
                left = [unpack(self.mul(local, pack(a))) for a in left]
                vectors.append((left, right))
            result = certified_low_rank_profile(self.diagonals[side], vectors,
                                                 cutoff=self.jet_precision, work_precision=work_precision)
            profiles.append(tuple(a-2*self.compiler.scalars[side] for a in result['orders']))
            certificates.append({k: v for k, v in result.items() if k != 'orders'})
        return tuple(profiles), certificates

    def original_sum_matrices(self, updates):
        result = [[row[:] for row in matrix] for matrix in self.compiler.base]
        for update in updates:
            current = self.original_matrices(update['vectors'], update['polynomial'], update.get('coefficient', 1))
            for side in (0, 1):
                base = self.compiler.base[side]
                for i, row in enumerate(current[side]):
                    for j, (low, high) in enumerate(row):
                        old_low, old_high = result[side][i][j]
                        base_low, base_high = base[i][j]
                        result[side][i][j] = (old_low ^ low ^ base_low, old_high ^ high ^ base_high)
        return result
