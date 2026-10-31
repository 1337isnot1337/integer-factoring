#!/usr/bin/env python3
"""Exact Fourier-positive relation kernel and a finite presentation-to-torsion pipeline."""
from fractions import Fraction as Q
import json

from arithmetic_two_group import project_to_two_part, two_group_basis


def determinant(matrix):
    work = [row[:] for row in matrix]
    previous, sign = 1, 1
    for k in range(len(work)-1):
        pivot = next((i for i in range(k, len(work)) if work[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            work[k], work[pivot] = work[pivot], work[k]
            sign = -sign
        value = work[k][k]
        for i in range(k+1, len(work)):
            for j in range(k+1, len(work)):
                numerator = value*work[i][j]-work[i][k]*work[k][j]
                assert numerator % previous == 0
                work[i][j] = numerator//previous
            work[i][k] = 0
        previous = value
    return sign*work[-1][-1] if work else 1


def main():
    size, dimension = 16, 4
    group = list(range(size))
    sample = [x for x in group if ((x & 1)*((x >> 1) & 1) ^
                                  ((x >> 2) & 1)*((x >> 3) & 1)) == 1]
    selected = [x for x in group if x not in sample][:5]
    character = lambda a, x: 1 if (a & x).bit_count() % 2 == 0 else -1
    density = Q(len(selected), size)
    bias = max(abs(Q(sum(character(a, x) for x in sample), len(sample))) for a in group[1:])
    assert len(sample) == 6 and bias == Q(1, 3)
    assert not set(sample) & set(selected)
    correlation = {g: Q(sum((x ^ g) in selected for x in selected), size) for g in group}
    transforms = []
    for a in group:
        fhat = Q(sum(character(a, x) for x in selected), size)
        chat = sum(correlation[g]*character(a, g) for g in group)/size
        assert chat == fhat*fhat >= 0
        transforms.append(chat)
    assert transforms[0] == density*density and sum(transforms) == density
    lower_bound = density*density-bias*(density-density*density)
    averages = {g: sum(correlation[g ^ u] for u in sample)/len(sample) for g in group}
    assert lower_bound > 0 and min(averages.values()) >= lower_bound
    witnesses = {}
    for target in group:
        witnesses[target] = next((u, a, b) for u in sample for a in selected for b in selected
                                 if target ^ u == a ^ b)
    word = lambda x: [(x >> i) & 1 for i in range(dimension)]
    smooth_word_bound = max(x.bit_count() for x in selected)
    sample_word_bound = max(x.bit_count() for x in sample)
    correction_bound = sample_word_bound+2*smooth_word_bound
    diagonal = correction_bound+1
    relations = []
    for i in range(dimension):
        target = (1 << i) if diagonal % 2 else 0
        u, a, b = witnesses[target]
        row = [diagonal*int(i == j)+word(u)[j]-word(a)[j]+word(b)[j] for j in range(dimension)]
        assert row[i] > sum(abs(row[j]) for j in range(dimension) if j != i)
        assert all(value % 2 == 0 for value in row)
        relations.append(row)
    annihilator = abs(determinant(relations))
    assert annihilator > size and annihilator % 2 == 0
    half_relation_span = {0}
    for row in relations:
        image = sum(((value//2) & 1) << i for i, value in enumerate(row))
        half_relation_span |= {value ^ image for value in half_relation_span}
    assert len(half_relation_span) == 2
    generators = [1 << i for i in range(dimension)]
    projected = project_to_two_part(generators, annihilator, 0, int.__xor__)
    basis = two_group_basis(projected, 0, int.__xor__, lambda x: x, 1)
    assert len(basis) == 4 and all(order == 1 for _, order in basis)
    span = {0}
    for element, _ in basis:
        span |= {value ^ element for value in span}
    assert span == set(group)
    print(json.dumps({'scope': 'Exact finite calibration; pair witnesses are obtained by exhaustive finite search',
                      'group': '(Z/2Z)^4', 'smooth_set_control': selected, 'sample': sample,
                      'sample_bias': str(bias), 'density': str(density),
                      'correlation_lower_bound': str(lower_bound),
                      'actual_minimum_sampled_correlation': str(min(averages.values())),
                      'translated_pair_witnesses': len(witnesses),
                      'original_set_missed_by_sample': True,
                      'relation_matrix': relations, 'relation_determinant': annihilator,
                      'quotient_index_exceeds_group_size': True,
                      'halved_relation_image_size': len(half_relation_span),
                      'extracted_two_torsion_rank': len(basis)}, indent=2))


if __name__ == '__main__':
    main()
