#!/usr/bin/env python3
"""Independent encoder, norm, unranking and finite-state crossing controls."""
from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

from bounded_spin_crossing import (paths,suffix_counts,unrank_path,binary_path,
                                   path_entry,path_norm,recursive_vector,
                                   polynomial_entry,literal_entry,integral_polynomial_entry)
from check_multipoint_hermite_crossings import rank

ROOT = Path(__file__).resolve().parents[1]


def literal_integral_entry(selected):
    vectors = [recursive_vector(path) for path in selected]
    length = len(selected[0])
    total = [0]*(length+1)
    for choices in product(*(tuple(vector.items()) for vector in vectors)):
        words = [word for word,_ in choices]
        polynomial = [1]
        for _,c in choices:
            polynomial = [a*c for a in polynomial]
        for i in range(length):
            x,y,u,v = [(word >> i)&1 for word in words]
            identity,crossing = int(x == u and y == v),int(x == v and y == u)
            out = [0]*(len(polynomial)+1)
            for j,c in enumerate(polynomial):
                out[j] += identity*c
                out[j+1] += (crossing-identity)*c
            polynomial = out
        total = [a+b for a,b in zip(total,polynomial)]
    denominator = path_norm(selected[2])*path_norm(selected[3])
    assert 12**length % denominator == 0
    return [c*(12**length//denominator) for c in total]


def check():
    coefficients,grams,raising,unranking = 0,0,0,0
    reduced_sites = 0
    dimensions = []
    for length in range(11):
        row = []
        for height in range(4):
            ps = paths(length,height)
            assert len(ps) == suffix_counts(length,height)[length][0]
            vs = [recursive_vector(path) for path in ps]
            row.append(len(ps))
            for index,(path,vector) in enumerate(zip(ps,vs)):
                assert unrank_path(length,index,height) == path
                unranking += 1
                for word in range(1 << length):
                    assert path_entry(path,word) == vector.get(word,0)
                    coefficients += 1
                e = defaultdict(int)
                for word,c in vector.items():
                    for j in range(length):
                        if word >> j & 1:
                            e[word^(1 << j)] += c
                assert not any(e.values())
                raising += 1
                assert sum(c*c for c in vector.values()) == path_norm(path)
                for other,w in zip(ps,vs):
                    dot = sum(c*w.get(word,0) for word,c in vector.items())
                    assert dot == (path_norm(path) if path == other else 0)
                    grams += 1
            if length >= 2 and length % 2 == 0 and height == 0:
                dimensions.append(len(ps))
                for v,w in product(vs,repeat=2):
                    gram = sum(c*w.get(word,0) for word,c in v.items())
                    for site in range(length):
                        reduced = [[0,0],[0,0]]
                        for word,c in v.items():
                            a = (word >> site)&1
                            for b in range(2):
                                other = (word & ~(1 << site)) | (b << site)
                                reduced[a][b] += c*w.get(other,0)
                        assert [[2*c for c in row] for row in reduced] == [[gram,0],[0,gram]]
                        reduced_sites += 1
    assert dimensions == [1,2,5,13,34]
    binary = []
    for width in range(1,9):
        chosen = [binary_path(width,label) for label in range(1 << width)]
        assert len(set(chosen)) == 1 << width
        assert all(sum(path) == 0 for path in chosen)
        binary.append({'width':width,'underlying_length':len(chosen[0]),'labels':len(chosen)})
    four = paths(4)
    table = {}
    literal = 0
    for labels in product(range(len(four)),repeat=4):
        selected = tuple(four[i] for i in labels)
        actual = polynomial_entry(selected)
        assert actual == literal_entry(selected)
        expected = int(labels[0] == labels[3] and labels[1] == labels[2])
        assert sum(actual) % 5 == expected
        integral = integral_polynomial_entry(selected)
        assert integral == literal_integral_entry(selected)
        assert sum(integral) == 12**4*expected
        assert 2*integral[1]+4*integral[0] == 0
        assert 2*sum(j*c for j,c in enumerate(integral))-4*sum(integral) == 0
        assert 2*integral[-2]+4*integral[-1] == 0
        assert [c % 5 for c in integral] == [pow(12,4,5)*c % 5 for c in actual]
        table[labels] = actual
        literal += 1
    rank_pairs = []
    for j in range(5):
        pairs = []
        for left,right in (((0,1),(2,3)),((0,2),(1,3))):
            matrix = []
            for a,b in product(range(2),repeat=2):
                row = []
                for c,d in product(range(2),repeat=2):
                    entry = dict(zip(left,(a,b)))|dict(zip(right,(c,d)))
                    row.append(table[tuple(entry[i] for i in range(4))][j])
                matrix.append(row)
            pairs.append(rank(matrix,5))
        rank_pairs.append(pairs)
    random = Random(314159)
    six = paths(6)
    endpoints = 0
    for i in range(12):
        labels = [random.randrange(len(six)) for _ in range(4)]
        if i % 3 == 0:
            labels[2],labels[3] = labels[1],labels[0]
        polynomial = polynomial_entry(tuple(six[j] for j in labels))
        assert sum(polynomial) % 5 == int(labels[0] == labels[3] and labels[1] == labels[2])
        endpoints += 1
    first = paths(2)[0]
    assert sum(polynomial_entry((first,)*4,normalize=False)) % 5 == 4
    assert sum(polynomial_entry((first,)*4)) % 5 == 1
    return {'status':'pass','prime':5,'height_cap':3,'encoder_word_coefficients':coefficients,
            'integer_gram_entries':grams,'highest_weight_vectors':raising,
            'integral_one_site_reduced_matrices':reduced_sites,
            'independent_integral_three_jet_controls':literal,
            'complete_path_unrank_controls':unranking,'return_to_zero_dimensions':dimensions,
            'exact_binary_width_controls':binary,'independent_literal_polynomial_entries':literal,
            'independent_integral_polynomial_entries':literal,'integral_crossing_scale_at_length_four':12**4,
            'four_bit_coefficient_rank_pairs':rank_pairs,'six_bit_endpoint_samples':endpoints,
            'unscaled_negative_control_endpoint':4,
            'scope':'Uniform coherent crossing construction with finite encoder and endpoint controls; no uniform rank saving or factoring improvement.',
            'source_sha256':{name:sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                              for name in ('bounded_spin_crossing.py','check_bounded_spin_crossing.py')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
