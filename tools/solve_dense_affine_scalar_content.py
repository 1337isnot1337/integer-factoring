#!/usr/bin/env python3
"""Solve scalar-zero constraints in the full finite rank-one frame span.

This is exact linear algebra over F4, not a pressure or optimality claim.
Each column is the complete coefficient vector of T_G,k minus T. The saved
rank and remainder describe precisely the declared finite catalogue.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from dense_coupled_reflection_frames import rank_one_terms, transformed_from_terms
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import add, bell_normalize, scale
from screen_affine_finite_places import evaluate


def pack(values):
    values = list(values)
    return (sum((value & 1) << i for i, value in enumerate(values)),
            sum((value >> 1) << i for i, value in enumerate(values)))


class Span:
    def __init__(self, points, base):
        self.points = points
        self.target = pack(evaluate(value, p) for p in points for value in base['values'])
        self.pivots = {}
        self.increments = []

    def reduce(self, value, combination=(0, 0)):
        while value[0] | value[1]:
            position = ((value[0] | value[1]) & -(value[0] | value[1])).bit_length() - 1
            if position not in self.pivots:
                break
            coefficient = ((value[0] >> position) & 1) + 2 * ((value[1] >> position) & 1)
            pivot, witness = self.pivots[position]
            value = add(value, scale(pivot, coefficient))
            combination = add(combination, scale(witness, coefficient))
        return value, combination

    def include(self, values, index):
        column = pack(evaluate(value, p) for p in self.points for value in values)
        column = add(column, self.target)
        reduced, witness = self.reduce(column, (1 << index, 0))
        if not (reduced[0] | reduced[1]):
            return
        position = ((reduced[0] | reduced[1]) & -(reduced[0] | reduced[1])).bit_length() - 1
        coefficient = ((reduced[0] >> position) & 1) + 2 * ((reduced[1] >> position) & 1)
        inverse = (0, 1, 3, 2)[coefficient]
        self.pivots[position] = scale(reduced, inverse), scale(witness, inverse)
        self.increments.append({'candidate': index, 'pivot': position})

    def report(self, candidates):
        remainder, word = self.reduce(self.target)
        consistent = not (remainder[0] | remainder[1])
        result = {'points': self.points, 'constraint_count': 256 * len(self.points),
                  'rank': len(self.pivots), 'nullity': len(candidates) - len(self.pivots),
                  'consistent': consistent, 'rank_increments': self.increments,
                  'target_remainder_bitplanes': [hex(x) for x in remainder]}
        if consistent:
            result['particular_word'] = [dict(candidates[i], coefficient=((word[0] >> i) & 1) + 2 * ((word[1] >> i) & 1))
                                         for i in range(len(candidates)) if ((word[0] | word[1]) >> i) & 1]
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue', type=Path, default=Path('evidence/dense-rank-one-coupled-frames.json'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--limit', type=int)
    args = parser.parse_args()
    raw = args.catalogue.read_bytes()
    catalogue = json.loads(raw)
    base = tensor(4)
    _, factors = bell_normalize(base)
    spans = [Span(points, base) for points in ([0], [2], [3], [2, 3], [0, 2, 3])]
    candidates = []
    began = time.monotonic()
    records = catalogue['records'][:args.limit]
    for record in records:
        terms = rank_one_terms(base, record['u'], record['v'], factors)
        for candidate in record['candidates']:
            changed = transformed_from_terms(base, terms, candidate['pole'])
            digest = sha256(json.dumps(changed['values'], separators=(',', ':')).encode()).hexdigest()
            if digest != candidate['seed_sha256']:
                raise ArithmeticError('Exact reconstruction does not match the input catalogue.')
            index = len(candidates)
            candidates.append({'generator_index': record['index'], 'pole': candidate['pole'], 'seed_sha256': digest})
            for span in spans:
                span.include(changed['values'], index)
        if len(candidates) % 100 < len(record['candidates']):
            print(json.dumps({'candidates': len(candidates), 'ranks': [len(s.pivots) for s in spans],
                              'elapsed_seconds': time.monotonic() - began}), flush=True)
    names = (Path(__file__).name, 'dense_coupled_reflection_frames.py', 'reflection_dual_rail_gaussian.py',
             'reflection_laurent_physical_frames.py', 'screen_affine_finite_places.py')
    report = {'schema_version': 1, 'scope': __doc__, 'catalogue_sha256': sha256(raw).hexdigest(),
              'generators': len(records), 'candidate_count': len(candidates),
              'complete_catalogue': len(records) == len(catalogue['records']),
              'endpoint': 1, 'alphabet': 4, 'degree': 24,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'results': [s.report(candidates) for s in spans], 'elapsed_seconds': time.monotonic() - began}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ('generators', 'candidate_count', 'complete_catalogue', 'elapsed_seconds')}), flush=True)
    for result in report['results']:
        print(json.dumps({k: result[k] for k in ('points', 'rank', 'nullity', 'consistent')}), flush=True)


if __name__ == '__main__':
    main()
