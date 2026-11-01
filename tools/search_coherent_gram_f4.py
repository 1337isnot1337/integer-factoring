#!/usr/bin/env python3
"""Four declared Gram-paired quartic probes, with complete exact profiles."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

from coherent_hot_compounds_f4 import profiles
from coherent_inverse_systems_f4 import inverse_system
from coherent_orthogonal_f4 import orthogonal_profiles, paired_retraction
from search_coherent_hot_retractions import pressure


def quartic(a, b):
    result = {word: 1 for word in (60, 90, 102, 153, 165, 195)}
    result.update({word: value for word, value in ((105, a), (150, b)) if value})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=int, choices=range(4), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    twisted = quartic(0, 0)
    twisted[60] = 2
    cases = ((quartic(2, 3), quartic(2, 3)),
             (quartic(2, 0), quartic(2, 0)),
             (quartic(1, 0), quartic(2, 0)),
             (quartic(0, 0), twisted))
    first, second = cases[args.case]
    p = inverse_system([first], 8)['embedding']
    q = inverse_system([second], 8)['embedding']
    maps = paired_retraction(p, q, 8)
    result = orthogonal_profiles(p, 8) if args.case < 2 else profiles(p, maps['retraction'], 8)
    names = (Path(__file__).name, 'coherent_orthogonal_f4.py', 'coherent_inverse_systems_f4.py',
             'coherent_hot_compounds_f4.py', 'coherent_hot_compounds_f4.cpp',
             'coherent_hot_compounds.py', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp', 'search_coherent_hot_retractions.py')
    report = dict(case=args.case, top_form=first, dual_top_form=second, maps=maps,
                  profile=result, pressure=pressure(result),
                  source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                 for name in names},
                  scope='One specified finite code pair; numerical pressure fit, no rate certificate or optimum.')
    args.output.write_text(json.dumps(report, separators=(',', ':'))+'\n')
    print(json.dumps(dict(case=args.case, q=result['q'], pressure=report['pressure'])))


if __name__ == '__main__':
    main()
