#!/usr/bin/env python3
"""Complete two-cut profiles of one declared paired divided-power code."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from coherent_divided_powers import graded_dimensions, paired_top
from coherent_hot_compounds import profiles
from coherent_inverse_systems import inverse_system
from search_coherent_hot_retractions import pressure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes', type=int, required=True)
    parser.add_argument('--power', type=int, default=2)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 2 <= args.modes <= 16 or args.modes % 2 or not 0 <= args.power <= args.modes//2:
        parser.error('The native probe supports positive even mode counts through sixteen.')
    top = paired_top(args.modes, args.power)
    code = inverse_system([sum(1 << word for word in top)], args.modes)
    assert Counter(zip(code['charges'], code['heights'])) == graded_dimensions(args.modes, args.power)
    charges = Counter(code['charges'])
    maximum = max(sum(a*b for i, a in charges.items() for j, b in charges.items() if i+j == r)
                  for r in range(4*args.power+1))
    print(json.dumps(dict(modes=args.modes, power=args.power, q=len(code['embedding']),
                          maximum_block=maximum)), flush=True)
    result = profiles(code['embedding'], code['retraction'], args.modes)
    names = (Path(__file__).name, 'coherent_divided_powers.py', 'coherent_inverse_systems.py',
             'coherent_hot_compounds.py', 'coherent_hot_compounds.cpp',
             'selected_gaussian_coefficients.cpp', 'check_binary_coefficient_entropy.py',
             'search_coherent_hot_retractions.py')
    report = dict(modes=args.modes, power=args.power, top=top, code=code,
                  profile=result, pressure=pressure(result),
                  source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                 for name in names},
                  scope='Complete natural profiles of one specified finite code; numerical pressure fit only.')
    args.output.write_text(json.dumps(report, separators=(',', ':'))+'\n')
    print(json.dumps(dict(pressure=report['pressure'], seconds=result['seconds'])))


if __name__ == '__main__':
    main()
