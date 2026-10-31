#!/usr/bin/env python3
"""Independent full-core entry and arbitrary four-leg weight checks at M4."""
import argparse
import json
from pathlib import Path

from check_rectangular_defect_core_basis import data_for
from check_rectangular_defect_coordinates import point_control, weighted_controls
from rectangular_defect_coordinates import coordinate_matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    points, profiles = [], []
    for place in ('infinity', 'zero'):
        data = data_for(place, 4, full_core=True)
        compiled = coordinate_matrix(data)
        points.extend(point_control(data, compiled, z) for z in (2, 3, 4, 5))
        profiles.append(weighted_controls(data, compiled, 5))
    result = {'status': 'independent original-tensor full-core coordinate and weight controls',
              'field': {'modulus': 'X^4+X+1', 'omega_image': 6},
              'entrywise_cases': points, 'weighted_profiles': profiles}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
