#!/usr/bin/env python3
"""Independent dense weighted-profile controls for refined core corrections."""
import argparse
import json
from pathlib import Path

from check_rectangular_f4_repair import build
from check_rectangular_refined_repairs import postcorrect_seed, single_shear
from check_rectangular_defect_core_basis import data_for
from check_rectangular_defect_coordinates import weighted_controls
from rectangular_defect_coordinates import coordinate_matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    base = build(tuple(range(4)), 4, 1, 2, check_types=False)
    records = []
    for row, column in ((3, 8), (8, 3), (10, 11), (1, 11)):
        spec = single_shear(row, column, 2)
        seed = postcorrect_seed(base, *spec, check_types=False)
        for place in ('infinity', 'zero'):
            data = data_for(place, 4)
            compiled = coordinate_matrix(data, repair_spec=spec)
            record = weighted_controls(data, compiled, 2, seed)
            record.update({'row': row, 'column': column, 'coefficient': 2})
            records.append(record)
    result = {'status': 'independent original-tensor weighted profiles of refined repairs', 'cases': records}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
