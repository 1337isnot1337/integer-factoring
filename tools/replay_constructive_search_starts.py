#!/usr/bin/env python3
"""Replay saved exploratory starts; these floating values are not certificates."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from check_coordinate_core_crossings import build as coordinate
from involutive_core_crossing import build as involutive
from native_rectangular_smith import RectangularOracle

ROOT = Path(__file__).resolve().parent.parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path,
                        default=ROOT/'evidence/constructive-frontier-search-starts.json')
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    records = []
    for record in data['starts']:
        if record['family'] == 'coordinate':
            seed, _ = coordinate(record['indices'], record['ambient'])
        elif record['family'] == 'involutive':
            seed = involutive(record['indices'], record['ambient'])
        else:
            raise ValueError('Unknown saved construction')
        digest = hashlib.sha256(json.dumps(seed['values'],separators=(',',':')).encode()).hexdigest()
        assert digest == record['coefficient_sha256']
        q = seed['q']
        oracle = RectangularOracle(seed['values'],q,q,
                                   denominator_degree=seed.get('denominator_degree',0),
                                   endpoint_check=False)
        x_type, y_type = record['type_x_v'],record['type_y_u']
        cost, _ = oracle.evaluate(record['params'],x_type,y_type,record['orientation'])
        entropy = sum(-p*math.log(p) for p in x_type+y_type if p)
        rho = cost/(2*entropy)
        difference = abs(rho-record['observed_rho'])
        assert difference < 1e-7, (record['name'],rho,record['observed_rho'])
        records.append(dict(name=record['name'],replayed_rho=rho,absolute_difference=difference))
    print(json.dumps(dict(status='floating replay only; no certified pressure gain',starts=records),indent=2))


if __name__ == '__main__':
    main()
