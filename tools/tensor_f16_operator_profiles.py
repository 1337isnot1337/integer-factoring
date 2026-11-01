#!/usr/bin/env python3
"""Profile an explicitly supplied list of homogeneous F16 differential operators."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

import tensor_hot_compounds_f16 as field
from tensor_f16_fast_kernel_profiles import profile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--operators', required=True,
                        help='JSON list of objects mapping occupation masks to coefficients.')
    parser.add_argument('--label', default='explicit homogeneous kernel')
    parser.add_argument('--height', type=int, required=True)
    parser.add_argument('--order-cap', type=int, default=64)
    parser.add_argument('--retraction', choices=('orthogonal', 'coordinate'), default='orthogonal')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    operators = [{int(mask): coefficient for mask, coefficient in row.items()}
                 for row in json.loads(args.operators)]
    code = field.differential_kernel(operators, 16, args.height)
    left = code['retraction'] if args.retraction == 'coordinate' else field.orthogonal_retraction(code)
    assert left is not None
    names = [Path(__file__).name, 'tensor_f16_fast_kernel_profiles.py',
             'tensor_hot_compounds_f16_fast.py', 'tensor_hot_compounds_f16_fast.cpp',
             'tensor_hot_compounds_f16.py', 'tensor_hot_compounds_f16.cpp',
             'tensor_orthogonal_scaling.py', 'coherent_hot_compounds_f4.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp',
             'search_coherent_hot_retractions.py']
    report = {'scope': __doc__, 'label': args.label, 'code': code, 'retraction': left,
              'retraction_kind': args.retraction, 'order_cap': args.order_cap,
              'completed': False, 'in_progress': None,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    started = time.monotonic()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        report['elapsed_seconds'] = time.monotonic() - started
        temporary = args.output.with_suffix(args.output.suffix + '.tmp')
        temporary.write_text(json.dumps(report, separators=(',', ':')) + '\n')
        temporary.replace(args.output)

    def callback(block):
        report['in_progress'] = block
        save()

    print(json.dumps({'starting': args.label, 'height': args.height, 'q': code['q']}), flush=True)
    data, fits = profile(code, left, args.retraction, args.order_cap, callback)
    report.update(completed=True, in_progress=None, profile=data, pressure_fits=fits)
    save()
    print(json.dumps({'completed': True, 'q': code['q'], 'fits': fits, 'seconds': data['seconds']}), flush=True)


if __name__ == '__main__':
    main()
