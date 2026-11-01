#!/usr/bin/env python3
"""Replay a registered list of independently paired coordinate alphabets."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from permuted_pascal_kernels import permuted_pascal_kernel
from rectangular_pascal_sectors import CoordinateTensorBuilder, RectangularPascalPressure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, started = Path(__file__).resolve().parents[1], time.monotonic()
    manifest = json.loads(args.manifest.read_text())
    seed = permuted_pascal_kernel(manifest['images'])
    builder = CoordinateTensorBuilder(seed)
    names = ['tools/screen_rectangular_pascal_sectors.py','tools/rectangular_pascal_sectors.py',
             'tools/permuted_pascal_kernels.py','tools/weighted_pascal_flags.py',
             'tools/cube_f4_bivariate_smith.py','tools/native_matrix_smith_f4.py',
             'tools/native_matrix_smith_f4.cpp','tools/native_rectangular_smith_f4.cpp',
             'tools/native_weighted_smith_f4.py','tools/pascal_coordinate_sectors.py',
             'tools/tensor_pascal_box_crossings.py']
    report = {'status':'running_incomplete', 'registration':manifest,
              'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'as_run_sources':{n:(root/n).read_text() for n in names}, 'completed':[],
              'scope':'Complete original-polynomial profiles; numerical inherited-gauge balanced '
                      'and supporting-line pressures only. No exact contraction certificate.'}
    def save():
        report['total_seconds'] = time.monotonic()-started
        args.output.parent.mkdir(parents=True,exist_ok=True)
        tmp = args.output.with_suffix('.tmp')
        tmp.write_text(json.dumps(report,indent=2)+'\n')
        tmp.replace(args.output)
    save()
    for case in manifest['cases']:
        start, attempts = time.monotonic(), []
        for precision in manifest['precision_policy']:
            try:
                model = RectangularPascalPressure(seed, case['alphabet_b'], case['alphabet_x'],
                                                 precision=precision, builder=builder)
                parameters = model.fitted_natural_start()
                model.evaluate(parameters)
                attempts.append({'precision':precision,'status':'pass'})
                break
            except ArithmeticError as error:
                attempts.append({'precision':precision,'status':'unresolved','error':str(error)})
        else:
            report['failure'] = {'case':case['label'],'attempts':attempts}
            save()
            raise ArithmeticError('Every registered precision failed; retain the incomplete receipt.')
        baseline = model.last.copy()
        unit = model.natural(1.)
        # Scale covariance permits the entire inherited supporting-line curve
        # to be evaluated from these two complete resolved profile lists.
        profiles = [np.asarray(c['orders'])/parameters[-1]
                    for c in baseline['complete_original_profiles']]
        physical = model.physical(unit)
        def objective(v):
            scale, angle = math.exp(v[0]), math.atan(math.exp(v[1]))
            weights = [math.cos(angle),math.sin(angle)]
            log_b = sum(w*float(logsumexp(-scale*p/w)) for w,p in zip(weights,profiles))
            return model.rate(log_b,[scale*p for p in physical])[0]
        supports = []
        for log_tangent in manifest['support_log_tangent_starts']:
            fitted = minimize(objective,[math.log(parameters[-1]),log_tangent],method='L-BFGS-B',
                              bounds=manifest['support_bounds'],options={'maxiter':300,'ftol':1e-13})
            supports.append({'value':float(fitted.fun),'log_scale':float(fitted.x[0]),
                             'log_tangent':float(fitted.x[1]), 'success':bool(fitted.success),
                             'iterations':int(fitted.nit),'message':str(fitted.message)})
        record = {'case':case,'attempts':attempts,'construction':model.construction,
                  'balanced':baseline,'support_fits':supports,'queries':model.queries,
                  'seconds':time.monotonic()-start}
        report['completed'].append(record)
        save()
        print(json.dumps({'case':case['label'],'rho':baseline['rho'],
                          'support_value':min(r['value'] for r in supports),
                          'seconds':record['seconds']}),flush=True)
        del model
    report['status'] = 'complete_registered_screen'
    save()
    print(json.dumps({'status':report['status'],'cases':len(report['completed']),
                      'seconds':report['total_seconds']}),flush=True)


if __name__ == '__main__':
    main()
