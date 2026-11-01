#!/usr/bin/env python3
"""Three registered free-height fits of the fixed locally twisted CSS source."""
import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np

from bounded_complete_support_fit import fit
from nonclifford_css_crossing import seed

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    source,planes,metadata = seed()
    randomizer = np.random.default_rng(20261101)
    graded = np.array([i//3+3*(i % 3) for i in range(9)],float)[1:]
    zero = np.zeros(8)
    starts = [
        ('near-zero-balanced',np.r_[randomizer.normal(0,.05,64),math.log(.1),0.],True),
        ('graded-balanced',np.r_[np.r_[zero,zero,graded,graded,graded,graded,zero,zero]+randomizer.normal(0,.05,64),math.log(.1),0.],True),
        ('near-zero-free-angle',np.r_[randomizer.normal(0,.05,64),math.log(.1),0.],False)]
    report = dict(scope='Finite fixed-code weighted-type fits; no global optimum or factoring improvement.',
                  status='running',metadata=metadata,seed=20261101,
                  limits=dict(seconds_per_trial=20,maximum_iterations=100,free_heights=64,
                              height_bounds=[-128,128],log_scale_bounds=[-9,math.log(4)],log_tangent_bounds=[-2,2]),
                  bare_optimization='skipped by proved q9 box floor',fits=[],
                  source_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in
                      ('tools/screen_nonclifford_css_crossing.py','tools/nonclifford_css_crossing.py',
                       'tools/bounded_complete_support_fit.py','tools/scaled_asymmetric_crossing_pressure.py',
                       'tools/asymmetric_crossing_pressure.py')})
    begun = time.monotonic()

    def save():
        report['seconds'] = time.monotonic()-begun
        args.output.write_text(json.dumps(report,indent=2)+'\n')

    save()
    if metadata['structural_gate']:
        for label,start,balanced in starts:
            attempt = fit(source,start,seconds=20,maxiter=100,balanced=balanced,height_bound=128.)
            attempt['label'] = label
            report['fits'].append(attempt)
            save()
            print(json.dumps({'label':label,'status':attempt['status'],'queries':attempt['queries'],
                              'best_support':attempt['best']['rate'] if attempt['best'] else None}),flush=True)
        report['status'] = 'completed_registered_attempts'
    else:
        report['status'] = 'stopped_at_structural_gate'
    candidates = [attempt['best'] for attempt in report['fits'] if attempt['best']]
    report['continuation_gate'] = any(candidate['rate'] < 1-1e-6 for candidate in candidates)
    save()


if __name__ == '__main__':
    main()
