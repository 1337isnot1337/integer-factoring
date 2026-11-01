#!/usr/bin/env python3
"""Registered complete principal-pair screen on 67 retained F4 pencil points."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import time

from check_scaled_asymmetric_pressure import inputs
from native_rectangular_smith_f4 import RectangularOracle
from pencil_principal_alphabets import pressure_key, restrict
from search_f4_recentered_crossings import root_for_profiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--directory',type=Path,required=True)
    args = parser.parse_args()
    root,started = Path(__file__).resolve().parents[1],time.monotonic()
    manifest = json.loads(args.manifest.read_text())
    original,cases = inputs(root)
    assert sha256((root/'evidence/two-place-free-gauge-asymmetric-manifest.json').read_bytes()).hexdigest() == manifest['source_manifest_sha256']
    assert len(cases)==manifest['source_cases']==67
    subsets = [[i for i in range(5) if mask >> i & 1] for mask in manifest['subset_masks']]
    assert len(subsets)==26 and sum(map(len,subsets))==75
    assert len(cases)*len(subsets)**2 == manifest['registered_cases']==45292
    args.directory.mkdir(parents=True,exist_ok=True)
    if (args.directory/'profiles.jsonl').exists():
        raise ValueError('Use a fresh output directory; do not overwrite a previous ledger.')
    names = ['tools/screen_pencil_principal_alphabets.py','tools/pencil_principal_alphabets.py',
             'tools/check_scaled_asymmetric_pressure.py','tools/two_place_f4_census.py',
             'tools/search_f4_recentered_crossings.py','tools/native_rectangular_smith_f4.py',
             'tools/native_rectangular_smith_f4.cpp','tools/native_weighted_smith_f4.py',
             'tools/joint_energy_conjugations.py']
    report = {'status':'running_incomplete','registration':manifest,'completed_cases':0,
              'resolved_orders':0,'unresolved_profiles':[],'scalar_curve_classes':0,
              'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'as_run_sources':{n:(root/n).read_text() for n in names},
              'scope':'Complete registered inherited-gauge principal restrictions, not all pencil '
                      'members at new gauges or a global pressure exclusion. Fits are numerical.'}
    cache,fit_records,best = {},[],None
    def save():
        report['seconds'] = time.monotonic()-started
        (args.directory/'summary.tmp').write_text(json.dumps(report,indent=2)+'\n')
        (args.directory/'summary.tmp').replace(args.directory/'summary.json')
    save()
    with (args.directory/'profiles.jsonl').open('w') as ledger:
        for case,values in cases:
            for ax in subsets:
                for ay in subsets:
                    label = {'source_case_index':case['index'],'alphabet_x':ax,'alphabet_y':ay}
                    seed,parameters,gx,gy = restrict(values,original['natural_heights_and_scale'],ax,ay)
                    try:
                        native = RectangularOracle(seed['values'],len(ax),len(ay)).prefixes(parameters,1)
                    except ArithmeticError as error:
                        report['unresolved_profiles'].append(dict(label,error=str(error)))
                        save()
                        raise
                    n = len(ax)*len(ay)
                    orders = list(map(round,native['orders']))
                    assert orders == native['orders'] and len(orders)==2*n
                    profiles = [orders[:n],orders[n:]]
                    key = pressure_key(profiles,gx,gy)
                    if key not in cache:
                        cache[key] = len(fit_records)
                        fit_records.append(root_for_profiles(profiles,gx,gy))
                    fit_index = cache[key]
                    result = dict(label,parameters=parameters,orders=profiles,physical_energies=[gx,gy],
                        scalar_fit_index=fit_index,
                        literal_values_sha256=sha256(json.dumps(seed['values'],separators=(',',':')).encode()).hexdigest())
                    ledger.write(json.dumps(result,separators=(',',':'))+'\n')
                    fit = fit_records[fit_index]
                    if best is None or fit['rho'] < best['rho']:
                        best = dict(label,**fit)
                    report['completed_cases'] += 1
                    report['resolved_orders'] += 2*n
            ledger.flush()
            report.update(scalar_curve_classes=len(fit_records),best=best)
            save()
            print(json.dumps({'source_case':case['index'],'completed':report['completed_cases'],
                              'fit_classes':len(fit_records),'best':best,'seconds':report['seconds']}),flush=True)
    assert report['completed_cases']==45292 and report['resolved_orders']==753750
    (args.directory/'scalar-fits.json').write_text(json.dumps(fit_records,indent=2)+'\n')
    report.update(status='complete_registered_screen',best=best,
                  profile_ledger_sha256=sha256((args.directory/'profiles.jsonl').read_bytes()).hexdigest(),
                  scalar_fits_sha256=sha256((args.directory/'scalar-fits.json').read_bytes()).hexdigest())
    save()


if __name__ == '__main__':
    main()
