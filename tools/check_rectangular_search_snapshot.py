#!/usr/bin/env python3
"""Replay frozen exploratory pressure points with full rank and precision guards."""
import argparse
import hashlib
import json
from pathlib import Path
from search_rectangular_core_pressure import CoreDefectPressure

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=Path('evidence/rectangular-search-snapshots.json'))
    parser.add_argument('--output',type=Path,default=Path('evidence/rectangular-search-snapshot-replay.json'))
    args=parser.parse_args();raw=args.input.read_bytes();source=json.loads(raw)
    names=['rectangular_defect_profiles.py','rectangular_defect_coordinates.py','native_fock_f4.cpp','native_matrix_smith_f4.py','native_matrix_smith_f4.cpp','native_rectangular_smith_f4.cpp','search_rectangular_core_pressure.py',Path(__file__).name]
    hashes={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}
    records=[]
    for case in source['cases']:
        oracle=CoreDefectPressure(case['M'],case['alpha'],case['beta'],case['quartic_terms'],full_core=case['full_core'],repair_kind=case.get('repair_kind','partial'))
        rho,gradient,detail=oracle.evaluate(case['best']['parameters'])
        error=abs(rho-case['best']['rho'])
        profile_error=max(abs(a-b) for actual,expected in zip(detail['active_orders_times_log_s'],case['best']['active_orders_times_log_s']) for a,b in zip(actual,expected))
        assert error<1e-9 and profile_error<1e-8
        records.append({'M':case['M'],'replayed_rho':rho,'rho_discrepancy':error,'maximum_profile_discrepancy':profile_error,'entropy':detail['entropy'],'gradient_maximum':max(abs(v) for v in gradient),'active_matrix_orders':[o.n for o in oracle.oracles],'profile_windows':[m['profile_precision'] for m in oracle.metadata],'builder_windows':[m['builder_precision'] for m in oracle.metadata],'known_below':[m['known_below'] for m in oracle.metadata],'complete_numerators':[m['exact_polynomial_matrix'] for m in oracle.metadata]})
        print(json.dumps(records[-1]),flush=True)
    result={'status':'pass','scope':'Independent rerun of frozen floating points and all active orders. Exact field arithmetic and proved ranks are used; these are not rational pressure certificates or claims of global optimality.','cases':records,'input_sha256':hashlib.sha256(raw).hexdigest(),'source_sha256':hashes}
    args.output.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
