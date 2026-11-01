#!/usr/bin/env python3
"""Replay complete coefficient fibers and exact generic polynomial row frames."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import cube_coefficient_fibers as fibers
import cube_bivariate_smith as smith

ROOT=Path(__file__).resolve().parents[1]
BASE=[0,2,4,6,8,10,10,12,12,14,14,16,18,20,22,24]

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def dominated(orders):
    return all(sum(orders[:k])<=sum(BASE[:k]) for k in range(1,17))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=ROOT/'evidence/cube-coefficient-floor4-input.json')
    parser.add_argument('--census',type=Path,default=ROOT/'evidence/cube-coefficient-floor4-census.json')
    parser.add_argument('--recompute-census',action='store_true')
    parser.add_argument('--output',type=Path,default=ROOT/'evidence/cube-coefficient-fibers-controls.json')
    args=parser.parse_args()
    retained=json.loads(args.input.read_text())
    for name,expected in retained['source_sha256'].items():
        assert digest(ROOT/'tools'/name)==expected,name
    directions,symbolic,orbits,group=fibers.basis('none')
    variables=len(directions)
    assert variables==12672 and len(orbits)==12800 and len(group)==1
    seed,left,right,diagonal,_=fibers.frames()
    reference=[exponent[0] for exponent in diagonal]
    assert sorted(reference)==BASE
    equations=fibers.three_cut_equations(symbolic,left,right,diagonal,variables)
    pivots_by_loss={}
    for loss,dimension in ((3,0),(4,22)):
        targets=tuple(max(0,order-loss) for order in reference)*3
        pivots=fibers.solve(equations,targets,variables)
        assert pivots is not None and len(pivots)==variables-dimension
        assert not any(row>>variables for row in pivots.values())
        pivots_by_loss[loss]=pivots
    reconstructed=fibers.floor_four_input(directions,seed,pivots_by_loss[4],variables)
    for key,value in reconstructed.items():assert retained[key]==value,key
    assert retained['reference_orders']==reference
    targets=[max(0,order-4) for order in reference]*3
    assert retained['row_targets']==targets
    assert reconstructed['specialization_image_dimension']==22
    assignments=[int(value,16) for value in reconstructed['assignments']]
    coefficient_controls=0
    for assignment in assignments:
        tensor=fibers.literal(seed,directions,assignment)
        smith.validate_tensor(tensor)
        coefficient_controls+=len(tensor)
        for i,target in enumerate(targets):
            for power in range(target):
                for equation in equations[i][power]:
                    assert ((equation&assignment).bit_count()&1)==((equation>>variables)&1)
    census_sources={name:digest(ROOT/'tools'/name) for name in
                    ('cube_binary_fiber_census.cpp','cube_bivariate_smith.py')}
    if args.recompute_census:
        census=smith.census(reconstructed,args.census)
        census.update(input_sha256=digest(args.input),source_sha256=census_sources)
        args.census.write_text(json.dumps(census,indent=2,sort_keys=True)+'\n')
    else:
        census=json.loads(args.census.read_text())
        assert census['input_sha256']==digest(args.input)
        assert census['source_sha256']==census_sources
    assert census['dimension']==22 and census['cases']==1<<22
    assert census['cutwise_weakly_majorized']==1<<22
    assert census['complete_binary_space'] and census['completed_requested_cases']
    assert census['exception_classes']==[]
    controls=census['profile_controls']
    assert len(controls)==95
    independent_orders=0
    for control in controls:
        value=control['assignment']
        assert 0<=value<1<<22
        vector=list(reconstructed['base'])
        for bit,direction in enumerate(reconstructed['directions']):
            if value>>bit&1:vector=[a^b for a,b in zip(vector,direction)]
        orders=[]
        for axis in range(3):
            current=[]
            for parity in range(2):
                offset=64*(2*axis+parity)
                matrix=[vector[offset+8*i:offset+8*i+8] for i in range(8)]
                current+=smith.specialized_profile_python(matrix)
            current.sort()
            assert dominated(current)
            orders+=current
        assert orders==control['orders']
        independent_orders+=len(orders)
    full_assignment=0
    for assignment in assignments:full_assignment^=assignment
    representatives=[0,assignments[0],assignments[7],assignments[-1],full_assignment]
    generic=[]
    for value in representatives:
        tensor=fibers.literal(seed,directions,value)
        for axis in range(3):
            record=smith.profile(tensor,axis)
            certificate=smith.verify_profile(tensor,axis,record)
            assert dominated(record['orders'])
            generic.append(dict(assignment=hex(value),axis=axis,orders=record['orders'],certificate=certificate))
    names=('check_cube_coefficient_fibers.py','cube_coefficient_fibers.py',
           'cube_bivariate_smith.py','cube_bivariate_smith.cpp','cube_binary_fiber_census.cpp',
           'cube_homogeneous_gaussian.py')
    result=dict(status='passed',factoring_exponent_improved=False,
                coefficient_variables=variables,loss_three_dimension=0,loss_four_dimension=22,
                specialization_image_dimension=22,complete_binary_census_cases=census['cases'],
                cutwise_weakly_majorized=census['cutwise_weakly_majorized'],
                independently_replayed_specialized_orders=independent_orders,
                endpoint_coefficient_controls=coefficient_controls,exact_generic_profile_controls=generic,
                exact_generic_orders=16*len(generic),
                polynomial_row_certificate_entries=sum(g['certificate']['polynomial_product_entries'] for g in generic),
                full_census_recomputed=args.recompute_census,
                input_sha256=digest(args.input),census_sha256=digest(args.census),
                source_sha256={name:digest(ROOT/'tools'/name) for name in names},
                scope='Fixed Gaussian row-floor fibers; uniqueness after loss three over every characteristic-two field; complete loss-four census only over F2.')
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({key:result[key] for key in
          ('status','coefficient_variables','loss_three_dimension','loss_four_dimension',
           'complete_binary_census_cases','independently_replayed_specialized_orders',
           'exact_generic_orders','polynomial_row_certificate_entries')}))

if __name__=='__main__':main()
