#!/usr/bin/env python3
"""Search complete-entry four-symbol crossing fibers with exact constraints."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path
import random
import time

import numpy as np
from scipy.optimize import brentq,minimize_scalar
from scipy.special import logsumexp

from full_alphabet_two_cut_lattices import FullAlphabetTwoCutLattice,CompressedFullAlphabetFiber
from reflection_dual_rail_gaussian import tensor


def root_for_profiles(profiles,energy_x,energy_y):
    first,second=map(lambda values:np.asarray(values,dtype=float),profiles)
    ex=np.asarray(energy_x,dtype=float);ey=np.asarray(energy_y,dtype=float)
    ex_min=min(ex);ey_min=min(ey);ex=ex-ex_min;ey=ey-ey_min
    def root(log_scale):
        scale=math.exp(log_scale)
        constant=logsumexp(-scale*first)+logsumexp(-scale*second)+scale*(ex_min+ey_min)
        if constant<1e-12:return 1.
        def pressure(r):return constant-r*(logsumexp(-scale*ex/r)+logsumexp(-scale*ey/r))
        upper=2.
        while pressure(upper)>0 and upper<64:upper*=2
        if pressure(upper)>0:return 32.
        return brentq(pressure,1e-8,upper,xtol=2e-11)/2
    grid=np.linspace(-8.,1.,25);values=[root(point) for point in grid]
    candidates=[(value,float(point)) for point,value in zip(grid,values)]
    for index in sorted(range(len(grid)),key=lambda i:values[i])[:3]:
        lo=grid[max(0,index-1)];hi=grid[min(len(grid)-1,index+1)]
        if lo==hi:continue
        result=minimize_scalar(root,bounds=(lo,hi),method='bounded',options={'xatol':1e-9})
        candidates.append((float(result.fun),float(result.x)))
    value,log_scale=min(candidates);scale=math.exp(log_scale)
    critical=logsumexp(-scale*first)+logsumexp(-scale*second)
    critical-=math.sqrt(2)*(logsumexp(-scale*np.asarray(energy_x)/math.sqrt(2))+
                            logsumexp(-scale*np.asarray(energy_y)/math.sqrt(2)))
    return {'rho':value,'scale':scale,'critical_pressure':float(critical)}


def increment(fiber,pivots,row,target):
    result=dict(pivots)
    for power,value,rhs in fiber.equations[row]:
        if power<target:continue
        if power>target:break
        while value:
            pivot=value.bit_length()-1
            if pivot not in result:result[pivot]=(value,rhs);break
            old,constant=result[pivot];value^=old;rhs^=constant
        if not value and rhs:return None
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--energies',default='0,12,4,16')
    parser.add_argument('--lower-window',type=int,default=6)
    parser.add_argument('--upper-window',type=int,default=8)
    parser.add_argument('--trials',type=int,default=20000)
    parser.add_argument('--seconds',type=float,default=3600)
    parser.add_argument('--samples',type=int,default=3)
    parser.add_argument('--seed',type=int,default=311739)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.monotonic();deadline=started+args.seconds
    energies=list(map(int,args.energies.split(',')))
    if len(energies)!=4 or any(value%2 for value in energies):raise ValueError('Use four even integer paired energies.')
    source=tensor(4);degree=max(a.bit_length()-1 for a,b in source['values'])
    half=[value//2 for value in energies];zero=[0]*4
    heights=[zero,zero,half,half,half,half,zero,zero]
    lattice=FullAlphabetTwoCutLattice(source['values'],4,degree,heights,raise_cap=args.upper_window)
    floors=tuple(max(0,value-args.lower_window) for value in lattice.orders)
    fiber=CompressedFullAlphabetFiber(lattice,floors)
    reference=lattice.profiles(lattice.seed(0));reference_cost=root_for_profiles(reference,energies,energies)
    print(json.dumps({'event':'constructed','seconds':time.monotonic()-started,'variables':lattice.variables,
        'fiber_variables':fiber.variables,'reference':reference_cost}),flush=True)
    rng=random.Random(args.seed);classes={};parents=[lattice.orders];seen={lattice.orders}
    retained=[];best=reference_cost['rho'];queries=full_replays=0;trials=0
    source_names=[Path(__file__).name,'full_alphabet_two_cut_lattices.py','reflection_dual_rail_gaussian.py',
                  'weyl_lattice_filter_synthesis.py','native_matrix_smith_f4.py','native_matrix_smith_f4.cpp']
    result={'status':'Exact complete-entry affine constraints with exploratory root optimization; no all-family claim.',
        'arguments':vars(args)|{'output':str(args.output)},'degree':degree,'physical_alphabet':4,
        'heights_infinity_then_zero':heights,'variables':lattice.variables,'fiber_variables':fiber.variables,
        'reference_profiles':reference,'reference_cost':reference_cost,'row_target_floors':floors,
        'row_target_ceilings':lattice.caps,'reference_row_targets':lattice.orders,
        'source_sha256':{name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in source_names}}
    def save(completed=False):
        result.update(completed_requested_trials=completed,completed_trials=trials,
            distinct_target_vectors=len(seen),complete_profile_queries=queries,
            independent_full_tensor_replays=full_replays,distinct_profile_classes=len(classes),
            classes=[{'profiles':key,'samples':value} for key,value in sorted(classes.items())],
            retained_witnesses=sorted(retained,key=lambda row:row['rho'])[:40],
            elapsed_seconds=time.monotonic()-started,best_recorded_rho=best,
            strict_crossing_candidate_found=any(row['critical_pressure']<-1e-8 for row in retained))
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    for trial in range(args.trials):
        if time.monotonic()>=deadline:break
        trials=trial+1
        if trial==0:targets=list(lattice.orders)
        else:
            targets=list(rng.choice(parents))
            for index in rng.sample(range(32),rng.randrange(2,17)):
                targets[index]=max(floors[index],targets[index]-rng.randrange(1,args.lower_window+1))
        pivots=fiber.solve(targets);assert pivots is not None
        order=list(range(32));rng.shuffle(order)
        if trial%4==1:order.sort(key=lambda i:targets[i]+rng.random()*args.lower_window)
        for row in order:
            while targets[row]<lattice.caps[row]:
                proposal=increment(fiber,pivots,row,targets[row])
                if proposal is None:break
                pivots=proposal;targets[row]+=1
        target=tuple(targets)
        if trial and target in seen:continue
        seen.add(target)
        if len(parents)<256:parents.append(target)
        elif rng.randrange(4)==0:parents[rng.randrange(1,len(parents))]=target
        for sample in range(args.samples):
            values=fiber.sample(pivots,rng);profiles=lattice.profiles(values);queries+=1
            if profiles in classes:classes[profiles]+=1;continue
            classes[profiles]=1;cost=root_for_profiles(profiles,energies,energies)
            witness={'trial':trial,'sample':sample,'targets':target,'profiles':profiles,'values':values,
                     'fiber_solution_dimension':fiber.variables-len(pivots),**cost}
            if len(retained)<40 or cost['rho']<max(row['rho'] for row in retained):
                retained.append(witness);retained=sorted(retained,key=lambda row:row['rho'])[:40]
            if cost['rho']<best-1e-10 or full_replays<4:
                assert lattice.check_witness(values,target)==profiles;full_replays+=1
                if cost['rho']<best-1e-10:
                    best=cost['rho'];print(json.dumps({'event':'new_best','trial':trial,**cost}),flush=True)
        if trial%25==0:
            save();print(json.dumps({'event':'progress','trial':trial,'distinct_targets':len(seen),
                'classes':len(classes),'best_rho':best,'seconds':time.monotonic()-started}),flush=True)
    save(trials==args.trials)
    print(json.dumps({'event':'finished','trials':trials,'classes':len(classes),'best_rho':best,
        'completed_requested_trials':result['completed_requested_trials']}),flush=True)


if __name__=='__main__':main()
