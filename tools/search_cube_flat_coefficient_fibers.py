#!/usr/bin/env python3
"""Bounded exact-frame searches in the complete degree-six cube coefficients.

The search is heuristic. Retained tensors have exact generic profiles;
reported decimal pressures are numerical and do not prove a score.
"""
import argparse
import hashlib
import json
import math
import random
import time
from pathlib import Path
import cube_flat_coefficient_fibers as fibers
from cube_flat_seeds import rank_two_seed, odd_edge_face_seed

def run(args):
    seed=(fibers.parse_tensor(json.loads(args.reference_file.read_text())) if args.reference_file
          else (rank_two_seed() if args.reference=='rank-two' else odd_edge_face_seed()))
    started=time.monotonic();rng=random.Random(args.seed)
    system=fibers.build_system(seed,args.field,args.window)
    variables=system['variables'];mask=(1<<variables)-1;reference=system['reference']
    best=fibers.pressure([r['orders'] for r in system['frames']])['rho']
    witnesses=[];history=[];classes=set();attempts=completed=0
    second_order=[]
    if args.second_order:
        for first in range(1,6):
            for second in range(1,7-first):
                third=7-first-second
                if third<1:continue
                q=[first,second,third];targets=sum(([0]+[v]*15 for v in q),[])
                pivots=fibers.solve(system,targets)
                item=dict(second_orders=q,feasible=pivots is not None,
                          dimension=None if pivots is None else variables-len(pivots))
                second_order.append(item);print(json.dumps(dict(event='second_order',**item)),flush=True)
    def save(finished=False):
        sources=fibers.source_hashes()
        sources[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        record=dict(schema='flat-coefficient-search-v1',status='completed bounded search' if finished else 'finite search prefix',
                    field=args.field,random_seed=args.seed,extra_order_window=args.window,requested_seconds=args.seconds,
                    requested_runs=args.runs,elapsed_seconds=time.monotonic()-started,completed_runs=completed,
                    increment_attempts=attempts,distinct_target_classes=len(classes),variables=variables,
                    reference=fibers.serialize_tensor(seed),reference_sha256=fibers.tensor_digest(seed),
                    reference_orders=[r['orders'] for r in system['frames']],best_numerical_rho=best,
                    second_order_feasibility=second_order,retained_improvements=witnesses,history=history,sources=sources,
                    scope='Fixed generic row frames; sampled coefficient assignments; no exhaustive search or factoring exponent claim.')
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    save()
    for run in range(args.runs):
        if time.monotonic()-started>=args.seconds:break
        loss=rng.choice((2,3,3,4,4,5,6));targets=[max(0,s-loss) for s in reference]
        if run%3==1:targets=[max(0,s-rng.randrange(7)) for s in reference]
        pivots=fibers.solve(system,targets);assert pivots is not None
        blocked=set()
        while len(blocked)<48 and time.monotonic()-started<args.seconds:
            at=fibers.pressure([targets[a:a+16] for a in (0,16,32)])['u']
            at=1.0 if at is None else at
            options=[]
            for row in range(48):
                if row in blocked:continue
                if targets[row]>=system['caps'][row]:blocked.add(row);continue
                values=targets[row//16*16:row//16*16+16];minimum=min(values)
                partition=sum(math.exp(-at*(v-minimum)) for v in values)
                priority=-at*(targets[row]-minimum)-math.log(partition)+rng.gauss(0,rng.choice((.1,.4,1.,2.)))
                options.append((priority,row))
            if not options:break
            row=max(options)[1];proposed=dict(pivots);attempts+=1
            if all(fibers.insert_equation(proposed,value,mask) for value in system['equations'][row][targets[row]]):
                targets[row]+=1;pivots=proposed
            else:blocked.add(row)
        key=tuple(targets);dimension=variables-len(pivots)
        if key not in classes:
            classes.add(key);free=[i for i in range(variables) if i not in pivots]
            free_value=0 if not free or run%2==0 else sum(1<<i for i in rng.sample(free,min(len(free),rng.choice((1,2,4,8)))))
            assignment=fibers.complete_assignment(pivots,free_value,variables)
            tensor=fibers.tensor_from_assignment(seed,system['directions'],assignment)
            orders=([r['orders'] for r in fibers.profiles(tensor)] if assignment
                    else [r['orders'] for r in system['frames']])
            numerical=fibers.pressure(orders);floor=fibers.pressure([targets[a:a+16] for a in (0,16,32)])
            assert numerical['rho']<=floor['rho']+1e-10
            history.append(dict(run=run,dimension=dimension,numerical_pressure=numerical['rho'],floor_pressure=floor['rho']))
            if numerical['rho']<best-1e-12:
                products=fibers.verify_row_floors(tensor,system,targets)
                best=numerical['rho'];witnesses.append(dict(run=run,assignment=hex(assignment),targets=targets,
                    tensor=fibers.serialize_tensor(tensor),tensor_sha256=fibers.tensor_digest(tensor),orders=orders,
                    numerical_pressure=numerical,polynomial_product_entries=products['polynomial_product_entries']))
                print(json.dumps(dict(event='improvement',run=run,rho=best)),flush=True)
        completed=run+1
        if run%10==0:
            save();print(json.dumps(dict(event='progress',runs=completed,classes=len(classes),best=best)),flush=True)
    save(True)
    print(json.dumps(dict(event='finished',runs=completed,classes=len(classes),best=best)),flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference',choices=('rank-two','odd-edge-face'),default='rank-two')
    parser.add_argument('--reference-file',type=Path)
    parser.add_argument('--field',type=int,choices=(2,4),default=2)
    parser.add_argument('--window',type=int,default=6)
    parser.add_argument('--seed',type=int,default=2026110113)
    parser.add_argument('--seconds',type=int,default=600)
    parser.add_argument('--runs',type=int,default=10000)
    parser.add_argument('--second-order',action='store_true')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.seconds<0 or args.runs<0 or args.window<5:parser.error('Nonnegative budgets and window at least five are required.')
    run(args)

if __name__=='__main__':main()
