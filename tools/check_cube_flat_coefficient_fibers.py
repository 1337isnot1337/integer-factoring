#!/usr/bin/env python3
"""Replay complete flat coefficient spaces, fixed-frame ranks, and witnesses."""
import argparse
import hashlib
import json
import random
from pathlib import Path
import cube_flat_coefficient_fibers as fibers
from cube_flat_seeds import rank_two_seed, odd_edge_face_seed

def replay_search(path):
    record=json.loads(path.read_text())
    assert record['schema']=='flat-coefficient-search-v1'
    for name,digest in record['sources'].items():
        assert Path(name).name==name
        assert hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()==digest
    seed=fibers.parse_tensor(record['reference'])
    assert fibers.tensor_digest(seed)==record['reference_sha256']
    system=fibers.build_system(seed,record['field'],record['extra_order_window'])
    assert record['variables']==system['variables']
    assert record['reference_orders']==[r['orders'] for r in system['frames']]
    feasible=0
    for item in record['second_order_feasibility']:
        targets=sum(([0]+[q]*15 for q in item['second_orders']),[])
        pivots=fibers.solve(system,targets)
        assert item['feasible']==(pivots is not None)
        assert item['dimension']==(None if pivots is None else system['variables']-len(pivots))
        feasible+=int(pivots is not None)
    best=fibers.pressure(record['reference_orders'])['rho'];products=profiles=0
    for witness in record['retained_improvements']:
        assignment=int(witness['assignment'],16)
        tensor=fibers.tensor_from_assignment(seed,system['directions'],assignment)
        assert tensor==fibers.parse_tensor(witness['tensor'])
        assert fibers.tensor_digest(tensor)==witness['tensor_sha256']
        verified=fibers.verify_row_floors(tensor,system,witness['targets'])
        assert verified['polynomial_product_entries']==witness['polynomial_product_entries']
        products+=verified['polynomial_product_entries']
        orders=[r['orders'] for r in fibers.profiles(tensor)];profiles+=3
        assert orders==witness['orders']
        numerical=fibers.pressure(orders)
        assert abs(numerical['rho']-witness['numerical_pressure']['rho'])<1e-12
        assert numerical['rho']<best;best=numerical['rho']
    assert abs(best-record['best_numerical_rho'])<1e-12
    return dict(record=str(path),scope='Exact retained witnesses and recorded feasibility systems; search trajectory counts are not replayed.',
                coefficient_variables=system['variables'],second_order_targets=len(record['second_order_feasibility']),
                feasible_second_order_targets=feasible,retained_improvements=len(record['retained_improvements']),
                polynomial_products=products,complete_generic_profiles=profiles,best_numerical_rho=best)

def run():
    rows=[];products=0;profile_checks=0
    references=(('rank_two_F2',rank_two_seed(),2,[0,0,0,72,284]),
                ('odd_edge_face_F2',odd_edge_face_seed(),2,[0,0,0,0,60]),
                ('rank_two_F4',rank_two_seed(),4,[0,0,0,144,568]))
    for name,seed,field,expected in references:
        system=fibers.build_system(seed,field);dimensions=[];retained=None
        for loss in range(5):
            target=[max(0,s-loss) for s in system['reference']]
            pivots=fibers.solve(system,target);assert pivots is not None
            dimensions.append(system['variables']-len(pivots))
            if loss==4:retained=target,pivots
        assert dimensions==expected,(name,dimensions)
        target,pivots=retained;free=[i for i in range(system['variables']) if i not in pivots]
        rng=random.Random(2026110105+field+int(name.startswith('odd')))
        free_value=sum(1<<i for i in rng.sample(free,min(4,len(free))))
        value=fibers.complete_assignment(pivots,free_value,system['variables'])
        tensor=fibers.tensor_from_assignment(seed,system['directions'],value)
        checked=fibers.verify_row_floors(tensor,system,target);products+=checked['polynomial_product_entries']
        actual=fibers.profiles(tensor);profile_checks+=3
        orders=[r['orders'] for r in actual]
        floor=fibers.pressure([target[a:a+16] for a in (0,16,32)])
        result=fibers.pressure(orders);assert result['rho']<=floor['rho']+1e-10
        assert fibers.parse_tensor(fibers.serialize_tensor(tensor))==tensor
        rows.append(dict(name=name,variables=system['variables'],dimensions_at_losses_0_through_4=dimensions,
                         reference_orders=[r['orders'] for r in system['frames']],
                         assignment=hex(value),tensor_sha256=fibers.tensor_digest(tensor),orders=orders,
                         numerical_pressure=result,polynomial_product_entries=checked['polynomial_product_entries']))
        print('PASS:',name,dimensions,flush=True)
    orders=rows[0]['reference_orders'];approx=fibers.pressure(orders)
    assert abs(approx['rho']-0.678171780875)<1e-10
    base=round(10000*approx['base'])
    certificate=fibers.certify_pressure(orders,679,1000,base,10000)
    assert certificate['passed'] and not certificate['upper_bound_below_cube_threshold']
    assert not fibers.certify_pressure(orders,577,1000,base,10000)['passed']
    return dict(scope='Complete degree-six endpoint coefficient spaces in three fixed generic row frames; no exhaustive search of the fibers.',
                sources=fibers.source_hashes(),controls=rows,independent_polynomial_products=products,
                independent_complete_generic_profiles=profile_checks,reference_pressure_certificate=certificate,
                achieved_factoring_score='1/2; unchanged')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('/tmp/cube-flat-coefficient-fibers.json'))
    parser.add_argument('--compare',type=Path)
    parser.add_argument('--search',type=Path,action='append')
    args=parser.parse_args();result=([replay_search(p) for p in args.search] if args.search else run())
    if args.compare:assert result==json.loads(args.compare.read_text())
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('PASS: complete flat coefficient spaces and exact polynomial witnesses')

if __name__=='__main__':main()
