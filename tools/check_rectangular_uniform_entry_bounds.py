#!/usr/bin/env python3
"""Certify uniform rectangular entry-slope bounds by exact formal-Z arithmetic.

The small variable z is specialized to omega in F4 while Z=z^M remains
formal. The proof of specialization safety is in the companion note.
Only standard-library exact arithmetic is used; no floating solver or
native Smith oracle is called.
"""
from __future__ import annotations
import argparse
import json
import sys
from hashlib import sha256
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from rectangular_geometric_core import concentrated,Rat,rat,mono,fpow,MUL,INV,Z0,degree,scalar,mul4
import rectangular_geometric_core as geometric
from rectangular_defect_profiles import Program

INF=1000000
TERMS=80


@dataclass(frozen=True)
class Jet:
 p:tuple
 base:int
 known:int
 def exact_zero(self):return not any(self.p)and self.known==INF
 def order(self):return self.base


ZERO=Jet((0,0),INF,INF)
ONE=Jet((1,0),0,INF)


def normalized(p,base,known):
 length=known-base
 if length<=0:return Jet((0,0),known,known)
 if length<max(p).bit_length():p=(p[0]&((1<<length)-1),p[1]&((1<<length)-1))
 v=p[0]|p[1]
 if not v:return ZERO if known==INF else Jet((0,0),known,known)
 shift=(v&-v).bit_length()-1
 return Jet((p[0]>>shift,p[1]>>shift),base+shift,known)
@lru_cache(maxsize=200000)
def add(a,b):
 if a.exact_zero():return b
 if b.exact_zero():return a
 base=min(a.base,b.base);known=min(a.known,b.known)
 return normalized((a.p[0]<<(a.base-base)^b.p[0]<<(b.base-base),a.p[1]<<(a.base-base)^b.p[1]<<(b.base-base)),base,known)
@lru_cache(maxsize=200000)
def multiply(a,b):
 if a.exact_zero()or b.exact_zero():return ZERO
 known=min(a.known+b.base,b.known+a.base,INF)
 return normalized(mul4(a.p,b.p),a.base+b.base,known)


def constant(c,exponent=0):return ZERO if not c else Jet((c&1,c>>1),exponent,INF)


@lru_cache(maxsize=20000)
def series(value):
 if not value:return ZERO
 vn=((value.n[0]|value.n[1])&-(value.n[0]|value.n[1])).bit_length()-1
 vd=((value.d[0]|value.d[1])&-(value.d[0]|value.d[1])).bit_length()-1
 numerator=(value.n[0]>>vn,value.n[1]>>vn);denominator=(value.d[0]>>vd,value.d[1]>>vd)
 def coefficient(p,k):return((p[0]>>k)&1)+2*((p[1]>>k)&1)
 inverse=[0]*TERMS;inverse[0]=INV[coefficient(denominator,0)]
 for k in range(1,TERMS):
  value=0
  for j in range(1,min(k,degree(denominator))+1):value^=MUL[coefficient(denominator,j)][inverse[k-j]]
  inverse[k]=MUL[inverse[0]][value]
 p=(sum((v&1)<<k for k,v in enumerate(inverse)),sum((v>>1)<<k for k,v in enumerate(inverse)))
 return normalized(mul4(numerator,p),vn-vd,vn-vd+TERMS)


@dataclass(frozen=True)
class Affine:
 a:int=0
 b:int=0
 def __add__(self,other):
  other=affine(other);return Affine(self.a+other.a,self.b+other.b)
 __radd__=__add__
 def __neg__(self):return Affine(-self.a,-self.b)
 def __sub__(self,other):return self+-affine(other)
 def __rsub__(self,other):return affine(other)+-self
 def __mul__(self,other):return Affine(self.a*other,self.b*other)
 __rmul__=__mul__
 def __le__(self,other):
  value=self-affine(other);return value.a<0 or(value.a==0 and value.b<=0)
 def __lt__(self,other):
  value=self-affine(other);return value.a<0 or(value.a==0 and value.b<0)
 def __str__(self):return 'M-1'if(self.a,self.b)==(1,-1)else f'({self.a}M+{self.b})'


def affine(value):return value if isinstance(value,Affine)else Affine(0,value)


class SymbolicProgram:
 word=Program.word
 number=Program.number
 unit=Program.unit
 build=Program.build
 def __init__(self,data,images):
  self.data={'M':Affine(1,0),'place':data['place'],'full_core':True}
  self.instructions=[('identity',),('zero',)];self.memo={};self.images={name:[(i,True,series(v))for i,v in enumerate(c)if v]+[(i,False,series(v))for i,v in enumerate(a)if v]for name,(c,a)in images.items()}
  self.orders=data['orders'];self.diagonals=[]
  for mask in range(1<<len(self.orders)):
   slope=sum(a for i,(a,b)in enumerate(self.orders)if mask>>i&1);offset=sum(b for i,(a,b)in enumerate(self.orders)if mask>>i&1)
   self.diagonals.append(constant(fpow(geometric.Z0,offset),slope))
 def node(self,key):
  if key not in self.memo:self.memo[key]=len(self.instructions);self.instructions.append(key)
  return self.memo[key]
 def add(self,a,b):
  if a==1:return b
  if b==1:return a
  if a==b:return 1
  return self.node(('add',min(a,b),max(a,b)))
 def scale(self,a,coefficient=1,z=0,u=0):
  if a==1 or coefficient==0:return 1
  z=affine(z);c=MUL[coefficient][MUL[fpow(geometric.Z0,z.b)][fpow(1^geometric.Z0,u)]]
  factor=constant(c,z.a)
  if factor==ONE:return a
  return self.node(('scale',a,factor))
 def car(self,a,name):return 1 if a==1 else self.node(('car',a,name))
 def diagonal(self,a):return self.node(('diagonal',a))


def basis_names(place,side):
 letters=(('out_v','out_y')if place=='infinity'else('left_y','left_x'))if side=='left'else(('in_u','in_x')if place=='infinity'else('right_v','right_u'))
 return[prefix+index for prefix in letters for index in('0','M-1')]


def compile_program(data):
 images=dict(data['images']);names={}
 for side in('left','right'):
  originals=basis_names(data['place'],side);names[side]=['basis_'+side+'_'+name for name in originals]
  for original,name in zip(originals,names[side]):
   for suffix in('_create','_annihilate'):
    c,a=images[original+suffix]
    if side=='right':
     c=[v*mono(a=-s,b=-t)for v,(s,t)in zip(c,data['orders'])]
     a=[v*mono(a=s,b=t)for v,(s,t)in zip(a,data['orders'])]
    images[name+suffix]=(c,a)
 program=SymbolicProgram(data,images)
 def vacuum(value,group):
  for name in group:value=program.add(value,program.number(value,name))
  return value
 right_vacuum=vacuum(0,names['right']);right_states={0:right_vacuum};outputs={}
 for state in range(16):
  if state:
   bit=state&-state;index=bit.bit_length()-1
   right_states[state]=program.car(right_states[state^bit],names['right'][index]+'_create')
  initial=program.diagonal(right_states[state]);value=program.build(initial=initial)
  left_states={0:value}
  for mask in range(16):
   if mask:
    bit=mask&-mask;index=bit.bit_length()-1
    left_states[mask]=program.car(left_states[mask^bit],names['left'][index]+'_annihilate')
   outputs[mask,state]=vacuum(left_states[mask],names['left'])
 return program,outputs


def evaluator(program,column):
 @lru_cache(maxsize=500000)
 def entry(register,row):
  instruction=program.instructions[register];kind=instruction[0]
  if kind=='zero':return ZERO
  if kind=='identity':return ONE if row==column else ZERO
  if kind=='add':return add(entry(instruction[1],row),entry(instruction[2],row))
  if kind=='scale':return multiply(instruction[2],entry(instruction[1],row))
  if kind=='diagonal':return multiply(program.diagonals[row],entry(instruction[1],row))
  if kind=='car':
   total=ZERO
   for mode,create,coefficient in program.images[instruction[2]]:
    bit=1<<mode
    if bool(row&bit)==create:total=add(total,multiply(coefficient,entry(instruction[1],row^bit)))
   return total
  raise ValueError(kind)
 return entry


def weight_feature(place, masks, coordinates):
    """Recover all 24 physical weights directly from CAR homogeneity."""
    mask, state = masks
    row, column = coordinates
    infinity = place == 'infinity'
    if infinity:
        row_charge = row.bit_count()+mask.bit_count()
        column_charge = column.bit_count()+state.bit_count()
        v, y, u, x = mask & 3, mask >> 2, state & 3, state >> 2
    else:
        def charge(word):
            return (word & 63).bit_count()-(word >> 6).bit_count()
        row_charge = charge(row)-(mask & 3).bit_count()+(mask >> 2).bit_count()
        column_charge = charge(column)+(state & 3).bit_count()-(state >> 2).bit_count()
        y, x, v, u = mask & 3, mask >> 2, state & 3, state >> 2
    assert y != 3 and u != 3
    result = [0]*24
    offset = 0 if infinity else 2
    result[offset], result[offset+1] = row_charge, column_charge
    start = 4+(0 if infinity else 10)
    for physical_state, position in ((v, 0), (y, 3), (u, 5), (x, 7)):
        if physical_state:
            result[start+position+physical_state-1] += 1
    return result


def finite_controls():
    """Compare the closed formulas with the original M-step concentration."""
    from check_rectangular_full_core_support import full_core_data

    def poly_value(polynomial, point):
        value = 0
        for j in range(polynomial.bit_length()-1, -1, -1):
            value = MUL[value][point]^((polynomial >> j) & 1)
        return value

    def f4_value(polynomial, point):
        return poly_value(polynomial[0], point)^MUL[2][poly_value(polynomial[1], point)]

    def rational_value(numerator, denominator, point, extension=False):
        evaluate = f4_value if extension else poly_value
        divisor = evaluate(denominator, point)
        if not divisor:
            return None
        return MUL[evaluate(numerator, point)][INV[divisor]]

    records = []
    for place in ('infinity', 'zero'):
        for M in (16, 32, 64):
            original = full_core_data(place, M)
            for c in (2, 3):
                formula = concentrated(place, c)
                assert original['active_orders'] == [a*M+b for a, b in formula['orders']]
                checks = 0
                poles = 0
                for name, (creation, annihilation) in formula['images'].items():
                    source = original['active_images'][name.replace('M-1', str(M-1))]
                    for key, values in (('creation', creation), ('annihilation', annihilation)):
                        for value, (numerator, denominator) in zip(values, source[key]):
                            actual = rational_value(numerator, denominator, c)
                            expected = rational_value(value.n, value.d, fpow(c, M), True)
                            if actual is None or expected is None:
                                assert actual is None and expected is None, (place, M, c, name, key, 'pole mismatch')
                                poles += 1
                            else:
                                assert actual == expected, (place, M, c, name, key)
                                checks += 1
                records.append({'place': place, 'M': M, 'small_variable': c,
                                'components_checked': checks, 'matching_poles_skipped': poles})
    return records


def verify_supported_duals(model, witnesses):
    """Check lower-bound duals using only the symbolically supported rows."""
    from fractions import Fraction
    from check_rectangular_leading_face_lp import physical, add_states, dot

    assert model['full_core']
    supported = {(w['place'], tuple(w['feature']), w['target_slope']) for w in witnesses}
    used = set(model['used_entry_constraints'])
    assert len(supported) == len(used) == 28
    for i in used:
        record = model['constraints'][i]
        assert (record['place'], tuple(record['feature']), record['fitted_slope']) in supported
    xs, xb, ys, yb = physical(True)
    d = len(xs[0])
    rows, rhs = [], []
    for record in model['constraints']:
        row = [-value for value in record['feature']]+[0, 0]
        row[d+(record['place'] == 'zero')] = 1
        rows.append(row)
        rhs.append(record['fitted_slope'])
    equalities = [[0]*(d+2) for _ in range(2)]
    equalities[0][0] = equalities[0][1] = 1
    equalities[1][2] = equalities[1][3] = 1

    def verify(certificate, objective, inequalities, bounds, equations, targets, constant):
        dual = [Fraction() for _ in inequalities]
        for i, value in certificate['dual_inequalities']:
            assert not dual[i]
            dual[i] = Fraction(value)
            assert dual[i] <= 0
            if i < len(rows) and dual[i]:
                assert i in used
        multipliers = [Fraction() for _ in equations]
        for i, value in certificate['dual_equalities']:
            assert not multipliers[i]
            multipliers[i] = Fraction(value)
        for j, coefficient in enumerate(objective):
            assert sum(row[j]*value for row, value in zip(inequalities, dual)) \
                 + sum(row[j]*value for row, value in zip(equations, multipliers)) == coefficient
        value = dot(bounds, dual)+dot(targets, multipliers)+constant
        assert value == Fraction(certificate['optimum'])
        return str(value)

    table = [[None]*3 for _ in range(4)]
    for case in model['cases']:
        sx, sy = case['selected_x'], case['selected_y']
        inequalities, bounds = add_states(rows, rhs, xs, xb, ys, yb, sx, sy)
        objective = [a+b for a, b in zip(xs[sx], ys[sy])]+[-1, -1]
        table[sx][sy] = verify(case, objective, inequalities, bounds, equalities,
                               [0, 0], xb[sx]+yb[sy])
    inequalities, bounds = add_states(rows, rhs, xs, xb, ys, yb, 1, 1)
    c_objective = [a+b for a, b in zip(xs[1], ys[1])]+[-1, -1]
    face_equations = equalities+[c_objective]
    face_targets = [0, 0, -xb[1]-yb[1]]
    gaps = []
    for case in model['C0_physical_gaps']:
        group, offsets = (xs, xb) if case['letter'] == 'X' else (ys, yb)
        state = case['state']
        objective = [a-b for a, b in zip(group[state], group[1])]+[0, 0]
        value = verify(case, objective, inequalities, bounds, face_equations,
                       face_targets, offsets[state]-offsets[1])
        gaps.append({'letter': case['letter'], 'state': state, 'minimum_gap': value})
    return {'leading_C_lower_bounds': table, 'C0_physical_gap_lower_bounds': gaps,
            'exact_dual_certificates': len(model['cases'])+len(gaps)}


def main():
    sys.setrecursionlimit(10000)
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path,
                        default=root/'evidence/rectangular-uniform-entry-inputs.json')
    parser.add_argument('--lp-input', type=Path,
                        default=root/'evidence/rectangular-full-leading-face-lp.json')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    inputs = json.loads(args.input.read_text())
    assert inputs['seed'] == {'alpha': 1, 'beta': 2, 'quartic_terms': [[1, 0]],
                              'repair_kind': 'partial', 'full_core': True}
    witnesses = inputs['witnesses']
    assert len(witnesses) == 28
    records, concentration = [], []
    for place in ('infinity', 'zero'):
        data = concentrated(place, 2)
        assert all(item['value'] not in (0, 1) for item in data['geometric_ratios'])
        program, outputs = compile_program(data)
        concentration.append({key: data[key] for key in
                              ('place', 'active', 'orders', 'geometric_ratios', 'pivot_checks')})
        for witness in (w for w in witnesses if w['place'] == place):
            assert witness['specialization'] == 2
            assert weight_feature(place, witness['physical_core_states'],
                                  witness['vacuum_coordinates']) == witness['feature']
            mask, state = witness['physical_core_states']
            row, column = witness['vacuum_coordinates']
            entry = evaluator(program, column)
            value = entry(outputs[mask, state], row)
            coefficient = (value.p[0] & 1)+2*(value.p[1] & 1)
            assert coefficient and value.known > witness['target_slope']
            assert value.base <= witness['target_slope']
            assert value.base == witness['expected_leading_slope']
            assert coefficient == witness['expected_leading_coefficient']
            records.append({'place': place, 'entry': witness['entry'],
                            'target_slope': witness['target_slope'],
                            'specialized_leading_slope': value.base,
                            'specialized_leading_coefficient': coefficient,
                            'coefficients_known_below': value.known,
                            'physical_core_states': witness['physical_core_states'],
                            'vacuum_coordinates': witness['vacuum_coordinates'],
                            'weight_feature_verified': True})
            entry.cache_clear()
    model = json.loads(args.lp_input.read_text())
    duals = verify_supported_duals(model, witnesses)
    controls = finite_controls()
    source_paths = [Path(__file__), Path(geometric.__file__),
                    root/'tools/rectangular_defect_profiles.py',
                    root/'tools/check_rectangular_f4_repair.py']
    output = {'status': 'All 28 uniform entry-slope upper bounds and 17 supported exact LP duals pass.',
              'scope': 'Fixed balanced-charge asymptotic rays for the alpha=1, beta=omega partial rectangular repair.',
              'formal_variable': 'Z=z^M; M ranges over sufficiently large powers of two.',
              'small_variable_specialization': 'z=omega; Z remains formal.',
              'witnesses_checked': len(records), 'witnesses': records,
              'concentration': concentration, 'dual_lower_bounds': duals,
              'finite_implementation_controls': controls,
              'finite_components_checked': sum(case['components_checked'] for case in controls),
              'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
              'lp_input_sha256': sha256(args.lp_input.read_bytes()).hexdigest(),
              'source_sha256': {path.name: sha256(path.read_bytes()).hexdigest() for path in source_paths}}
    rendered = json.dumps(output, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
