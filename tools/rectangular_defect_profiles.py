#!/usr/bin/env python3
"""Compute the proved fixed-gauge rectangular defect profiles over F4.

Only unit denominators (polynomials with constant term one) are cleared.
Polynomial arithmetic uses an adaptive exact coefficient window; a common
precision bound is checked before the Smith calculation. Generic rank is supplied by the
all-length endpoint theorem, never inferred from truncated zeros.
"""
from __future__ import annotations

import argparse
import ctypes
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tempfile
import time

from check_rectangular_f4_repair import partial_repair,full_inverse_repair,matrix_multiply
from native_matrix_smith_f4 import MatrixOracle

ROOT=Path(__file__).resolve().parents[1]
MASK=(1<<64)-1

def mul2(a,b):
    out=0
    while b:
        bit=b&-b;out^=a<<(bit.bit_length()-1);b^=bit
    return out

def div2(a,b):
    if not b:raise ZeroDivisionError
    quotient=0
    while a and a.bit_length()>=b.bit_length():
        shift=a.bit_length()-b.bit_length();quotient^=1<<shift;a^=b<<shift
    return quotient,a

def quotient2(a,b):
    quotient,remainder=div2(a,b)
    if remainder:raise ArithmeticError('Inexact unit-denominator division.')
    return quotient

def gcd2(a,b):
    while b:a,b=b,div2(a,b)[1]
    return a

def lcm2(a,b):return mul2(quotient2(a,gcd2(a,b)),b)

def power_two_ceiling(n):return 1<<max(0,(n-1).bit_length())

def precision_sizes(program,base,M,profile_precision=None):
    window=profile_precision or max(512,power_two_ceiling(16*M))
    if window<128 or window>8192 or window&(window-1):raise ValueError('Invalid profile precision.')
    coefficient_bits=max(max(pair).bit_length() for pair in program.coefficients)
    builder=max(1024,power_two_ceiling(max(coefficient_bits,window+max(0,-base)+64)))
    if builder>32768:raise ArithmeticError('Requested polynomial builder precision is too large.')
    return window,builder

def upow(e):
    if e<0:raise ValueError('Negative polynomial unit power.')
    out=1
    while e:
        bit=e&-e;out^=out<<bit;e^=bit
    return out

def mul4(a,b):
    p=mul2(a[0],b[0]);q=mul2(a[1],b[1]);r=mul2(a[0]^a[1],b[0]^b[1])
    return p^q,r^p

def denominator_form(d):
    if not d:raise ValueError('Zero denominator.')
    z=(d&-d).bit_length()-1;unit=d>>z
    return z,unit

@dataclass(frozen=True)
class Expr:
    register:int
    base:int=0
    denominator:int=1

class Program:
    def __init__(self,data):
        self.data=data;self.k=len(data['active_indices'])
        self.instructions=[[0,0,0,0,0]];self.coefficients=[];self.coefficient_ids={}
        self.cars=[];self.car_meta={};self.memo={};self.degree_bounds=[0]
        self.coefficient((0,0));self.coefficient((1,0))
        for name,image in data['active_images'].items():
            terms=[]
            for creation,label in ((1,'creation'),(0,'annihilation')):
                for mode,(num,den) in enumerate(image[label]):
                    if num:
                        a,b=denominator_form(den);terms.append((mode,creation,num,-a,b))
            if not terms:raise ValueError('Zero CAR generator.')
            base=min(t[3] for t in terms);den=1
            for term in terms:den=lcm2(den,term[4])
            packed=[]
            for mode,creation,num,a,b in terms:
                p=mul2(num,quotient2(den,b))<<(a-base)
                packed.append((mode,creation,self.coefficient((p,0))))
            self.car_meta[name]=(len(self.cars),base,den);self.cars.append(packed)

    def coefficient(self,p):
        if max(p).bit_length()>32768:raise ArithmeticError('Program coefficient exceeds 32768 bits.')
        if p not in self.coefficient_ids:
            self.coefficient_ids[p]=len(self.coefficients);self.coefficients.append(p)
        return self.coefficient_ids[p]

    def add(self,a,b):
        base=min(a.base,b.base);den=lcm2(a.denominator,b.denominator)
        ca=quotient2(den,a.denominator)<<(a.base-base)
        cb=quotient2(den,b.denominator)<<(b.base-base)
        key=('add',a,b)
        if key not in self.memo:
            reg=len(self.instructions)
            self.instructions.append([2,a.register,b.register,self.coefficient((ca,0)),self.coefficient((cb,0))])
            self.degree_bounds.append(max(self.degree_bounds[a.register]+ca.bit_length()-1,self.degree_bounds[b.register]+cb.bit_length()-1))
            self.memo[key]=Expr(reg,base,den)
        return self.memo[key]

    def scale(self,a,coefficient=1,z=0,u=0):
        if coefficient not in (0,1,2,3):raise ValueError('Expected an F4 coefficient.')
        denominator=a.denominator;factor=1
        if u<0:denominator=mul2(denominator,upow(-u))
        elif u>0:
            factor=upow(u);common=gcd2(factor,denominator)
            factor=quotient2(factor,common);denominator=quotient2(denominator,common)
        out=Expr(a.register,a.base+z,denominator)
        if coefficient==1 and factor==1:return out
        scalar=(factor if coefficient&1 else 0,factor if coefficient&2 else 0)
        key=('scale',a.register,scalar)
        if key not in self.memo:
            reg=len(self.instructions)
            self.instructions.append([2,a.register,a.register,self.coefficient(scalar),0])
            self.degree_bounds.append(self.degree_bounds[a.register]+max(scalar).bit_length()-1)
            self.memo[key]=reg
        return Expr(self.memo[key],out.base,out.denominator)

    def car(self,a,name):
        op,base,den=self.car_meta[name];key=('car',a.register,op)
        if key not in self.memo:
            reg=len(self.instructions);self.instructions.append([1,a.register,0,op,0]);self.memo[key]=reg
            self.degree_bounds.append(self.degree_bounds[a.register]+max(max(self.coefficients[term[2]]).bit_length()-1 for term in self.cars[op]))
        return Expr(self.memo[key],a.base+base,mul2(a.denominator,den))

    def word(self,a,names):
        for name in reversed(names):a=self.car(a,name)
        return a

    def polynomial(self,a,factor):
        if factor==1:return a
        reg=len(self.instructions)
        self.instructions.append([2,a.register,a.register,self.coefficient((factor,0)),0])
        self.degree_bounds.append(self.degree_bounds[a.register]+factor.bit_length()-1)
        return Expr(reg,a.base,a.denominator)

    def diagonal(self,a,orders):
        low=min(orders);reg=len(self.instructions)
        self.instructions.append([3,a.register,0,0,0])
        self.degree_bounds.append(self.degree_bounds[a.register]+max(orders)-low)
        return Expr(reg,a.base+low,a.denominator)

    def number(self,a,name):
        return self.word(a,[name+'_create',name+'_annihilate'])

    def unit(self,a,names,row,column):
        for bit,name in enumerate(names):
            r=(row>>bit)&1;s=(column>>bit)&1
            if r==s:
                occupied=self.number(a,name);a=occupied if r else self.add(a,occupied)
            else:a=self.car(a,name+('_create' if r else '_annihilate'))
        return a

    def build(self,alpha=1,beta=2,quartic_terms=None,initial=None,repair_spec=None,input_spec=None,repair_kind='partial',identity_spec=None):
        data=self.data;M=data['M'];T=M-1;infinity=data['place']=='infinity'
        if quartic_terms is None:quartic_terms=[(1,0)]
        endpoint=0
        for coefficient,exponent in quartic_terms:
            if coefficient not in (0,1,2,3) or not isinstance(exponent,int):raise ValueError('Invalid quartic Laurent term.')
            endpoint^=coefficient
        if endpoint!=1:raise ValueError('The quartic numerator must equal one at w=1.')
        u0=('in_' if infinity else 'right_')+'u0';uT=('in_' if infinity else 'right_')+f'u{T}'
        y0=('out_' if infinity else 'left_')+'y0';yT=('out_' if infinity else 'left_')+f'y{T}'
        v0=('out_' if infinity else 'right_')+'v0';vT=('out_' if infinity else 'right_')+f'v{T}'
        if repair_kind not in ('partial','inverse'):raise ValueError('Unknown residual repair.')
        input_matrix=None
        if identity_spec is not None:
            if repair_spec is not None or input_spec is not None:raise ValueError('Identity transport and independent flag repairs are separate constructors.')
            if not data.get('full_core',False):raise ValueError('Identity transports require the full-core construction.')
            from check_rectangular_identity_transports import validate_identity_transport,validate_typed_transport
            if len(identity_spec)==2:correction,input_matrix=validate_identity_transport(alpha,beta,*identity_spec)
            elif len(identity_spec)==3:correction,input_matrix=validate_typed_transport(alpha,beta,*identity_spec)
            else:raise ValueError('Supply C,W and an optional certified typed endpoint.')
        else:correction=(partial_repair if repair_kind=='partial' else full_inverse_repair)(alpha,beta)
        if repair_spec is not None:
            from check_rectangular_refined_repairs import validate_postcorrection
            validate_postcorrection(*repair_spec)
            correction=matrix_multiply(repair_spec[0],correction)
        if input_spec is not None:
            if not data.get('full_core',False):raise ValueError('Input repairs require the full-core construction.')
            from check_rectangular_refined_repairs import validate_precorrection,validate_postcorrection
            validate_precorrection(*input_spec)
            input_matrix=input_spec[0]
            if repair_spec is not None:
                validate_postcorrection(repair_spec[0],input_spec[1],input_spec[2],False)
        core=(1,0,2,3);heightX=(0,0,1,1);heightY=(0,0,1)
        def middle(a):
            a=self.add(a,self.number(self.number(a,uT),u0))
            if alpha:
                b=self.car(self.number(a,uT),u0+'_create')
                a=self.add(a,self.scale(b,alpha,z=M-1 if infinity else 0))
            if infinity:word=[vT+'_create',yT+'_create',u0+'_annihilate','in_x0_annihilate']
            else:word=['left_x0_create',yT+'_create',u0+'_annihilate',vT+'_annihilate']
            quartic=self.word(a,word);parts=[]
            for coefficient,exponent in quartic_terms:
                if coefficient:parts.append(self.scale(quartic,coefficient,z=2-exponent if infinity else exponent,u=-2))
            for part in parts:a=self.add(a,part)
            if beta:
                b=self.car(self.number(a,yT),y0+'_annihilate')
                a=self.add(a,self.scale(b,beta,z=M if infinity else -1))
            return a
        def physical(v,y):
            if infinity:return (1-M)*bool(core[v]&1)+(M-1)*bool(core[v]&2)-M*(core[y].bit_count())
            return bool(core[y]&1)+(2*M-1)*bool(core[y]&2)
        if initial is None:initial=Expr(0)
        def output(source):
            base=middle(source) if infinity else None;parts=[]
            for row,line in enumerate(correction):
                v,y=row%4,row//4
                for column,coefficient in enumerate(line):
                    if not coefficient:continue
                    vp,yp=column%4,column//4
                    if infinity:
                        a=self.unit(base,[v0,vT],core[v],core[vp])
                        a=self.unit(a,[y0,yT],core[y],core[yp])
                    else:
                        # Right multiplication transposes this matrix unit.
                        a=self.unit(source,[v0,vT],core[vp],core[v])
                        a=middle(a);a=self.unit(a,[y0,yT],core[y],core[yp])
                    dh=heightX[vp]+heightY[yp]-heightX[v]-heightY[y]
                    dz=physical(v,y)-physical(vp,yp)-(dh if infinity else 0)
                    parts.append(self.scale(a,coefficient,z=dz,u=dh))
            out=parts[0]
            for part in parts[1:]:out=self.add(out,part)
            return out
        out=output(initial)
        if input_matrix is not None:
            x0=('in_' if infinity else 'left_')+'x0';xT=('in_' if infinity else 'left_')+f'x{T}'
            def input_physical(x,u):
                if infinity:return M*core[x].bit_count()+(1-M)*bool(core[u]&1)+(M-1)*bool(core[u]&2)
                return (1-2*M)*bool(core[x]&1)-bool(core[x]&2)
            for row,line in enumerate(input_matrix):
                xa,ua=row%4,row//4
                for column,coefficient in enumerate(line):
                    coefficient^=int(row==column)
                    if not coefficient:continue
                    xb,ub=column%4,column//4
                    a=self.unit(initial,[u0,uT],core[ua],core[ub])
                    if infinity:
                        a=self.unit(a,[x0,xT],core[xa],core[xb]);a=output(a)
                    else:
                        a=output(a);a=self.unit(a,[x0,xT],core[xb],core[xa])
                    dh=heightX[xb]+heightY[ub]-heightX[xa]-heightY[ua]
                    dz=input_physical(xb,ub)-input_physical(xa,ua)-(dh if infinity else 0)
                    out=self.add(out,self.scale(a,coefficient,z=dz,u=dh))
        return out

def library(precision=1024):
    source=Path(__file__).with_name('native_fock_f4.cpp');included=source.with_name('native_rectangular_smith_f4.cpp')
    compiler=shlex.split(os.environ.get('CXX','c++'))
    if not compiler or shutil.which(compiler[0]) is None:raise RuntimeError('C++17 is required.')
    flags=['-O3','-std=c++17','-shared','-fPIC',f'-DF4_SERIES_WORDS={precision//64}']
    try:cpu=Path('/proc/cpuinfo').read_text()
    except OSError:cpu=''
    if platform.machine().lower() in ('x86_64','amd64') and 'pclmulqdq' in cpu:flags+=['-mpclmul','-msse2']
    digest=hashlib.sha256(source.read_bytes()+included.read_bytes()+repr((compiler,flags)).encode()).hexdigest()[:24]
    cache=Path(tempfile.gettempdir())/('integer-factoring-fock-f4-'+str(os.getuid()));cache.mkdir(mode=0o700,parents=True,exist_ok=True)
    target=cache/(digest+'.so')
    if not target.exists():
        fd,name=tempfile.mkstemp(prefix=digest+'-',suffix='.so',dir=cache);os.close(fd)
        try:
            run=subprocess.run(compiler+flags+[str(source),'-o',name],capture_output=True,text=True)
            if run.returncode:raise RuntimeError(run.stderr[-10000:])
            os.replace(name,target)
        finally:
            if os.path.exists(name):os.unlink(name)
    lib=ctypes.CDLL(str(target));I=ctypes.POINTER(ctypes.c_int);U=ctypes.POINTER(ctypes.c_uint64);L=ctypes.POINTER(ctypes.c_int64)
    lib.fock_f4_evaluate_v2.argtypes=[ctypes.c_int,ctypes.c_int,I,ctypes.c_int,U,ctypes.c_int,I,I,ctypes.c_int,I,ctypes.c_int,ctypes.c_int,U,L,L]
    lib.fock_f4_evaluate_v2.restype=ctypes.c_int;return lib

def active_matrix(data,alpha=1,beta=2,quartic_terms=None,repair_spec=None,profile_precision=None,input_spec=None,repair_kind='partial',identity_spec=None):
    program=Program(data);last=program.build(alpha,beta,quartic_terms,repair_spec=repair_spec,input_spec=input_spec,repair_kind=repair_kind,identity_spec=identity_spec);n=1<<program.k
    diagonal=[sum(s for j,s in enumerate(data['active_orders']) if mask>>j&1) for mask in range(n)]
    window,builder=precision_sizes(program,last.base+min(diagonal),data['M'],profile_precision);words=builder//64;output_words=window//64
    instructions=(ctypes.c_int*(5*len(program.instructions)))(*(v for row in program.instructions for v in row))
    packed=(ctypes.c_uint64*(2*words*len(program.coefficients)))(*(value>>(64*w)&MASK for pair in program.coefficients for value in pair for w in range(words)))
    terms=[];offsets=[0]
    for car in program.cars:terms.extend(v for term in car for v in term);offsets.append(len(terms)//3)
    cc=(ctypes.c_int*len(terms))(*terms);co=(ctypes.c_int*len(offsets))(*offsets)
    dd=(ctypes.c_int*n)(*diagonal);output=(ctypes.c_uint64*(2*output_words*n*n))();shifts=(ctypes.c_int64*(n*n))();meta=(ctypes.c_int64*4)()
    start=time.monotonic();status=library(builder).fock_f4_evaluate_v2(program.k,len(program.instructions),instructions,len(program.coefficients),packed,len(program.cars),co,cc,last.register,dd,last.base,output_words,output,shifts,meta)
    if status:raise ArithmeticError(f'Fock program failed or lacked common precision (status {status}, metadata {list(meta)}).')
    matrix=[];shift_matrix=[]
    for i in range(n):
        row=[];shift_row=[]
        for j in range(n):
            index=i*n+j;row.append(tuple(sum(int(output[(2*index+plane)*output_words+w])<<(64*w) for w in range(output_words)) for plane in range(2)))
            shift_row.append(int(shifts[index]))
        matrix.append(row);shift_matrix.append(shift_row)
    exact=program.degree_bounds[last.register]<builder and int(meta[3])<window
    return matrix,shift_matrix,{'registers':len(program.instructions),'coefficient_count':len(program.coefficients),'unit_denominator_degree':last.denominator.bit_length()-1,'unit_denominator_hex':hex(last.denominator),'program_degree_bound':program.degree_bounds[last.register],'maximum_normalized_degree':int(meta[3]),'exact_polynomial_matrix':exact,'minimum_valuation':int(meta[0]),'known_below':int(meta[1]),'builder_precision':int(meta[2]),'profile_precision':window,'matrix_seconds':time.monotonic()-start}

def profile(data,alpha=1,beta=2,quartic_terms=None):
    matrix,shifts,metadata=active_matrix(data,alpha,beta,quartic_terms);n=len(matrix);rank=3*n//4
    start=time.monotonic();result=MatrixOracle(matrix,rank,shifts,precision=metadata['profile_precision']).profile()
    if not metadata['exact_polynomial_matrix'] and max(result['orders'])>=metadata['known_below']:raise ArithmeticError('Profile exceeds known polynomial precision.')
    orders=[]
    for value in result['orders']:
        if value!=round(value):raise ArithmeticError('Noninteger unweighted profile.')
        orders.append(round(value))
    return {'place':data['place'],'M':data['M'],'alpha':alpha,'beta':beta,'quartic_terms':quartic_terms or [(1,0)],'matrix_order':n,'proved_generic_rank':rank,'orders':orders,'spectator_orders':data['spectator_orders'],'profile_seconds':time.monotonic()-start,**metadata}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('data',type=Path);parser.add_argument('--alpha',type=int,default=1);parser.add_argument('--beta',type=int,default=2);parser.add_argument('--quartic-power',type=int,default=0);parser.add_argument('--output',type=Path)
    args=parser.parse_args();data=json.loads(args.data.read_text());result=profile(data,args.alpha,args.beta,[(1,args.quartic_power)])
    result['data_sha256']=hashlib.sha256(args.data.read_bytes()).hexdigest();rendered=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(rendered)
    print(rendered)

if __name__=='__main__':main()
