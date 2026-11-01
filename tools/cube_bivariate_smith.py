#!/usr/bin/env python3
"""Exact generic profiles and polynomial row frames for binary cube tensors.

The supported entries are homogeneous of degree 3*occupation and have even
physical parity. Two 8-by-8 blocks suffice. Their minors have total degree
at most 96, so stride 256 separates both coordinates throughout Bareiss
products, whose degree is at most 192. No finite-field specialization is
used by the generic profile routine.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tempfile

import cube_homogeneous_gaussian as cube


def native_program(name):
    source=Path(__file__).with_name(name+'.cpp')
    compiler=shlex.split(os.environ.get('CXX','c++'))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError('The cube profile tools require a C++17 compiler.')
    flags=['-O3','-std=c++17']
    try: cpu=Path('/proc/cpuinfo').read_text()
    except OSError: cpu=''
    if platform.machine().lower() in ('x86_64','amd64','i386','i686') and 'pclmulqdq' in cpu:
        flags+=['-mpclmul','-msse2']
    digest=hashlib.sha256(source.read_bytes()+repr((compiler,flags,platform.machine())).encode()).hexdigest()[:24]
    cache=Path(tempfile.gettempdir())/('integer-factoring-cube-polynomial-'+str(os.getuid()))
    cache.mkdir(mode=0o700,parents=True,exist_ok=True)
    target=cache/(name+'-'+digest)
    if not target.exists():
        fd,temporary=tempfile.mkstemp(prefix=digest+'-',dir=cache)
        os.close(fd)
        try:
            result=subprocess.run(compiler+flags+[str(source),'-o',temporary],capture_output=True,text=True)
            if result.returncode:
                raise RuntimeError(result.stderr[-12000:])
            os.replace(temporary,target)
        finally:
            if os.path.exists(temporary):os.unlink(temporary)
    return target


def encode(polynomial,axis):
    other=next(i for i in range(3) if i!=axis)
    value=0
    for exponent in polynomial:
        value^=1<<(exponent[axis]+256*exponent[other])
    return value


def validate_tensor(tensor):
    if len(tensor)!=256:
        raise ValueError('Expected the complete eight-leg binary tensor.')
    for mask,polynomial in enumerate(tensor):
        if mask.bit_count()%2 and polynomial:
            raise ValueError('The parity-block interface requires even physical parity.')
        if any(len(e)!=3 or min(e)<0 or sum(e)!=3*mask.bit_count() for e in polynomial):
            raise ValueError('Entry degree must be exactly three times occupation.')
        endpoint=all((mask>>i&1)==(mask>>(i^7)&1) for i in range(4))
        if len(polynomial)%2!=endpoint:
            raise ValueError('The literal antipodal crossing endpoint failed.')


def profile(tensor,axis):
    validate_tensor(tensor)
    if axis not in range(3):raise ValueError('Unknown cube axis.')
    matrix=cube.axis_matrix(tensor,axis)
    blocks=[];orders=[]
    for parity in range(2):
        labels=[i for i in range(16) if i.bit_count()%2==parity]
        block=[[encode(matrix[i][j],axis) for j in labels] for i in labels]
        data='8\n'+'\n'.join(format(value,'x') for row in block for value in row)+'\n'
        result=subprocess.run([str(native_program('cube_bivariate_smith'))],input=data,
                              text=True,capture_output=True)
        if result.returncode:
            raise ArithmeticError(result.stderr)
        record=json.loads(result.stdout)
        record['labels']=labels
        record['frames']=[[int(value,16) for value in row] for row in record['frames']]
        blocks.append(record);orders+=record['orders']
    return dict(orders=sorted(orders),blocks=blocks)


def polynomial_multiply(a,b):
    if a.bit_count()<b.bit_count():a,b=b,a
    result=0
    while b:
        bit=b&-b;b^=bit
        result^=a<<(bit.bit_length()-1)
    return result


def polynomial_divide_exact(a,b):
    if not b:raise ZeroDivisionError
    result=0;degree=b.bit_length()
    while a.bit_length()>=degree:
        shift=a.bit_length()-degree
        result^=1<<shift;a^=b<<shift
    if a:raise ArithmeticError('Polynomial division was not exact.')
    return result


def polynomial_determinant(matrix):
    matrix=[list(row) for row in matrix]
    n=len(matrix);previous=1
    if not n:return 1
    for k in range(n-1):
        candidates=[i for i in range(k,n) if matrix[i][k]]
        if not candidates:return 0
        at=candidates[0];matrix[k],matrix[at]=matrix[at],matrix[k]
        pivot=matrix[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=(polynomial_multiply(pivot,matrix[i][j])^
                           polynomial_multiply(matrix[i][k],matrix[k][j]))
                matrix[i][j]=polynomial_divide_exact(numerator,previous)
            matrix[i][k]=0
        previous=pivot
    return matrix[-1][-1]


def coordinate_terms(value):
    while value:
        bit=value&-value;value^=bit
        exponent=bit.bit_length()-1
        yield exponent%256,exponent//256


def x_coefficient(value,power):
    result=0
    for x,y in coordinate_terms(value):
        if x==power:result^=1<<y
    return result


def verify_profile(tensor,axis,record):
    """Check the exact row certificate without repeating native pivot choices."""
    validate_tensor(tensor)
    matrix=cube.axis_matrix(tensor,axis)
    entries=0;determinants=[];orders=[]
    assert len(record['blocks'])==2
    for parity,block in enumerate(record['blocks']):
        labels=block['labels'];frame=block['frames'];s=block['orders']
        assert labels==[i for i in range(16) if i.bit_count()%2==parity]
        source=[[encode(matrix[i][j],axis) for j in labels] for i in labels]
        assert len(frame)==8 and all(len(row)==8 for row in frame)
        assert len(s)==8 and s==sorted(s) and min(s)>=0
        assert all(x+y<=96 for row in frame for value in row for x,y in coordinate_terms(value))
        assert all(x+y<=24 for row in source for value in row for x,y in coordinate_terms(value))
        product=[[0]*8 for _ in range(8)]
        for i in range(8):
            for j in range(8):
                for k in range(8):
                    product[i][j]^=polynomial_multiply(frame[i][k],source[k][j])
                assert all(x>=s[i] for x,y in coordinate_terms(product[i][j]))
                entries+=1
        constant=[[x_coefficient(value,0) for value in row] for row in frame]
        quotient=[[x_coefficient(value,s[i]) for value in row] for i,row in enumerate(product)]
        first=polynomial_determinant(constant)
        second=polynomial_determinant(quotient)
        assert first and second
        determinants.append([format(first,'x'),format(second,'x')]);orders+=s
    assert sorted(orders)==record['orders']
    return dict(polynomial_product_entries=entries,unit_determinants=determinants,
                complete_orders=len(orders))


def specialized_profile_python(matrix,precision=128):
    """Independent binary series elimination for finite census controls."""
    matrix=[list(row) for row in matrix];n=len(matrix);mask=(1<<precision)-1
    orders=[]
    def order(value):return (value&-value).bit_length()-1 if value else precision
    def divide_unit(numerator,denominator):
        assert denominator&1
        result=0;numerator&=mask
        while numerator:
            bit=numerator&-numerator;result^=bit
            numerator^=denominator<<(bit.bit_length()-1);numerator&=mask
        return result
    for k in range(n):
        value,i,j=min((order(matrix[i][j]),i,j) for i in range(k,n) for j in range(k,n))
        if value==precision:raise ArithmeticError('Unresolved promised pivot.')
        matrix[k],matrix[i]=matrix[i],matrix[k]
        for row in matrix:row[k],row[j]=row[j],row[k]
        orders.append(value);unit=matrix[k][k]>>value
        for i in range(k+1,n):
            if not matrix[i][k]:continue
            quotient=divide_unit(matrix[i][k]>>value,unit)
            for j in range(k+1,n):
                matrix[i][j]^=polynomial_multiply(quotient,matrix[k][j])&mask
            matrix[i][k]=0
    return orders


def census(input_record,output,max_cases=None,seconds=1800):
    dimension=input_record['specialization_image_dimension']
    assert dimension==input_record['free_variables']==len(input_record['directions'])
    assert 0<dimension<=24
    arrays=[input_record['base']]+input_record['directions']
    assert all(len(row)==384 and all(0<=value<1<<25 for value in row) for row in arrays)
    with tempfile.TemporaryDirectory(prefix='integer-factoring-cube-census-') as temporary:
        source=Path(temporary)/'input.txt'
        source.write_text(str(dimension)+'\n'+'\n'.join(' '.join(map(str,row)) for row in arrays)+'\n')
        command=[str(native_program('cube_binary_fiber_census')),str(source),str(output),
                 str((1<<dimension) if max_cases is None else max_cases),str(seconds)]
        result=subprocess.run(command,text=True,capture_output=True)
        if result.returncode:raise ArithmeticError(result.stderr)
    return json.loads(Path(output).read_text())
