"""Acquire direct synchronized row polynomials for normal-basis bitwise AND.

This is a representation diagnostic, not a factoring or contraction API.
The acquired rows are ADDED. There is no free sum-to-product compilation.
"""


def binary_rank(values):
    pivots={}
    for value in values:
        while value:
            lead=value.bit_length()-1
            if lead in pivots:value^=pivots[lead]
            else:
                pivots[lead]=value
                break
    return len(pivots)


def trace(field,width,value):
    answer=0
    for i in range(width):answer^=field.power(value,1 << i)
    if answer not in (0,1):raise ValueError('Expected an absolute binary trace')
    return answer


def normal_basis_data(field,width,beta):
    basis=[field.power(beta,1 << i) for i in range(width)]
    if binary_rank(basis)!=width:raise ValueError('Generator is not normal')
    equations=[sum(trace(field,width,field.mul(1 << k,basis[j])) << k
                   for k in range(width)) | (int(j==0) << width)
               for j in range(width)]
    for k in range(width):
        pivot=next(j for j in range(k,width) if (equations[j] >> k)&1)
        equations[k],equations[pivot]=equations[pivot],equations[k]
        for j in range(width):
            if j!=k and (equations[j] >> k)&1:equations[j]^=equations[k]
    gamma=sum(((equations[k] >> width)&1) << k for k in range(width))
    dual=[field.power(gamma,1 << i) for i in range(width)]
    assert all(trace(field,width,field.mul(x,y))==int(i==j)
               for i,x in enumerate(basis) for j,y in enumerate(dual))
    coefficients=[[trace(field,width,field.mul(beta,field.mul(dual[i],dual[j])))
                   for j in range(width)] for i in range(width)]
    offsets={k:next((i for i in range(width) if coefficients[i][(i+k)%width]),None)
             for k in range(width)}
    assert all(i is not None for i in offsets.values())
    return dict(beta=beta,gamma=gamma,basis=basis,dual_basis=dual,
                coefficients=coefficients,offset_witnesses=offsets)


def synchronized_and_rows(width,t,coefficients):
    """Matching degree 2^floor(t/2)+1, allowing positive exponent wraparound."""
    if width<1 or t<1 or width%t:raise ValueError('Row modulus must divide width')
    columns=width//t
    rows=[[] for _ in range(t)]
    for i in range(width):
        for j in range(width):
            if not coefficients[i][j]:continue
            offset=(i-j)%t
            if offset<=t//2:
                a,b=offset,0
                row=j%t
            else:
                a,b=0,t-offset
                row=i%t
            assert (i-row-a)%t==(j-row-b)%t==0
            xcolumn=((i-row-a)//t)%columns
            ycolumn=((j-row-b)//t)%columns
            counts=[0]*(2*columns)
            counts[xcolumn]=1 << a
            counts[columns+ycolumn]=1 << b
            rows[row].append((coefficients[i][j],tuple(counts)))
    degree=max([0]+[sum(counts) for terms in rows for c,counts in terms if c])
    minimum=(1 << (t//2))+1
    assert degree==minimum
    return rows,dict(row_modulus=t,relative_columns=columns,
                     maximum_formal_degree=degree,proved_minimum_degree=minimum)


def evaluate_rows(field,width,t,rows,x,y):
    columns=width//t
    answer=0
    for row,terms in enumerate(rows):
        for c,counts in terms:
            a=sum(counts[j] << (row+t*j) for j in range(columns))
            b=sum(counts[columns+j] << (row+t*j) for j in range(columns))
            answer^=field.mul(c,field.mul(field.power(x,a),field.power(y,b)))
    return answer
