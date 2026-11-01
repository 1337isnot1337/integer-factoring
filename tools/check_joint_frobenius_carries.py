"""Independent finite controls for the joint carry theorem and row synchronization."""
import itertools,json,random
from pathlib import Path
from hashlib import sha256
from joint_frobenius_carries import coefficient,joint_trace,synchronization
from check_frobenius_rational_trace import field_of_degree
from check_carry_fourier_elimination import BinaryField


def ordinary(field,variables,t,factors):
    out={(0,)*variables:1}
    for row,terms in factors:
        nxt={}
        for old,a in out.items():
            for b,h in terms:
                D=len(h)//variables;e=tuple(old[v]+sum(h[v*D+j]<<(row+t*j)for j in range(D))for v in range(variables))
                nxt[e]=nxt.get(e,0)^field.mul(a,b)
        out={e:c for e,c in nxt.items()if c}
    return out


def evaluate(field,variables,t,factor,point):
    row,terms=factor;out=0
    for c,h in terms:
        D=len(h)//variables;value=c
        for v,u in enumerate(point):value=field.mul(value,field.power(u,sum(h[v*D+j]<<(row+t*j)for j in range(D))))
        out^=value
    return out


def main():
    rng=random.Random(917733);counts={'coefficient_queries':0,'joint_trace_queries':0,'direct_assignments':0,'zero_denominators':0,'synchronization_graphs':0}
    field=field_of_degree(3)
    for variables in (1,2,3):
        for t in (1,2,3):
            for trial in range(8):
                D=1+(trial%2);factors=[]
                for _ in range(1+trial%3):
                    terms=[]
                    for _ in range(3):
                        h=[0]*(variables*D)
                        for _ in range(rng.randrange(3)):h[rng.randrange(len(h))]+=1
                        terms.append((rng.randrange(field.size),tuple(h)))
                    factors.append((rng.randrange(t),terms))
                expanded=ordinary(field,variables,t,factors)
                targets=list(expanded)[:8]+[tuple(rng.randrange(10)for _ in range(variables))for _ in range(3)]
                for target in targets:
                    assert coefficient(field,variables,t,factors,target)==expanded.get(target,0),(variables,t,factors,target)
                    counts['coefficient_queries']+=1
    for width,t,extension in [(1,1,1),(2,2,2),(2,1,2),(3,3,3),(2,2,4)]:
        field=BinaryField(3)if extension==1 else field_of_degree(extension);D=width//t;variables=2
        domain=[u for u in range(field.size)if field.power(u,1<<width)==u]
        for trial in range(4):
            factors=[]
            for _ in range(2):
                h=[0]*(variables*D);h[0]=1;h[D+(trial%D)]=1
                terms=[(rng.randrange(field.size),(0,)*(variables*D)),(1,tuple(h))]
                factors.append((trial%t,terms))
            numerators=factors[:1];denominators=factors[1:]
            if trial==3:denominators=[(0,[(0,(0,)*(variables*D))])]
            for nonzero in ((False,False),(True,True),(False,True),(True,False)):
                expected=0
                for point in itertools.product(*([u for u in domain if u or not nz]for nz in nonzero)):
                    value=1
                    for f in numerators:value=field.mul(value,evaluate(field,variables,t,f,point))
                    for f in denominators:
                        z=evaluate(field,variables,t,f,point)
                        if not z:counts['zero_denominators']+=1
                        value=field.mul(value,field.power(z,field.size-2)if z else 0)
                    expected^=value;counts['direct_assignments']+=1
                actual=joint_trace(field,extension,width,t,variables,numerators,denominators,nonzero)
                assert actual==expected,(width,t,extension,trial,nonzero,actual,expected)
                counts['joint_trace_queries']+=1
    for trial in range(24):
        width=1+trial%3;variables=2;t=width;field=BinaryField(3)if width==1 else field_of_degree(width)
        factors=[]
        for index in range(3):
            terms=[]
            for _ in range(4):
                h=[0]*variables
                for _ in range(rng.randrange(3)):h[rng.randrange(variables)]+=1
                terms.append((rng.randrange(field.size),tuple(h)))
            factors.append((rng.randrange(t),terms))
        mask=(bool(trial&1),bool(trial&2));expected=0
        for point in itertools.product(*(range(1 if nz else 0,field.size)for nz in mask)):
            value=1
            for f in factors[:2]:value=field.mul(value,evaluate(field,variables,t,f,point))
            z=evaluate(field,variables,t,factors[2],point)
            if not z:counts['zero_denominators']+=1
            value=field.mul(value,field.power(z,field.size-2)if z else 0)
            expected^=value;counts['direct_assignments']+=1
        assert joint_trace(field,width,width,t,variables,factors[:2],factors[2:],mask)==expected
        counts['joint_trace_queries']+=1
    for width in range(1,7):
        for trial in range(12):
            variables=3;nf=3;edges=[(rng.randrange(nf),rng.randrange(variables),rng.randrange(width))for _ in range(trial)]
            t,shifts,rows,residuals=synchronization(width,variables,nf,edges)
            possible=[]
            for divisor in range(1,width+1):
                if width%divisor:continue
                for guess in itertools.product(range(divisor),repeat=variables):
                    sets=[{(e+guess[v])%divisor for f,v,e in edges if f==i}for i in range(nf)]
                    if all(len(s)<=1 for s in sets):possible.append(divisor);break
            assert t==max(possible)
            counts['synchronization_graphs']+=1
    field=field_of_degree(2)
    assert joint_trace(field,2,2,2,2,[],[],False)==0
    assert joint_trace(field,2,2,2,2,[],[],True)==1
    names=['joint_frobenius_carries.py','check_joint_frobenius_carries.py','check_frobenius_rational_trace.py','check_carry_fourier_elimination.py']
    result={'completed':True,'counts':counts,'source_sha256':{s:sha256(Path('tools',s).read_bytes()).hexdigest()for s in names},'scope':'Exact finite controls for synchronized multivariate products, including coupled terms and zero denominators. No representation of the complete factoring network or improved factoring bound is claimed.'}
    Path('evidence/joint-frobenius-carries.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
