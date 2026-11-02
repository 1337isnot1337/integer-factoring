from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import sys
sys.path.insert(0,str(Path.cwd()/'tools'))
from cyclotomic_character_decoder import (binary_cyclotomic, polynomial_divmod,
    polynomial_remainder, modular_power, modular_product)
from sparse_faithful_character import primitive_component, is_prime
from odd_character_residues import primitive_generator, root_of_unity_powers

records=[]
for level in (47,503):
    u=(level-1)//2
    remaining=binary_cyclotomic(u)
    fields=[]
    while remaining!=1:
        field=primitive_component(remaining)
        fields.append(field)
        remaining,rem=polynomial_divmod(remaining,field)
        assert rem==0
    assert len(fields)==(2 if level==47 else 5)
    fullmod,powers=root_of_unity_powers(level-1)
    for k in (1,3):
        assert gcd(k,level-1)==1
        g=pow(primitive_generator(level),k,level)
        logs={pow(g,j,level):j for j in range(level-1)}
        for left,right in ((2,3),(3,5),(3,97),(7,101),(101,113),(3,1009),(1009,1013)):
            assert is_prime(left) and is_prime(right) and left!=right
            target=left*right
            # Independent integer cyclotomic fixture, reduced after summation.
            exact=[sum(powers[logs[d%level]][j] for d in (1,left,right,target))
                   for j in range(len(fullmod)-1)]
            full=sum((v&1)<<j for j,v in enumerate(exact))
            for field in fields:
                observed=polynomial_remainder(full,field)
                c=modular_power(2,logs[target%level]%u,field)
                s=observed^1^c
                roots=[]
                for j in range(u):
                    z=modular_power(2,j,field)
                    if modular_product(z,z,field)^modular_product(s,z,field)^c==0:
                        roots.append(j)
                assert len(roots) in (1,2)
                repeated=len(roots)==1
                if repeated:roots*=2
                squares=sorted(pow(g,2*j,level) for j in roots)
                assert squares==sorted((left*left%level,right*right%level))
                records.append([level,field,k,target,repeated,squares])
result=dict(scope='Scoped ordinary inspection of PR101 head35b0d74; independent integer-cyclotomic fixtures across every field at q47/q503 and two supplied primitive generators. No coefficient engine or arbitrary-input theorem.',controls=len(records),equal_roots=sum(r[4] for r in records),record_sha256=sha256(json.dumps(records).encode()).hexdigest(),records=records,source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
Path('/tmp/sol2-review-semiprime-field.json').write_text(json.dumps(result,indent=2)+'\n')
print({k:result[k] for k in ('controls','equal_roots','record_sha256','source_sha256')})
