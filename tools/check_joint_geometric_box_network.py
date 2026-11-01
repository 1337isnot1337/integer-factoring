"""Direct finite-field oracle for global, coupled geometric-box networks."""
import itertools,json,time
from pathlib import Path
from hashlib import sha256
from joint_geometric_box_network import network_sum
from check_frobenius_rational_trace import field_of_degree


def main():
    records=[];start=time.monotonic();evaluations=0
    for width,left,right,edges in [(1,1,1,[(0,0,1)]),(2,1,1,[(0,0,1)]),(3,1,1,[(0,0,3)]),(2,2,1,[(0,0,1),(1,0,2)]),(2,1,2,[(0,0,2),(0,1,3)]),(2,2,2,[(0,0,1),(0,1,2),(1,0,3),(1,1,4)])]:
        f=field_of_degree(2*width);b=1<<width
        U=[u for u in range(1,f.size)if f.power(u,b)==u];V=[v for v in range(1,f.size)if f.power(v,b+1)==1]
        omega=next(a for a in range(f.size)if f.power(a,b)!=a);bar=f.power(omega,b)
        image={1}
        for z in [0]+U:image.add(f.mul(z^omega,f.power(z^bar,f.size-2)))
        assert image==set(V)
        left_charges=[b-2]*left;right_charges=[b]*right
        expected=0
        for us in itertools.product(U,repeat=left):
            for vs in itertools.product(V,repeat=right):
                term=1
                for u,a in zip(us,left_charges):term=f.mul(term,f.power(u,a))
                for v,a in zip(vs,right_charges):term=f.mul(term,f.power(v,a))
                for i,j,tag in edges:
                    x=f.mul(tag,f.mul(us[i],vs[j]));value=0;power=1
                    for _ in range(b):value^=power;power=f.mul(power,x)
                    term=f.mul(term,value)
                expected^=term;evaluations+=1
        actual=network_sum(f,width,left,right,edges,left_charges,right_charges);assert actual==expected,(width,left,right,actual,expected)
        records.append({'width':width,'left':left,'right':right,'edges':edges,'sum':actual,'left_charges':left_charges,'right_charges':right_charges,'direct_assignments':len(U)**left*len(V)**right})
        print(json.dumps(records[-1]),flush=True)
    large_controls=[]
    for width in (4,8):
        f=field_of_degree(2*width);b=1<<width
        actual=network_sum(f,width,1,1,[(0,0,3)],[b-2],[b])
        assert actual==3
        large_controls.append({'width':width,'sum':actual,'oracle':'Unique surviving edge digit1; no field assignments enumerated'})
    names=['joint_frobenius_carries.py','joint_geometric_box_network.py','check_joint_geometric_box_network.py','check_frobenius_rational_trace.py','check_carry_fourier_elimination.py']
    out={'completed':True,'records':records,'direct_assignments':evaluations,'symbolic_large_width_controls':large_controls,'source_sha256':{s:sha256(Path('tools',s).read_bytes()).hexdigest()for s in names},'scope':'Whole bipartite geometric-box products. This is not a representation of the complete factoring query or a stronger factoring bound.'}
    Path('evidence/joint-geometric-box-networks.json').write_text(json.dumps(out,indent=2)+'\n');print('passed',time.monotonic()-start,flush=True)
if __name__=='__main__':main()
