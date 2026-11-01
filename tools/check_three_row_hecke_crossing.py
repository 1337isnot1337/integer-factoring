#!/usr/bin/env python3
"""Independent confluent-jet replay of finite three-row Hecke crossings."""
import argparse,hashlib,json,time
from pathlib import Path
from three_row_hecke_paths import paths,M,q,F,local_motif_certificate,check_relations
from three_row_hecke_crossing import coefficient_tensor,example_paths,actions
from binary_extension_smith import ExtensionMatrixOracle


def trimmed(p):
    p=list(p)
    while p and not p[-1]:p.pop()
    return p


def remainder(a,b):
    a=trimmed(a);b=trimmed(b);inv=F.power(b[-1],6)
    while len(a)>=len(b):
        shift=len(a)-len(b);c=M[a[-1]][inv]
        for j,v in enumerate(b):a[shift+j]^=M[c][v]
        a=trimmed(a)
    return a


def common_factor(T):
    g=[]
    for p in T.reshape(len(T),-1).T:
        b=trimmed(p.tolist())
        if b:
            a=g
            while b:a,b=b,remainder(a,b)
            inv=F.power(a[-1],6);g=[M[inv][x] for x in a]
        if g==[1]:break
    return g


def jet_replay(T,ps,scalar,progress=None):
    import numpy as np
    length=len(ps[0]);m=length//3;D=length*length;J=(D+1+7)//8;N=len(ps)
    qs=paths((2*m,2*m,2*m),4);ix={p:i for i,p in enumerate(qs)};boundary=np.array([ix[p+r] for p in ps for r in ps]);ops=actions(qs)
    MM=np.array(M,dtype=np.uint8);word=[i for j in range(length) for i in range(length+j-1,j-1,-1)];records=[]
    for point in range(8):
        states=np.zeros((len(qs),N*N,J),dtype=np.uint8);states[boundary,np.arange(N*N),0]=1
        for site in word:
            diag,src,dst=ops[site]
            E=MM[diag[:,None,None],states]
            E[dst]^=MM[diag[src,None,None],states[src]]
            following=MM[q,states]^MM[point,E];following[:,:,1:]^=E[:,:,:-1];states=following
        actual=MM[F.power(scalar,6),states[boundary].transpose(2,1,0).reshape(J,N,N,N,N)]
        expected=np.zeros_like(actual)
        for j in range(J):
            for k in range(j,D+1):
                if j&k==j:expected[j]^=MM[F.power(point,k-j),T[k]]
        assert np.array_equal(actual,expected),(N,point)
        records.append({'point':point,'jet_count':J,'scalar_identities':J*N**4,'jet_sha256':hashlib.sha256(actual.tobytes()).hexdigest()})
        if progress:progress({'point':point,'jets':J,'identities':J*N**4})
    assert 8*J>D
    return {'degree_ceiling':D,'distinct_points':8,'jets_per_point':J,'total_multiplicity':8*J,'coefficient_determination_complete':True,'scalar_identities':8*J*N**4,'records':records}


def main():
    import numpy as np
    ap=argparse.ArgumentParser();ap.add_argument('--alphabet',type=int,choices=(5,9),default=5);ap.add_argument('--output',type=Path);args=ap.parse_args();start=time.monotonic()
    motifs=local_motif_certificate();relations=check_relations((4,4,4),4)
    ps=example_paths(args.alphabet);T,construction=coefficient_tensor(ps,progress=lambda v:print(json.dumps({'construction':v}),flush=True));construction.pop('seconds')
    assert common_factor(T)==[1]
    jets=jet_replay(T,ps,construction['braid_scalar'],progress=lambda v:print(json.dumps({'jet_replay':v}),flush=True))
    N=len(ps);D=len(T)-1;profiles=[]
    for cut,C in enumerate((T.reshape(D+1,N*N,N*N),T.transpose(0,1,3,2,4).reshape(D+1,N*N,N*N))):
        matrix=[[C[:,i,j].tolist() for j in range(N*N)] for i in range(N*N)]
        oracle=ExtensionMatrixOracle(matrix,11,N*N);zero=oracle.profile();infinity=oracle.profile(True);regular=[D+v for v in infinity['orders']];assert min(regular)>=0
        profiles.append({'cut':cut,'zero':zero,'infinity':infinity,'regular_infinity_orders':regular})
    floors=[]
    for a,b in ((0,1),(1,0)):
        zero=profiles[a]['zero']['orders'];inf=profiles[b]['regular_infinity_orders'];assert max(zero)+max(inf)<=D
        floors.append({'zero_cut':a,'infinity_cut':b,'maximum_order_sum':max(zero)+max(inf),'degree':D,'exact_bare_pressure_infimum':1})
    names=[Path(__file__).name,'three_row_hecke_paths.py','three_row_hecke_crossing.py','binary_extension_smith.py','check_fourier_multiplication_gates.py','native_matrix_smith_f4.py','native_matrix_smith_f4.cpp','native_rectangular_smith_f4.cpp']
    report={'status':'passed','local_motifs':motifs,'global_twelve_site_relations':relations,'construction':construction,'common_scalar_polynomial':[1],'confluent_jet_replay':jets,'profiles':profiles,'bare_profile_exclusions':floors,'source_sha256':{n:hashlib.sha256(Path('tools',n).read_bytes()).hexdigest() for n in names},'scope':'Uniform local Hecke relations, finite checked block-swap endpoints, complete local profiles and exact bare zero/infinity failure. No uniform endpoint theorem or improved factoring bound.'}
    output=args.output or Path('evidence/three-row-hecke-q'+str(N)+'-replay.json');output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','alphabet':N,'degree':D,'endpoint_entries':N**4,'confluent_coefficient_identities':jets['scalar_identities'],'bare_profile_infimum':1,'seconds':time.monotonic()-start}),flush=True)

if __name__=='__main__':main()
