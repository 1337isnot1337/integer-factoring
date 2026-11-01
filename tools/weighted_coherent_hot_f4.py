#!/usr/bin/env python3
"""Two-cut F4 coherent Smith oracle inside a deletion-operator weight cone.

Floating weighted pivots and pressure are discovery diagnostics. A bound needs
independent exact rational replay and a certified strict pressure inequality.
"""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.special import logsumexp,softmax
from scipy.optimize import brentq
from coherent_hot_compounds_f4 import profiles,validate_maps,dot
from native_matrix_smith_f4 import MatrixOracle


def deletion_edges(vectors,duals,modes):
    edges=[]
    for target,v in enumerate(vectors):
        for bit in range(modes):
            d={w^(1<<bit):c for w,c in v.items() if w>>bit&1}
            for source,l in enumerate(duals):
                if dot(l,d):edges.append((source,target,bit))
    return edges


class CoherentWeightedOracle:
    def __init__(self,p,l,modes):
        self.heights,self.charges=validate_maps(p,l,modes)
        self.q=len(p);self.offsets=np.arange(5)*self.q;self.height_count=4*self.q
        self.edges_p=deletion_edges(p,l,modes);self.edges_l=deletion_edges(l,p,modes)
        self.raw=profiles(p,l,modes,keep_blocks=True);self.blocks=[]
        charges=sorted({a+b for a in self.charges for b in self.charges})
        for cut in self.raw['cuts']:
            blocks=[]
            for charge,b in zip(charges,cut['blocks']):
                pairs=np.array([(x,y) for x in range(self.q) for y in range(self.q) if self.charges[x]+self.charges[y]==charge],dtype=int)
                precision=b['precision_attempts'][-1][0]
                blocks.append({'pairs':pairs,'matrix':b['matrix'],'precision':precision,'oracle':MatrixOracle(b['matrix'],len(pairs),precision=precision)})
            self.blocks.append(blocks)

    def cone(self):
        if hasattr(self,"_cone"):return self._cone
        # Both natural cuts, with mirrored physical weights. Every deletion
        # edge gives its transpose creation edge in the corresponding factor.
        union=sorted(set(self.edges_p)|set(self.edges_l));ri=[];ci=[];va=[];bounds=[]
        for first,second,aa,bb in [(0,1,self.edges_l,union),(2,3,union,self.edges_p)]:
            for s,t,a in aa:
                for u,v,b in bb:
                    row=len(bounds)
                    for column,value in [(first*self.q+t,1),(first*self.q+s,-1),(second*self.q+v,1),(second*self.q+u,-1)]:
                        ri.append(row);ci.append(column);va.append(value)
                    bounds.append(-2*(a&b))
        self._cone=csr_matrix((va,(ri,ci)),shape=(len(bounds),self.height_count)),np.array(bounds,dtype=float)
        return self._cone

    def log_profiles(self,g,scale):
        values=[];gradients=[];derivatives=[];all_orders=[]
        for blocks in self.blocks:
            orders=[];labels=[];bases=[]
            for b in blocks:
                pairs=b['pairs'];rw=g[pairs[:,0]]+g[self.q+pairs[:,1]];cw=g[2*self.q+pairs[:,1]]+g[3*self.q+pairs[:,0]]
                while True:
                    try:r=b['oracle'].profile(rw,cw,scale);break
                    except ArithmeticError:
                        if b['precision']>=8192:raise
                        b['precision']*=2;b['oracle']=MatrixOracle(b['matrix'],len(pairs),precision=b['precision'])
                orders.extend(r['orders']);bases.extend(np.diff(np.r_[0.,r['base_prefix']]))
                rr=pairs[np.asarray(r['rows'])];cc=pairs[np.asarray(r['cols'])]
                labels.extend(np.column_stack((rr[:,0],rr[:,1],cc[:,1],cc[:,0])))
            orders=np.array(orders);labels=np.array(labels);bases=np.array(bases)
            assert len(orders)==self.q*self.q
            value=float(logsumexp(-orders));prob=softmax(-orders);gradient=np.zeros(4*self.q)
            for leg in range(4):gradient[leg*self.q:(leg+1)*self.q]=-np.bincount(labels[:,leg],weights=prob,minlength=self.q)
            values.append(value);gradients.append(gradient);derivatives.append(-float(prob@bases));all_orders.append(sorted(orders.tolist()))
        return np.mean(values),np.mean(gradients,axis=0),float(np.mean(derivatives)),all_orders

    def evaluate(self,g,scale):
        g=np.asarray(g,dtype=float)
        if g.shape!=(4*self.q,) or not np.all(np.isfinite(g)) or not np.isfinite(scale) or scale<=0:raise ValueError('Invalid weights or scale.')
        A,b=self.cone()
        if np.min(A@g-b*scale)<-1e-10:raise ValueError('Weights violate the sufficient coherent integrality cone.')
        raw,dg,ds,orders=self.log_profiles(g,scale)
        h=np.array(self.heights,dtype=float);cx=g[:self.q]+g[3*self.q:]+scale*h;cy=g[self.q:2*self.q]+g[2*self.q:3*self.q]+scale*h
        equation=lambda rho:raw-rho*(float(logsumexp(-cx/rho))+float(logsumexp(-cy/rho)))
        rho=brentq(equation,1e-8,4,xtol=2e-13)
        px=softmax(-cx/rho);py=softmax(-cy/rho);entropy=-sum(float(x)*np.log(float(x)) for x in np.r_[px,py] if x>0)
        if entropy<1e-10:raise ArithmeticError('Degenerate physical type.')
        gradient=(dg+np.r_[px,py,py,px])/entropy
        derivative=(ds+float(px@h+py@h))/entropy
        return rho,gradient,derivative,{'rho':rho,'scale_t':scale,'heights':g.tolist(),'orders_t_scaled':orders,'type_x':px.tolist(),'type_y':py.tolist(),'entropy_sum':entropy}
