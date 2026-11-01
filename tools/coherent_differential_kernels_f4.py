#!/usr/bin/env python3
"""Homogeneous F4 differential kernels with explicit deletion-closed retracts.

This is the field-extension version of coherent_differential_kernels.py.
Sparse vectors and operators map occupation masks to nonzero F4 coefficients.
"""
from collections import defaultdict
from coherent_hot_compounds_f4 import mul,axpy,basis,height,validate_maps


def apply_operator(vector,operator):
    out={}
    for word,c in vector.items():
        for removed,a in operator.items():
            if word&removed==removed:
                target=word^removed;out[target]=out.get(target,0)^mul(c,a)
                if not out[target]:del out[target]
    return out


def differential_kernel(operators,modes,height_cap):
    if not 1<=modes<=16 or not operators or type(height_cap)is not int or height_cap<0:raise ValueError('Invalid kernel dimensions.')
    for op in operators:
        if not op or any(type(w)is not int or not 0<=w<1<<modes or c not in (1,2,3) for w,c in op.items()):raise ValueError('Invalid F4 differential operator.')
        if len({height(w) for w in op})!=1 or len({w.bit_count() for w in op})!=1:raise ValueError('Operator must be homogeneous in height and occupation.')
    cells=defaultdict(list)
    for word in range(1<<modes):
        if height(word)<=height_cap:cells[height(word),word.bit_count()].append(word)
    generators=[];ambient=0;image_ranks=[]
    for cell,words in sorted(cells.items()):
        pivots={};ambient+=len(words)
        for word in words:
            pre={word:1};image={}
            for oi,op in enumerate(operators):
                for target,c in apply_operator(pre,op).items():image[(oi,target)]=c
            while image:
                k=max(image);factor=image[k]
                if k in pivots:
                    old_image,old_pre=pivots[k];image=axpy(image,factor,old_image);pre=axpy(pre,factor,old_pre)
                else:
                    inv=mul(factor,factor);image={k:mul(inv,v) for k,v in image.items()};pre={k:mul(inv,v) for k,v in pre.items()};pivots[k]=(image,pre);break
            if not image:generators.append(pre)
        image_ranks.append({'height':cell[0],'charge':cell[1],'ambient':len(words),'image_rank':len(pivots)})
    pivots=basis(generators)
    for pivot in sorted(pivots):
        for other in pivots:
            if other!=pivot and pivot in pivots[other]:pivots[other]=axpy(pivots[other],pivots[other][pivot],pivots[pivot])
    words=sorted(pivots,key=lambda w:(height(w),w));p=[pivots[w] for w in words];l=[{w:1} for w in words]
    hs,charges=validate_maps(p,l,modes)
    for v in p:
        for op in operators:assert not apply_operator(v,op)
    assert len(p)==ambient-sum(c['image_rank'] for c in image_ranks)
    return {'modes':modes,'height_cap':height_cap,'operators':operators,'ambient_dimension':ambient,'q':len(p),'embedding':p,'retraction':l,'pivot_words':words,'heights':hs,'charges':charges,'cell_image_ranks':image_ranks}
