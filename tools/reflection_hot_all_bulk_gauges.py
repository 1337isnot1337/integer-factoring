#!/usr/bin/env python3
"""All nonnegative physical-leg bulk orders through two regular cross pivots."""
from __future__ import annotations

from functools import lru_cache

import numpy as np

from native_matrix_smith_f4 import MatrixOracle
from reflection_hot_alternate_coordinates import gaussian_coordinates as alternate_coordinates
from reflection_hot_compensated_gauges import exterior_product,retained_internal_factor
from reflection_hot_block_profiles import normalize_factor
from reflection_hot_leg_bulk_gauges import LegBulkHotStatePressure,BOUNDARIES
from reflection_hot_state_coordinates import gaussian_coordinates as original_coordinates
from reflection_hot_state_gauges import defect_block


def vacuum_block_totals(labels,indices,m):
    return tuple(sum(sum(1 for bit,original in enumerate(indices)
                         if label>>bit&1 and (original<m)==(block==0))for label in labels)
                 for block in range(2))


@lru_cache(maxsize=32)
def branch_coordinates(interval,place,ramification,core,pivot,swapped,profile_precision):
    number=0 if place=='zero'else 1;left,right=BOUNDARIES[number]
    if pivot not in ('I','T'):raise ValueError('Choose the I or T cross pivot.')
    while True:
        source=(original_coordinates(interval,place,ramification,core,profile_precision,True)
                if pivot=='I'else alternate_coordinates(interval,place,ramification,core,profile_precision))
        if source['metadata']['exact_polynomial_matrix']:break
        profile_precision*=2
        if profile_precision>8192:
            raise ArithmeticError('The complete local-unit-cleared source exceeds the native precision limit.')
    replacements={}
    if swapped:
        for first,second in (left,right):
            replacements['xyuv'[first]]='xyuv'[second];replacements['xyuv'[second]]='xyuv'[first]
    def rename(name):
        side,label=name.split('_',1);letter=label[0]
        return side+'_'+replacements.get(letter,letter)+label[1:]
    answer=dict(source)
    for side in ('left','right'):
        basis=source[side];out=dict(basis)
        out['names']=[rename(name)for name in basis['names']]
        out['labels']={replacements.get(letter,letter):values for letter,values in basis['labels'].items()}
        out['regular_images']={rename(name):value for name,value in basis['regular_images'].items()}
        answer[side]=out
    actual_left=left[::-1]if swapped else left;actual_right=right[::-1]if swapped else right
    legs={'left':actual_left[::-1]if pivot=='I'else actual_left,'right':actual_right}
    answer['data']={**source['data'],'bulk_branch':{'pivot':pivot,'swapped':swapped},'bulk_block_legs':legs}
    for side in ('left','right'):
        basis=answer[side];data=answer['data']
        if vacuum_block_totals(basis['vacuum_rows'],data['active_indices'],interval//2)!=vacuum_block_totals(
                basis['vacuum_columns'],data['active_indices'],interval//2):
            raise ArithmeticError('The saved vacuum minor is not balanced in the two diagonal-block occupations.')
    return answer


class AllBulkHotStatePressure(LegBulkHotStatePressure):
    def __init__(self,interval,factor,ramification=1,core=1,profile_precision=512):
        self.interval,self.ramification,self.core,self.q=interval,ramification,core,1<<core
        self.factor=normalize_factor(factor,core);self.profile_precision=profile_precision
        self.local_factors=[exterior_product(self.factor,retained_internal_factor(interval,place,ramification,core))
                            for place in ('zero','infinity')]
        self.records=[None,None];self.oracles=[None,None];self.queries=0;self.branch_cache={}

    def validate_bulk_order(self,place,bulk):
        left,right=BOUNDARIES[place]
        swapped=bool(bulk[right[0]]>bulk[right[1]])
        actual_left=left[::-1]if swapped else left
        pivot='I'if bulk[actual_left[0]]>=bulk[actual_left[1]]else'T'
        key=place,pivot,swapped
        if key not in self.branch_cache:
            coordinates=branch_coordinates(self.interval,('zero','infinity')[place],self.ramification,
                                           self.core,pivot,swapped,self.profile_precision)
            sides=[];metadata=[]
            for sector in range(2):
                matrix,shifts,md=defect_block(coordinates,self.local_factors[place],sector)
                sides.append(MatrixOracle(matrix,len(matrix),shifts,md['profile_precision']));metadata.append(md)
            self.branch_cache[key]=((coordinates,metadata),sides)
        self.records[place],self.oracles[place]=self.branch_cache[key]

    def block_legs(self,place,side):
        return self.records[place][0]['data']['bulk_block_legs'][side]

    def coordinate_counts(self,coordinates,side,place):
        basis=coordinates[side];data=coordinates['data'];m=self.interval//2
        vacuum=basis['vacuum_rows'if side=='left'else'vacuum_columns']
        size=len(next(iter(basis['labels'].values())));result=np.zeros((size,4))
        for letter,values in basis['labels'].items():result[:,'xyuv'.index(letter)]=[int(v).bit_count()for v in values]
        low,high=self.block_legs(place,side)
        for index in range(size):
            state=vacuum[index%len(vacuum)]
            for bit,original in enumerate(data['active_indices']):
                if state>>bit&1:result[index,low if original<m else high]+=1
        return result
