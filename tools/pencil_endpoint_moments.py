#!/usr/bin/env python3
"""Second-moment objective on the saturated four-letter pencil face."""
import numpy as np

from native_rectangular_smith_f4 import RectangularOracle
from pencil_principal_alphabets import restrict
from search_f4_recentered_crossings import five_letter_source


class EndpointMoments:
    def __init__(self):
        values,heights,_,_ = five_letter_source()
        natural = [v for group in heights for v in group[1:]]+[1]
        self.seed,self.natural,_,_ = restrict(values,natural,[0,2,3,4],[0,2,3,4])
        self.oracle = RectangularOracle(self.seed['values'],4,4)
        self.stride,self.n = 12,16

    def evaluate(self, variables, strong=0):
        if strong not in (0,1):
            raise ValueError('Choose infinity or zero as the strong channel.')
        variables = np.asarray(variables,dtype=float)
        if variables.shape!=(12,) or not np.all(np.isfinite(variables)):
            raise ValueError('Supply twelve finite strong-place symbol heights.')
        parameters = np.asarray(self.natural,dtype=float)
        parameters[strong*12:(strong+1)*12] = variables
        raw = self.oracle.prefixes(parameters,1)
        profiles = np.asarray(raw['orders']).reshape(2,16)
        r,s = profiles[strong],profiles[1-strong]
        gx,gy = np.zeros(4),np.zeros(4)
        for place in range(2):
            for leg in range(4):
                target = gx if leg in (0,3) else gy
                target[1:] += parameters[place*12+3*leg:place*12+3*(leg+1)]
        rc,xc,yc = r-r.mean(),gx-gx.mean(),gy-gy.mean()
        value = float(rc@rc/16-xc@xc/4-yc@yc/4)
        derivative = np.zeros(12)
        for i,weight in enumerate(2*rc/16):
            p = strong*16+i
            x,second = divmod(raw['rows'][p],4)
            third,v = divmod(raw['cols'][p],4)
            labels = (x,third,second,v) if strong==0 else (x,second,third,v)
            for leg,symbol in enumerate(labels):
                if symbol:
                    derivative[3*leg+symbol-1] += weight
        for leg in range(4):
            derivative[3*leg:3*(leg+1)] -= .5*(xc[1:] if leg in (0,3) else yc[1:])
        gap = float(gx.mean()+gy.mean()-r.mean()-min(s))
        multiplicity = int(sum(s==min(s)))
        if abs(gap)>1e-9 or multiplicity!=1:
            raise ArithmeticError('The registered saturated simple-minimum locus changed.')
        return value,derivative,{'parameters':parameters.tolist(),'strong_place':strong,
            'profiles':profiles.tolist(),'physical_energies':[gx.tolist(),gy.tolist()],
            'first_order_gap':gap,'weak_minimum_multiplicity':multiplicity,
            'variance_difference':value,'weak_order_gap':float(sorted(s)[1]-min(s))}
