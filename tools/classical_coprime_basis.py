"""Exact reference for Bernstein's constructive coprime-basis algorithms.

Algorithm numbers refer to D. J. Bernstein, Factoring into coprimes in
essentially linear time, Journal of Algorithms 54 (2005), 1--30.
This Python implementation charges native arithmetic operations and bits;
the published fast-arithmetic analysis is separate from Python timings.
"""
from math import gcd

from oriented_norm_basis_lift import verify_coprime_norms


class CoprimeArithmetic:
    def __init__(self):
        self.statistics = {name:0 for name in (
            'gcd_calls','gcd_operand_bits','multiplications','multiplication_operand_bits',
            'exact_divisions','division_operand_bits','maximum_operand_bits',
            'pair_calls','maximum_pair_depth','support_splits','comparison_splits',
            'basis_extensions','basis_merges','coprime_merge_shortcuts',
            'family_factor_calls','single_factor_calls','power_reductions')}

    def _charge(self, kind, first, second):
        bits = first.bit_length()+second.bit_length()
        count, weight = {'gcd':('gcd_calls','gcd_operand_bits'),
                         'mul':('multiplications','multiplication_operand_bits'),
                         'div':('exact_divisions','division_operand_bits')}[kind]
        self.statistics[count] += 1
        self.statistics[weight] += bits
        self.statistics['maximum_operand_bits'] = max(self.statistics['maximum_operand_bits'],
                                                     first.bit_length(),second.bit_length())

    def common(self, first, second):
        self._charge('gcd',first,second)
        return gcd(first,second)

    def multiply(self, first, second):
        self._charge('mul',first,second)
        return first*second

    def divide(self, first, second):
        self._charge('div',first,second)
        value, remainder = divmod(first,second)
        assert remainder==0
        return value

    def product(self, values):
        layer = list(values)
        if not layer:
            return 1
        while len(layer)>1:
            layer = [self.multiply(layer[i],layer[i+1]) if i+1<len(layer) else layer[i]
                     for i in range(0,len(layer),2)]
        return layer[0]

    def support_split(self, value, support):
        """Return the gcd and the full supported/unsupported parts (11.3)."""
        self.statistics['support_splits'] += 1
        original = inside = self.common(value,support)
        outside = self.divide(value,inside)
        while True:
            overlap = self.common(inside,outside)
            if overlap==1:
                return original,inside,outside
            inside = self.multiply(inside,overlap)
            outside = self.divide(outside,overlap)

    def comparison_split(self, first, second):
        """Partition first by whether its prime exponent exceeds second (11.4)."""
        self.statistics['comparison_splits'] += 1
        original = lower = self.common(first,second)
        greater = self.divide(first,lower)
        while True:
            overlap = self.common(greater,lower)
            if overlap==1:
                return original,greater,lower
            greater = self.multiply(greater,overlap)
            lower = self.divide(lower,overlap)

    def pair_basis(self, first, second, depth=0):
        """Natural two-integer basis by left-shift exponent separation (13.2)."""
        self.statistics['pair_calls'] += 1
        self.statistics['maximum_pair_depth'] = max(self.statistics['maximum_pair_depth'],depth)
        if second==1:
            return [first] if first>1 else []
        _,inside,outside = self.support_split(first,second)
        answer = [outside] if outside>1 else []
        common,greater,lower = self.comparison_split(inside,second)
        first_lower = accumulated = lower
        shift = 0
        while True:
            common,greater,lower = self.comparison_split(greater,self.multiply(common,common))
            divisor = self.common(lower,second)
            accumulated = self.multiply(accumulated,divisor)
            power = divisor
            for _ in range(shift):
                power = self.multiply(power,power)
            answer.extend(self.pair_basis(self.divide(lower,power),divisor,depth+1))
            if greater==1:
                break
            shift += 1
        answer.extend(self.pair_basis(self.divide(second,accumulated),first_lower,depth+1))
        return answer

    def split_over(self, value, basis):
        """Associate each coprime basis entry with its supported part (15.3)."""
        if not basis:
            return []
        _,inside,_ = self.support_split(value,self.product(basis))
        if len(basis)==1:
            return [(basis[0],inside)]
        middle = len(basis)//2
        return self.split_over(inside,basis[:middle])+self.split_over(inside,basis[middle:])

    def extend_basis(self, basis, value):
        """Refine a coprime family together with one new integer (16.2)."""
        self.statistics['basis_extensions'] += 1
        if not basis:
            return [value] if value>1 else []
        if value==1:
            return list(basis)
        _,inside,outside = self.support_split(value,self.product(basis))
        answer = [outside] if outside>1 else []
        for old,part in self.split_over(inside,basis):
            answer.extend(self.pair_basis(old,part))
        return answer

    def merge_bases(self, first, second):
        """Use logarithmically many binary-index products to merge bases (17.3)."""
        self.statistics['basis_merges'] += 1
        if not first or not second:
            return list(first or second)
        if self.common(self.product(first),self.product(second))==1:
            self.statistics['coprime_merge_shortcuts'] += 1
            return list(first)+list(second)
        result = list(first)
        for bit in range(max(1,(len(second)-1).bit_length())):
            for digit in (0,1):
                packed = self.product(value for i,value in enumerate(second) if (i>>bit)&1==digit)
                result = self.extend_basis(result,packed)
        return result

    def natural_basis(self, values):
        """Balanced recursion over the supplied list (18.1)."""
        if len(values)<=1:
            return [values[0]] if values and values[0]>1 else []
        middle = len(values)//2
        return self.merge_bases(self.natural_basis(values[:middle]),
                                self.natural_basis(values[middle:]))

    def reduce_power(self, base, value):
        """Extract the largest power by recursive exponent doubling (19.2)."""
        self.statistics['power_reductions'] += 1
        if self.common(base,value)!=base:
            return 0,value
        exponent,remainder = self.reduce_power(self.multiply(base,base),self.divide(value,base))
        if self.common(base,remainder)==base:
            return 2*exponent+2,self.divide(remainder,base)
        return 2*exponent+1,remainder

    def factor_one(self, value, basis):
        """Return sparse powers over a coprime basis, or reject it (20.1)."""
        self.statistics['single_factor_calls'] += 1
        if value==1:
            return []
        if not basis:
            raise ValueError('The supplied coprime basis does not factor the integer.')
        if len(basis)==1:
            exponent,remainder = self.reduce_power(basis[0],value)
            if remainder!=1:
                raise ValueError('The supplied coprime basis does not factor the integer.')
            return [(basis[0],exponent)] if exponent else []
        middle = len(basis)//2
        _,left,right = self.support_split(value,self.product(basis[:middle]))
        return self.factor_one(left,basis[:middle])+self.factor_one(right,basis[middle:])

    def factor_family(self, values, basis):
        """Filter relevant basis entries before dividing the input family (21.2)."""
        self.statistics['family_factor_calls'] += 1
        if not values:
            return {}
        _,relevant,_ = self.support_split(self.product(basis),self.product(values))
        selected = [base for base,part in self.split_over(relevant,basis) if part==base]
        if len(values)==1:
            return {values[0]:self.factor_one(values[0],selected)}
        middle = len(values)//2
        result = self.factor_family(values[:middle],selected)
        result.update(self.factor_family(values[middle:],selected))
        return result


def factor_into_coprimes(values, *, verify=True):
    """Compute a natural coprime basis and every original sparse exponent word."""
    values = list(values)
    if any(not isinstance(value,int) or value<1 for value in values):
        raise ValueError('Supply positive integers.')
    distinct = sorted(set(values))
    arithmetic = CoprimeArithmetic()
    basis = sorted(arithmetic.natural_basis(distinct))
    factored = arithmetic.factor_family(distinct,basis)
    indices = {value:i for i,value in enumerate(basis)}
    words = [sorted((indices[base],exponent) for base,exponent in factored[value]) for value in values]
    result = {'input_values':values,'basis':basis,'words':words,
              'statistics':dict(arithmetic.statistics),
              'input_bits':sum(value.bit_length() for value in values),
              'distinct_inputs':len(distinct),
              'factoring_oracle_calls':0,'primality_tests':0}
    if verify:
        check = verify_coprime_norms(basis)
        for value,word in zip(values,words):
            reconstructed = 1
            for column,exponent in word:
                assert 0<exponent<=value.bit_length()
                reconstructed *= basis[column]**exponent
                assert reconstructed<=value
            assert reconstructed==value
        result['verification'] = {'complete_norm_words':len(values),'coprime_basis':check}
    return result


def classical_oriented_coprime_basis(discriminant, keys, *, verify=True):
    """Construct all norm data, then attach orientations without prime factoring."""
    from oriented_norm_basis_lift import lift_coprime_norm_basis
    keys = list(keys)
    norm_record = factor_into_coprimes([norm for norm,_ in keys],verify=verify)
    result = lift_coprime_norm_basis(discriminant,keys,norm_record['basis'],norm_record['words'],verify=verify)
    result['norm_construction'] = {'algorithm':'Bernstein 2005 natural coprime basis and sparse factoring',
                                   'statistics':norm_record['statistics'],
                                   'distinct_norms':norm_record['distinct_inputs']}
    result['scope'] = 'Complete supplied-key norm construction and orientation; no relation-supply bound.'
    return result
