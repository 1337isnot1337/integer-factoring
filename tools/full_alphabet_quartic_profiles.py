#!/usr/bin/env python3
"""Bounded-support CAR program for the full-alphabet quartic family.

The endpoint proof is in full-alphabet-quartic-repair.md. This program
uses exactly the original canonical quartic and the explicit 16-state
correction. Its generic cut rank is full, supplied by the typed endpoint.
"""
from check_full_alphabet_quartic_repair import residual, inverse, partial_repair
from rectangular_defect_profiles import Program, Expr


class FullAlphabetProgram(Program):
    def __init__(self, data, gamma):
        if gamma not in (0, 2, 3):
            raise ValueError('The full-alphabet F4 coefficient must be 0, 2, or 3.')
        super().__init__(data)
        self.gamma = gamma

    def build(self, alpha=1, beta=2, quartic_terms=None, initial=None,
              repair_spec=None, input_spec=None, repair_kind='partial',
              identity_spec=None, extra_quartics=None):
        if any(value is not None for value in
               (repair_spec, input_spec, identity_spec, extra_quartics)):
            raise ValueError('Full-alphabet coordinates currently use their explicit repair only.')
        if quartic_terms is None:
            quartic_terms = [(1, 0)]
        endpoint = 0
        for coefficient, exponent in quartic_terms:
            if coefficient not in (0, 1, 2, 3) or not isinstance(exponent, int):
                raise ValueError('Invalid quartic Laurent coefficient.')
            endpoint ^= coefficient
        if endpoint != 1:
            raise ValueError('The scalar quartic numerator must equal one at w=1.')
        if repair_kind == 'partial':
            correction = partial_repair(self.gamma)[0]
        elif repair_kind == 'inverse':
            correction = inverse(residual(self.gamma))[0]
        else:
            raise ValueError('Use partial or inverse repair.')
        M = self.data['M']
        top = M-1
        infinity = self.data['place'] == 'infinity'
        v0 = ('out_' if infinity else 'right_')+'v0'
        vT = ('out_' if infinity else 'right_')+f'v{top}'
        y0 = ('out_' if infinity else 'left_')+'y0'
        yT = ('out_' if infinity else 'left_')+f'y{top}'
        u0 = ('in_' if infinity else 'right_')+'u0'
        core, heights = (1, 0, 2, 3), (0, 0, 1, 1)

        def middle(source):
            if not self.gamma:
                return source
            if infinity:
                word = [vT+'_create', yT+'_create', u0+'_annihilate',
                        'in_x0_annihilate']
            else:
                word = ['left_x0_create', yT+'_create', u0+'_annihilate',
                        vT+'_annihilate']
            quartic = self.word(source, word)
            result = source
            from check_rectangular_f4_repair import fmul
            for coefficient, exponent in quartic_terms:
                coefficient = fmul(self.gamma, coefficient)
                if coefficient:
                    result = self.add(result, self.scale(
                        quartic, coefficient,
                        z=2-exponent if infinity else exponent, u=-2))
            return result

        def physical(v, y):
            if infinity:
                return ((1-M)*bool(core[v]&1)+(M-1)*bool(core[v]&2)
                        -M*core[y].bit_count())
            return bool(core[y]&1)+(2*M-1)*bool(core[y]&2)

        if initial is None:
            initial = Expr(0)
        base = middle(initial) if infinity else None
        parts = []
        for row, line in enumerate(correction):
            v, y = row % 4, row // 4
            for column, coefficient in enumerate(line):
                if not coefficient:
                    continue
                source_v, source_y = column % 4, column // 4
                if infinity:
                    value = self.unit(base, [v0, vT], core[v], core[source_v])
                    value = self.unit(value, [y0, yT], core[y], core[source_y])
                else:
                    value = self.unit(initial, [v0, vT], core[source_v], core[v])
                    value = middle(value)
                    value = self.unit(value, [y0, yT], core[y], core[source_y])
                height = heights[source_v]+heights[source_y]-heights[v]-heights[y]
                power = physical(v, y)-physical(source_v, source_y)
                if infinity:
                    power -= height
                parts.append(self.scale(value, coefficient, z=power, u=height))
        result = parts[0]
        for part in parts[1:]:
            result = self.add(result, part)
        return result
