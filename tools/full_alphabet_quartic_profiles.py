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
              identity_spec=None, extra_quartics=None, extra_monomials=None,
              residual_monomials=None):
        if any(value is not None for value in
               (repair_spec, input_spec, identity_spec)):
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
        M = self.data['M']
        top = M-1
        if residual_monomials is not None:
            from check_full_alphabet_monomial_residual import prepare_residual_repair
            endpoint_data = prepare_residual_repair(self.gamma, residual_monomials, top, repair_kind)
            correction = endpoint_data['repair_matrix']
        elif repair_kind == 'partial':
            correction = partial_repair(self.gamma)[0]
        elif repair_kind == 'inverse':
            correction = inverse(residual(self.gamma))[0]
        else:
            raise ValueError('Use partial or inverse repair.')
        infinity = self.data['place'] == 'infinity'
        from check_rectangular_extra_quartics import validate_extra_quartics
        from check_rectangular_core_monomials import monomial_program_terms
        extras = validate_extra_quartics(top, extra_quartics or [])
        extra_operations = monomial_program_terms(M, self.data['place'], extras+list(extra_monomials or []))
        extra_operations += monomial_program_terms(M, self.data['place'], residual_monomials or [], allow_marginal=True)
        for word, _, _, _ in extra_operations:
            if any(name not in self.car_meta for name in word):
                raise ValueError('This monomial needs the full-core support; a required CAR image is absent.')
        v0 = ('out_' if infinity else 'right_')+'v0'
        vT = ('out_' if infinity else 'right_')+f'v{top}'
        y0 = ('out_' if infinity else 'left_')+'y0'
        yT = ('out_' if infinity else 'left_')+f'y{top}'
        u0 = ('in_' if infinity else 'right_')+'u0'
        core, heights = (1, 0, 2, 3), (0, 0, 1, 1)

        def middle(source):
            if not self.gamma and not extra_operations:
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
            for extra_word, gauge, terms, pole in extra_operations:
                monomial = self.word(source, extra_word)
                for coefficient, exponent in terms:
                    if coefficient:
                        result = self.add(result, self.scale(
                            monomial, coefficient,
                            z=gauge+pole-exponent if infinity else gauge+exponent, u=-pole))
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
