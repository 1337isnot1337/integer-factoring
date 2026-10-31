#!/usr/bin/env python3
"""Regular physical hot-core coordinates, exported in exact parity blocks."""
from __future__ import annotations

import ctypes
from functools import lru_cache
import hashlib
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tempfile
import time

from check_rectangular_defect_support import Rat, Z, car_checks
from check_rectangular_defect_core_basis import residue_operator, left_car, right_car, independent_rows
from rectangular_defect_profiles import Program, Expr, precision_sizes, retained_known_below, MASK, lcm2, quotient2
from reflection_hot_state_support import core_data


@lru_cache(maxsize=32)
def coordinate_basis(interval, place, side, ramification=1, core=1, internal_regular=False):
    data = core_data(interval, place, ramification, core, internal_regular)
    return coordinate_basis_for_data(data, side)


def coordinate_basis_for_data(data, side):
    """Construct regular physical coordinates from a proved CAR support record."""
    core = data['core']
    names = data['core_names'][side]
    images = {}
    for name in names:
        for suffix in ('create', 'annihilate'):
            source = data['active_images'][name + '_' + suffix]
            creation = [Rat(*entry) for entry in source['creation']]
            annihilation = [Rat(*entry) for entry in source['annihilation']]
            if side == 'right':
                creation = [value * Z ** (-order) for value, order in zip(creation, data['active_orders'])]
                annihilation = [value * Z ** order for value, order in zip(annihilation, data['active_orders'])]
            assert all(value.order() >= 0 for value in creation + annihilation if value)
            images[name + '_' + suffix] = (creation, annihilation)
    car_checks(images, names)
    active, rank = len(data['active_indices']), len(names)
    size, multiplicity = 1 << active, 1 << (active - rank)
    operators = {name: residue_operator(pair) for name, pair in images.items()}
    projector = [1 << i for i in range(size)]
    for name in names:
        occupied = left_car(left_car(projector, operators[name + '_annihilate']),
                            operators[name + '_create'])
        projector = [first ^ second for first, second in zip(projector, occupied)]
    rows, columns = independent_rows(projector)
    assert len(rows) == len(columns) == multiplicity
    minor = [sum(((projector[i] >> j) & 1) << k for k, j in enumerate(columns)) for i in rows]
    assert len(independent_rows(minor)[0]) == multiplicity
    left, right = [], [0] * size
    for mask in range(1 << rank):
        row_map, column_map = projector, projector
        for i, name in enumerate(names):
            if mask >> i & 1:
                row_map = right_car(row_map, operators[name + '_annihilate'])
                column_map = left_car(column_map, operators[name + '_create'])
        left.extend(row_map[i] for i in rows)
        for i, value in enumerate(column_map):
            right[i] |= sum(((value >> j) & 1) << (mask * multiplicity + k)
                            for k, j in enumerate(columns))
    assert len(independent_rows(left)[0]) == len(independent_rows(right)[0]) == size
    # Every transported CAR generator is odd; the vacuum projector is even.
    for i, value in enumerate(projector):
        while value:
            bit = value & -value
            value ^= bit
            assert (i.bit_count() ^ (bit.bit_length() - 1).bit_count()) & 1 == 0
    letters = [names[0].split('_')[1][0], names[core].split('_')[1][0]]
    labels = {letter: [(mask >> (group * core)) & ((1 << core) - 1)
                       for mask in range(1 << rank) for _ in range(multiplicity)]
              for group, letter in enumerate(letters)}
    parities = lambda vacuum: [(mask.bit_count() + value.bit_count()) & 1
                               for mask in range(1 << rank) for value in vacuum]
    return {'names': names, 'regular_images': {
                name: {'creation': [value.pair() for value in pair[0]],
                       'annihilation': [value.pair() for value in pair[1]]}
                for name, pair in images.items()},
            'multiplicity': multiplicity, 'vacuum_rows': rows, 'vacuum_columns': columns,
            'row_parities': parities(rows), 'column_parities': parities(columns),
            'labels': labels, 'unit_basis_rank': size}


@lru_cache(None)
def library(precision):
    source = Path(__file__).with_name('native_fock_f4_parity.cpp')
    includes = [source.with_name(name) for name in ('native_fock_f4.cpp', 'native_rectangular_smith_f4.cpp')]
    compiler = shlex.split(os.environ.get('CXX', 'c++'))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError('The optional parity-coordinate engine requires C++17.')
    flags = ['-O3', '-std=c++17', '-shared', '-fPIC', f'-DF4_SERIES_WORDS={precision // 64}']
    try:
        cpu = Path('/proc/cpuinfo').read_text()
    except OSError:
        cpu = ''
    if platform.machine().lower() in ('x86_64', 'amd64') and 'pclmulqdq' in cpu:
        flags += ['-mpclmul', '-msse2']
    digest = hashlib.sha256(source.read_bytes() + b''.join(path.read_bytes() for path in includes)
                           + repr((compiler, flags)).encode()).hexdigest()[:24]
    cache = Path(tempfile.gettempdir()) / ('integer-factoring-fock-parity-' + str(os.getuid()))
    cache.mkdir(mode=0o700, parents=True, exist_ok=True)
    target = cache / (digest + '.so')
    if not target.exists():
        descriptor, temporary = tempfile.mkstemp(prefix=digest + '-', suffix='.so', dir=cache)
        os.close(descriptor)
        try:
            result = subprocess.run(compiler + flags + [str(source), '-o', temporary],
                                    capture_output=True, text=True)
            if result.returncode:
                raise RuntimeError(result.stderr[-12000:])
            os.replace(temporary, target)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    value = ctypes.CDLL(str(target))
    integer, unsigned, signed = ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_int64)
    value.fock_f4_coordinates_parity.argtypes = [ctypes.c_int, ctypes.c_int, integer, ctypes.c_int, unsigned,
        ctypes.c_int, integer, integer, integer, ctypes.c_int, integer, ctypes.c_int, integer, integer,
        ctypes.c_int, ctypes.c_int, integer, integer, integer, integer, unsigned, signed, signed]
    value.fock_f4_coordinates_parity.restype = ctypes.c_int
    return value


@lru_cache(maxsize=8)
def gaussian_coordinates(interval, place, ramification=1, core=1, profile_precision=512, internal_regular=False):
    data = core_data(interval, place, ramification, core, internal_regular)
    if len(data['active_indices']) > 11:
        raise ValueError('The parity-coordinate exporter currently supports at most eleven active modes.')
    left, right = [coordinate_basis(interval, place, side, ramification, core, internal_regular) for side in ('left', 'right')]
    if place == 'infinity':
        zero = gaussian_coordinates(interval, 'zero', ramification, core, profile_precision, internal_regular)
        # The compression theorem gives identical gauged kernels in the two
        # stated boundary orders. Check the coordinate bases as well before
        # reusing the one numerical matrix with its relabelled physical legs.
        for current, original in ((left, zero['left']), (right, zero['right'])):
            assert current['vacuum_rows'] == original['vacuum_rows']
            assert current['vacuum_columns'] == original['vacuum_columns']
            for a, b in zip(current['names'], original['names']):
                for kind in ('create', 'annihilate'):
                    assert current['regular_images'][a + '_' + kind] == original['regular_images'][b + '_' + kind]
        return {**zero, 'data': data, 'left': left, 'right': right,
                'metadata': {**zero['metadata'], 'reused_identical_zero_matrix': True}}
    return gaussian_coordinates_for_data(data, left, right, profile_precision)


def gaussian_coordinates_for_data(data, left, right, profile_precision=512):
    """Export complete parity matrices for an explicit regular coordinate basis."""
    interval = data['M']
    extended = dict(data)
    extended['active_images'] = {}
    names = []
    for side, basis in (('left', left), ('right', right)):
        prefix = 'basis_' + side + '_'
        names.append([prefix + name for name in basis['names']])
        extended['active_images'].update({prefix + name: value for name, value in basis['regular_images'].items()})
    program = Program(extended)
    size = 1 << program.k
    diagonal = [sum(order for j, order in enumerate(data['active_orders']) if mask >> j & 1)
                for mask in range(size)]
    def vacuum(value, pairs):
        for name in pairs:
            value = program.add(value, program.number(value, name))
        return value
    left_names, right_names = names
    right_vectors = {0: vacuum(Expr(0), right_names)}
    outputs = {}
    for state in range(1 << len(right_names)):
        if state:
            bit = state & -state
            right_vectors[state] = program.car(right_vectors[state ^ bit], right_names[bit.bit_length() - 1] + '_create')
        value = program.diagonal(right_vectors[state], diagonal)
        annihilated = {0: value}
        for mask in range(1 << len(left_names)):
            if mask:
                bit = mask & -mask
                annihilated[mask] = program.car(annihilated[mask ^ bit], left_names[bit.bit_length() - 1] + '_annihilate')
            outputs[mask, state] = vacuum(annihilated[mask], left_names)
    base = min(value.base for value in outputs.values())
    denominator = 1
    for value in outputs.values():
        denominator = lcm2(denominator, value.denominator)
    output_ids = []
    for mask in range(1 << len(left_names)):
        for state in range(1 << len(right_names)):
            value = outputs[mask, state]
            factor = quotient2(denominator, value.denominator) << (value.base - base)
            output_ids.append(program.polynomial(value, factor).register)
    window, builder = precision_sizes(program, base, interval, profile_precision)
    words, output_words = builder // 64, window // 64
    instructions = (ctypes.c_int * (5 * len(program.instructions)))(*(value for row in program.instructions for value in row))
    packed = (ctypes.c_uint64 * (2 * words * len(program.coefficients)))(
        *(value >> (64 * w) & MASK for pair in program.coefficients for value in pair for w in range(words)))
    terms, offsets = [], [0]
    for car in program.cars:
        terms.extend(value for term in car for value in term)
        offsets.append(len(terms) // 3)
    cc, co = (ctypes.c_int * len(terms))(*terms), (ctypes.c_int * len(offsets))(*offsets)
    dd = (ctypes.c_int * size)(*(value - min(diagonal) for value in diagonal))
    ii = (ctypes.c_int * len(left['vacuum_rows']))(*left['vacuum_rows'])
    jj = (ctypes.c_int * len(right['vacuum_columns']))(*right['vacuum_columns'])
    rr = (ctypes.c_int * len(output_ids))(*output_ids)
    row_parities, column_parities = left['row_parities'], right['column_parities']
    row_groups = [[i for i, parity in enumerate(row_parities) if parity == sector] for sector in range(2)]
    column_groups = [[i for i, parity in enumerate(column_parities) if parity == sector] for sector in range(2)]
    half = size // 2
    assert all(len(group) == half for group in row_groups + column_groups)
    row_positions, column_positions = [0] * size, [0] * size
    for groups, positions in ((row_groups, row_positions), (column_groups, column_positions)):
        for group in groups:
            for position, index in enumerate(group):
                positions[index] = position
    rp, cp = (ctypes.c_int * size)(*row_parities), (ctypes.c_int * size)(*column_parities)
    ri, ci = (ctypes.c_int * size)(*row_positions), (ctypes.c_int * size)(*column_positions)
    entries = 2 * half * half
    output = (ctypes.c_uint64 * (2 * output_words * entries))()
    shifts = (ctypes.c_int64 * entries)()
    meta = (ctypes.c_int64 * 4)()
    started = time.monotonic()
    status = library(builder).fock_f4_coordinates_parity(program.k, len(program.instructions), instructions,
        len(program.coefficients), packed, len(program.cars), co, cc, dd,
        left['multiplicity'], ii, right['multiplicity'], jj, rr, base, output_words,
        rp, ri, cp, ci, output, shifts, meta)
    if status:
        raise ArithmeticError(f'Parity-coordinate construction failed: {status}, metadata={list(meta)}')
    matrices, shift_matrices = [], []
    for sector in range(2):
        matrix, shift_matrix = [], []
        for i in range(half):
            row, shift_row = [], []
            for j in range(half):
                index = sector * half * half + i * half + j
                pair = tuple(sum(int(output[(2 * index + plane) * output_words + w]) << (64 * w)
                                 for w in range(output_words)) for plane in range(2))
                assert not pair[1]
                row.append(pair)
                shift_row.append(int(shifts[index]))
            matrix.append(row)
            shift_matrix.append(shift_row)
        matrices.append(matrix)
        shift_matrices.append(shift_matrix)
    bound = max(program.degree_bounds[index] for index in output_ids)
    metadata = {'matrix_order': size, 'parity_block_order': half, 'active_modes': program.k,
                'registers': len(program.instructions), 'profile_precision': window, 'builder_precision': builder,
                'unit_denominator_degree': denominator.bit_length() - 1, 'program_degree_bound': bound,
                'maximum_entry_coefficient_degree': int(meta[3]),
                'minimum_entry_laurent_shift': int(meta[0]),
                'builder_known_below': int(meta[1]),
                'exact_polynomial_matrix': bound < builder and int(meta[3]) < window,
                'known_below': retained_known_below(meta, window), 'matrix_seconds': time.monotonic() - started}
    return {'data': data, 'metadata': metadata, 'left': left, 'right': right,
            'row_groups': row_groups, 'column_groups': column_groups,
            'row_positions': row_positions, 'column_positions': column_positions,
            'matrices': matrices, 'shifts': shift_matrices}
