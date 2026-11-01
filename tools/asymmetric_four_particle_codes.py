"""Exact two-alphabet four-particle restrictions of the strict Gaussian."""

from balanced_four_particle_codes import canonical_code, compress_square, natural_energies
from native_matrix_smith_f4 import MatrixOracle
from saturated_dual_kernel_directions import divide, polynomial_gcd
from selected_gaussian_coefficients import selected_tensor


def construct(code_x, code_y):
    code_x, code_y = canonical_code(code_x), canonical_code(code_y)
    original = selected_tensor(tuple(range(8)), 8, code_x, code_y)
    compressed = []
    for low, high in original['values']:
        assert not high
        compressed.append(compress_square(low))
    common = 0
    for value in compressed:
        common = polynomial_gcd(common, value)
    assert common and common.bit_count() % 2 == 1
    values = []
    for index, value in enumerate(compressed):
        quotient, remainder = divide(value, common)
        assert not remainder
        x, y, u, v = ((index // 5 ** leg) % 5 for leg in range(4))
        assert quotient.bit_count() % 2 == int(x == v and y == u)
        values.append((quotient, 0))
    return {'values': tuple(values), 'code_x': code_x, 'code_y': code_y,
            'dimensions': [5] * 4,
            'energies_x': natural_energies(code_x), 'energies_y': natural_energies(code_y),
            'entry_gcd_in_z': common, 'entry_gcd_degree': common.bit_length() - 1,
            'degree': max(value.bit_length() - 1 for value, _ in values),
            'raw_degree_in_z': max(value.bit_length() - 1 for value in compressed),
            'strict_endpoint': True,
            'native_complete_degree_bound_in_w': original['coefficient_degree_bound'],
            'native_exact_precision': original['exact_precision']}


def complete_profiles(values, energies_x, energies_y):
    matrices = [[[(0, 0)] * 25 for _ in range(25)] for _ in range(2)]
    for index, pair in enumerate(values):
        x, y, u, v = ((index // 5 ** leg) % 5 for leg in range(4))
        matrices[0][5 * x + u][5 * y + v] = pair
        matrices[1][5 * x + y][5 * u + v] = pair
    hx, hy = (tuple(value // 2 for value in energy) for energy in (energies_x, energies_y))
    infinity = MatrixOracle(matrices[0], 25).profile(
        row_weights=[hy[u] for x in range(5) for u in range(5)],
        column_weights=[hx[v] for y in range(5) for v in range(5)], infinity=True)['orders']
    zero = MatrixOracle(matrices[1], 25).profile(
        row_weights=[hx[x] + hy[y] for x in range(5) for y in range(5)])['orders']
    assert all(abs(value - round(value)) < 1e-8 for value in infinity + zero)
    return tuple(map(round, infinity)), tuple(map(round, zero))
