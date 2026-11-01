#!/usr/bin/env python3
"""Eight-k-mode transport for arbitrary fixed physical reflection maps.

Both internal Gaussian factors are retained.  The output is the complete
creation/annihilation operator interface, not an asserted Smith profile.
The independent low-band calculation uses the bivariate specialization
T=t**(2*M) -> 0; physical normalization weights are handled separately.
"""

from check_cancellation_module_automata import pmul


def validate(M, k):
    assert isinstance(k, int) and k > 0 and not k & (k - 1)
    assert isinstance(M, int) and M >= 4 * k and not M & (M - 1)


def assemble(k, f, h, x, T, T_over_x, first_sum, second_sum, multiply):
    """The four physical operator groups, before any physical map.

    An operator is (creation coefficients, annihilation coefficients).
    Each physical mode has its creator and annihilator in that order.
    Active blocks: x0,x1,x(N-2),x(N-1),y0,y1,y(N-2),y(N-1).
    """
    size = 8 * k
    pair = lambda: ([0] * size, [0] * size)
    groups = [[], [], [], []]
    for a in range(k):
        cx, ax, cy, ay = pair(), pair(), pair(), pair()
        cv, av, cu, au = pair(), pair(), pair(), pair()
        cx[0][a] = ax[1][a] = cv[0][a] = av[1][a] = 1
        cy[0][4 * k + a] = ay[1][4 * k + a] = 1
        cu[0][4 * k + a] = au[1][4 * k + a] = 1
        reverse = k - 1 - a
        ax[0][3 * k + reverse] = 1 ^ T
        cv[1][k + reverse] = cv[1][3 * k + reverse] = 1 ^ T
        for c in range(k):
            fca, fac, hca = f(c, a), f(a, c), h(c, a)
            cx[0][4 * k + c] ^= fca
            cx[0][5 * k + c] ^= multiply(x, fca)
            Th = multiply(T, hca)
            for block in (4, 6, 7):
                ax[0][block * k + c] ^= Th
            for block, value in ((0, fac), (1, multiply(first_sum, fac)), (3, fac)):
                ay[1][block * k + c] ^= value
            for block in (0, 2, 3):
                ay[0][block * k + c] ^= hca
            for block in (4, 5, 7):
                cv[1][block * k + c] ^= Th
            cv[1][6 * k + c] ^= multiply(Th, second_sum)
            av[1][4 * k + c] ^= fca
            av[1][5 * k + c] ^= multiply(multiply(x, first_sum), fca)
            av[1][7 * k + c] ^= multiply(T_over_x, fca)
            for block in (0, 1):
                cu[0][block * k + c] ^= fac
            for block, value in ((0, hca), (1, hca),
                                 (2, multiply(second_sum, hca)), (3, hca)):
                cu[1][block * k + c] ^= value
        for group, operators in zip(groups, ((cx, ax), (cy, ay), (cv, av), (cu, au))):
            group.append(operators)
    return groups


def polynomial_images(M, k=2):
    """Complete binary polynomial images at one finite physical interval."""
    validate(M, k)
    N = M // k
    first = sum(1 << (2 * k * j) for j in range(N - 1))
    second = sum(1 << (2 * k * j) for j in range(N - 2))
    f = lambda a, b: 3 << (2 * (a & ~b))
    h = lambda a, b: 3 << (2 * (a & b))
    images = assemble(k, f, h, 1 << (2 * k), 1 << (2 * M),
                      1 << (2 * M - 2 * k), first, second, pmul)
    retained = list(range(2 * k)) + list(range(M - 2 * k, M))
    orders = [2 * j + 1 for j in retained]
    orders += [2 * M + 2 * j + 1 for j in retained]
    spectators = [2 * j + 1 for j in range(M) if j not in retained]
    spectators += [2 * M + 2 * j + 1 for j in range(M) if j not in retained]
    return images, orders, spectators


def concentration_numerator(M, k=2):
    """C multiplied by its common local-unit denominator.

    C retains tail columns 1,N-2,N-1 in each residue.  Its inverse,
    and the J-weighted inverse, are regular because their constant
    determinants are units.  No endpoint regularity is asserted.
    """
    validate(M, k)
    N = M // k
    denominator = 1 ^ (1 << (2 * k * (N - 3)))
    C = [[0] * M for _ in range(M)]
    for c in range(k):
        C[c][c] = denominator
        first, penultimate, last = c + k, c + k * (N - 2), c + k * (N - 1)
        for B in range(1, N):
            row = c + k * B
            C[row][first] = denominator << (2 * k * (B - 1))
            if B < N - 1:
                C[row][penultimate] = denominator
        C[last][last] = denominator
        for j in range(2, N - 2):
            column = c + k * j
            C[column][column] = denominator
            C[first][column] = 1 ^ (1 << (2 * k * (N - 2 - j)))
            C[penultimate][column] = ((1 << (2 * k * (N - 2 - j)))
                                      ^ (1 << (2 * k * (N - 3))))
    return C, denominator


def finite_field_images(k, t, T, multiply, power):
    """Bivariate images at an exact finite-field specialization.

    The caller must use characteristic two and t*(1+t**(2*k)) != 0.
    T=0 is the fixed-code, infinite-tail leading band.
    """
    x = power(t, 2 * k)
    inverse_x = power(x, -1)
    denominator_inverse = power(1 ^ x, -1)
    first = multiply(1 ^ multiply(T, inverse_x), denominator_inverse)
    second = multiply(1 ^ multiply(T, multiply(inverse_x, inverse_x)), denominator_inverse)
    f = lambda a, b: multiply(1 ^ t, power(t, 2 * (a & ~b)))
    h = lambda a, b: multiply(1 ^ t, power(t, 2 * (a & b)))
    return assemble(k, f, h, x, T, multiply(T, inverse_x), first, second, multiply)
