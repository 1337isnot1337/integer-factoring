"""Exact crossing pencils from affine parity selectors and dual linear sections."""

from code_kernel_crossing import xor_shift


def slope_family(width, free_mask):
    dimension = 1 << width
    nonunits = [p for p in range(1,dimension) if p.bit_count()>1]
    assert 0 <= free_mask < 1 << (width * len(nonunits))
    slopes = [0]*dimension
    for j in range(width):
        slopes[1 << j] = 1 << j
    for index,p in enumerate(nonunits):
        v = (free_mask >> (width*index)) & (dimension-1)
        slopes[p] = v
        for j in range(width):
            if p >> j & 1:
                slopes[1 << j] ^= v
    return slopes


def constant_representatives(slopes,width):
    dimension = 1 << width
    basis = {}
    for j in range(width):
        column = sum(((slopes[p] >> j) & 1) << (p-1) for p in range(1,dimension))
        while column:
            pivot = column.bit_length()-1
            if pivot in basis:
                column ^= basis[pivot]
            else:
                basis[pivot] = column
                break
    assert len(basis)==width
    free = [j for j in range(dimension-1) if j not in basis]
    for value in range(1 << len(free)):
        yield sum(((value >> i) & 1) << j for i,j in enumerate(free))


def general_kernels(width,slopes,constants):
    dimension = 1 << width
    result = [[0]*dimension for _ in range(width+1)]
    for s in range(dimension):
        coefficients = [1]+[0]*width
        for p in range(1,dimension):
            selected = ((constants >> (p-1)) & 1) ^ ((slopes[p] & s).bit_count() & 1)
            if not selected:
                continue
            assert coefficients[width] == xor_shift(coefficients[width],p,width)
            for j in range(width,0,-1):
                coefficients[j] ^= coefficients[j-1] ^ xor_shift(coefficients[j-1],p,width)
        for j,kernel in enumerate(coefficients):
            kernel = xor_shift(kernel,s,width)
            for b in range(j+1):
                if not b & ~j:
                    result[width-j+b][s] ^= kernel
    return result


def validate_slopes(slopes,width):
    for j in range(width):
        result = 0
        for p in range(1,1 << width):
            if slopes[p] >> j & 1:
                result ^= p
        assert result == 1 << j
