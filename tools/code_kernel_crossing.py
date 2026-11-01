"""Convolution-kernel normal form for binary-code compressed crossings."""

from functools import lru_cache


def xor_shift(mask, vector, width):
    """Permute a 2^width-bit truth vector by x -> x XOR vector."""
    dimension = 1 << width
    for j in range(width):
        if vector >> j & 1:
            step = 1 << j
            lower = ((1 << dimension)-1) // ((1 << step)+1)
            mask = ((mask & lower) << step) | ((mask >> step) & lower)
    return mask


@lru_cache(maxsize=None)
def convolution_rank(width, kernel):
    basis = {}
    for x in range(1 << width):
        vector = xor_shift(kernel, x, width)
        while vector:
            pivot = vector.bit_length()-1
            if pivot in basis:
                vector ^= basis[pivot]
            else:
                basis[pivot] = vector
                break
    return len(basis)


def base_kernels(width):
    dimension = 1 << width
    result = [[0]*dimension for _ in range(width+1)]
    for s in range(dimension):
        for w in range(dimension):
            if s & ~w:
                continue
            shift = s.bit_count()
            exponent = (w ^ s).bit_count()
            for degree in range(shift, shift+exponent+1):
                if ((degree-shift) & ~exponent) == 0:
                    result[degree][s] |= 1 << w
    return result


def toggle(kernels, parity, width):
    assert 0 < parity < 1 << width
    dimension = 1 << width
    answer = [[0]*dimension for _ in range(width+1)]
    assert all(row == xor_shift(row, parity, width) for row in kernels[0])
    for degree in range(width+1):
        next_rows = kernels[degree+1] if degree < width else [0]*dimension
        for s in range(dimension):
            answer[degree][s] = (
                xor_shift(kernels[degree][s], parity, width)
                ^ next_rows[s] ^ xor_shift(next_rows[s], parity, width)
            )
    return answer


def transpose(rows, width):
    dimension = 1 << width
    return [sum(((rows[s] >> w) & 1) << s for s in range(dimension))
            for w in range(dimension)]


def rank_profile(kernels, width):
    answer = []
    for rows in kernels:
        r = sum(convolution_rank(width, row) for row in rows)
        c = sum(convolution_rank(width, col) for col in transpose(rows,width))
        answer.append((r,c))
    return answer


def matrices_from_kernels(kernels, width):
    dimension = 1 << width
    result_rows, result_cols = [], []
    for kernel in kernels:
        rows, cols = [0]*dimension**2, [0]*dimension**2
        for x in range(dimension):
            for y in range(dimension):
                s = x ^ y
                for u in range(dimension):
                    w = x ^ u
                    if kernel[s] >> w & 1:
                        v = x ^ y ^ u
                        rows[x*dimension+y] |= 1 << (u*dimension+v)
                        cols[x*dimension+u] |= 1 << (y*dimension+v)
        result_rows.append(rows)
        result_cols.append(cols)
    return result_rows,result_cols
