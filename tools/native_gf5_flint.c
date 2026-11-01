#include <stdint.h>
#include <flint/nmod_mat.h>

static void load(nmod_mat_t matrix, const uint8_t *data)
{
    for (slong i = 0; i < matrix->r; i++)
        for (slong j = 0; j < matrix->c; j++)
            nmod_mat_entry(matrix, i, j) = data[i * matrix->c + j];
}

int64_t gf5_rank(const uint8_t *data, int64_t rows, int64_t columns)
{
    nmod_mat_t matrix;
    nmod_mat_init(matrix, rows, columns, 5);
    load(matrix, data);
    slong result = nmod_mat_rank(matrix);
    nmod_mat_clear(matrix);
    return result;
}

int64_t gf5_rref(uint8_t *out, const uint8_t *data, int64_t rows, int64_t columns)
{
    nmod_mat_t matrix;
    nmod_mat_init(matrix, rows, columns, 5);
    load(matrix, data);
    slong result = nmod_mat_rref(matrix);
    for (slong i = 0; i < rows; i++)
        for (slong j = 0; j < columns; j++)
            out[i * columns + j] = nmod_mat_entry(matrix, i, j);
    nmod_mat_clear(matrix);
    return result;
}

void gf5_multiply(uint8_t *out, const uint8_t *left, const uint8_t *right,
                  int64_t rows, int64_t inner, int64_t columns)
{
    nmod_mat_t a, b, c;
    nmod_mat_init(a, rows, inner, 5);
    nmod_mat_init(b, inner, columns, 5);
    nmod_mat_init(c, rows, columns, 5);
    load(a, left);
    load(b, right);
    nmod_mat_mul(c, a, b);
    for (slong i = 0; i < rows; i++)
        for (slong j = 0; j < columns; j++)
            out[i * columns + j] = nmod_mat_entry(c, i, j);
    nmod_mat_clear(a);
    nmod_mat_clear(b);
    nmod_mat_clear(c);
}
