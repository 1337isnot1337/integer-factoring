# The coarse class sum and the coefficient-entry class sum differ

## Result

The twelve-site coarse `(T,R)` certificate and the coefficient table use
two different class-sum normalizations. The coarse operator is the
Kronecker sum

```text
T_coarse = sum_i M_i kron M_i,
```

while the coefficient-entry functional on a quadruple `(a,b,c,d)` gives

```text
T_entry[(a,b),(c,d)] = sum_i M_i[a,d] M_i[b,c].
```

The two are related by the path-pair swap `R` on the column pair,

```text
T_entry = T_coarse R.
```

At ten sites both matrices act on the `34^2 = 1156` path-pair module and
fit in memory. Their exact generalized nullities are

```text
eigenvalue      0      1      2      3
T_coarse      441    225    490      0
T_entry       441    225    335    155.
```

The coarse `T=2` block splits into `335 + 155` under `T_entry`. The
published ten-site primary dimensions are

```text
34, 1, 225, 300, 35, 406, 155,
```

so the `335` coarse block is the sum `34+1+300`, and the `155` primary
is the new `T_entry` eigenvalue `3` block. This is not a relabeling of
one block: the swap changes the operator, not just its name.

## Consequence for twelve sites

The exact coarse `T`-primary dimensions `1485,3281,3155` belong to
`T_coarse`. The separator joint characters and the row trace groups use
the coefficient-table `T`, which is the entry-convention operator. Their
group sums are `0,0,1` modulo five, while the coarse dimensions are
`0,1,0`. The ten-site control explains why the two need not agree: the
path-pair swap can split a coarse generalized eigenspace into several
entry-convention eigenvalues.

Therefore the exact coarse dimensions should not be combined with the
joint `T,S` characters until the twelve-site `T_entry` spectrum is
computed. That computation is the next target. The all-three-group row
split is independent of this normalization question; it uses the
coefficient-table characters and the row trace directly.

## Reproduction

```sh
python3 tools/check_bounded_spin_T_normalization.py
make verify-bounded-spin-T-normalization
```

The checker builds both ten-site matrices, verifies the entry formula on
sample entries, computes the generalized nullities, and checks the split
of the coarse `T=2` block. It does not build a twelve-site matrix or claim
a cut rank. No factoring exponent changes; the complete bound remains
`1/2`.
