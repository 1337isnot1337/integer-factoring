# Twelve-site T-primary cut dimensions modulo five

## Result

Let `T` be the transposition class sum in the twelve-site coefficient
algebra and let `C_row`, `C_col` be the two balanced cut modules.  On both
modules `T` satisfies

```math
T^5=4T+2T^3,
\qquad
T(T-1)^2(T+1)^2=0.
```

The generalized eigenspace dimensions for eigenvalues `0,1,4`, reduced
modulo five, are

| cut | `d_0` | `d_1` | `d_4` | `m_1` | `m_4` |
|---|---:|---:|---:|---:|---:|
| row `(x,y)|(u,v)` | `2` | `3` | `4` | `2` | `4` |
| column `(x,u)|(y,v)` | `1` | `3` | `2` | `1` | `2` |

Here `m_lambda` is the trace of the nilpotent part of `T` on the
`lambda`-primary block, also modulo five.  The total cut dimensions are
therefore `4` and `1` modulo five, respectively.  In particular the
`T=1` primary has the same residue `3` in both cuts, while the `T=0` and
`T=4` primaries differ; the two cut modules are not isomorphic as
`T`-modules.

These are congruences, not exact dimensions.  They do not determine the
individual primitive-central-idempotent ranks, the two-cut pairing, or
the score.

## Method

On the `T=0` primary, `T` is nilpotent and `T^4-2T^2-4` is invertible
with eigenvalue `-4=1`, so the relation forces `T=0`.  On the `T=1` and
`T=4` primaries, respectively,

```math
T=1+N,\quad N^2=0,
\qquad
T=4(1+M),\quad M^2=0,
```

because the corresponding minimal-polynomial factors are squared.  Thus
for `k>=1`,

```math
\operatorname{Tr}(T^k)
 =d_1+d_4 4^k+k\,m_1+k\,4^k m_4
 \pmod 5,
```

and `d_0+d_1+d_4=\operatorname{Tr}(1)`.  The traces
`\operatorname{Tr}(T^k)` for `k=0,\ldots,5` are computed from the exact
coefficient table and the row/column cut trace vectors.  The four
equations `k=1,\ldots,4` have a unique solution modulo five; `k=5` is an
independent check.

As a cross-check, the eight idempotent trace residues from the
[cut-trace certificate](../evidence/bounded-spin-twelve-cut-trace-form.json)
are grouped by their joint `T` character from the
[`R,T` coefficient certificate](../evidence/bounded-spin-twelve-RT-coefficient.json).
The group sums reproduce exactly the three `d`-values in the table.  The
checker is
[`check_bounded_spin_twelve_T_primary_dimensions.py`](../tools/check_bounded_spin_twelve_T_primary_dimensions.py);
its receipt is
[`bounded-spin-twelve-T-primary-dimensions.json`](../evidence/bounded-spin-twelve-T-primary-dimensions.json).

## Scope

This is a structural checksum on the two cut modules, not a rank
computation.  The exact dimensions may differ from the displayed values
by multiples of five.  The `R`-primary split shows that the ordinary
idempotent ranks are mixed, so the `T`-primary residues here are not a
substitute for the missing exact two-cut ranks.  No score change is
claimed.
