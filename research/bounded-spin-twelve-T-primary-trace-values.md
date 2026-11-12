# Twelve-site T-primary trace weights modulo five

## Result

Let `T` be the transposition class sum in the twelve-site coefficient
algebra, and let `tau_row`, `tau_col` be the two cut trace functionals
from the diagonal trace forms.  These are **linear functionals on the
coefficient algebra**, not traces of operators on the cut modules:
`tau(C_0)` is `4` in the row orientation and `1` in the column
orientation, while `C_0` is not the cut-module identity.

In the coefficient algebra `T` satisfies

```math
T^5=4T+2T^3,
\qquad
T(T-1)^2(T+1)^2=0.
```

Group the eight primitive central idempotents by the eigenvalue of `T`:

```text
T=0: e_0,e_6
T=1: e_1,e_2,e_5,e_7
T=4: e_3,e_4
```

All three groups have scalar `T` action except that the `T=1` and `T=4`
dual blocks contribute square-zero nilpotent corrections.  Write

```math
a_\lambda=\sum_{e_i:\ T e_i=\lambda e_i}\tau(e_i),
\qquad
m_\lambda=\tau(N_\lambda),\quad N_\lambda^2=0
```

for each trace functional `tau`.  Then, for `k>=1`,

```math
\tau(T^k)
 =a_1+a_4 4^k+k\,m_1+k\,4^k m_4
 \pmod 5,
```

and `a_0+a_1+a_4=\tau(C_0)`.  The exact values are

| cut functional | `a_0` | `a_1` | `a_4` | `m_1` | `m_4` |
|---|---:|---:|---:|---:|---:|
| row `(x,y)|(u,v)` | `2` | `3` | `4` | `2` | `4` |
| column `(x,u)|(y,v)` | `1` | `3` | `2` | `1` | `2` |

These are **trace-functional weights**, not generalized eigenspace
dimensions.  The sums `a_0+a_1+a_4` are `4` and `1` modulo five, not
`1`; this is consistent with `C_0` not being the cut-module identity.
The two trace functionals differ, but that gives no module isomorphism or
non-isomorphism statement.

## Method and check

The power traces `tau(T^k)` for `k=0,\ldots,5` are computed from the exact
coefficient table and the stored row/column trace vectors.  The four
equations `k=1,\ldots,4` have a unique solution modulo five, and `k=5` is
an independent power-relation check.

The idempotent trace residues grouped by joint `T` character reproduce
the same three sums.  This is a consistency check of the recurrence
against the idempotent residues, not an independent module-dimension
certificate.

## Scope

This records structural trace-functional data on the two cut trace
orientations.  It does not give exact dimensions, individual central
idempotent ranks, the two-cut pairing, or a score change.  The retraction
of the `R`-twist interpretation is in
[the column-trace correction](bounded-spin-twelve-column-trace-correction.md);
this note makes the corresponding correction for the `T`-power table.

## Reproduction

```sh
python3 tools/check_bounded_spin_twelve_T_primary_trace_values.py \
  --output evidence/bounded-spin-twelve-T-primary-trace-values.json
make verify-bounded-spin-twelve-T-primary-trace-values
```

The receipt is
[`bounded-spin-twelve-T-primary-trace-values.json`](../evidence/bounded-spin-twelve-T-primary-trace-values.json).
