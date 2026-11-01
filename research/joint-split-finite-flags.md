# Finite-field joint splits with free low-rank factors

This checkpoint queries the rank conditions left open by the
[protected-core guards](joint-split-jet-guards.md). It allows the leading
and first residual matrices to vary. Over F2, F4 and F8, the recorded Z3
runs return only physical diagonal phases of the coordinate crossing.
Together with the existing guards and an exact baseline inequality, these
runs give a **solver-trusted balanced natural floor `rho >= 5/7`** for the
entire declared protected split family over those three fields. This is
not an all-field theorem. The F16 classification remains unresolved after
a timeout, and asymmetric charges or different physical heights are outside
the floor.

The [eight query receipts](../evidence/joint-split-flag-queries.json) record
source hashes, formula hashes, field arithmetic controls, positive controls,
solver version, timeout and outcomes. The optional
[driver](../tools/probe_joint_split_rank_flags.py) requires `z3-solver`;
the retained runs used version 5.3.3. The ordinary Python
[receipt checker](../tools/check_joint_split_flag_queries.py) does not
independently certify UNSAT and does not rerun the solver.

## The complete declared finite query

Use the literal five-letter source `A(z)` and target physical heights
`h=(0,1,3,4)` from the preceding note. For every joint target label `j`,
allow all entries of the embedding column `P_j` and retraction row `L_j`
in the source energy fiber of height `h_x+h_y`. Enforce `LP=I_16` and fix
all coefficients of the two-letter physical subtensor: the rows and
columns indexed by `{0,1,4,5}` in `(x,y)` order.

There are 88 homogeneous map entries before normalization. The source
vacuum energy fiber is one dimensional, so `LP=I` makes its two entries
nonzero inverses. Scaling the entire embedding and inversely scaling the
entire retraction sets both vacuum entries to one and leaves `LAP`
unchanged. The query therefore has **86 field variables**, with no
additional fixed map entries. Protecting the tensor core is weaker than
fixing the individual low-energy maps.

Construct every coefficient of `C=LAP`, then both actual flattenings.
Let `M0,M1` be either normalized cut's first two matrices. The fixed core
makes `M0[0,0]=1`. Leading rank one is exactly

```
M0[i,j] = M0[i,0] M0[0,j].
```

After eliminating this unit pivot, its first residual coefficient is

```
R[i,j] = M1[i,j] + M0[i,0] M1[0,j] + M1[i,0] M0[0,j]
         + M0[i,0] M1[0,0] M0[0,j]       (i,j > 0).
```

The core gives `R[1,1]=1`; residual rank one is exactly
`R[i,j]=R[i,1]R[1,j]`. These equations are imposed at both cuts. They
are equivalent to one zero order and one order one at each cut in this
core family, rather than fixing the two jet matrices.

## Two smaller queries replace a difficult orbit disequality

Four literal coefficients of `C` recover the independent physical phase
anchors: coefficient `z^8` at positions `(8,1),(2,4),(12,5),(3,5)`.
All four are one in the baseline. Write them as `a,b,c,d` and define

```
e_(x,y) = s_x t_y,       s=(1,1,a,c), t=(1,1,b,d).
```

The first query asks whether any anchor can be zero. The second fixes
all four anchors to one and asks whether any coefficient of `C` differs
from baseline. When all anchors are nonzero, replace

```
P <- P D,       L <- D^-1 L,       D[j,j]=e_j.
```

This sends `C` to `D^-1 C D` and normalizes its four anchors to one.
It preserves `LP=I`, energy homogeneity and the core. Since `D` factors
into independent physical characters, its action on each reshuffled cut
is by invertible constant diagonal row and column matrices. Hence it
preserves the complete Smith profiles, including the rank conditions.
Thus two UNSAT results classify the entire finite family up to physical
phases; no list of individual phase assignments needs enumeration.

| Field | Any zero anchor? | Nonbaseline with anchors one? | Outcome |
| --- | --- | --- | --- |
| F2 | UNSAT | UNSAT | Only baseline phases |
| F4 | UNSAT | UNSAT | Only baseline phases |
| F8 | UNSAT | UNSAT | Only baseline phases |
| F16 | UNSAT | UNKNOWN: 180-second timeout | Classification open |

These are solver results for finite fields defined by the recorded
irreducible binary polynomials. No characteristic-zero lifting or
algebraic-closure conclusion follows. Earlier unsplit orbit queries timed
out over F8 and F16; those timeouts are not negative results.

## Checks and the resulting balanced floor

Every field multiplication table is independently compared with polynomial
reduction in `Field`. Each solver run also fixes the maps to known physical
phase examples, checks feasibility, and checks that the phase-exclusion
predicate rejects them. The separate Python receipt checker reconstructs
ten such tensors from literal source entries over the four fields. It
checks 48,640 coefficients, 2,560 endpoint entries and 60 exact matrix
ranks. In particular it checks the block-Toeplitz identity `rank T2=3`
independently of the driver equations, along with the Schur residual.

At both cuts the baseline orders are

```
0,1,3,4,5,6,7,8,8,9,10,11,12,13,15,16.
```

Let `B(x)` be their exponential-sum polynomial, let
`Q(x)=1+x+x^3+x^4+12x^10` be the earlier optimistic cut bound, and let
`F(x)=1+x+x^3+x^4`. Exact integer Bernstein coefficients prove

```
(B(x^5) Q(x^5))^7 > F(x^7)^20,        0 < x < 1.
```

After removing its leading monomial, the certificate has degree 905 and
minimum integer Bernstein coefficient 14. The receipt checker rebuilds
this certificate. Thus even retaining one complete baseline profile and
only the optimistic other-cut bound gives `rho >= 5/7` at all balanced
natural scales. For F2/F4/F8, any protected split either violates a rank
condition and is covered by the earlier exact guard, or satisfies them
and has baseline profiles by the solver classification. This closes the
declared finite-field family **conditional on the recorded UNSAT results**.

Reproduce the receipt, positive controls and integer inequality with:

```
make verify-joint-split-flags
```

Rerun an actual solver query in an environment containing the optional
solver and the repository's numerical dependencies, for example:

```
python tools/probe_joint_split_rank_flags.py --field-degree 3 --tensor-core --query zero-anchor --timeout-ms 180000 --output /tmp/joint-f8-zero.json
python tools/probe_joint_split_rank_flags.py --field-degree 3 --tensor-core --query normalized --timeout-ms 180000 --output /tmp/joint-f8-normalized.json
```

`--field-degree` accepts 1,2,3,4. Output includes the exact SMT-LIB formula
beside its JSON receipt. The finite solver work and this focused replay
postdate the clean full default run at `decce64`; no later full-suite pass
is claimed. No Lean work or stronger factoring exponent is established.
