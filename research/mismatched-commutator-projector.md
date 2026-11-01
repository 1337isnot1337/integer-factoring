# A constructive commutator projector for weights 1 and 3

The characteristic-two mismatch produces an explicit pairing from a short
noncommutative word.  This is an exact identity at every energy, not an
extrapolation from a finite matrix.  It is a building block, not a complete
crossing decomposition or an improved factoring bound.

## 1. Operators and the two moving modes

Use the [mismatched Weyl curve](mismatched-mode-weyl-crossings.md), with first two
weights 1 and 3 and all other weights greater than 3.  Write

\[
B(z)=P(1+z)=\sum_i B_i z^i,\qquad D=B_1,\qquad
H=DB_3+B_3D.                                      \tag{1}
\]

On a row state `(alpha,y)`, `D` is diagonal with value
`(alpha_0 mod 2)(y_0 mod 2)`.  The coefficient `B_3` has two diagonal terms
and two off-diagonal moves.  The diagonal terms commute with `D`.

Writing an input row state as `(beta,u)`, the forward move is

\[
U:(\beta,u)\longmapsto
(\beta+(3,-1),\ u+(3,-1)),
\]

with coefficient `[beta_1 odd][u_0=0 mod 4]`, provided `u_1>=1`.
The reverse move is

\[
V:(\beta,u)\longmapsto
(\beta-(3,-1),\ u-(3,-1)),
\]

with coefficient `[beta_0=3 mod 4][u_1 even]`, provided `u_0>=3`.
These formulas follow directly from the two cross terms in `C_3D_3`:
ordinary binomial coefficients on the first label and divided-power
creation coefficients on the second.  All unmentioned coordinates stay
fixed.  Both moves preserve each label's energy.

The commutator in (1) retains a move exactly when the input and output
values of `D` differ.  In particular, it vanishes when the two mode-zero
occupations have different parities.

## 2. The exact projector identity

Partition the following two subsets of the label set:

\[
F_0=\{\alpha:\alpha_0=0\pmod4,\ \alpha_1\text{ odd}\},\qquad
F_1=\{\alpha:\alpha_0=3\pmod4,\ \alpha_1\text{ even}\}.
\]

The shift `tau(alpha)=alpha+(3,-1)` bijects `F_0` with `F_1`.  The inverse
is well defined because every element of `F_1` has `alpha_0>=3`.
Let `m=|F_0|=|F_1|`.

Two forward moves vanish: after a forward move, the second mode-zero
occupation is 3 modulo 4, so it cannot move forward again.  Two reverse
moves vanish for the analogous first-coordinate reason.  Therefore the
only two-step paths are a move followed by its reverse.  Tracking the
four binomial parity conditions gives

\[
P:=H^2=\operatorname{diag}
 \bigl([\alpha,y\in F_0]+[\alpha,y\in F_1]\bigr).  \tag{2}
\]

This accounts for the boundary conditions as well: the odd mode-one
occupations in `F_0` are positive, and the mode-zero occupations in `F_1`
are at least 3.  Hence

\[
H^4=H^2,\qquad P^2=P.                           \tag{3}
\]

The cube `H^3=PH` swaps `F_0 x F_0` with `F_1 x F_1`, applying `tau` or
`tau^{-1}` simultaneously to both labels; it is zero elsewhere.  The
unpaired one-way moves of `H` disappear under this projection.

Thus the short word

\[
R=P+H^3=H^2+H^3                               \tag{4}
\]

is an embedded two-level pairing.  Under the explicit label bijection
`F_0 union F_1 = {0,1} x F_0`, it is the binary row pairing tensored with
the crossing on the `m` residual labels.  It consequently has exact cuts

\[
\operatorname{rank}_{\rm row}R=m^2,\qquad
\operatorname{rank}_{\rm col}R=4m^2.           \tag{5}
\]

The projector `P` itself has both cut ranks `2m^2`.

## 3. Verified finite instance and scope

At energy 63 with weights `(1,3,9,27)`, there are 132 labels and `m=19`.
Direct extraction of `B_1,B_3` from all 632732 supported polynomial entries
gives

| Operator | Row rank | Column rank |
|---|---:|---:|
| `H` | 2380 | Not needed |
| `P=H^2` | 722 | 722 |
| `R=H^2+H^3` | 361 | 1444 |

The independent local checker tests every pair of mode energies from 0
through 24, allowing different energies in the two labels.  This includes
the spectator-mode cases that an equal-local-energy test would miss.
An optional full-shell replay checks the source coefficients and the exact
projector and swap matrices against (2)--(4):

```sh
python tools/check_mismatched_commutator.py
python tools/check_mismatched_commutator.py --prefix PATH_TO_H63_DUMPS
```

At energy 63, this operator is provably outside the polynomial algebra of
`A=P(0)`.  Take the row pair
`alpha=y=(57,2,0,0)` and column pair `beta=u=(3,20,0,0)`.
The row pair is outside both active sectors, so `(RA)` has entry zero.
The column pair is in `F_1 x F_1`, and its partner is
`beta'=u'=(0,21,0,0)`.  The two contributing entries of `A` are

\[
\binom{60}{57}\binom{22}{2}=0,\qquad
\binom{57}{57}\binom{23}{2}=1\quad\text{in characteristic two}.
\]

Thus `(AR)` has entry one and `[A,R]` is nonzero.  Every polynomial in `A`
commutes with `A`, so this is a concrete exclusion of `R` from that algebra
at the stated shell.  No such exclusion is asserted for every small energy.

Once restricted to this one paired sector, (4) is still an ordinary binary
pairing, so it does not itself
escape known binary pairing barriers.  A useful next construction would
have to combine overlapping mode-exchange projectors or other coefficient
words into a complete crossing with a better aggregate rank profile.
