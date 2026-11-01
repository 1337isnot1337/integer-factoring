# Optimize a Smith channel basis at each lattice state

## Result and scope

The [state-adaptive channel construction](state-adaptive-smith-channels.md)
can also choose a new basis inside each equal-order Smith block. For a
given positive state potential, the best such basis is a minimum-weight
basis of a vector matroid. The greedy algorithm finds it exactly.

This yields a stronger companion-crossing control. At the rational scale
`s=5`, ten full-tail lattice states at precision two certify

```math
B_{\mathrm{adaptive}}\le\frac{24073}{1875}
=12.838933\ldots
<12.9
\le\inf_{s\ge1}B_{\mathrm{fixed}}(s).                 \tag{1}
```

`B_fixed` refers to the specific fixed-channel representation in
[exact rational module orbits](exact-rational-module-orbits.md). Its
infinite-precision orbit is closed, and its boundary rate agrees with its
free word rate by [minimal boundary completeness](minimal-boundary-module-completeness.md).
Thus neither choosing another scale nor increasing that representation's
module precision recovers (1).

The earlier unrotated adaptive control proved a strict comparison only
at the same scale. The present result compares against every scale of
the preceding fixed-channel certificate. It still does not improve the
factoring exponent: its rate product exceeds both the elementary scalar
certificate `256/27` and the required binary threshold `2^(2 sqrt(2))`.
No optimum over other local representations or actual coefficient ranks
is asserted.

## 1. Equal-order channels form a vector matroid

Work over `F_2[[tau]]`, and retain the notation of the state-adaptive
construction. At a state `M`, flatten the full physical action `B_M K`.
Minimum-valuation Schur elimination writes it as

```math
A_M=\sum_\nu\tau^\nu L_{M,\nu}R_{M,\nu}.             \tag{2}
```

For one order `nu`, let the block have `r` channels. The `r` columns of
`L` and the `r` rows of `R` are regular. The columns of `L` are linearly
independent modulo `tau`: their entries on the ordered pivot rows form
a triangular matrix with unit diagonal. Every nonzero binary combination
of those columns therefore contains a unit.

For `v` in `F_2^r` other than zero, reshape `L v` into auxiliary rows and
form their exact row lattice. At fixed precision `P`, enlarge it by
`tau^P O^D`; call the resulting state `M(v)`. The selected construction
uses a finite family of states closed under all these targets.

Let `v_1,...,v_r` be any basis of `F_2^r`, and let `S` have these vectors
as its rows. Then

```math
LR=(LS^T)((S^{-1})^T R).                                \tag{3}
```

The new right factors are regular scalar functions of the opposite
physical half. The new left factors have regular coordinates in the
states `M(v_i)`. Hence (3) is another valid local state-adaptive identity.
It leaves the crossing, the endpoint order, and the coefficient-rank
argument unchanged.

Fix a positive state potential `h`. At zero this block contributes

```math
s^{-\nu}\sum_{i=1}^r h_{M(v_i)}
```

to the potential image of state `M`. At infinity its common prefactor is
`s^(L-nu)`. In either case, choosing the best basis means minimizing
`sum_i h_(M(v_i))` over all binary bases.

The ground set is the nonzero vectors of `F_2^r`, with independence given
by linear independence and weight `h_(M(v))`. Sort these vectors by
weight, then accept a vector exactly when it increases the span. The
usual exchange argument proves that this greedy choice minimizes the
basis weight. The implementation uses exact rational weights and fixed
tie-breaking by the binary vector mask.

There are at most `2^r-1` candidates. In the binary physical interface,
the flattened matrix has four columns, so `r<=4` and at most fifteen
candidates suffice. This avoids enumerating the 20,160 ordered bases of
`F_2^4` for each state and potential.

## 2. Positive potentials select an exact stationary certificate

For each state and order block, choose the greedy basis for `h`. Let
`T_h(s)` be the resulting transition matrix. Equivalently, define

```math
(F_s(h))_M
=\sum_{\text{blocks }(M,\nu)}
  s^{\epsilon L-\nu}
  \min_{\{v_1,\ldots,v_r\}\text{ basis}}
        \sum_i h_{M(v_i)},                              \tag{4}
```

where `epsilon=0` at zero and `epsilon=1` at infinity. Then
`T_h(s) h=F_s(h)` exactly.

Consequently a rational inequality `F_s(h)<=beta h` provides a fixed
finite transition matrix with `T_h(s)h<=beta h`. The factorization
selected at a state is used every time that state appears. Selection
never depends on a physical input word, so the separated-path proof
of coefficient ranks continues to apply.

There is also an exact finite-policy characterization. Let `A` contain
all transition matrices obtained by independently choosing one allowed
basis for every state and block. Then

```math
\inf_{h>0}\max_M\frac{(F_s(h))_M}{h_M}
=\min_{T\in A}\rho(T).                                 \tag{5}
```

For one direction, the selected matrix `T_h` has spectral radius at
most `max F_s(h)/h`. For the other, `F_s(h)<=Th` for every allowed
matrix `T`. Taking the infimum over positive `h` and using the
Collatz--Wielandt formula gives the reverse inequality. The formula for
reducible nonnegative matrices follows, for example, by applying the
positive resolvent potential `(beta I-T)^(-1) 1` for `beta>rho(T)` and
letting `beta` decrease to the spectral radius.

Equation (5) concerns the spectral radius on the complete declared state
set. It does not assert an optimum for one initial state, for arbitrary
changes of lattice bases, for nonconstant channel rotations, or for all
possible time-dependent factorizations.

## 3. Ten-state exact witness

Use the same binary companion crossing as the preceding adaptive note:

```math
K(z)=\begin{pmatrix}z(C+E)&C\\zE&0\end{pmatrix},
\qquad f=(1,0),\quad e=(1,1)^T,
```

where `C=[x=u][y=v]` and `E=[x=v][y=u]`. Apply the same constant
auxiliary conjugation, so both boundaries become the first coordinate.

At precision `P=2`, all normalized full-tail lattices in `O^2` are the
following ten states, in the order used by the checker:

```math
\begin{gathered}
I,\quad
\begin{pmatrix}1&0\\0&\tau\end{pmatrix},\quad
\begin{pmatrix}1&1\\0&\tau\end{pmatrix},\quad
\begin{pmatrix}\tau&0\\0&1\end{pmatrix},\quad
\begin{pmatrix}1&0\\0&\tau^2\end{pmatrix},\\
\begin{pmatrix}1&1\\0&\tau^2\end{pmatrix},\quad
\begin{pmatrix}1&\tau\\0&\tau^2\end{pmatrix},\quad
\begin{pmatrix}1&1+\tau\\0&\tau^2\end{pmatrix},\quad
\begin{pmatrix}\tau^2&0\\0&1\end{pmatrix},\quad
\begin{pmatrix}\tau&1\\0&\tau\end{pmatrix}.
\end{gathered}
```

The start state is the fifth displayed lattice. The terminal orders are
`(0,0,0,1,0,0,0,0,2,1)` at both endpoints. The state family is closed
under every nonzero combination in every equal-order block.

At `s=5`, use the integer potentials

```math
\begin{aligned}
h_\infty={}&(10000,5320,10000,10000,5320,
                 10000,5320,10000,10000,10000),\\
h_0={}&(10000,9365,6316,6168,9365,
             6211,9365,6316,5423,6062).
\end{aligned}
```

Exact greedy selection and rational arithmetic verify

```math
T_\infty h_\infty\le\frac{266}{25}h_\infty,
\qquad
T_0h_0\le\frac{181}{150}h_0.                             \tag{6}
```

In fact the largest zero-side row ratio is `7619/6316`, strictly below
the slightly simpler rate used in (6). The boundary prefactors in the
earlier all-length theorem are `C_infinity=1` and `C_0=9365/6211`.
Thus the rank-sum rate product is bounded by `24073/1875`.

The evidence file records every chosen basis, every regular coordinate
factor, and both transfer matrices. The checker reconstructs all 1,280
entries of the selected local identities. It also checks the chosen
basis weights against all 8,690 independent unordered competing bases
across the declared blocks.

## 4. A global lower bound for the old fixed-channel certificate

The previous exact factorization has an infinity orbit of size two and
a zero orbit of size seven. The latter has a four-state exact lumping.
The checker rebuilds those orbits and verifies the lumping as Laurent
polynomial identities in `s`. The matrices are

```math
U(s)=\begin{pmatrix}0&4s\\3&s\end{pmatrix},\qquad
V(s)=\begin{pmatrix}
s^{-1}&1&2s^{-1}&0\\
s^{-1}&0&3s^{-1}&0\\
0&1&2s^{-1}&s^{-1}\\
0&1+s^{-1}&2s^{-2}&0
\end{pmatrix}.                                          \tag{7}
```

The quotient and original zero matrix have the same Perron root:
the quotient's positive Perron vector lifts to a positive eigenvector
of the original matrix. This also follows directly from the connected
transition graphs. Hence
`B_fixed(s)=rho(U(s)) rho(V(s))`.

Put `h=rho(U(s))`. Its quadratic equation gives

```math
h^2-sh-12s=0,\qquad s=\frac{h^2}{h+12}.
```

For `s>=1`, we have `h>=4`. With `y=s^-1`, the characteristic polynomial
of `V` is

```math
P(x,y)=x^4-3yx^3+(-4y+2y^2-2y^3)x^2
       +(-3y^3+2y^4)x+y^3+3y^4.                        \tag{8}
```

Take `b=129/10`. Clearing the positive denominator after substituting
`x=b/h` and `y=(h+12)/h^2` gives

```math
Q(h)=10000h^9P\left(\frac{b}{h},\frac{h+12}{h^2}\right)
```

with

```math
\begin{aligned}
Q(h)={}&5349888000-3345753600h-1676246400h^2
       +247845600h^3\\
      &-704173440h^4+135976611h^5-6646400h^6.             \tag{9}
\end{aligned}
```

The polynomial is strictly negative for all `h>=4`. An exact compact
certificate uses `h=4+t/(1-t)` and the degree-six polynomial

```math
G(t)=-(1-t)^6 Q\left(4+\frac{t}{1-t}\right).
```

Dyadic de Casteljau subdivision gives positive Bernstein coefficients
on each interval in the following cover. The table gives the least
coefficient on that interval; the evidence file contains all seven.

| Interval | Least Bernstein coefficient of `G` |
|---|---:|
| `[0,1/2]` | `141335830625/64` |
| `[1/2,3/4]` | `202499261963/4096` |
| `[3/4,13/16]` | `111301720877375/16777216` |
| `[13/16,7/8]` | `45184387297/524288` |
| `[7/8,1]` | `89533909679/262144` |

Therefore `det((b/h)I-V(s))<0`. Since the characteristic polynomial is
positive for sufficiently large positive arguments, it has a positive
real root greater than `b/h`. Every eigenvalue has modulus at most the
Perron root, so `rho(V(s))>b/h`. Multiplying by `h` proves

```math
B_{\mathrm{fixed}}(s)>\frac{129}{10}\qquad(s\ge1).      \tag{10}
```

This proof uses no floating-point optimizer, mesh assumption, or claim
that the characteristic polynomial has only one positive root.

## 5. Reproduce the certificate

```sh
python3 tools/check_matroid_smith_channels.py \
  --output evidence/matroid-smith-channels.json
```

The checker and engine use only the Python standard library. The evidence
includes the full policy, state data, exact factors, integer potentials,
global polynomial proof, and source hashes. Literal polynomial products
independently check 4,368 endpoint entries and eighteen coefficient-rank
inequalities through length three. Those finite checks validate the
implemented indexing; the local identities and positive potentials prove
the bounds at every length.

Files:

- [Engine](../tools/matroid_smith_channels.py).
- [Exact checker](../tools/check_matroid_smith_channels.py).
- [Certificate](../evidence/matroid-smith-channels.json).

The remaining useful question is whether state-dependent refactorization
and basis selection can improve a competitive local family. This control
establishes additional certificate strength, but its absolute pressure
does not yet approach the half-rate threshold.
