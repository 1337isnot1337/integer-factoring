# Moment identities and a cumulative bosonic crossing profile

## Result and scope

The full rational Appendix B profile has an all-weight first-moment
relation. Its second moment is a scalar operator on each fixed-energy
state space. The scalar pairing has an explicit Virasoro formula and an
exact rank:

```math
\sum_s sQ_s(a,b)=0,\qquad
\sum_s\binom{s}{2}Q_s(a,b)=2E\,\Gamma_2(a,b)\,I_{V_E}.       \tag{1}
```

For the source dual substitution `Theta`, restricted to parameters of
weight `h`, the rank of `Gamma_2(Theta a,b)` is

```math
r_h=p(h)-p(h-2),\qquad p(0)=1,\quad p(n)=0\ (n<0).         \tag{2}
```

These facts give a new complete cumulative profile. Put

```math
B_s=\sum_{t\le s}tQ_t=-\sum_{t\ge s+1}tQ_t.
```

The two expressions give both channel factorizations, and

```math
\sum_{s=-E-h}^{E+h-1}B_s=-4E\,\Gamma_2\otimes I_{V_E}.     \tag{3}
```

After choosing and normalizing a nonsingular block of `Gamma_2`, this is
an exact crossing identity. For `E=h`, its parameter dimension can be
`r_h`; its declared row and column channels have top energies
`2h+s` and `2h-s-1`, respectively. The identical middle pieces
`B_(-1)=B_0` have both channels in energies `1,...,2h-1`.

This is an ordinary proof with exact finite checks. It does not establish
a stronger factoring exponent. In particular, the leading channel
contraction still approaches the existing half-rate frontier. The rank
formula concerns the scalar second-moment pairing, not the ranks of every
complete cumulative piece. The checker measures those complete ranks
separately at small weights.

The polynomial ring, coefficient operators, source duals, and exchange
identity are Anthropic's Appendix B construction at revision
[`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The new moment identities, exact second-moment rank argument, and
cumulative profile below used OpenAI assistant support. No new Lean proof
is claimed.

## 1. Definitions and the fixed-energy Laurent polynomial

Work in `V=Q[t_1,t_2,...]`, with `t_m` of weight `m`. Let `V_E` be its
weight-`E` subspace. The source linear coefficient operator is

```math
t_m[r]=\binom{r-1}{m-1}\nabla_{-r},
\qquad
\nabla_{-k}=M_{t_k},\quad \nabla_k=k\partial_{t_k}\ (k>0),
\quad\nabla_0=0.                                      \tag{4}
```

Monomial parameter operators are normally ordered products of these
linear fields. Normal ordering places every multiplication before every
derivative. Let `b` be homogeneous of weight `h`, and let `a` have weight
at most `h`. Write `a=sum_w a_w` for its homogeneous decomposition and

```math
a^{\langle i\rangle}=\sum_w\binom{h+w}{i}a_w.
```

The source profile is

```math
Q_s(a,b)=\sum_{i=0}^{2h}(-1)^i
 a^{\langle i\rangle}[s-h+i]\,b[h-s-i].                \tag{5}
```

On `V_E`, exchanging the two coefficient operators in every summand
preserves `Q_s`; also `Q_s=0` for `|s|>E+h`. These are the source
[Lemma Q](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Series/ReadOff.lean)
properties. The coefficient and dual definitions are in
[AppHalfPiecesStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).

For monomials `a,b` of weights `w,h`, the source constructs one Laurent
polynomial by multiplying their two-field matrix coefficient by
`(z_1-z_2)^(w+h)`. Setting `z_1=z,z_2=1` gives

```math
F_{a,b;E}(z)=\sum_{s=-E-h}^{E+h}Q_s(a,b)|_{V_E}\,z^s.
                                                               \tag{6}
```

The construction and coefficient read-off are in
[Expansion.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Series/Expansion.lean),
[OneLaurent.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Series/OneLaurent.lean),
and the preceding read-off file. For an inhomogeneous `a`, add (6) over
its homogeneous components, each with its own regularizing exponent.

The constant value is `C=F(1)=Gamma_0(a,b) I_(V_E)`, where `Gamma_0`
is the full contraction pairing, equivalently the source vacuum pairing.

## 2. The first and second moments

For two linear parameter parts `m,n`, the contraction of an annihilator
from the first field with a creator from the second is

```math
\gamma_{mn}(z_1-z_2)^{-m-n},\qquad
\gamma_{mn}=(-1)^{m-1}
 \frac{(m+n-1)!}{(m-1)!(n-1)!}.                       \tag{7}
```

One can verify (7) directly from (4), using the binomial series for
`(1-z_2/z_1)^(-m-n)`. Expanding a product of fields by commuting
annihilators past creators gives a sum over partial matchings between the
parameter parts. Each contracted pair contributes (7), and the unmatched
parts remain normally ordered. This is also the Wick expansion used by
the source's Laurent construction.

Let `W` be the total weight of the unmatched parts. After multiplying by
the regularizing power, each such term has the form

```math
\Gamma\,(z_1-z_2)^W N(z_1,z_2).                       \tag{8}
```

Here `N` is the fixed-input, fixed-output coefficient of the remaining
normal product. It is a Laurent polynomial: a surviving creator has
mode at most the output energy, and a surviving annihilator has mode at
most the input energy. Thus Taylor expansion at `z_1=z_2=1` is valid
entry by entry.

If `W=0`, the term is a scalar identity. If `W=1`, the only unmatched
field is `t_1`. A coefficient of this single field between states of
the same energy is its zero mode, which is zero by `nabla_0=0`.
Consequently every nonconstant contribution in (8) has a factor
`(z_1-z_2)^2`. This proves

```math
F(z)-C\in(z-1)^2\operatorname{End}(V_E)[z,z^{-1}],
\qquad \sum_s sQ_s=0.                                \tag{9}
```

For `W=2`, a single `t_2` field again has zero energy-preserving mode.
The other possibility is two `t_1` fields. At coincident arguments the
energy-preserving part of their normal product is

```math
2\sum_{k\ge1} k t_k\partial_{t_k}=2H,
```

which equals `2E I` on `V_E`. Let `partial_1^[j]=partial_(t_1)^j/j!`
be the Hasse derivative. The sum of contraction coefficients leaving
two `t_1` fields is

```math
\Gamma_2(a,b)=
 \Gamma_0(\partial_1^{[2]}a,b)
 +\Gamma_0(\partial_1a,\partial_1b)
 +\Gamma_0(a,\partial_1^{[2]}b).                       \tag{10}
```

The full contraction pairing in (10) is defined for arbitrary weights;
it is zero if no full matching exists. Repeated parts are accounted for
by the Hasse derivatives. Terms of unmatched weight greater than two do
not contribute to the second Taylor coefficient, and the weight-one
terms vanish identically. Hence

```math
[u^2]F(1+u)=\sum_s\binom{s}{2}Q_s=2E\Gamma_2(a,b)I.
                                                               \tag{11}
```

Generalized binomial coefficients in this formula have their usual
integer values even when `s` is negative. All identities above hold
integrally for integral parameters. Reduction modulo a prime is valid,
but division by the integer factors in later normalizations requires
that those factors be nonzero in the chosen field.

## 3. An exact Hermite normal form, including its contact terms

Put `m=E+h`. Equation (9) gives a unique polynomial `H(z)` of degree at
most `2m-2` such that

```math
z^mF(z)=Cz^m+(z-1)^2H(z).                            \tag{12}
```

Equivalently, write

```math
\frac{F(z)-C}{(z-1)^2}=\sum_{s=-m}^{m-2}H_s z^s.
```

Finite polynomial division gives both formulas

```math
\begin{aligned}
H_s&=\sum_{t\le s}(s-t+1)Q_t-(s+1)_+C,\\
H_s&=\sum_{t\ge s+2}(t-s-1)Q_t-(-s-1)_+C.
\end{aligned}                                       \tag{13}
```

The notation `x_+` means `max(x,0)`. The agreement of the two expressions
uses both `sum Q_t=C` and `sum tQ_t=0`. In particular,

```math
H_{-1}=\sum_{t<0}(-t)Q_t=\sum_{t>0}tQ_t,
\qquad \sum_s H_s=2E\Gamma_2\otimes I.               \tag{14}
```

The subtraction of `C` in (13) matters for rank arguments. Only `s=-1`
has no such term in either expression. A short factorization of the
prefix sum does not, by itself, factor its difference with `C` through
the same channel. Also, evaluating an infinite intermediate-state
expansion at `z=1` and interchanging its sums is invalid; that operation
would incorrectly cancel the full contraction term itself. The finite
formulas (12)--(14) avoid that interchange.

## 4. The second-moment pairing is an explicit rank-deficient operator

Use the Hall inner product

```math
\langle t^\mu,t^\nu\rangle=z_\mu[\mu=\nu],\qquad
z_\mu=\prod_{k\ge1}k^{m_k(\mu)}m_k(\mu)!.
```

The source dual substitution is

```math
\Theta(t_m)=(-1)^{m-1}\sum_{j=1}^m\binom mjt_j.
```

Define

```math
D=\sum_{a\ge2}a t_{a-1}\partial_{t_a},\qquad
U=\sum_{a\ge1}a t_{a+1}\partial_{t_a},\qquad
H=\sum_{a\ge1}a t_a\partial_{t_a}.                    \tag{15}
```

These are the derivations used in
[the homogeneous-dual construction](homogeneous-duals.md); they satisfy
`[D,U]=2H` and `D*=U`. For `k>=1` put

```math
L_{-k}=\frac12\sum_{a=1}^{k-1}t_at_{k-a}
       +\sum_{a\ge1}a t_{a+k}\partial_{t_a},
\qquad L_0=H,\quad L_1=D.                            \tag{16}
```

For a homogeneous polynomial of positive weight `h`, define

```math
K_h=\sum_{k\ge0}\frac{(-1)^k}{k!}L_{-k}D^k.
                                                               \tag{17}
```

The sum terminates before `k=h`. It is an explicit rational operator
obtained from sparse differentiations and multiplications.

### The pairing formula

For all `a,b in V_h`,

```math
\boxed{\Gamma_2(\Theta a,b)=\langle a,K_hb\rangle.}     \tag{18}
```

Here is a direct derivation. A binomial finite-difference identity in
(7) gives

```math
\Gamma_0(\Theta t_m,t_n)=m\binom nm.
```

Thus the full contraction is the Hall pairing with
`Y=exp(D)`, where `Y(t_n)=sum_(m=1)^n binom(n,m)t_m`.
Introduce a variable `u` which marks unmatched `t_1` fields. An unmatched
field on the `Theta a` side is represented, in the second argument of
the Hall pairing, by multiplication by

```math
P=\sum_{m\ge1}(-1)^{m-1}t_m.
```

An unmatched field on the `b` side is represented by `partial_(t_1)`.
Consequently the pairing leaving two unmatched fields is represented by

```math
[u^2]e^{uP}Y e^{u\partial_{t_1}}.
```

The equality

```math
Y\partial_{t_1}Y^{-1}
 =P^*=\sum_{m\ge1}(-1)^{m-1}m\partial_{t_m}
```

turns this expression into

```math
\left(\tfrac12P^2+PP^*+\tfrac12(P^*)^2\right)e^D.
```

This is the normally ordered square of the current at `z=-1`, divided
by two. Its homogeneous energy-shift components are `(-1)^nL_n`.
Projection from `V_h` back to `V_h` selects `n=-k` from the term
`D^k/k!`, giving (17). All manipulations can be performed at a fixed
weight cutoff; no divergent product of unbounded series is required.

### Two derivatives annihilate its image

The explicit formulas (15)--(16), or their commutators on generators,
give

```math
[D,L_n]=(1-n)L_{n+1}.                                \tag{19}
```

Only `n<=1` is needed here. Put
`S_r=sum_(k>=0)(-1)^k L_(r-k)D^k/k!`. Moving the leftmost `D` to the
right and shifting one summation index yields

```math
DS_r=(1-r)S_{r+1}.
```

Since `K=S_0`, this proves

```math
D^2K_h=0.                                           \tag{20}
```

### Exact rank and spectrum

Let `Prim_j=ker(D:V_j->V_(j-1))`. The positive Hall form and `D*=U`
give the orthogonal primitive decomposition from the homogeneous-dual
note. In particular, for `h>=3`,

```math
\ker D^2|_{V_h}=Prim_h\oplus UPrim_{h-1},
\qquad \dim\ker D^2=p(h)-p(h-2).                     \tag{21}
```

The vacuum exception is relevant at `h=1,2`: `U(1)=0`.
For `v in Prim_h`, equation (17) gives `K_hv=hv`. For
`v in Prim_(h-1)`, use `DUv=2(h-1)v` and `D^2Uv=0` to obtain

```math
K_hUv=(2-h)Uv.                                      \tag{22}
```

When `h>=3`, both eigenvalues are nonzero. Thus (20) and (22) prove
`im K_h=ker D^2`, and (2) follows. Moreover,

```math
\begin{array}{c|c}
\text{eigenvalue of }K_h&\text{multiplicity}\hline
h&p(h)-p(h-1)\\
2-h&p(h-1)-p(h-2)\\
0&p(h-2)
\end{array}                                         \tag{23}
```

These are actual eigenspaces, since the image is the sum of the first
two and is disjoint from the kernel. At `h=1`, `K_1=I`; at `h=2`,
its matrix on `(t_1^2,t_2)` is `diag(2,0)`. These cases also satisfy
(2).

For `h>=3`, an explicit projector onto the image is

```math
\Pi_h=\frac{K_h(K_h-2I)}{h(h-2)}.                     \tag{24}
```

This is not a claim that `K_h` is a polynomial in `UD`, nor that its
kernel is `U^2V_(h-2)`. For example, `K_3(t_3)=3t_1^3/2`, which
exhibits a component changing the number of parts.

## 5. A full cumulative crossing profile with two short cuts

Let `m=E+h`. Define the finite family

```math
B_s=\sum_{t=-m}^s tQ_t,
\qquad -m\le s\le m-1.                              \tag{25}
```

Equation (9) gives the second expression in the introduction. Summing
(25) and changing the order of these finite sums gives

```math
\sum_{s=-m}^{m-1}B_s
 =\sum_{t=-m}^{m}(m-t)tQ_t
 =-\sum_t t^2Q_t
 =-4E\Gamma_2\otimes I.                              \tag{26}
```

The last step uses (9) and (11). Unlike a general Hermite coefficient,
each `B_s` therefore has a prefix and a suffix expression with no
full-crossing contact term.

For explicit channel formulas, take `E=h` and a homogeneous component
`a_w`; write `L=h+w`. For the source duals at positive weight,
`L>=2`. Let `q` be an intermediate energy, and interpret every binomial
coefficient outside `0<=k<=n` as zero. The row-prefix expression
uses the coefficient

```math
\theta_{s,q}(a_w)
 =(-1)^k\left[q\binom{L-1}{k-1}
                    -L\binom{L-2}{k-2}\right]a_w,
\qquad k=q-s.                                       \tag{27}
```

Its intermediate energies satisfy `0<=q<=2h+s`. The column-suffix
expression uses

```math
\psi_{s,q}(a_w)
 =-(-1)^k\left[(2h-q)\binom{L-1}{k}
                    -L\binom{L-2}{k-1}\right]a_w,
\qquad k=2h-q-s-1,                                  \tag{28}
```

with `0<=q<=2h-s-1`. If `k<0`, set these expressions to zero.
Both formulas follow from finite alternating binomial sums. For example,
the left side of (27) before summing is

```math
\sum_{i\ge q-s}(q-i)(-1)^i\binom Li a_w.
```

For fixed `q`, the other operator factor is independent of `i` and `t`.
Thus all equal-energy channels combine into one copy of `V_q`; one
does not multiply its dimension by the number of coefficient indices.

In coordinates, for source parameters `a_x,b_v` and weight-`h` state
labels `y,u`, the row expression is

```math
B_s(x,y,u,v)=\sum_{q=0}^{2h+s}\sum_{\tau\vdash q}
 [t^\tau]\theta_{s,q}(a_x)[q-h]t^y\,
 [t^u]b_v[h-q]t^\tau.                               \tag{29}
```

The column expression exchanges the order, with `psi_(s,q)(a_x)`
acting from `V_q` to `V_h` and `b_v` acting from `V_h` to `V_q`.
These are factorizations of the complete pieces, not lower bounds
deduced from individual blocks.

Since the `t=0` summand of (25) vanishes, `B_(-1)=B_0`. Use `s=-1`
in the row factorization and `s=0` in the column factorization. At
energy `q=0` the corresponding coefficient also vanishes, so both
channels have dimension at most

```math
\sum_{q=1}^{2h-1}p(q).                              \tag{30}
```

The central Hermite coefficient in (14) is `-B_(-1)`, giving the
same two factorizations.

### Turning the scalar pairing into a crossing

For each fixed `h>=1`, form the matrix of `Gamma_2` on the original
parameter labels. By (2), a nonsingular `r_h` by `r_h` minor exists.
Choose its row and column labels and apply its inverse as a basis change
on the first parameter slot. Choose any `r_h` independent states of
energy `h`, with the corresponding coordinate functionals. Equation
(26), followed by multiplication by `-1/(4h)`, is the dimension-`r_h`
crossing. Restriction and invertible changes in individual tensor slots
do not increase either channel rank in (29).

This is an explicit finite construction: exact elimination locates and
normalizes the minor. The statement alone does not give a favorable
uniform coefficient-height bound for that elimination, and it does not
replace those costs in a factoring-machine transfer theorem. Any fixed
member can instead be treated as an explicit finite seed, with all its
constants included.

## 6. Rank implications and the higher-moment obligation

The exact scalar rank (2) removes most parameter directions at a fixed
weight in the ratio sense, but not in the leading exponential entropy.
For fixed energy loss, partition asymptotics give
`log r_h=log p(h)+O(log h)`. The middle channel in (30) still has
`log dimension=sqrt(2) log p(h)+O(log h)`. Thus these declared
channels do not prove a score above `1/2`.

The explicit normal form changes what can be investigated next. One can
measure the complete ranks of the normalized cumulative pieces, search
their full coefficient span, or seek a higher-order analogue in which a
linear fraction of the intermediate energy is removed while the scalar
pairing retains almost all label entropy. Such an analogue needs an
additional theorem. Starting at unmatched weight three, normally
ordered zero modes need not be scalar on `V_E`; they cannot simply be
replaced by `E`. Likewise, subtracting lower moments creates contact
terms as in (13). A valid higher-order transfer must account for those
operators and both complete cuts.

## 7. Reproduction and finite evidence

Run from the repository root:

```bash
python3 tools/check_bosonic_hermite_profile.py \
  --output /tmp/bosonic-hermite-profile.json
```

The committed evidence is
[bosonic-hermite-profile.json](../evidence/bosonic-hermite-profile.json).
The checker independently constructs the source's integral coefficient
operators. It verifies both moments and the full Hermite reconstruction
through `h=4`, both independently collapsed cumulative cuts through
`h=3`, the Hall pairing representation through `h=9`, and
`D^2K=0` plus the predicted ranks and eigenspace dimensions through
`h=12` (77 original parameter labels). The larger rank checks are
over the stated prime; the rank theorem is the all-weight proof above.

For every complete normalized crossing piece through `h=4`, the checker
computes both ranks by exact rational elimination, rather than inferring
one from a factor dimension or from a single modular rank. The evidence
records the selected minor and state labels. Negative controls omit a
term from the Virasoro telescoping sum and test false scalar-moment
claims. No finite check is used as an asymptotic proof.
