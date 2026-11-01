# Fixed-precision count lifting and local nilpotent corrections

## Results and scope

Two separate constructions transfer modular counting to a fixed prime
power. The first applies to arbitrary linear-volume deterministic gate
networks. The second preserves the leading crossing-width estimate of a
fixed local polynomial matrix-product family.

* For fixed prime `p` and precision `a>=1`, counts modulo `p^a` reduce to
  counts modulo `p` on `O_{p,a}(n)`-volume networks with unique auxiliary
  settings. Every fixed logarithmic saving exponent is preserved.
* Every fixed finite crossing seed over `F_p` lifts asymptotically to
  `Z/p^a Z` with the same leading rank-sum exponent. A fixed local
  polynomial matrix-product family, certified by two Smith profiles as in
  [smith-profile-crossing-transfer.md](smith-profile-crossing-transfer.md),
  has such a lift with the same leading `rho` as well.
  The supplied row and column factorization widths increase by only
  `O_a(log(t+2))` at word length `t`.

The second statement concerns explicit factorizations over a finite ring,
not ordinary field rank. It does not claim that arbitrary lifts of large
low-rank matrices retain their ranks. Both results require fixed
precision; their constants depend on `p` and `a`.

These are resource-preserving reductions. They do not supply a new
crossing family with `rho<1/sqrt(2)` or improve the established factoring
exponent of one half.

## 1. Recover count digits by counting sorted tuples

Let `C` be the number of accepted witnesses of a finite verifier. Write

```math
C=\sum_{j\ge0} c_jp^j,\qquad 0\le c_j<p.
```

Then

```math
c_j={C\choose p^j}\pmod p.                              \tag{1}
```

An elementary proof uses formal polynomials over `F_p`:

```math
(1+X)^C=\prod_{j\ge0}(1+X^{p^j})^{c_j}.
```

In coefficient degree `p^j`, all lower-index factors together have
degree at most `p^j-1`; higher-index factors cannot contribute. The
coefficient is exactly `c_j`. Thus (1) does not require an oracle for
binomial coefficients of an unknown count.

Fix a total order on the original complete witness encodings. The number
of strictly increasing `r`-tuples of accepted witnesses is
`binom(C,r)`. For `r=p^j`, count such tuples modulo `p`, for
`0<=j<a`, and reconstruct `sum_j c_j p^j` modulo `p^a`.

The construction is uniform. It does not enumerate the accepted
witnesses or know `C` in advance.

### Exact network budget

Use the deterministic-message gate model in
[weighted-network-top-compiler.md](weighted-network-top-compiler.md).
Suppose its total setting bits plus ordered incident-wire width are at
most `Cn`, its maximum gate load is at most `Cn^(3/4)`, and local queries
take polynomial time. A gate setting is encoded by its supplied unique
index; pad each index to its setting-bit width. The concatenation has
length at most `Cn` and gives a unique order even if some encodings are
unused.

Take `r` disjoint copies of the network. Expose the binary setting index
of every original gate to the comparisons. A gate in one copy has at
most two adjacent comparisons, so it transmits at most twice its setting
bit count in these new wires. Its maximum load stays within a constant
multiple of the original maximum.

Compare each consecutive pair of full encodings by a deterministic
three-state lexicographic scan. Its state is equal, smaller, or larger;
the initial state is equal and the final state is required to be smaller.
Every scan state is uniquely determined by the two input strings. A
separate constant-size gate per bit gives total setting and wire cost
`O(rn)` and constant comparison-gate load. Therefore the whole sorted
tuple verifier has

```math
W'=O(rn),\qquad D'=O(n^{3/4}),                         \tag{2}
```

with exactly `binom(C,r)` accepted complete settings. No multiplicity
comes from comparison witnesses.

For fixed `p,a`, every `r<=p^(a-1)` is a constant, so the weighted
compiler applies with an enlarged fixed constant. If the base evaluator
takes `poly(n) exp(O(F(n)))` on this class, the lift takes
`poly(n) exp(O(F(C_{p,a}n)))`. In particular, it preserves every fixed
target `F(n)=n/(log n)^alpha` and every fixed target `F(n)=n^beta` with
`0<beta<1`. This statement concerns resource closure for this network
class, not a new complexity-class identity.

## 2. A full-root detector over `Z/4Z`

For odd `N>1`, there are exactly `2^omega(N)` roots of one modulo `N`.
Count all roots with canonical representatives `0<=z<N`, and subtract
the two known roots `1` and `N-1`. The remaining count satisfies

```math
2^{\omega(N)}-2\equiv
\begin{cases}
0&\omega(N)=1,\\
2&\omega(N)\ge2
\end{cases}\pmod4.                                      \tag{3}
```

After the usual deterministic perfect-power preprocessing, this is a
valid composite detector. Prefix counts satisfy the exact child-sum
identity in `Z/4Z`, so a nonzero parent has a nonzero child. Subtract each
trivial root only from its own prefix. A recovered nontrivial root gives
the proper factor `gcd(N,z-1)`.

The quadratic verifier and its unique quotient are described in
[quadratic-congruence-factor-recovery.md](quadratic-congruence-factor-recovery.md).
Equation (3) removes the need to select one representative from each
negation pair. The range `0<=z<N` remains part of the verifier. It is
not permission to replace the integer root equation by a congruence
modulo four.

If only a characteristic-two evaluator is available, Section 1 computes
the necessary `Z/4Z` count using the original verifier and its sorted
two-copy version. This has a constant volume cost. A literal lift of
every field tensor identity is unnecessary for that reduction.

## 3. Every fixed finite seed has an asymptotic prime-power lift

Suppose a fixed alphabet size `q>=2` and a fixed list of nonzero field
tensors satisfy

```math
E_q=\sum_{i=1}^m A_i\quad\hbox{over }\mathbb F_p.
```

Supply row and column factorizations of lengths `r_i,c_i>=1`, and put

```math
S=\sum_{i=1}^m\sqrt{r_i c_i},\qquad b=q^2.
```

Choose arbitrary lifts `widetilde A_i` over `R=Z/p^a Z`. Their endpoint
error is divisible by `p`; select a tensor `D` such that

```math
E_q=\sum_i\widetilde A_i+pD.                           \tag{A1}
```

The supplied field row factorization of each `A_i` also lifts to

```math
\widetilde A_i=\widetilde L_i\widetilde R_i+pD_i^{\rm row}.
                                                               \tag{A2}
```

There is an independent analogous column factorization. Every residual
tensor has a full factorization of length at most `b` in either cut.
No compatibility between its two full factorizations is needed.

Take the `t`-fold tensor power of (A1). Group words by their counts of
each ordinary piece and the endpoint error piece `D`. Terms with at
least `a` copies of `pD` vanish. There are only polynomially many grouped
pieces, since `m` and `a` are fixed.

Consider one group with `j<a` endpoint-error sites and ordinary
multiplicities `k_1,...,k_m`, summing to `t-j`. In each row factorization
use (A2). At most `a-1-j` additional residual sites can survive. Replacing
one ordinary factor of length `r_i` by a full residual factor multiplies
its channel bound by at most `b`, since `r_i>=1`. Therefore the row
factorization has length at most

```math
{t!\over j!\prod_i k_i!}\,b^j\prod_i r_i^{k_i}
 \sum_{h=0}^{a-1-j}{t-j\choose h}b^h,                 \tag{A3}
```

with binomial coefficients interpreted as zero outside their ranges.
The column bound is identical with `c_i` in place of `r_i`. The factor
`p^(j+h)` is absorbed into one side of every factorization.

Take square roots of the products of these two bounds and sum over
ordinary multiplicities. The multinomial theorem gives the explicit
rank-sum estimate

```math
\begin{aligned}
\mathcal S_{R,t}
&\le S^t\sum_{\ell=0}^{\min(a-1,t)}
 {t\choose\ell}\left[b\left(1+{1\over S}\right)\right]^\ell\\
&=S^t\,O_{p,a,q,A}( (t+1)^{a-1}).                    \tag{A4}
\end{aligned}
```

For the first line, set `ell=j+h` in the double sum from (A3), and use
`binom(t,j) binom(t-j,h)=binom(t,ell) binom(ell,j)`.
This also bounds every individual square-root rank product. At binary
width `lambda=t log2(q)`, the leading parameter is unchanged:

```math
\rho={\log S\over 2\log q},\qquad
\log_2 r_R+\log_2 c_R
 \le4\rho\lambda+O_a(\log(t+2)).                      \tag{A5}
```

All local data have fixed finite size. Enumerating words, their grouped
coefficient tables, and the factor channels in (A3) takes `2^O(t)` actual
bit operations. Thus this is an effective asymptotic lift, even when no
one-copy lift preserves both field ranks. The polynomial matrix-product
argument below also covers the stronger Smith-profile bounds that need
not follow from the raw rank sum of the local coefficient pieces.

## 4. Lift a fixed local polynomial crossing exactly

Here is the more local result. Let the fixed data over `F_p` be an
external alphabet of size `q`, auxiliary dimension `D`, polynomial
matrix `K(z;x,y,u,v)` of degree at most `L`, and vectors `e,f` satisfying

```math
f^T e=1,\qquad K(1;x,y,u,v)e=E_q(x,y,u,v)e,             \tag{4}
```

where `E_q(x,y,u,v)=[x=v][y=u]`. Put `R=Z/p^a Z`. Choose arbitrary
coefficient lifts of `K` and `e`. Some coordinate of `e` is a unit in
`R`, because its reduction is nonzero. Choose a row vector `h^T` with
`h^T e=1` by inverting that coordinate. Modify a lift of `f` by

```math
\widetilde f^{,T}\longmapsto
\widetilde f^{,T}+(1-\widetilde f^{,T}e)h^T.
```

For each external label quadruple, modify the constant coefficient of
the lifted local polynomial by

```math
\widetilde K(z)\longmapsto
\widetilde K(z)+(E_q e-\widetilde K(1)e)h^T.             \tag{5}
```

Both corrections are divisible by `p`. They preserve the original
reduction, keep degree at most `L`, and make (4) hold exactly over `R`.
Consequently

```math
\widetilde f^{,T}\widetilde K(z)_1\cdots
  \widetilde K(z)_t e
  =\sum_{\ell=0}^{Lt} z^\ell\widetilde A_{t,\ell}
```

is an exact crossing decomposition over `R`: its coefficient tensors
sum to `E_(q^t)`. The auxiliary matrices retain their original position
order.

## 5. Nilpotent errors in the two local factorizations

Use the two rectangular flattenings from the Smith-profile note, with
both auxiliary indices on the same side of each cut. Put `b=q^2`.
The fixed field profiles are

```math
d_1,\ldots,d_b\le L,\qquad \nu_1,\ldots,\nu_b\ge0.
```

At infinity, write the field factorization as a sum of `b` outer
products with powers `w^(-d_i)` and outside formal series of nonnegative
order. Lift the outside series coefficientwise to `R`. The difference
from `widetilde K(w^-1)` is divisible by `p`, and multiplication by
`w^L` makes that difference a power-series matrix. Thus it has the form

```math
\widetilde M(w^{-1})
 =\sum_{i=1}^b w^{-d_i}\widetilde U_i(w)\widetilde V_i(w)
   +p\,w^{-L}H_\infty(w).                              \tag{6}
```

The residual `H_infinity` has an explicit factorization through at most
`b` channels: keep its entire rectangular matrix on the auxiliary side
and put an identity matrix on the other side. This is a full local
factorization, of fixed dimension, not a rank assertion about a large
word tensor.

Likewise, at zero the lifted field Smith factorization gives

```math
\widetilde M(z)
 =\sum_{i=1}^b z^{\nu_i}\widetilde U'_i(z)\widetilde V'_i(z)
   +p H_0(z),                                           \tag{7}
```

with nonnegative series orders and at most `b` residual channels.
Division by `p` in (6)--(7) means selecting any coefficient preimage
under multiplication by `p` for a known divisible coefficient. It does
not invert a nonunit.

In a product of `t` local matrices, a term with `j` residual sites has
a factor `p^j`. All terms with `j>=a` vanish. There are only
`sum_(j<a) binom(t,j)` residual-position sets. Noncommuting auxiliary
matrices stay in their original order throughout this expansion.

## 6. Explicit width bounds

For fixed real `s>=1`, put

```math
U(s)=\sum_i s^{d_i},\qquad V(s)=\sum_i s^{-\nu_i},
\qquad B(s)=U(s)V(s),
```

and define

```math
E_\infty(t)=\sum_{j=0}^{\min(a-1,t)}
 {t\choose j}\left(\frac{b s^L}{U(s)}\right)^j,
\qquad
E_0(t)=\sum_{j=0}^{\min(a-1,t)}
 {t\choose j}\left(\frac{b}{V(s)}\right)^j.              \tag{8}
```

There are explicit row and column factorizations of every lifted
coefficient tensor with channel counts bounded by

```math
\begin{aligned}
r_R(\widetilde A_{t,\ell})
 &\le (Lt-\ell+1)s^{-\ell}U(s)^t E_\infty(t),\\
c_R(\widetilde A_{t,\ell})
 &\le (\ell+1)s^\ell V(s)^t E_0(t).
\end{aligned}                                                    \tag{9}
```

Here `r_R,c_R` denote the lengths of supplied ring factorizations;
minimal inner rank is no larger.

To prove the first bound, choose `j<a` residual positions, their `b^j`
channels, and ordinary channels on the remaining sites. A residual
position contributes degree `L`; the other positions contribute their
`d_i`. The total must be at least `ell`. The Chernoff count is at most
`s^(-ell+jL) U(s)^(t-j)`. The two full-word outside factors split the
remaining series degree in at most `Lt-ell+1` ways. Their auxiliary
factors are ordered matrix products, so every such split remains one
outer product across the cut. Summing over residual positions gives
the first bound. At zero, residual positions have order zero, and the
same argument with total order at most `ell` proves the second.

Thus

```math
r_R(\widetilde A_{t,\ell})c_R(\widetilde A_{t,\ell})
 \le (Lt+1)^2 B(s)^t E_\infty(t)E_0(t).                 \tag{10}
```

For fixed precision `a`, both polynomials in (8) have degree at most
`a-1`. With `lambda=t log2(q)` and
`rho=log(B(s))/(4 log(q))`, taking binary logarithms in (10) gives

```math
\log_2 r_R+\log_2 c_R
 \le 4\rho\lambda+O_a(\log(t+2)).                     \tag{11}
```

Summing the square-root rank products over the `Lt+1` pieces also costs
only another polynomial factor. The leading `rho` is unchanged.

The fixed local Smith algorithms supply the needed series to precision
`O(t)`: all extracted coefficient indices and total outside degrees are
at most `Lt`. Coefficient construction and enumeration of the explicit
channel factorizations take `2^O(t)` operations in the fixed finite
ring. Their actual bit cost has the same form. This closes the
effectiveness requirement for the lifted algebraic family.

Any crossing evaluator whose proof uses these supplied factorizations
and ring addition and multiplication can use the lift. An evaluator
that invokes field division must first be checked or replaced at that
step; a nonunit is not an admissible pivot. The fixed-precision network
reduction in Section 1 is available independently of this algebraic
interface.

## 7. Limits and exact checks

The result does not lift an arbitrary global matrix of field rank `r`
to ring rank `poly(t)r`. For example, the zero matrix over `F_p` can
lift to `p I_m`, which requires `m` nonzero diagonal channels over
`Z/p^2 Z`. The proof succeeds because every error occurs in a fixed
local tensor before the `t`-site product is expanded.

For variable precision, retain (8) explicitly. Replacing it by
`poly(t)` is valid only when `a` is fixed. If `a` grows with `t`, the
number of surviving error-location sets can be exponential in `t`.
Likewise, Section 1 uses up to `p^(a-1)` copies and is not a
constant-volume reduction at growing precision.

The accompanying checker verifies (1), sorted-tuple multiplicities and
unique comparison traces, full-root recovery modulo four, endpoint
corrections, and literal one- and two-defect expansions of small lifted
polynomial crossings. These finite checks support the constructions;
the all-length and resource statements follow from the arguments above.

```sh
python3 tools/check_fixed_precision_count_lifting.py \
  --output evidence/fixed-precision-count-lifting.json
```

The binomial digit identity and nilpotent product expansion are standard
algebraic facts. The contributions here are their explicit linear-volume
network accounting and their application to local Smith-profile crossing
certificates. This work used substantial OpenAI assistant support.
