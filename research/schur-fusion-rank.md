# Asymptotically full Schur fusion rank

## Result

Let `Lambda_h` be the degree-`h` integral symmetric-function lattice,
reduced modulo two. Its dimension is the ordinary partition number
`p(h)`. For the source modes, put

```math
\mu_h:\Lambda_{\mathbf F_2,h}\otimes\Lambda_{\mathbf F_2,h}
        \longrightarrow\Lambda_{\mathbf F_2,2h},
\qquad \mu_h(a,b)=a[h]b.
```

There is an absolute constant `c>0` such that

```math
\operatorname{rank}\mu_h
\ge\bigl(1-O(e^{-c h^{1/4}})\bigr)p(2h).
\tag{1}
```

This resolves the full-space image question in
[schur-fusion-gram.md](schur-fusion-gram.md). At powers of two
`h=2^r`, the central coefficient of the full homogeneous Schur
crossing family has equal row and column ranks satisfying

```math
r_h=c_h=(1-o(1))p(2h).
\tag{2}
```

Consequently

```math
\frac{\log(r_hc_h)}{4\log p(h)}\longrightarrow\frac1{\sqrt2}.
\tag{3}
```

The actual central coefficient therefore enforces the half-rate width
on this full-domain construction. Improving only its rank estimate
cannot produce a better exponent.

The result concerns the full homogeneous Schur domain and its stated
coefficient pieces. It does not rule out further label restrictions,
new dual choices, or recombinations of different pieces. In particular,
the claim is not a lower bound on factoring or on arbitrary crossings.

The proof strengthens the
[regular quotient result](regular-fusion-rank.md). It uses the integral
monomial symmetric basis to retain information that the power-sum basis
loses modulo two. A partition can then be split through the binary
digits of its part multiplicities.

## 1. An integral filtration by power-sum length

Write `p_nu` for a power-sum monomial and `m_nu` for a monomial
symmetric function. Let `ell(nu)` be its number of parts. Define

```math
F_l\Lambda_{\mathbf Z,n}
=\Lambda_{\mathbf Z,n}\cap
 \operatorname{span}_{\mathbf Q}\{p_\nu:\ell(\nu)\le l\}.
\tag{4}
```

The exact integral basis of this sublattice is

```math
F_l\Lambda_{\mathbf Z,n}
=\operatorname{span}_{\mathbf Z}
   \{m_\nu:|\nu|=n,\ \ell(\nu)\le l\}.
\tag{5}
```

To see this, inclusion-exclusion on repeated variable indices gives

```math
m_\nu=
\frac{p_\nu}{\prod_j a_j!}
 +\text{terms of strictly smaller power-sum length},
\tag{6}
```

where `a_j` is the multiplicity of part `j` in `nu`.
More explicitly, if `nu` has `l` listed parts, the lower terms are
obtained by summing over set partitions `pi` of those `l` positions:

```math
m_\nu=\frac1{\prod_j a_j!}
\sum_{\pi}
\left(\prod_{B\in\pi}(-1)^{|B|-1}(|B|-1)!\right)
\prod_{B\in\pi}p_{\sum_{i\in B}\nu_i}.
\tag{7}
```

The discrete set partition gives (6); every other partition merges
some positions. Distinct `m_nu` of the same length have distinct top
power-sum monomials. Expanding an integral symmetric function in its
integral `m`-basis now proves (5), by considering its largest length.
In particular this filtration is saturated and can be reduced modulo
two without making a division by a factorial.

## 2. Fusion has ordinary product as its leading term

For partitions `A,B` of `h`, normal ordering gives

```math
\mu_h(m_A,m_B)-m_A m_B
\in F_{\ell(A)+\ell(B)-2}\Lambda_{\mathbf Z,2h}.
\tag{8}
```

First prove the claim over the rationals in power sums.
In a term with no annihilation, every created mode index is at least
its parameter index. The total creation weight is `h`, so all these
indices must be unchanged. The sum of these terms is exactly the
ordinary product. A term with `k>=1` annihilations has two fewer
particles per annihilation. Expanding (6) therefore gives the length
bound in (8).

Both fusion and ordinary multiplication preserve the integral lattice,
by [integral-schur-crossing.md](integral-schur-crossing.md).
Intersecting the rational length bound with the lattice and using (5)
proves the integral statement (8).

If `A,B` have part multiplicities `a_j,b_j`, their leading product is

```math
m_A m_B=
\left(\prod_j\binom{a_j+b_j}{a_j}\right)m_{A\cup B}
 +\text{terms of smaller length}.
\tag{9}
```

The coefficient counts the choices of which equal parts came from the
first factor. All other terms identify variable positions and reduce
the length. The coefficient in (9) is odd exactly when each binary
addition `a_j+b_j` has no carry. This parity identity follows directly
by comparing coefficients in
`(1+X)^m=prod_r(1+X^(2^r))^(bit_r(m))` over `F_2`.

Thus a carry-free equal-weight split of a partition `nu` gives

```math
\mu_h(m_A,m_B)=m_\nu+\text{terms of smaller length}
\quad\text{over }\mathbf F_2.
\tag{10}
```

Choose one split for each such `nu`. Ordering these distinct top
monomials by length produces a triangular minor with diagonal one.
If `b_2(2h)` counts the carry-free splittable partitions, then

```math
\operatorname{rank}\mu_h\ge b_2(2h).
\tag{11}
```

No containment of the ordinary product image in the fusion image
is asserted. Equation (10) supplies the required independence directly.

## 3. Multiplicity digits are independent weighted atoms

For an ordinary partition, write the multiplicity of part `j` as

```math
a_j=\sum_{r\ge0}\varepsilon_{j,r}2^r,
\qquad \varepsilon_{j,r}\in\{0,1\}.
```

Call each selected pair `(j,r)` an atom of weight `j*2^r`.
The tags distinguish atoms with the same numerical weight.
A subset of these atoms with total weight `h` gives exactly a
carry-free split into two partitions of weight `h`.

Fix a large even `N`, put `A=pi^2/6`, and set `t=sqrt(A/N)`.
The ordinary partition Boltzmann distribution makes the multiplicity
of each part an independent geometric random variable. Its binary
digits are independent Bernoulli variables, with

```math
\Pr(\varepsilon_{j,r}=1)
=\frac1{1+e^{t j2^r}}.
\tag{12}
```

This follows from the exact product identity

```math
\prod_{r\ge0}(1+x^{j2^r})=\frac1{1-x^j}.
\tag{13}
```

The selected atoms `(j,0)` with `j` odd have precisely the independent
odd occupations used in Sections 2--3 of
[regular-fusion-rank.md](regular-fusion-rank.md).
Use `K=floor(N^(1/4))` there. Its adjacent pairs and small selected
odd part produce a subset-sum interval. Its dense doubling intervals
extend the width beyond `N^(3/4)/100`, except with probability
`O(exp(-c_0 N^(1/4)))`.

The change from `pi^2/24` to `pi^2/6` in `t` does not affect that
argument: every needed occupation probability is still at least `1/3`
up to weight `2N^(3/8)`, for sufficiently large `N`.

All additional atoms can be adjoined with the same interval lemma.
Their weights need not be odd or distinct. It remains to ensure that
none exceeds `R=N^(3/4)/1000`.

A numerical weight `s` has `v_2(s)+1<=1+log2(s)` atom tags.
Using (12), a union bound gives

```math
\begin{aligned}
\Pr(\exists\text{ selected atom of weight }>R)
&\le\sum_{s>R}(1+\log_2s)e^{-ts}\\
&\le e^{-tR/2}\sum_{s\ge1}(1+\log_2s)e^{-ts/2}\\
&=O\!\left(t^{-2}e^{-tR/2}\right)
 =O(e^{-c_1N^{1/4}}).
\end{aligned}
\tag{14}
```

For the last estimate, use `1+log2(s)<=s+1` and absorb the
polynomial factor into the exponential. On the complementary event,
every remaining atom fits into the existing interval. The final
interval is `[L,H-L]`, where `H` is the total weight and
`L=O(K^2)=O(sqrt(N))`. It therefore contains `N/2` when `H=N`
and `N` is sufficiently large.

The unconditioned failure probability is `O(exp(-c_2 N^(1/4)))`.

## 4. Exact-weight conditioning costs only a polynomial factor

Let

```math
P(t)=\prod_{j\ge1}(1-e^{-tj})^{-1}.
```

Conditioning the atom model on total weight `H=N` gives the uniform
distribution on ordinary partitions of `N`. Its conditioning mass is

```math
\Pr(H=N)=p(N)e^{-tN}/P(t).
```

The decreasing function `-log(1-e^(-x))` has integral `A=pi^2/6`.
The right-endpoint Riemann sum gives the exact upper bound

```math
\log P(t)\le A/t.
\tag{15}
```

The Hardy--Ramanujan asymptotic is

```math
p(N)\sim\frac{e^{2\sqrt{AN}}}{4\sqrt3\,N}.
\tag{16}
```

This classical estimate is the partition input already used by the
pinned source and the existing rank notes. One primary exposition is
Section 4 of [Mark Wildon, *Counting partitions on the abacus*](https://www.ma.rhul.ac.uk/~uvah099/Maths/abacus.pdf),
which derives the constant from the standard exponential asymptotic.

Since `tN=A/t=sqrt(AN)`, (15)--(16) give

```math
\Pr(H=N)\ge c_3/N
\tag{17}
```

for all sufficiently large `N`. This lower bound is sufficient;
no local limit theorem is needed. Divide the exceptional probability
from Section 3 by (17), and reduce the positive exponential constant.
We obtain

```math
b_2(N)=\bigl(1-O(e^{-c_4N^{1/4}})\bigr)p(N)
\quad\text{as even }N\longrightarrow\infty.
\tag{18}
```

Together, (11) and (18) prove (1).

## 5. Consequences for the homogeneous crossing

At `h=2^r`, the homogeneous Gram operator `B_h` is an involution
modulo two. The exact central formula from
[schur-fusion-gram.md](schur-fusion-gram.md) gives, with
`f_h=rank mu_h` and `Q=p(2h)`,

```math
\max(0,2f_h-Q-1)\le r_h=c_h\le f_h+1\le Q+1.
\tag{19}
```

Equation (1) now proves (2). Applying (16) at `h` and `2h`
proves (3). The possible radical of the fusion image in characteristic
two is fully charged by `2f_h-Q`; no positivity assumption is used.

There is also a rational consequence along the same weights. The
integral Gram determinant is odd, so the homogeneous rational duals
have coefficients in `Z_(2)`, the rationals with odd denominators.
Their central matrices reduce to those in (19). Rank over `Q`
is at least rank modulo two. Conversely, the existing source channel
counts bound either central rank by `poly(h)p(2h)`. Hence the full
rational homogeneous family has the same leading central width
coefficient along powers of two. This conclusion still concerns the
homogeneous dual choice; it does not identify all alternative gauges.

## 6. Verification scope

Run:

    python3 tools/check_schur_fusion_rank.py \
      --output evidence/schur-fusion-rank.json

The finite checker uses set-partition inclusion-exclusion for (7),
independent finite-variable products for (9), and the source recursion
for (8). It counts carry-free splittable partitions through weight `32`.
The existing independently checked weight-`8` Schur calculation has
fusion rank `229` in dimension `231` and central rank `228` in both cuts.
This last record is reused explicitly, with its source hash; the new
checker does not present it as a fresh large-matrix calculation.

The new checks include `67` unique leading expansions, `134`
independent finite-variable evaluations, `1,616` monomial product
coefficients, `144` carry-free graded products, and `671` integral
fusion entries. A negative control distinguishes ordinary subset
splitting from carry-free splitting: `(1,1)` splits ordinarily,
but its multiplicity has only one binary atom of weight two.

The asymptotic proof comes from the integral filtration, interval
construction, and conditioning estimates. Finite ranks are supporting
checks, not an extrapolation to (1).
