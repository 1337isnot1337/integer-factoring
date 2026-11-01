# Singular endpoint loops: a construction and an all-field cost floor

The three-place regularity assumption in
[global rational rigidity](joint-evaluation-algebra.md) has a real boundary.
Allowing a gauge to be singular at the crossing point `t=1` produces
genuine new strict polynomial crossings with the original complete
zero-place Smith profile and leading infinity coefficient.

For the five-letter reference, the nilpotent sector of the first
endpoint-jet commutant reduces every such rational gauge to three
effective scalar parameters. Every **nontrivial** tensor in this family
has typed Smith parameter `rho>=3/4` at the inherited heights, over
every characteristic-two field. Since `3/4>1/sqrt(2)`, this family
cannot improve the half-rate criterion. Here `rho` is a cost parameter,
not the factoring saving exponent; the complete score remains `1/2`.

## 1. A nontrivial nilpotent endpoint direction

Use the exact reference `A(t)=T(t)S`, with `A(1)=I`, and set `s=t+1`.
Write its Hasse expansion

```math
A(1+s)=I+\sum_{j=1}^{18}s^j D_j.
```

The energy-block centralizer of `D_1` has dimension five, one greater
than the full stabilizer of the polynomial family. It is commutative.
Its nilpotent sector `R_1` has dimension three and satisfies `R_1^2=0`.
The intersection `Z` with the full commutant has dimension two.
The two complementary full-stabilizer idempotents span the quotient
by `R_1`, so the same nilpotent-sector description holds after every
coefficient-field extension.

A concrete generator for `R_1/Z` is the following binary matrix `N`:
on each joint-state pair

```text
(4,20), (7,11), (8,16)
```

put a two-square all-ones block, and put zero elsewhere. Indices are
`x+5y`. The two states in each pair have equal integer energy. Each
all-ones block squares to zero in characteristic two, so `N^2=0`.

Exact multiplication verifies, for all polynomial coefficients `A_j`,

```math
R_iA_jR_i=0,\qquad R_iA_jR_k+R_kA_jR_i=0\quad(i<k),     \tag{1}
```

where `R_i` are a binary basis of `R_1`. There are 114 identities.
They imply `P(t)A(t)P(t)=0` for arbitrary rational coefficients
`P(t) in R_1 tensor K(t)`. Also

```math
[N,D_0]=\cdots=[N,D_3]=0,\qquad [N,D_4]\ne0.            \tag{2}
```

Thus the matrix polynomial `[N,A(t)]` vanishes to exactly order four
at one. In particular, for arbitrary `a,b,c` in any characteristic-two
field,

```math
h(t)=\frac a{t+1}+\frac b{(t+1)^2}+\frac c{(t+1)^3},
\qquad U(t)=I+h(t)N                                  \tag{3}
```

has inverse itself and gives

```math
A_U=A+h[N,A].                                        \tag{4}
```

Equation (2) makes (4) polynomial and preserves `A_U(1)=I`.
The gauge is regular and invertible at zero, so the whole zero-place
weighted Smith profile is unchanged. It equals `I` at infinity, so
the leading weighted infinity coefficient remains the original one.
For a nonzero parameter triple, (4) is different from the reference.

## 2. Every rational gauge in this sector reduces to (3)

This is also a completeness statement within the specified sector.
Let `U=I+P(t)`, with rational coefficients in `R_1`, be regular at zero,
equal to the identity at infinity, and give a strict polynomial seed.
By (1) its action is linear. Splitting

```math
P(t)=Z(t)+h(t)N,\qquad Z(t)\in\mathcal Z\otimes K(t)
```

shows that the tensor depends only on `h`. At every finite point
`d!=0,1`, the single-evaluation algebra theorem implies
`[N,A(d)]!=0`: otherwise this energy-preserving `N` would commute with
the whole coefficient algebra, contrary to its choice outside `Z`.
This remains true after passing to an algebraic closure.

Consequently polynomiality of (4) forbids a pole of `h` at every such
point. Regularity at zero removes that point, and `U(infinity)=I`
gives `h(infinity)=0`. Thus `h` is a principal part supported only at
one. By (2), a pole of order four would alter the crossing endpoint,
and any larger pole would destroy polynomiality. Its order is at most
three, proving (3). The two rational `Z` directions preserve the tensor.

This reduction holds for every field in scope and every rational
degree; it is not an inference from a finite pole-order search.

## 3. Two symbolic minors give a uniform partition lower bound

At infinity put `C(z)=z^8 A(z^-1)D_H(z)` and realign into the opposite
cut. The regular matrix of the transformed seed is

```math
C_U(z)=C(z)+g(z)[N,C(z)],\qquad g(z)=h(z^{-1}).         \tag{5}
```

For a nonzero triple in (3), `ell=ord_z(g)` is one, two or three.
Build the realigned matrix `C+u[N,C]` with an independent formal
parameter `u`. Two literal principal minors are enough:

- Rows and columns `(0,3)` give exactly `z(1+z)`, independent of `u`.
- Rows and columns `(0,2,3,7,8,11,12,13,17)` give
  `d_0(z)+u^2d_2(z)+u^4d_4(z)+u^6d_6(z)`, with respective valuations
  `52,34,18,4` and binary nonzero leading coefficients.

After substituting `u=g(z)`, the last term has uniquely smallest
valuation `4+6ell<=22`. Hence this minor cannot suffer a parameter
cancellation over any field in scope. The checker computes both full
determinants by subset dynamic programming, retaining every polynomial
coefficient in the evidence.

Let `eta_1<=...<=eta_25` be the regular infinity Smith orders of (5).
Its unchanged rank-one leading coefficient gives `eta_1=0`. The two
minors give the determinantal prefix bounds

```math
\eta_1+\eta_2\le1,\qquad \sum_{i=1}^9\eta_i\le22.      \tag{6}
```

Since the orders are increasing, (6) makes the first nine weakly
prefix-dominated by `(0,1,3,3,3,3,3,3,3)`. For example, for `2<=k<=9`
the average of the first `k-2` terms among `eta_3,...,eta_9` is no
greater than their total average, so

```math
\sum_{i=1}^k\eta_i
\le\eta_2+\frac{k-2}{7}(22-\eta_2)
\le1+3(k-2).
```

The standard hinge-sum argument for decreasing convex exponentials
therefore gives, for every `tau>0`,

```math
\sum_i e^{-\tau\eta_i}\ge1+e^{-\tau}+7e^{-3\tau}.       \tag{7}
```

Only nine orders are needed; no all-field classification of the
remaining sixteen orders is asserted.

## 4. An exact polynomial proves rho at least three quarters

Write the unchanged zero orders as `nu_i` and define

```math
P_0(y)=\sum_i y^{\nu_i-8}.
```

The checker recomputes the complete reference zero profile; its
smallest order is eight. The unregularized infinity orders are
`eta_i-8`, so their shift cancels the zero-place shift. By (7), the
product of partition sums is at least

```math
(1+e^{-\tau}+7e^{-3\tau})P_0(e^{-\tau}).                \tag{8}
```

The inherited physical energies are `E=(0,10,6,2,8)` on each side.
At typed exponent `r=2rho=3/2`, the pressure is

```math
\log(Z_\infty Z_0)
       -3\log\!\left(\sum_i e^{-2\tau E_i/3}\right).
```

Set `x=e^(-tau/3)`, so `0<x<1`. A sufficient strict positivity
inequality is the entirely integral polynomial inequality

```math
(1+x^3+7x^9)P_0(x^3)
       >(1+x^4+x^{12}+x^{16}+x^{20})^3.                \tag{9}
```

The difference is exactly `x^3 Q(x)`, where `Q` has degree 66.
If `Q(x)=sum a_i x^i`, its scaled Bernstein coefficients are

```math
b_k=\sum_{i=0}^k a_i\binom{66-i}{k-i}\qquad(0\le k\le66).
```

All 67 integers are positive; their minimum is two. Indeed
`Q(x)=sum b_k x^k(1-x)^(66-k)`, so `Q` is positive on `[0,1]`.
The checker regenerates the degree-69 difference, divides off `x^3`,
and computes every `b_k` in integers. Thus (9) holds at every scale.

Pressure decreases with `r`: the derivative of
`r log sum exp(-tau E_i/r)` is the corresponding nonnegative entropy.
Strict positive pressure at `r=3/2` therefore forces the required
typed parameter to exceed `3/4` at every finite scale. Its infimum is
at least `3/4`, excluding the target `1/sqrt(2)` uniformly.

## 5. Replay and precise limits

The [checker](../tools/check_joint_endpoint_loops.py) recomputes the
five-dimensional endpoint commutant, its three-dimensional nilpotent
sector, all 114 identities, (2), both complete symbolic minors and the
positive Bernstein certificate. It also independently checks every
strict endpoint and both complete weighted Smith profiles for all
64 triples over `F4`: 130 fraction-free profiles including the two
reference computations. There are four infinity profile classes;
every profile is weakly prefix-dominated by the reference at the
natural heights. This proves a separate all-scale comparison over
`F4`. The all-field floor uses (6)--(9), not that finite catalogue.

```sh
make verify-joint-endpoint-loops
python3 tools/check_joint_endpoint_loops.py --generate
```

The [evidence](../evidence/joint-endpoint-loops.json) stores the
generator, symbolic determinants, all 67 positive Bernstein integers,
the complete finite profile catalogue and source hashes. All proof
checks use only Python's standard library. Arbitrary gauges leaving
this nilpotent sector, other physical energies and stronger
coefficient recombinations remain outside the theorem. This work
used substantial OpenAI assistant support and retains the Gaussian
reference attribution. No new Lean, independent review, general
factoring lower bound or improved saving exponent is claimed.
