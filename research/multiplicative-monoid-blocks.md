# Explicit multiplication blocks and Fourier--Mellin closure

## Result and scope

The residue multiplication gate has a constructive decomposition that
retains the large characteristic-two nilpotent part of the unit group.
For every odd modulus `M`, its multiplication algebra splits into
explicit truncated polynomial blocks over a fully charged coefficient
algebra. For the adaptive modulus `M=2^m+1`, the resulting number of
bilinear products is

```math
R\le M\,4^{\omega(M)}=M^{1+O(1/\log D)}.                 \tag{1}
```

Here `D|m` is the power-of-two digit parameter, and `omega(M)` counts
distinct prime factors. The product count is **over the explicit
coefficient algebra**, not an uncharged rank bound over the original
binary field.

The same coordinates also give exact formulas for additive Fourier
transforms, addition, multiplication, unit scaling, and joint copying.
They work at every prime power dividing `M`, hence over `M` by CRT.
The Fourier transform reverses one valuation index and applies an
invertible polynomial on the characteristic-two unit index. It does
not require a semisimple transform of an even-order group.

This supplies a closed arithmetic tensor family. It does **not** yet
give a faster evaluator for that family. The bounded input moments,
correlations created by copying, and sums over free character labels
remain charged operations. No factoring exponent above `1/2` follows
from this note.

This construction follows the
[adaptive transform](adaptive-transform-constructor.md),
[operand-length refinement](operand-length-transform-constructor.md), and
[weighted Fourier contraction](weighted-dyadic-carry-contraction.md).
The coefficient extension is compatible with the
[sparse recovery theorem](frobenius-compatible-sparse-recovery.md) and
[constant-field weights](constant-field-carry-weights.md).

## 1. An explicit splitting algebra

Write

```math
M=\prod_i p_i^{e_i},\qquad
p_i-1=d_i h_i,\qquad d_i=2^{s_i},\qquad h_i\text{ odd},
\qquad H_i=h_i p_i^{e_i-1}.
```

Let `L=lcm_i H_i`; thus `L` is odd and
`L<=prod_i phi(p_i^e_i)=phi(M)<=M`. In the adaptive application the
existing field is `F=F_(2^(2m))`, with its explicit primitive `M`th
root. Form

```math
\mathcal B=F[Z]/(\Phi_L(Z)\bmod2),\qquad \theta=Z\bmod\Phi_L.
                                                               \tag{2}
```

For `L=1` this is just `F`, with `theta=1`. This algebra can be a
product of fields. None of the following proofs assumes it is a field.
It is free of rank `phi(L)<=M` over `F`.

Every `theta^j-1`, for `0<j<L`, is a unit. Indeed, `X^L-1` is
squarefree in characteristic two. Its reductions of the cyclotomic
factors are therefore pairwise coprime; every root of `Phi_L mod2`
has exact order `L`. Consequently `Phi_L` and `X^j-1` are coprime.
The extended Euclidean algorithm supplies the inverse. In particular,
the odd character transforms below have unit Vandermonde differences.

All ingredients are deterministic and have bit cost polynomial in
`M` and `m`. Trial division factors `M` and the relevant group orders.
Monic polynomial division constructs the cyclotomic polynomial from
`X^L-1`. Primitive unit generators and discrete logarithm tables can
be found by exhaustive enumeration. There is no uncharged large
finite-field compositum search.

For completeness, the unit group of an odd prime power is cyclic:
a generator modulo `p` exists because a polynomial of degree `r`
has at most `r` roots in `F_p`; one of its lifts `g` or `g+p` has
`g^(p-1)!=1 mod p^2`, and the binomial expansion then gives order
`(p-1)p^(e-1)` modulo `p^e`. The exhaustive search just described
therefore terminates.

If sparse recovery uses `E=F_(2^ell)`, use
`B tensor_(F_2) E`. Its bit dimension is `2m phi(L) ell`.
The inclusion of the previous coefficient algebra `F tensor E`
is injective because (2) is free. Its nonzero coefficients remain
nonzero. All established scalar and character unit identities remain
valid after this extension.

## 2. The local blocks, including the modular radical

Fix `q=p^e`, put `p-1=dh` as above, and set `H=h p^(e-1)`.
Choose a primitive unit `g mod q`, and write

```math
u=g^H,\qquad v=g^d.
```

Then `u` has order `d`, `v` has odd order `H`, and each unit has a
unique expression `u^a v^b`. Let `chi` range over the `H` characters
of the odd factor. If the `p`-primary part of `chi` has exact order
`p^r`, where `0<=r<=e-1`, define

```math
\ell_\chi=e-r.
```

This parameter is an order exponent. It is not being called the
usual `p`-adic conductor exponent, which is `r+1` when `r>0`.
There are `h phi(p^r)` characters with this value of `r`.

Let `B[Z/q]` mean the algebra of the **multiplicative monoid**:
its basis elements `[a]` satisfy `[a][b]=[ab mod q]`.
There is an explicit algebra isomorphism

```math
\mathcal B[\mathbb Z/q]
\simeq \mathcal B\ \times
\prod_\chi \mathcal B[X,Y]/(X^d,Y^{\ell_\chi}).       \tag{3}
```

The first factor is the augmentation `sum_a f(a)`. To give the other
factors on a basis element, write `a=p^t a_0!=0`, with `a_0` a unit,
and let `eta(a_0)` be its exponent in the factor generated by `u`.
The `chi` component is

```math
[a]\longmapsto
\begin{cases}
\chi(a_0)(1+X)^{\eta(a_0)}Y^t,&t<\ell_\chi,\\
0,&t\ge\ell_\chi.
\end{cases}                                                   \tag{4}
```

Every non-augmentation component of `[0]` is zero.

### Proof of the isomorphism

The idempotent `[0]` splits off the first factor. In the complementary
algebra the basis is `[a]-[0]` for `a!=0`, and products that reach
the zero residue become zero.

The odd character idempotents split this complementary algebra.
There is no division by an even integer: `H` is odd and is the
scalar one in characteristic two. On valuation layer `t`, the
stabilizer of a unit is the kernel of
`U(p^e)->U(p^(e-t))`, of odd order `p^t`. A character survives
precisely when it is trivial on that stabilizer, namely
`t<e-r=ell_chi`.

The `2`-primary unit factor acts freely on every surviving layer.
Its algebra is

```math
\mathcal B[U]/(U^d-1)
=\mathcal B[X]/X^d,\qquad U=1+X.                   \tag{5}
```

Multiplication by the residue `p` moves to the next valuation
layer. It gives `Y`, with exactly `ell_chi` surviving layers.
This proves (3), and also proves that (4) does not depend on the
chosen lift of `a_0`.

One can verify bijectivity without an arbitrary matrix inversion over
a ring with zero divisors. On layer `t`, apply the odd character
transform of order `H/p^t`, followed by the binary binomial matrix
`binom(a,k) mod2`. Both are explicitly invertible. For the latter,
Lucas's identity gives `binom(a,k)=1` precisely when every set bit
of `k` is a set bit of `a`; the Boolean zeta matrix squares to the
identity in characteristic two. Recover the coefficient of `[0]`
from the augmentation.

The dimensions agree exactly:

```math
1+dh\sum_{r=0}^{e-1}\varphi(p^r)(e-r)=p^e.            \tag{6}
```

For example, the nonzero blocks for `q=9` have shapes `2 by2`,
`2 by1`, and `2 by1`; the zero block adds one dimension. Replacing
these blocks by scalar odd characters would lose their nilpotent
coordinates.

CRT identifies the multiplicative monoid modulo `M` with the
product of these local monoids. Tensor (3) over all prime powers.
The resulting blocks are truncated polynomial rings with at most
`2 omega(M)` variables, and their dimensions sum to `M`.

## 3. Constructive bilinear product counts

For a block

```math
\mathcal B[Z_1,\ldots,Z_a]/(Z_1^{c_1},\ldots,Z_a^{c_a}),
```

evaluate each input polynomial on a rectangular grid with
`2c_i-1` distinct scalar nodes in coordinate `i`. Multiply the two
values at each grid point, interpolate the untruncated product,
and discard coefficients outside the box. This uses exactly

```math
\prod_{i=1}^a(2c_i-1)                                  \tag{7}
```

bilinear products in `B`. All evaluation and interpolation maps are
explicit linear maps. In the adaptive application `m>=2` and
`|F|=2^(2m)>2M`, so the required nodes can be taken in `F` itself.
Their nonzero differences remain units in `B tensor E`.

The exact product count supplied by this construction for `p^e` is

```math
\begin{aligned}
R_{p^e}
 &=1+h(2d-1)\sum_{r=0}^{e-1}\varphi(p^r)(2(e-r)-1)\\
 &=1+h(2d-1)
 \left(2\frac{p^e-1}{p-1}-p^{e-1}\right).
\end{aligned}                                                   \tag{8}
```

For a general odd modulus the count is
`R=prod_(p^e||M) R_(p^e)`. In particular,

```math
R\le M4^{\omega(M)},\qquad
R_{p^e}\le2(1+1/p)p^e.                                \tag{9}
```

The looser first inequality also follows directly from (7) and
the `2 omega(M)` variable bound. Some exact examples are:

| Modulus | Dimension | Bilinear products over `B` |
| --- | ---: | ---: |
| 3 | 3 | 4 |
| 7 | 7 | 10 |
| 9 | 9 | 16 |
| 15 | 15 | 32 |
| 25 | 25 | 50 |
| 27 | 27 | 52 |
| 65 | 65 | 176 |
| 257 | 257 | 512 |

For the source family `M=2^m+1`, every prime divisor `p` satisfies

```math
v_2(\operatorname{ord}_p2)=v_2(m)+1,\qquad
2^{v_2(m)+1}\mid p-1.                                \tag{10}
```

To see this, `2^m=-1 mod p` gives order dividing `2m` but not `m`.
Its odd part divides the odd part of `m`, so its `2`-adic valuation
is exactly the one stated. Since the power of two `D` divides `m`,
every `d_i>=2D` and every `p_i>=2D+1`. Hence

```math
(2D+1)^{\omega(M)}\le M,
\qquad
R\le M^{\,1+2/\log_2(2D+1)}.                         \tag{11}
```

The sharper local bound in (9) can replace the numerator `2` by
`log2(2(1+1/(2D+1)))`. This is a near-linear product count for one
gate, with all modular radicals present. The identity element also
gives a flattening of dimension `M`, so a bilinear decomposition
over any field component needs at least `M` terms.

This is not a near-linear bit-operation claim. A `B tensor E`
element has `O(mM ell)` bits, and its operations must be charged.
The explicit arrays, scalar arithmetic, and their construction all
have cost `2^O(m) poly(ell)` in the adaptive family. Factoring the
small modulus by trial division is included in that cost; it is
not an oracle for factoring the original `n`-bit input.

## 4. Coordinates for complete arithmetic gates

The useful linear coordinates are the actual zero coefficient
`z=f(0)` and the block moments

```math
t_{\chi,v}(X)
=\sum_{\substack{a\bmod p^e\\v_p(a)=v}}
 f(a)\chi(a/p^v)U^{\eta(a/p^v)},
\qquad U=1+X,\quad 0\le v<\ell_\chi.                 \tag{12}
```

These differ from the augmentation coordinates of (3) by the
explicit identity

```math
\operatorname{aug}(f)=z+\sum_{v=0}^{e-1}t_{1,v}(0).    \tag{13}
```

There are exactly `p^e` coordinates. Write `t_(chi,v,k)` for the
coefficient of `X^k` in (12). On each valuation layer the inverse
of (12) is the odd inverse character transform followed by

```math
f(u^a v^b)
=\sum_\chi\chi(v^b)^{-1}
       \sum_{k:\ a\subseteq k}t_{\chi,k},             \tag{14}
```

with the appropriate layer understood. Here `a subseteq k` means
containment of their binary sets of bits. The odd inverse-transform
normalization is one in characteristic two.

### Multiplication, pointwise products, and copying

Multiplicative convolution becomes ordinary block multiplication:

```math
T_\chi(f*_\times g)=T_\chi(f)T_\chi(g)
       \pmod{X^d,Y^{\ell_\chi}},\qquad
T_\chi(f)=\sum_v t_{\chi,v}(X)Y^v.                   \tag{15}
```

The augmentation multiplies, and (13) recovers the zero coefficient.
Thus nonunit products that become zero are included exactly.

Pointwise multiplication of functions has a different exact rule:

```math
t_{\chi,v,k}(fg)
=\sum_{\substack{\psi\phi=\chi\\
                  v<\ell_\psi,\ell_\phi}}
  \ \sum_{i\mathbin{\&}j=k}
 t_{\psi,v,i}(f)t_{\phi,v,j}(g),\qquad z(fg)=z(f)z(g).
                                                               \tag{16}
```

The symbol `&` is bitwise AND. To prove (16), use (14) and sum
over the unit exponents. The coefficient in question is the parity
of the interval of bit sets
`k subseteq a subseteq (i AND j)`. It is one exactly when
`k=i AND j`. Odd character orthogonality supplies `psi phi=chi`.
Only characters surviving the common valuation layer are present.

Copying is a joint, two-output tensor; replacing it by independent
marginals would lose the arithmetic correlations. Its exact moment
formula is particularly short:

```math
\Delta(t)_{(\chi,v,i),(\psi,w,j)}
=\begin{cases}
t_{\chi\psi,v,i\mathbin{|}j},&v=w,\\
0,&v\ne w.
\end{cases}                                                   \tag{17}
```

The symbol `|` is bitwise OR. The zero-zero coordinate is `z`,
and mixed zero/nonzero coordinates vanish. Equation (17) follows
coefficientwise from
`binom(a,i)binom(a,j)=binom(a,i OR j) mod2`. It retains the full
joint copy tensor, including all character and nilpotent indices.

Unit scaling `a -> ca` multiplies the `chi` block by
`chi(c) U^eta(c)` and leaves the valuation index unchanged.
These formulas already cover residue multiplication and the
copy/scale operations in the arithmetic circuit.

## 5. The additive Fourier transform at every prime power

Let `zeta` be a primitive `p^e`th root in the coefficient algebra,
with unit root differences, and put

```math
(\mathcal F f)(s)=\sum_{a\bmod p^e}\zeta^{as}f(a),
\qquad
\iota(X)=\frac{X}{1+X}\pmod{X^d}.                    \tag{18}
```

For a character whose `p`-primary part has order `p^r`, define
the Gauss polynomial

```math
G_\chi(X)=
\sum_{a\in(\mathbb Z/p^{r+1})^\times}
  \chi(a)U^{\eta(a)}\zeta_{p^{r+1}}^a,
\qquad
\zeta_{p^{r+1}}=\zeta^{p^{e-r-1}}.                  \tag{19}
```

The character and the even unit exponent descend to this smaller
unit group, so the expression is well defined.

For every **nontrivial** odd character, additive Fourier has the
exact block formula

```math
t_{\chi,w}(\mathcal F f)
=G_\chi(X)\,
 \iota\!\left(t_{\chi^{-1},\ell_\chi-1-w}(f)
       \right),\qquad 0\le w<\ell_\chi.             \tag{20}
```

The `Y` coefficient index is reversed. This is a linear reversal,
not inversion of a nilpotent variable. For the trivial character,
whose block length is `e`, the formula is

```math
\begin{aligned}
t_{1,w}(\mathcal F f)
={}&G_1(X)\,
 \iota\!\left(t_{1,e-1-w}(f)\right)\\
&+X^{d-1}\left(z(f)+
       \sum_{v=e-w}^{e-1}t_{1,v}(f)(0)\right),\\
z(\mathcal F f)
={}&z(f)+\sum_{v=0}^{e-1}t_{1,v}(f)(0).
\end{aligned}                                                   \tag{21}
```

Empty sums in this formula are zero.

### Exact Gauss support and proof

For an input valuation `v` and an output valuation `w`, the
relevant sum runs over units with additive character of conductor
`p^(e-v-w)`, unless `v+w>=e`, when that character is trivial.

The weight `chi(a) U^eta(a)` is constant on residue classes modulo
`p^(r+1)`. If the additive conductor is larger, summing over its
last `p` lifts gives zero. If the additive conductor is smaller
and `chi` is nontrivial on the corresponding kernel, multiply
the summation variable by a kernel unit: the sum changes by a
nontrivial odd root of unity and hence is zero. The difference of
that root from one is a unit by Section 1.

Therefore a nontrivial character has support exactly on the
antidiagonal `v+w=e-r-1`. Reduction to `p^(r+1)` then gives
(19); its omitted fibers have odd cardinality and hence scalar
one. The substitution by the input unit contributes
`chi(a_0)^(-1) U^(-eta(a_0))`, producing the character inversion
and `iota` in (20).

For the trivial character, an additive conductor at least `p^2`
still gives zero, conductor `p` gives `G_1`, and the trivial
additive character gives

```math
\sum_{a\in U(p^a)}U^{\eta(a)}
=\sum_{j=0}^{d-1}U^j=X^{d-1}.                       \tag{22}
```

This accounts for the terms above the antidiagonal and for the
zero residue. It proves (21).

Every `G_chi` is a unit, including the one for the trivial
character. More precisely,

```math
G_\chi(X)\,\iota(G_{\chi^{-1}}(X))
=U^{d/2}+[\chi=1]X^{d-1}.                           \tag{23}
```

For nontrivial characters this follows by applying (20) twice:
`F^2 f(a)=f(-a)` because the odd integer `p^e` is the scalar one.
Multiplication by `-1=u^(d/2)` gives the first term. For the trivial
character, (19) is the prime-modulus sum. Restrict that prime DFT
to nonzero rows and columns; its square is reflection plus the
all-ones matrix. The latter acts as multiplication by `X^(d-1)`
on the trivial moment block, which proves the extra term.

The right side of (23) has constant coefficient one and squares
to one. It gives an explicit inverse for `G_chi`. Thus (20)--(21)
are invertible operations over `B`, including when `B` has zero
divisors.

Addition convolution is now completely explicit:

```math
f*_+g=\mathcal F^{-1}
       \big((\mathcal F f)(\mathcal F g)\big).        \tag{24}
```

Use (20)--(21), the pointwise rule (16), and inverse Fourier.
The inverse is Fourier followed by unit scaling by `-1`.

For a composite `M`, use the CRT idempotents `epsilon_i` and
the local roots `zeta_M^(epsilon_i)`. Then the additive DFT is
the tensor product of the local DFTs. Equations (15)--(24),
together with the CRT basis, therefore cover every source modulus,
including repeated prime factors.

## 6. What closes, and what still costs time

The local arithmetic formulas close on a finite explicit family:
odd characters, valuation layers, characteristic-two polynomial
indices, Gauss units, coefficient reversal, and AND/OR constraints.
There are exactly `M` basis coordinates. An arbitrary vector of
those `M` coordinates can be transformed or used in any of the
binary operations by `poly(M)` coefficient operations. The
bilinear count (8) is available for the multiplication operation.
Every scalar operation has the charged bit cost from Section 3.

This closure also describes the character conservation rules. At
a prime modulus, addition has `chi_out=chi_1 chi_2`, multiplication
has `chi_1=chi_2=chi_out`, and copying has
`chi_in=chi_out,1 chi_out,2`. These are linear relations on the
odd character indices. They do not determine all internal indices
in a general network. The remaining sums carry Gauss and nilpotent
coefficients; enumerating their solutions has not been made free.

The input moments are also explicit, but can be dense. A bounded
digit interval with scalar weight `lambda^a` contributes the sum
of (4) over its permitted digits. Direct construction of all `M`
coordinates costs `O(bM)` coefficient operations for a `b`-digit
box. Here `b<M`, so this is `2^O(m)`, including the CRT correlations
between different prime powers. A CRT projection of an integer
interval is not being replaced by independent residue intervals.
Prefix restrictions and a shorter final digit box can be treated
by the same finite sum.

More generally, a joint boundary with `a` residue variables has
`M^a` coordinates after this change of basis. The formulas do not
turn it into `a` independent vectors. In particular, a Fourier
butterfly has two coupled outputs, and (17) retains those copy
correlations. A complete table construction must still charge
every written boundary coefficient and every sum needed to
produce it.

The coefficient extension and local preparation have additive
logarithmic cost `O(m)` in the bit-operation exponent. For the
adaptive parameters `m=O(n/P^18)` and `P^40<=n`, this fits the
previous `2^O(n/P)` construction budget. It does not reduce the
number or dimensions of the joint boundary tables. Nor does a
near-linear decomposition of a single multiplication tensor by
itself prove a smaller recursion rate for those tables.

Finally, the nilpotent coordinates cannot be dropped even on
units. The nonzero vector `[1]+[u]` has zero scalar moment for
every odd character, but its full moment in each block is `X`.
All source prime factors have `d>=2D`, so these omitted blocks
would be large. The exact controls below include this failure of
the scalar-character shortcut.

The next concrete question is whether the coupled character-flow
and polynomial-index family admits a recursive evaluator with
smaller **joint** interface cost. Equations (15)--(24) provide the
local tensors and their construction costs for that question.
No such all-level cost improvement is asserted here.

## 7. Exact evidence and sources

Run:

```sh
python3 tools/check_multiplicative_monoid_blocks.py \
  --output evidence/multiplicative-monoid-blocks.json
```

The checker constructs explicit binary fields and unit generators.
It checks all basis products for the small monoids, full basis
ranks, dense multiplication convolutions, literal multivariate
interpolation, and the exact product-count formulas. It also
checks complete Fourier transforms and coupled addition and
multiplication at primes and prime powers, the Gauss-unit
identities, and the full joint copy coefficients. These are
finite controls for the ordinary proofs above, not a proof of
a stronger factoring bound.

The multiplication-network source remains
`IsaacRudich/SubexponentialDeterministicFactoring`
at `c91cc3f4947caca8977e81a5d14244f73b5ae2b1`.
General context for multiplicative semigroup algebras of finite
rings is Benjamin Steinberg,
[*Factoring the Dedekind--Frobenius determinant of a semigroup*](https://arxiv.org/abs/2108.02166),
arXiv:2108.02166, especially its finite-ring and commutative-monoid
sections. That paper explicitly distinguishes the modular group
algebra issue from the complex-character setting. The concrete
characteristic-two blocks, effective costs, and coupled formulas
used here are proved in this note; no semisimple even-order
character decomposition is imported from a characteristic-zero
argument.
