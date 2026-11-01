# The regular homogeneous central piece has the half-rate exponent

## Result and scope

Let `q(n)` count partitions of `n` into distinct odd parts. Work over
`F_2` in the regular Fock quotient

```math
W=\mathbf F_2[p_1,p_3,p_5,\ldots]/(p_1^2,p_3^2,p_5^2,\ldots).
```

Its weight-`h` basis consists of the `q(h)` regular monomials.
For the source mode convention, define

```math
\mu_h:W_h\otimes W_h\longrightarrow W_{2h},
\qquad \mu_h(a,b)=a[h]b.
```

There is an absolute constant `c>0` such that

```math
\operatorname{rank}\mu_h
\ge\bigl(1-O(e^{-c h^{1/4}})\bigr)q(2h)
\qquad(h\longrightarrow\infty).
\tag{1}
```

At weights `h=2^r`, use the homogeneous unit duals on the full
`q(h)`-label domain. The central crossing coefficient then has
the same rank in both cuts, and that rank is

```math
\bigl(1-o(1)\bigr)q(2h).
\tag{2}
```

Since `log q(2h)/log q(h)` tends to `sqrt(2)`, the actual central
rank forces `rho>=1/sqrt(2)` for this full-domain family.
Thus its half-rate bound cannot improve merely by replacing the
channel dimensions with the exact ranks of its existing pieces.

This theorem concerns the full regular homogeneous family at powers
of two. It does not cover arbitrary binary prefixes, further linear
label restrictions, or recombinations of different coefficient pieces.
It is not a lower bound on factoring algorithms.

The proof has two parts. A length filtration reduces the fusion rank
to an elementary subset-sum property of regular partitions. Almost
every large regular partition has that property. The latter statement
is proved below using independent random occupations and conditioning.

## 1. Splittable partitions give independent fusion vectors

Call a regular partition `nu` of `2h` splittable if some subset of
its parts has sum `h`. Its complementary subset also has sum `h`.
Write the two distinct-odd partitions as `A,B`; their supports are
disjoint. Normal ordering gives

```math
\mu_h(p_A,p_B)=p_{A\cup B}
 +\text{terms with fewer than }\ell(A)+\ell(B)\text{ parts}.
\tag{3}
```

Indeed, if no annihilation occurs, every creation index is at least
its parameter index. Their total weights are equal to `h`, so every
index must remain unchanged. This gives the first term with coefficient
one. If there are `k>=1` annihilations, the output has exactly
`ell(A)+ell(B)-2k` parts before zero terms are removed in the quotient.

Choose one split for each splittable `nu`. Order the corresponding
output monomials by their number of parts. The selected fusion columns
have a triangular minor with diagonal one, by (3). Therefore

```math
\operatorname{rank}\mu_h\ge b(2h),
\tag{4}
```

where `b(N)` counts splittable distinct-odd partitions of the even
integer `N`. Ordinary multiplication need not have its image contained
in the fusion image. The triangular-minor argument is the needed
rank statement and makes no such containment claim.

## 2. A deterministic subset-sum interval lemma

Suppose a set of positive odd integers contains disjoint pairs
`(a_i,a_i+2)`, for `i=1,...,l`, and another odd element `b`.
Assume `b<=2l-1`. Put `B=sum_i a_i`.
Choose exactly one element from each pair. The resulting sums are

```math
B,B+2,\ldots,B+2l.
```

Choosing whether to include `b` gives the opposite parity as well.
Hence every integer in

```math
[B+b,B+2l]
\tag{5}
```

is a subset sum. Its endpoints add to the total sum of the pair
elements and `b`.

More generally, if the subset sums of a finite set contain every
integer in `[L,U]`, adjoining a new element `a<=U-L+1` gives the
full interval `[L,U+a]`. This follows by taking the union of
`[L,U]` and its translate `[L+a,U+a]`.
The lower endpoint stays fixed. If the initial interval is symmetric
about half the set's total sum, that property also stays true.

These observations let a short set of adjacent pairs generate an
interval, after which other selected parts enlarge it.

## 3. Independent occupations contain such an interval

For a large integer `N`, put

```math
A=\pi^2/24,\qquad t=\sqrt{A/N},\qquad K=\lfloor N^{1/4}\rfloor.
```

Select each positive odd integer `j` independently, with probability

```math
\Pr(j\text{ selected})=\frac1{1+e^{tj}}.
\tag{6}
```

Only finitely many are selected with probability one, since the sum
of the probabilities is finite. For every odd `j<=2N^(3/8)`, the
probability in (6) is at least `1/3` once `N` is sufficiently large.

The following events hold together with failure probability
`O(exp(-c_0 K))`, for an absolute `c_0>0`.

1. Among disjoint adjacent-odd pairs in `[K/2,K)`, at least `K/200`
   pairs have both elements selected.
2. At least one odd `b<=K/1000` is selected.
3. At least `K/20000` odd integers in `[K/500,K/250]` are selected.
4. In every interval `[L,2L)`, where `L=2^j K` and
   `K<=L<=N^(3/8)`, at least `L/20` odd integers are selected.

Integer roundings do not affect these fixed positive proportions.
For (1), there are `K/8+O(1)` disjoint pairs, each selected with
probability at least `1/9`. The pair indicators are independent.
For (2), use the product of the failure probabilities.
For (3) and (4), the independent Bernoulli counts have means at least
a fixed positive multiple of the interval length. Each stated count
is below half that mean for all sufficiently large `N`.
The elementary multiplicative Chernoff bound
`Pr(X<EX/2)<=exp(-EX/8)` proves the claims. Summing the bounds
in (4) over the doubling intervals gives `O(exp(-c_0 K))`, since
their lengths increase geometrically. Independence between these
four events is not needed.

Use the pairs in (1) and the element from (2) in (5).
The initial interval has width at least
`2K/200-K/1000=9K/1000`, up to bounded rounding.
Every part in (3) is at most `K/250`, so it can be added by the
interval lemma. Their total contribution is at least

```math
\frac{K}{20000}\frac{K}{500}=\frac{K^2}{10^7}.
\tag{7}
```

For sufficiently large `N`, the new width exceeds `4K`.
Add every remaining selected part less than `K`. Now process the
doubling intervals in (4) in order. At the start of `[L,2L]`, the
current width is at least `2L`; every selected part in that interval
can therefore be added. Its selected parts increase the width by
at least `L^2/20`, which is at least `4L` for `L>=80`.
This proves the induction needed for the next interval.

There is a final such `L` with
`N^(3/8)/2<L<=N^(3/8)`. After it is processed, the width is at
least `N^(3/4)/80`. The half-open intervals have no repeated
boundary elements; the weaker bound `N^(3/4)/100` suffices below.

With additional failure probability `O(exp(-c_1 N^(1/4)))`, no
selected part is larger than `N^(3/4)/1000`. Indeed a union bound
and (6) give

```math
\Pr(\exists j>N^{3/4}/1000\text{ selected})
\le\sum_{j>N^{3/4}/1000}e^{-tj}
=O\!\left(\sqrt N\,e^{-c_2N^{1/4}}\right).
\tag{8}
```

All remaining selected parts can now be added to the interval.
The lower endpoint remains the `B+b` from (5), which is `O(K^2)`.
If the final total sum is `N`, the upper endpoint is `N-(B+b)`.
For large `N`, `B+b=O(sqrt(N))<N/2`. Thus `N/2` is a subset
sum whenever `N` is even and the total is `N`.

This proves a bound on the unconditioned failure probability.
The next step accounts for the cost of conditioning on the exact weight.

## 4. Conditioning produces the uniform partition distribution

Let `H` be the total selected weight and put

```math
Z(t)=\prod_{j\text{ odd}}(1+e^{-tj}).
```

Each distinct-odd partition `nu` has probability
`e^(-t|nu|)/Z(t)`. Consequently, conditioned on `H=N`, the
distribution is uniform among the `q(N)` regular partitions.

We use the classical self-conjugate partition asymptotic

```math
q(N)\sim
\frac{\exp(\pi\sqrt{N/6})}
     {2^{7/4}3^{1/4}N^{3/4}}.
\tag{9}
```

Distinct-odd partitions correspond to self-conjugate partitions through
their diagonal hook lengths. The asymptotic (9) is proved in Section
5.1 of [Mark Wildon, *Counting partitions on the abacus*](https://www.ma.rhul.ac.uk/~uvah099/Maths/abacus.pdf).
Only this stated coefficient estimate is imported here.

For the decreasing function `f(x)=log(1+e^(-x))`, comparison with
its integral on intervals of length `2t` gives

```math
\log Z(t)=\frac1{2t}\int_0^\infty f(x)\,dx+O(1)
          =\frac{A}{t}+O(1).
\tag{10}
```

The integral equals `pi^2/12`, by the convergent exponential series.
Using `t=sqrt(A/N)`, equations (9)--(10) yield

```math
\Pr(H=N)=q(N)e^{-tN}/Z(t)=\Theta(N^{-3/4}).
\tag{11}
```

Divide the unconditioned exceptional probability from Section 3
by (11). A polynomial factor can be absorbed by decreasing the
positive constant in the exponential. Therefore

```math
b(N)=\bigl(1-O(e^{-cN^{1/4}})\bigr)q(N)
\quad\text{as even }N\longrightarrow\infty.
\tag{12}
```

Equations (4) and (12) prove (1).

## 5. From fusion rank to the central crossing rank

The parameter field and its adjoint identities descend to the regular
quotient. For completeness, an even parameter has zero creation and
annihilation coefficients on odd states, by the lowest binary digit
in the source binomial coefficients. A squared odd parameter also
acts as zero: normal ordering and Frobenius leave only creation
squares or ordinary derivative squares, both zero in characteristic
two. Thus the parameter ideal, as well as the state ideal, vanishes.

The divided translations preserve this ideal, as proved in Section 6
of [modular-crossing.md](modular-crossing.md). The integral square
identity in [prime-power-gram-involution.md](prime-power-gram-involution.md)
therefore gives `B_h^2=I` on `W_h` at `h=2^r`.
The Hall form is diagonal with entries `z_mu`; on distinct odd
monomials every `z_mu` is one in `F_2`. It is nondegenerate.

The central Gram calculation of
[schur-fusion-gram.md](schur-fusion-gram.md) now applies with this
orthonormal basis and `q` in place of the ordinary partition count.
Writing `Q=q(2h)` and `f=rank mu_h`, it gives

```math
\max(0,2f-Q-1)\le r_h=c_h\le f+1\le Q+1.
\tag{13}
```

The lower bound uses nondegeneracy, not positivity: the radical
of a subspace of dimension `f` in a nondegenerate `Q`-dimensional
space has dimension at most `Q-f`. The central coefficient differs
from its fusion Gram matrix by a rank-one vacuum term and invertible
changes of the two label pairs.

Substitute (1) into (13) to obtain (2). Finally, (9) gives

```math
\frac{\log(r_hc_h)}{4\log q(h)}\longrightarrow\frac1{\sqrt2}.
```

This is an actual rank obstruction for the stated full-domain central
piece. It does not replace the separate analysis needed after a further
restriction of the external labels.

## 6. Reproducible checks

Run:

    python3 tools/check_regular_fusion_rank.py --weight-thirty-two \
      --output evidence/regular-fusion-rank.json

The default verification target omits only the largest fusion matrix.
The checker tests the deterministic interval lemma, exact
partition counts and subset-sum decisions, the length filtration,
regular Gram involutions, and the finite central rank bounds.
It verifies `78` initial intervals, `2,253` interval extensions, and
`2,592` parameter-ideal entries. The largest fusion check independently
compares all `279,841` entries of the central Gram formula with the
source operator construction.

| `h` | Labels `q(h)` | Ambient `q(2h)` | Fusion rank | Central rank, both cuts |
| ---: | ---: | ---: | ---: | ---: |
| 4 | 1 | 2 | 1 | 1 |
| 8 | 2 | 5 | 3 | 2 |
| 16 | 5 | 23 | 16 | 12 |
| 32 | 23 | 276 | 268 | 261 |

The probability argument and its conditioning step are ordinary proofs;
finite partition counts are not used to extrapolate (12).
