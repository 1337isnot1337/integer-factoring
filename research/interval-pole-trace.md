# A global pole trace with closed bilinear reduction rules

This continues the global elimination question in
[Discussion #53](https://github.com/1337isnot1337/integer-factoring/discussions/53).
It gives an exact pole filter, a family closed under both interval boundary
translations and dyadic reduction, and an explicit size bound for the
enumerative reference. It does **not** give a faster factoring algorithm.
The complete score remains `1/2`; there is no new Lean work.

## 1. Move the pole into a separate parameter

For positive integer intervals `I=[L,U)`, `J=[V,W)`, put

```math
C_{I,J}(Z)=\sum_{a\in I,b\in J}\frac1{Z-ab}.
```

The residue at `Z=N` is the exact number `D` of factor pairs in the
rectangle. The weighted version has the stronger identity

```math
\mathop{\rm Res}_{Z=N}
 \sum_{a\in I,b\in J}\frac1{(Z-ab)(X-a)}
 =\sum_{ab=N}\frac1{X-a}=G'(X)/G(X),                 \tag{1}
```

where `G` is the [interval orbit gcd](interval-orbit-gcd.md). Thus this
recovers the same global output as the rational-antiderivative projection,
with its sign reversed, by a residue in a new parameter. Unlike generic
Hermite reduction over `Q(N)`, the residue is taken at the specified input
`N` from the beginning. It does not discard the specialized common roots.

Writing `C=R'/R` for `R(Z)=prod_(a,b)(Z-ab)` is only an identity.
Neither a short circuit for `R` nor an inexpensive evaluator for its
normalized logarithmic derivative is supplied here.

## 2. Only polynomial precision is needed for a normalized trace

Work in the unramified quadratic 2-adic extension with
`w^2+w+1=0`. Let `p>=1`, choose `h` so that

```math
|N-ab|<2^h\quad(a\in I,b\in J),\qquad k=p+h,
```

and set `t=2^k`. Empty rectangles are immediate. For a nonempty rectangle,
one can use the bit length of the maximum absolute difference at the two
product corners for `h`. For endpoints and `N` of `O(n)` bits, `h=O(n)`.
Then

```math
t w\,C_{I,J}(N+t w)\equiv D\pmod {2^p}.             \tag{2}
```

Here is a termwise proof, including normalization. Write `delta=N-ab`.
If `delta=0`, the scaled term is exactly one. Otherwise, set
`v=v_2(delta)<h`. The denominator
`delta+t*w` has valuation `v`, and the scaled fraction has valuation
`k-v>p`. More explicitly, its exact coordinates in `Q(w)` are

```math
\frac{t w}{\delta+t w}
 =\frac{t^2+t\delta w}{\delta^2-\delta t+t^2}.        \tag{3}
```

The norm denominator has valuation `2v`; after exact cancellation both
coordinates have odd denominators and vanish modulo `2^p`.
Consequently every summand is integral, so the congruence survives the
whole sum without an error depending on the number of terms.

The same argument works with the unit weight `1/(w+T-a)` in any truncated
jet ring. It gives `G'(w+T)/G(w+T)` modulo `(2^p,T^(d+1))`.
The denominator is a unit because the residue of `w` lies outside `F_2`.
This is a **global projection identity**, not another evaluator for one
input product at a specified center.

Since `0<=D<=min(|I|,|J|)`, choosing
`p=max(1,bit_length(min(|I|,|J|)))` recovers the integer count exactly.
Both `p` and `k` are `O(n)`. An effective normalized-trace algorithm at
the desired subexponential cost would therefore give factoring by the
existing polynomial number of rectangle and bisection queries. A
polynomial-time trace algorithm would imply polynomial-time factoring;
no such algorithm is asserted.

Reducing the huge product `R` modulo `2^p`, differentiating, and dividing
does **not** implement (2). That product can have enormous content, and
its reduced derivative and denominator need not be invertible. Exact
normalization before reduction is part of the missing evaluation cost.

## 3. Close both boundary translations by enlarging the state

For `q(u,v)=A*u*v+B*u+C*v+D`, define

```math
\mathcal C_q(Z;I,J)=\sum_{u\in I,v\in J}\frac1{Z-q(u,v)}.
```

Replacing `u` by `u+1` sends the coefficient tuple to

```math
(A,B,C,D)\longmapsto(A,B,C+A,D+B),                 \tag{4}
```

and replacing `v` by `v+1` sends it to
`(A,B+A,C,D+C)`. These are exact identities in the whole parameter `Z`.
For example, if `q_+` is (4), cancellation of overlapping interval terms
gives the closed boundary rule

```math
\mathcal C_{q_+}(Z;[L,U),J)-\mathcal C_q(Z;[L,U),J)
=\mathcal C_q(Z;\{U\},J)-\mathcal C_q(Z;\{L\},J).   \tag{5}
```

There is an identical rule for the second interval. These are the
two interval-product boundary cancellations, expressed after passing to
logarithmic derivatives and a bilinear denominator. In fact, with
`L_I=P_I'/P_I`, the original trace satisfies
`C_(I,J)(Z)=sum_(a in I) L_J(Z/a)/a=sum_(b in J) L_I(Z/b)/b`.
For a singleton `u`,
the remaining denominator is affine in `v`; its logarithmic-derivative
boundary shift is the ordinary interval-product rule. The enlargement
keeps both translations inside one four-coefficient family. It does not
bound the number of distinct states reached by iterating the rules.

## 4. An exact projection-preserving dyadic rule

Assume `N` odd; even inputs already have smallest prime factor two.
Only odd `a,b` contribute to the pole, so write
`a=1+2u`, `b=1+2v`. Dividing `ab-N` by two gives the initial state

```math
q(u,v)=2uv+u+v+(1-N)/2.                           \tag{6}
```

The retained family has `A` positive and even, `B,C` positive and odd.
For each `e in {0,1}`, exactly one parity
`f=(B*e+D) mod 2` can contain a zero. Substitute
`u=e+2u'`, `v=f+2v'`, and divide by two. The child coefficients are

```math
A'=2A,\quad B'=B+A f,\quad C'=C+A e,\quad
D'=(Aef+Be+Cf+D)/2.                               \tag{7}
```

The interval endpoints transform by ceiling division. Direct substitution
proves `q(e+2u',f+2v')=2q'(u',v')`. The other two parity classes have no
zero. Thus the residue counts add over these two children. At the
resolvent level, a retained child contributes
`(1/2)*C_child(Z/2)`; its residue at zero is unchanged. The scaled filter
has the same identity with `k` replaced by `k-1`.

For the weighted projection (1), retain the affine map from the current
`u` to the original `a`. Under a section it changes by affine composition;
its coefficients have `O(n)` bits. The supplied fast reference implements
the scalar count, while the checker independently tests weighted traces.

Because the original factors are positive, `q` is strictly increasing in
each coordinate. If its two product corners have the same strict sign,
the state contains no zero. If the minimum or maximum corner is zero,
it is the unique zero. A singleton interval also reduces to one exact
linear divisibility test.

## 5. Charge the representation, including its width

Let `H` be the smaller number of odd integers in the two original
intervals. The binary recursion stops after at most `ceil(log_2 H)`
levels for `H>=1`, since each child interval has at most the ceiling of
half the parent length. It visits fewer than `4H` nodes. Empty rectangles
need one node. The stored coefficient tuple, affine map, and endpoints
have `O(n)` bits: after fixing `s` low bits they describe

```math
ab-N=2^s\left(2^suv+b_0u+a_0v+(a_0b_0-N)/2^s\right).
```

Depth-first traversal therefore costs `O((H+1)*poly(n))` bit operations
and polynomial working space for a scalar count. This is an enumerative
baseline, equivalent in purpose to the existing inverse-section search.
A four-coefficient state is small; there can still be exponentially
many states. No better bound follows from the closure identities alone.

The positive target is now specific: reduce or merge this global trace
family while preserving (2), (5), and (7), and charge the total number of
states and every normalized division. To improve the current score, the
complete query cost must beat `2^(n/(log n)^(1/2-o(1)))`. Another cheap
point evaluation, or a compact name for the state family, does not do so.

## 6. Checks, including artificial common factors

Run `make verify-interval-pole-trace`. The
[checker](../tools/check_interval_pole_trace.py) uses exact fractions in
`Q(w)` for (3), dense rational gcds for (1), literal pair sets for (7),
and direct rational sums for both boundary rules. Its
[receipt](../evidence/interval-pole-trace.json) records the actual counts.
These are finite controls supporting the ordinary arguments above.

For the [artificial modular-factor family](interval-packet-circuits.md),
take `A_0=2^r+1`, `B_0=2^(r+1)+1`, `N=A_0*B_0`,
`I=[1,A_0)`, `J=[1,B_0)`. Every product is strictly below `N`, so the
exact trace residue is zero. The raw modular input products nevertheless
share `(X-A_0)` to order at least `2^(r-1)-p+1` when this is positive.
The exact enclosure rule rejects the entire rectangle in one node.
Including the upper endpoints instead gives exactly one pair at the
maximum corner, also in one node.

The receipt checks these two global cases at `r=12,80,200`, and checks
exact filter coordinates at selected large pairs. It does not enumerate
the enormous rectangles. At `r=200,p=32`, the input has 402 bits and the
filter displacement has 434 bits despite the enormous raw fake factor.
This removes that false-factor issue from the definition of the target
operator; an efficient general evaluation of the operator remains open.
