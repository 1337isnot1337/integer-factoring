# A universal floor for the multipoint local-jet Smith criterion

## Result and scope

The [multipoint Hermite construction](multipoint-hermite-crossings.md)
uses both local Smith profiles to bound each retained jet's two ranks.
That sufficient pressure criterion **cannot improve the factoring
exponent one half**, for any coefficient field, physical alphabet,
polynomial degree, number of interpolation points, regular local frame,
degree allocation, or independent cut scales permitted by that criterion.

More precisely, let `q>=2` be the physical alphabet size and put `n=q^2`.
The criterion requires

```math
\rho>R_n:=\frac{\log n+H_n}{2\log n}>\frac34,
\qquad
H_n=\log(2n)-\frac{n+1}{2n}\log(n+1).                \tag{1}
```

Since `3/4 > 1/sqrt(2)`, it cannot supply a post-half crossing.
The corresponding saving exponent is strictly below
`log_2(4/3)=0.415037...`. The established factoring score remains **1/2**.
For example, the explicit floor `R_4=0.887198735...` for binary physical
letters is stronger than the earlier finite `3/4` certificates.

This is an obstruction to the displayed **sufficient Smith-profile
bound**, not to factoring, multipoint interpolation itself, or actual
ranks of tensor-power jets. It does not cover estimates that exploit
additional coefficient cancellation, auxiliary word modules, or physical
type projections with different complete costs. It uses no Lean work.
The recent [universal binary two-place theorem](https://github.com/1337isnot1337/integer-factoring/blob/0ac85be63d6eff434eddd92a32eb5c7e4ee75cb2/research/universal-binary-smith-floor.md)
has a different scope: it allows arbitrary physical heights, while this
argument covers the unweighted multipoint criterion over every alphabet.

## 1. Determinants give a shared global budget

Let `K(z)` be a polynomial tensor of degree at most `D`, with
`K(1;x,y,u,v)=[x=v][y=u]`. The two cut matrices have order `n`.
They are generically invertible, since each endpoint matrix is a
permutation matrix. Choose distinct interpolation points `a_j != 1`,
possibly including infinity. Use `K(a_j+t)` at a finite point and
`t^D K(1/t)` at infinity.

Write their nonnegative local Smith orders as

```math
\nu_{j,1}\le\cdots\le\nu_{j,n},\qquad
\mu_{j,1}\le\cdots\le\mu_{j,n}.
```

The determinant of either cut is a nonzero polynomial of degree at most
`nD`. Its order at a finite point is the sum of the corresponding Smith
orders. At infinity the order is `nD-deg(det K)`. Distinct finite
root factors divide the determinant with their multiplicities. Therefore,
with or without infinity among the selected points,

```math
\sum_j\sum_i\nu_{j,i}\le nD,\qquad
\sum_j\sum_i\mu_{j,i}\le nD.                         \tag{2}
```

Crucially, the minimum orders at any one place agree:

```math
\nu_{j,1}=\mu_{j,1}=:m_j.                            \tag{3}
```

Both minima are the minimum valuation among all tensor entries; changing
the flattening only permutes those entries. Independent regular local
row and column operations preserve these Smith orders. They cannot
invalidate either (2) or (3).

Define the local mean budget

```math
B_j=\frac{\sum_i\nu_{j,i}+\sum_i\mu_{j,i}}{2n}.
```

Equation (2) says `sum_j B_j <= D`.

## 2. One local profile pair must pay its entropy

Take any local pair `nu,mu` with common minimum `m`, and any slope
`d >= B=(sum nu+sum mu)/(2n)`. Subtract `m` from both profiles and
write `e=d-m>=0`. Let `A,C` be the sums of these shifted profiles,
ordered so that `A<=C`. Then

```math
A\le ne,\qquad C\le2ne.                            \tag{4}
```

For the smaller-sum profile, the uniform distribution on its `n`
channels has mean order at most `e` and entropy `log n`.
For the other profile use the distribution

```math
p=\frac12\delta_{\rm min}+\frac12\operatorname{Uniform}_n
 =\left(\frac{n+1}{2n},\frac1{2n},\ldots,\frac1{2n}\right).
                                                               \tag{5}
```

Its mean shifted order is `C/(2n)<=e`, and its entropy is exactly
`H_n` from (1).

For any probability vector `p_i` whose mean order is at most `d`,
weighted arithmetic-geometric mean gives, for every `s>=1`,

```math
\left(\sum_i s^{-\nu_i}\right)s^d
\ge\prod_i\left(\frac{s^{-\nu_i}}{p_i}\right)^{p_i}s^d
\ge e^{H(p)}.                                      \tag{6}
```

Apply (6) separately to the uniform distribution and (5), after restoring
the common minimum. Thus every pair of independent scales `s,r>=1`
satisfies

```math
\boxed{\log\sum_i s^{-\nu_i}+\log\sum_i r^{-\mu_i}
              +d\log(sr)\ge\log n+H_n.}            \tag{7}
```

The common-minimum condition is essential. If one profile is all ones
and the other all zeros, their mean budget is `1/2`, but letting the
first scale tend to infinity makes the left side tend to minus infinity.
Such a pair cannot be the two flattenings at the same tensor place.

The constant in (7) is sharp under just these abstract profile constraints.
Take one profile `(0,2n,...,2n)` and the other all zero, with slope
`d=n-1`. Choose the dominant cut's scale with `s^(2n)=n+1` and the
other scale one. Equality holds. This is an abstract profile example,
not a polynomial crossing attaining the global bound.

## 3. Every degree allocation hits a contributing place

For `D>0`, the criterion allocates nonnegative rational slopes with
`sum_j d_j=D`. Let `J` be the nonempty set of positive slopes. If every
contributing place had `d_j<B_j`, then

```math
D=\sum_{j\in J}d_j<\sum_{j\in J}B_j
 \le\sum_j B_j\le D,
```

a contradiction. Some positive-slope place therefore satisfies (7).
Since the criterion requires every contributing pressure to be strictly
less than `4 rho log q = 2 rho log n`, equation (1) follows. Omitting a
place does not evade this argument. If `D=0`, the tensor is the constant
crossing and every local pressure is `2 log n`, requiring `rho>1` under
the strict version of the sufficient criterion; it supplies no improvement.

Finally, strict concavity of entropy applied to (5) gives

```math
H_n>\tfrac12 H(\delta_{\rm min})
       +\tfrac12 H(\operatorname{Uniform}_n)
     =\tfrac12\log n,\qquad n>1.
```

Hence `R_n>3/4`. The explicit formula shows `R_n -> 3/4` as `n` grows.
The floor `3/4` is therefore uniform across all alphabet sizes even
though its displayed strict margin depends on the alphabet.

## 4. Exact controls and the completed cubic catalogue

The proof above is uniform ordinary mathematics. Its checker uses exact
rational pressure comparisons after raising to `2n`:

```math
\left[U(s)V(r)(sr)^B\right]^{2n}
\ge\frac{(2n^2)^{2n}}{(n+1)^{n+1}}.                 \tag{8}
```

It checks arbitrary finite profile pairs, sharp abstract corners,
literal original polynomial tensors in several characteristics and
alphabets, their complete local Toeplitz profiles, determinant budgets,
and nonuniform allocations including omitted points. An unequal-minimum
negative control fails (8), as it should.

The constructive cubic search finished before the uniform argument was
identified. It covers every binary-alphabet cubic over `F_3` with pure
values at `0,2,infinity` and strict crossing at one. Normalize the first
pure value using the [six stabilizer representatives](multipoint-rank-one-seeds.md),
and allow all 513 pure/zero tensors at each other point. This gives
**1,579,014 representatives** and **200** profile signatures. Its
independent replay retains every complete witness, both local cuts,
the infinity degree shift, and full two-copy Hermite reconstruction.
The native catalogue also has exact prefix certificates at `rho=3/4`;
the uniform theorem now gives the stronger floor (1).

```sh
g++ -O3 -std=c++17 -Wall -Wextra \
  tools/search_multipoint_pure_cubics.cpp -o /tmp/multipoint-pure-cubics
/tmp/multipoint-pure-cubics /tmp/multipoint-pure-cubics-replay.json
cmp evidence/multipoint-pure-cubics.json /tmp/multipoint-pure-cubics-replay.json
python3 tools/check_multipoint_smith_floor.py \
  --output /tmp/multipoint-smith-floor-replay.json
```

Further searches using only this sufficient criterion can be retired.
A constructive improvement needs a sharper actual-jet rank argument,
a physical type construction with complete costs, or a different crossing
interface. Finite profile gains within (7) cannot change the score.

This proof and its computations used substantial OpenAI Codex assistance.
