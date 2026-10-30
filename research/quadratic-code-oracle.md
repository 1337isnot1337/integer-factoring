# Polynomial-time coefficient queries on quadratic Boolean codes

## Result and scope

The characteristic-two 23-term crossing profile can be restricted to
quadratic Boolean evaluation codes without enumerating their exponentially
long words. Every coefficient entry of the resulting crossing is computable
in polynomial time in its binary label width. The coefficient set has
`2^O(sqrt(lambda))` positions, and full tables and exact factors remain
constructible in `2^O(lambda)` bit time.

This is a stronger local query interface for the quadratic-code candidate
identified in [PR #24, affine-code-barrier.md](https://github.com/1337isnot1337/integer-factoring/blob/05d9c1b8aed4f8494ce230f0869f9048ef879474/research/affine-code-barrier.md).
It supplies no uniform contraction estimate below `1/sqrt(2)`. The
factoring score therefore remains `1/2`.

A second result identifies a restriction on the support of this candidate.
If a Boolean quadratic has alternating rank at least six, its quadratic
annihilators are exactly zero and its complement. This forces the support
of every local term, on an appropriate high-rank quadratic subspace, to
depend only on the affine hull of a quartet's image. The numerical
multiplicities in that image can still differ from the affine-code case;
the existing affine-code rank obstruction cannot simply be imported.

All statements here are ordinary mathematics and exact computation, without
new Lean work. The original crossing construction retains its Appendix B
attribution. The 23-term polynomial is pinned to
`1337isnot1337/integer-factoring@05d9c1b8aed4f8494ce230f0869f9048ef879474`,
specifically [hamming-tensor-power.md](https://github.com/1337isnot1337/integer-factoring/blob/05d9c1b8aed4f8494ce230f0869f9048ef879474/research/hamming-tensor-power.md)
and [hamming_tensor_power.py](https://github.com/1337isnot1337/integer-factoring/blob/05d9c1b8aed4f8494ce230f0869f9048ef879474/tools/hamming_tensor_power.py).

## 1. Labels and the exact coefficient oracle

Work in the Boolean ring

```math
B_r=\mathbf F_2[t_1,\ldots,t_r]/(t_i^2-t_i).
```

Its degree-at-most-two subspace has basis
`1,t_1,...,t_r,t_1t_2,...,t_(r-1)t_r` and dimension

```math
q(r)=1+r+\binom r2.
```

For every positive binary width `lambda`, choose the least `r>=1` with
`q(r)>=lambda`, and use the first `lambda` basis monomials in the displayed
fixed order. Each label is one binary coefficient vector. Evaluation on
all `m=2^r` points gives a binary word. Distinct labels give distinct words:
square-free Boolean monomials are linearly independent as functions, as
follows by induction on the variables. Thus restriction preserves the exact
scale-one crossing. Also `r=O(sqrt(lambda))`.

The pinned source writes the profile of four words as

```math
P(z)=\sum_{j=1}^{23} L_j(z)A_j(0000;z)
                    \prod_{t\in\mathbf F_2^r}
                    A_j(f_x(t),f_y(t),f_u(t),f_v(t);z).       (1)
```

Every nonzero entry in this constant table is `z^a(1+z)^b` for
nonnegative constant integers `a,b`. This follows directly from the
six local edge formulas: their nonzero values are `1,z,1+z`, and the
matching factors multiply them. It is also checked independently by exact
polynomial division in the companion checker.

Let `H(p)` count points with quartet value `p in F_2^4`. A summand in (1)
vanishes if `H(p)>0` at one of its zero local entries. Otherwise it is

```math
z^{a_j}(1+z)^{b_j},\qquad
a_j=a_{j,0}+\sum_p H(p)a_{j,p},\quad
b_j=b_{j,0}+\sum_p H(p)b_{j,p}.                          (2)
```

Thus the complete profile is a sum of at most 23 binomial polynomials.
The exponents have `O(r)` bits. At degree `k`, Lucas's binary parity
identity gives the exact entry

```math
[z^k]P(z)=\sum_j [k\ge a_j]
              [(k-a_j)\mathbin{\&}\neg b_j=0]\quad\text{in }\mathbf F_2.  (3)
```

No expansion of `(1+z)^b` is needed by this oracle. Expanding such
polynomials in a finite checker is a separate operation.

## 2. Joint histograms from quadratic Gauss sums

For a Boolean quadratic `f`, set

```math
S(f)=\sum_{t\in\mathbf F_2^r}(-1)^{f(t)}.
```

Walsh inversion gives all sixteen histogram counts:

```math
H(p)=\frac1{16}\sum_{a\in\mathbf F_2^4}
 (-1)^{a\cdot p}S(a_xf_x+a_yf_y+a_uf_u+a_vf_v).         (4)
```

These are ordinary integer sums and exact divisions. They are used only
to obtain the nonnegative histogram counts, not as divisions in `F_2`.

Each Gauss sum has a polynomial-time elimination algorithm. If the current
quadratic contains the term `t_it_j`, write it as

```math
f=t_it_j+t_iL+t_jM+R,
```

where `L,M` are affine in the other variables and `R` is quadratic. Direct
summation over the two variables gives

```math
\sum_{t_i,t_j\in\mathbf F_2}(-1)^f
  =2(-1)^{R+LM}.                                      (5)
```

The product `LM` remains quadratic after reducing `t_k^2=t_k`. Remove
these two variables and repeat. When no quadratic term remains, the
remaining sum is zero if a nonconstant linear term remains, and is the
signed power of two determined by the constant term otherwise.

There are at most `r/2` eliminations. A quadratic has `O(r^2)` coefficients,
and multiplying two affine polynomials costs `O(r^2)` coefficient
operations. All final sums and counters have `O(r)` bits. Equations
(2)--(5) consequently give polynomial bit time and polynomial workspace
in `lambda` for a coefficient entry, with an explicit uniform algorithm.
The input contains four `lambda`-bit labels and an `O(r)`-bit degree.

## 3. Piece count and construction cost

The constant source table has local degree at most four and fixed
low/leading degree at most twelve. Consequently

```math
\deg P\le4m+12,
\qquad M\le4\,2^r+13=2^{O(\sqrt\lambda)}.              (6)
```

Use all positions from zero through this bound as pieces; zero pieces
may be retained as zero tensors. Their sum is the exact crossing because
`P(1)=E`. Degree enumeration is an ordinary counter with `O(r)` bits.

There are `2^(4lambda)` complete entries per coefficient tensor.
Generating all `M` tensors with (3), and obtaining both cut ranks and
exact factors by Gaussian elimination over `F_2`, costs `2^O(lambda)`
bit time. This is an effective crossing candidate in the extended piece
budget. It does not give polynomial-time factor queries: only complete
piece entries have the polynomial query bound established here.

The logarithmic piece count `O(sqrt(lambda))` is smaller than
`lambda/(log lambda)^delta` for every fixed positive `delta`. Hence piece
enumeration would preserve any fixed logarithmic exponent through the
existing extended transfer if an appropriate uniform width profile were
proved. That width profile is the remaining mathematical obligation.

## 4. A sharp rank threshold for quadratic annihilators

The alternating rank of a Boolean quadratic is the rank over `F_2` of
its symmetric zero-diagonal matrix of degree-two coefficients. Affine
terms do not change it. The rank is even.

**Theorem.** If `f in B_r` has degree at most two and alternating rank
at least six, then

```math
\{g\in B_r:\deg g\le2,\ fg=0\}=\{0,1+f\}.             (7)
```

**Proof.** Use the filtration by Boolean polynomial degree. Its associated
graded ring is the square-zero exterior algebra in characteristic two.
Write the degree-two part of `f`, after an invertible linear coordinate
change, as

```math
\omega=a_1b_1+\cdots+a_sb_s,\qquad s\ge3,
```

with the remaining variables spanning its radical. Such a coordinate
change can add lower-degree Boolean terms; this does not affect the
associated graded calculation.

First, the kernel of multiplication by `omega` from degree two to degree
four is exactly the span of `omega`. To prove this, inspect a degree-two
element `eta`. A monomial using variables from two different paired
blocks has a unique coefficient in `omega eta` after multiplying it by
a third full paired block. Therefore its coefficient in `eta` is zero.
A monomial involving radical variables is handled by a full paired block
disjoint from it. Only the terms `c_j a_jb_j` can remain. The coefficient
of `a_ib_i a_jb_j` is `c_i+c_j`, so all `c_j` agree. Conversely
`omega^2=0` in the graded ring.

Second, multiplication by `omega` is injective from degree one to degree
three. A variable's coefficient is recovered by adjoining a full paired
block disjoint from that variable; two paired blocks already suffice.

If `fg=0`, the degree-four part gives `g_2=c f_2` for `c in F_2`.
Hence `g=cf+ell` for an affine `ell`. Every element of a Boolean ring is
idempotent, so

```math
0=fg=f(c+ell).
```

Its degree-three part and the injectivity just proved force the linear
part of `ell` to vanish. Since `f` is nonzero, its constant part must be
`c`. Thus `g=c(1+f)`. Both claimed choices do annihilate `f`.

The rank threshold cannot be lowered to four in this statement. For
`f=t_1t_2+t_3t_4`, the quadratic annihilator has dimension five over
`F_2`, instead of one. This finite negative control is checked by exact
evaluation and elimination.

## 5. What this proves about local supports

Let `L` be a linear space of Boolean quadratics containing the constants,
and suppose every nonconstant element of `L` has alternating rank at
least six. Choose four labels in `L`, let `S` be their joint image in
`F_2^4`, and let `A` be the affine hull of `S`.

For any two affine functions `u,v` on `F_2^4`,

```math
uv\text{ vanishes on }S
\quad\Longleftrightarrow\quad
uv\text{ vanishes on }A.                              (8)
```

The reverse direction is immediate. For the forward direction, compose
`u,v` with the four labels. If either composition is zero, that affine
equation holds on `A`. If one composition is one, the other must be
zero. Otherwise (7) says the second composition is zero or the complement
of the first. The corresponding affine equation again holds on `A`,
and forces the product to vanish there.

Every one of the 23 local supports is the common zero set of products
of affine functions. Here is an explicit complete table; each displayed
expression is required to vanish. The bit order is `x,y,u,v`, and the
term numbers are zero-based in the pinned `local_terms()` output.

| Terms | Defining expressions |
| --- | --- |
| 0, 5 | `y+v` |
| 1, 6 | `vu`, `y(y+v)` |
| 2, 7 | `v(y+v)`, `u(u+v)` |
| 3, 8 | `vy`, `u(u+v)` |
| 4, 20 | `y+u`, `vu`, `v(x+v)` |
| 9, 21 | `y+u`, `u(u+v)`, `x(x+v)` |
| 10, 19 | `y+u`, `vu`, `vx` |
| 11, 12 | `v(u+v)`, `y(y+v)` |
| 13, 14 | `u+v` |
| 15, 18 | `y+u`, `vx`, `v(u+v)` |
| 16, 17 | `y+u`, `v(u+v)`, `v(x+v)` |
| 22 | `y+u`, `x+v` |

Linear equations are products with the constant one. Equation (8) implies
that a local term is supported on a quartet precisely when its support
contains `A`. This agrees with the support rule for affine evaluation
codes. It does not make the multiplicities `H(p)` uniform, so neither
the polynomial profile nor its ranks have been identified with those of
the affine code.

Such high-rank subspaces need not be small. Among alternating matrices
on `r` variables, at most `2^(4r)` have rank at most four: each is the
sum of at most two decomposable alternating forms, specified by four
vectors. Starting with zero, greedily choose a linear subspace avoiding
every nonzero rank-at-most-four form. If the current dimension is `s`,
at most `2^(s+4r)` next vectors are forbidden. There is another choice
whenever `s+4r<binom(r,2)`. This gives a homogeneous quadratic subspace
of dimension at least `binom(r,2)-4r`; adding constants preserves (7).
The greedy search and rank tests take `2^O(r^2)` time. The statement is
useful only when this lower dimension is positive.

## 6. Finite evidence

The checker verifies the Gauss-sum oracle against literal truth tables,
the joint histograms against direct word histograms, and the coefficient
profiles against the independent pinned source formula. It also verifies
the annihilator dimensions for canonical quadratic forms with affine
perturbations, including a rank-four negative control, and every entry of
the local-support table.

A raw rank experiment uses all 32 labels in the span
`1,t_1,t_2,t_3,t_1t_2`, evaluated at the eight points of `F_2^3`.
The maximum literal coefficient rank product is `384400`. At degree 22,
both cut ranks are 620. Its finite normalized value
`log2(384400)/(4*5)=0.9276124405...` exceeds the target threshold.
This rules out this literal finite profile as evidence of a new rate.
It is not a lower bound for recombined coefficients, other quadratic
codes, or an all-size family.

Adding the next monomial `t_1t_3` gives 64 labels on the same eight
evaluation points. The maximum product is `6635776`; at degree 22 both
cut ranks are 2576. The normalized value is `0.9442430731...`, which
also exceeds the target. These two finite profiles are recorded in
[`quadratic-ranks-3-5.json`](../evidence/quadratic-ranks-3-5.json) and
[`quadratic-ranks-3-6.json`](../evidence/quadratic-ranks-3-6.json).

Run the exact checks with
`python3 tools/check_quadratic_crossing.py` and the finite rank experiment
with `python3 tools/probe_quadratic_ranks.py --variables 3 --width 5`
(or `--width 6`). The three source helper modules are copied unchanged
from the pinned commit; their Git blob identities and SHA-256 hashes
are recorded in
[`quadratic-source-provenance.json`](../evidence/quadratic-source-provenance.json).
The [oracle certificate](../evidence/quadratic-crossing.json) records
the actual check counts. Substantial OpenAI assistant support was used.
