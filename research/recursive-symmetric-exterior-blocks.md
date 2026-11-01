# Polynomial-size blocks for symmetric exterior multiplication

## Result

The earlier [symplectic block reduction](cube-symmetric-exterior-blocks.md)
can be repeated. For `r>=1` square-zero generators, multiplication by
an arbitrary symmetric element has a simultaneous constant block
form with exactly `ceil(r/2)` nonzero-multiplicity block types. Every
block has order at most

```math
 \boxed{(r+2)^2/4.}
 \tag{1}
```

The original matrix has order `2^r`. The reduced matrices and their
integer multiplicities have polynomial description size in `r`.
This works over every commutative coefficient ring of characteristic
two, for all symmetric coefficients at once. It preserves full ranks,
Smith profiles over a discrete valuation ring, and every repeated
block; there is no discarded tail.

For example, 63 generators give an operator of order `2^63`, but
only 32 block types, the largest of order 512. At 127 generators,
there are 64 types and the largest has order 2,048. These statements
refer to a structured matrix representation, not to a new factoring
runtime. The complete factoring score remains `1/2`.

Files: [block constructor](../tools/recursive_symmetric_exterior_blocks.py),
[exact checker](../tools/check_recursive_symmetric_exterior_blocks.py),
and [receipt](../evidence/recursive-symmetric-exterior-blocks.json).

## 1. Repeat the simultaneous change of generators

Let

```math
 A_r=R[\varepsilon_1,\ldots,\varepsilon_r]/(\varepsilon_i^2),
 \quad \operatorname{char}R=2,
 \qquad f=\sum_{j=0}^r a_j e_j(\varepsilon),
 \tag{2}
```

where `e_j` is the elementary symmetric form. The preceding block
proof is a simultaneous constant change of basis for every `e_j`.
It therefore remains valid when the coefficients lie in a larger
commutative ring carrying previously exposed square-zero variables.

For odd `r=2s+1`, a block with `k` active symplectic pairs has
multiplicity `binom(s,k)2^(s-k)` and multiplier

```math
 \sum_{j=0}^k(a_{2j}+u a_{2j+1})e_j(z_1,\ldots,z_k),
 \quad u^2=z_i^2=0.
 \tag{3}
```

For even `r=2s`, its multiplicity is `binom(s-1,k)2^(s-1-k)` and
its multiplier is

```math
 \sum_{j=0}^k(a_{2j}+v a_{2j+1}+vw a_{2j+2})e_j(z),
 \quad v^2=w^2=z_i^2=0.
 \tag{4}
```

The remaining `z` dependence in both cases is still symmetric. Apply
the same theorem to those `k` variables, with coefficient ring enlarged
by `u` or by `v,w`. Repeating reaches zero symmetric variables after
at most `floor(log2(r+1))` nonempty stages. All changes are invertible
constant binary matrices on the represented module.

## 2. A block type is a word of ones and twos

Write `O` for a use of (3) and `E` for a use of (4). Equivalently,
write a type as `delta=(delta_0,...,delta_(ell-1))` with each digit
one or two. It exposes

```math
 b(\delta)=\sum_j\delta_j
```

square-zero variables, so its block order is `2^b`.
At level `j`, an `O` variable has weight `2^j`. An `E` pair has
these local possibilities:

| Local monomial | Coefficient-index contribution |
|---|---:|
| `1` | `0` |
| `v` | `2^j` |
| `w` | forbidden |
| `vw` | `2^(j+1)` |

For any monomial in the exposed variables, add these contributions
across levels. If a forbidden local pattern occurs, its coefficient
is zero. Otherwise its coefficient is exactly `a_I`, where `I` is
the summed index. This gives the full multiplier of the final block.
There are no additional signs or unspecified change-of-basis entries.

The largest possible index is

```math
 R(\delta)=\sum_{j=0}^{\ell-1}\delta_j2^j.
 \tag{5}
```

The [constructor](../tools/recursive_symmetric_exterior_blocks.py)
returns this explicit coefficient array together with its multiplicity.
Different types may give identical matrices for a special coefficient
choice. No irreducibility or minimality assertion is needed.

## 3. Exact reachability and the number of types

A nonempty word is reachable from `r` precisely when

```math
 R(\delta)\le r,\qquad R(\delta)\equiv r\pmod2.
 \tag{6}
```

Indeed, the first digit is forced by the parity of `r`. After it is
removed, the remaining generator count can be any integer from zero
through `(r-delta_0)/2`. The least count that can realize the tail is
its own value of (5). This proves (6) by induction. Conversely, choose
that least count at each step; every required binomial multiplicity
is positive.

Every positive integer has a unique finite expansion in the digits
one and two used in (5). Its first digit is one if the integer is odd
and two if it is even; subtract that digit and divide by two to recover
the tail. Thus the possible values in (6) are exactly
`1,3,...,r` when `r` is odd, or `2,4,...,r` when it is even.
There are exactly `ceil(r/2)` types.

## 4. The quadratic block-size bound is sharp

Let a word have length `ell`, with `k` of its digits equal to two.
Then `b=ell+k`. Its extra contributions above all-one digits occur
at `k` distinct powers of two, whose sum is at least `2^k-1`.
Consequently

```math
 R(\delta)+2
 =2^\ell+1+\sum_{j:\delta_j=2}2^j
 \ge2^\ell+2^k
 \ge2\,2^{(\ell+k)/2}.
 \tag{7}
```

Squaring gives `2^b <= (R+2)^2/4`, and (6) proves (1).
Equality is attained by an all-`E` word: then
`R=2^(ell+1)-2` and the block order is `2^(2ell)`.
This is a sharp worst-case bound for the block family.

## 5. Explicit multiplicities and polynomial construction cost

Let `M_r(w)` be the multiplicity of a word. Start with
`M_0(empty)=1` and every other empty-word multiplicity zero. The
complete recurrences are

```math
 M_{2s+1}(Ow)=\sum_{k=0}^s
        \binom sk2^{s-k}M_k(w),
 \tag{8}
```

```math
 M_{2s}(Ew)=\sum_{k=0}^{s-1}
        \binom{s-1}k2^{s-1-k}M_k(w).
 \tag{9}
```

Every word with the wrong first symbol has multiplicity zero.
Equal words arising through different active-pair choices are merged
by addition of their integer multiplicities. Inductively,

```math
 \sum_w M_r(w)2^{b(w)}=2^r.
 \tag{10}
```

These multiplicities count direct summands. They are not reduced
modulo the characteristic. Each uses at most `r+1` binary digits.

Memoizing (8)--(9) for generator counts through `r` uses polynomially
many integer operations. There are only `O(r)` final matrices, each
with `O(r^2)` rows and `O(r^2)` coefficients in its multiplier.
Constructing them from (5) has polynomial cost. For coefficients
given to series precision `p`, ordinary cubic elimination on all
blocks has the coarse bound `O(r^7)` series operations, or
`poly(r,p)` bit operations over a fixed finite field. This replaces
the earlier exponential matrix-size factor for this symmetric problem.

The explicit global change of basis would still have exponentially
many entries. It is unnecessary for profile computation: the proof
establishes conjugacy, while (8)--(10) retain the full dimensions and
multiplicities in binary.

## 6. Larger multiaxis screens

This reduction makes higher-axis symmetric specialization screens
practical. For an integer budget `D_0=128`, take

```math
 a_j(x)=x^{\lceil D_0/(j+1)\rceil},\qquad
 D=\max_{1\le k\le d}k\lceil D_0/k\rceil.
 \tag{11}
```

These are the one-variable crossblocks obtained by setting every
other axis variable to one in the Cayley construction

```math
 c_S=\prod_{a\in S}x_a^{\lceil D_0/|S|\rceil}
 \quad(S\ne\varnothing),\qquad
 c_\varnothing=\sum_{S\ne\varnothing}c_S.
 \tag{12}
```

The endpoint remains the antipodal matching. The kernel need not
be homogeneous, so the charged degree is the actual cap `D` in
(11), rather than `D_0`. This avoids a least-common-multiple degree
explosion while retaining a literal polynomial construction.

For its specialized multiplier, the square is `x^(2D_0)` times the
identity. Its determinant is `x^(D_0 2^(d-1))`. Therefore every Smith
order is at most `2D_0`, their sum is exactly `D_0 2^(d-1)`, and the
list pairs around `D_0`. A 512-coefficient series window resolves
every promised pivot. The exact receipt includes the complete
compressed specialized profiles at 16, 32, and 64 axes.

A specialization can raise determinantal valuations. These profiles
are therefore optimistic screens for the generic counting pressure,
as proved in the preceding note. They must not be presented as an
upper bound on generic multivariable contraction. The initial higher-axis
screens did not produce a useful pressure candidate. Further coefficient
searches can use the new block engine without claiming that one
specialized success would complete a factoring proof.

## 7. Independent replay

The checker compares the recursive profiles against 40 profiles from
the unrepeated symplectic construction and against 24 literal dense
multiplication matrices. It checks every type count, full dimension
sum, and quadratic size bound through 256 generators. Those finite
checks supplement the all-size proof above.

The optional large run reconstructs the actual blocks for (11),
resolves every exact pivot, and checks the total and paired-order
identities at 16, 32, and 64 axes.

```sh
python3 tools/check_recursive_symmetric_exterior_blocks.py --large \
  --output /tmp/recursive-symmetric-exterior-blocks.json
cmp /tmp/recursive-symmetric-exterior-blocks.json \
  evidence/recursive-symmetric-exterior-blocks.json
```

Without `--large`, the checker runs the independent algebraic controls
without the native C++ engine. The new result is an ordinary block
conjugacy theorem and a larger exact search capability. No new Lean
formalization or improved factoring exponent is claimed.
