# A lattice terminal removes the divisor-bound overhead

## Result and scope

The [positive ratio lookahead](circuit-box-ratio-lookahead.md) gives a
complete adaptive search in `2^(n/3+O(n/log n))` bit time. Its extra
term comes from a uniform divisor bound used to control the number of
live classes. A different terminal operation removes that term.

Let `N` be odd, let `I,J` be intervals of positive integers, and put

\[
 P=\max_{a\in I,b\in J}ab.
\]

For `R=2^h` with `R^3>=P`, all exact factor pairs in a prescribed odd
ratio class

\[
 ab=N,\qquad a\equiv cb\pmod R
 \tag{1}
\]

can be recovered in polynomial bit time. The algorithm examines at
most four low residue roots and at most five lattice lines per root.
It solves an integer quadratic on each line. No divisor oracle,
factorization of an auxiliary integer, or large interval enumeration
is needed.

If the endpoints have `O(n)` bits, `n=bits(N)`, and `P<=C N` for a
fixed constant `C`, streaming all ratio classes gives

\[
 \boxed{2^{n/3}\mathrm{poly}(n)\text{ bit time and polynomial bit space}.}
 \tag{2}
\]

The same bound computes the exact rectangle factor count and its least
row. An ordered dyadic rectangle search therefore returns the smallest
prime factor, including a correct answer on prime inputs. The
implementation does not first collect all factors or make uncharged
prefix queries.

This is a fully charged improvement within the inverse-window
representation. It is a classical cube-root factoring scale. It is
weaker than the established general factoring bounds discussed in the
preceding note and weaker than this repository's `1/2` theorem. It does
not improve the score. In particular, the claim here concerns the
explicit terminal reduction and the removal of the divisor-bound
overhead, not historical novelty of a cube-root factoring algorithm.

## 1. One ratio class lies in four affine lattices

Use half-open integer intervals throughout. For a solution of (1),
write

\[
 r=a\bmod R,\qquad s=b\bmod R,
 \qquad a=r+Rx,\quad b=s+Ry.
\]

Both residues are odd. The ratio condition and `ab=N` imply

\[
 r^2\equiv cN\pmod R,\qquad s=Nr^{-1}\bmod R.
 \tag{3}
\]

There are at most four odd square roots in (3). Successive binary
lifting computes their complete list in polynomial bit time, also
covering `h=1,2`. For each root, form the exact integer intervals

\[
 U_r=\{x:r+Rx\in I\},\qquad V_s=\{y:s+Ry\in J\}.
 \tag{4}
\]

Their endpoints are obtained by ceiling division. Set

\[
 t=\frac{N-rs}{R}s^{-1}\bmod R.
 \tag{5}
\]

The quotient in (5) is an integer. Expanding the product modulo `R^2`
gives the necessary condition

\[
 x+cy\equiv t\pmod R.
 \tag{6}
\]

Thus every genuine factor pair for this root lies in

\[
 (t,0)+\Lambda,\qquad
 \Lambda=\{(x,y)\in\mathbb Z^2:x+cy\equiv0\pmod R\}.
 \tag{7}
\]

The vectors `(R,0),(-c,1)` form a basis of `Lambda`, of determinant
`R`. Equation (6) is only a necessary condition for exact equality.
The later quadratic test enforces `ab=N` exactly.

Let `U=|U_r|`, `V=|V_s|`. If either is zero, this root contributes
nothing. If `U=1`, the unique possible row `a` is known, so one integer
division tests whether `b=N/a` is integral and belongs to the column
interval and the requested ratio class. The case `V=1` is symmetric.
These cases need no lattice reduction. For the rest of the argument,
assume `U,V>=2`.

## 2. A short direction with an explicit bit-cost potential

Use the weighted inner product

\[
 \langle v,w\rangle=V^2v_xw_x+U^2v_yw_y,
 \qquad g(v)=\langle v,v\rangle.
 \tag{8}
\]

The weights account for the possible imbalance between the two
conditioned interval widths. Starting with the basis in (7), perform
this exact two-dimensional reduction:

1. Replace `w` by `w-mu v`, where `mu` is a nearest integer to
   `\langle v,w\rangle/g(v)`.
2. If `4g(w)>=3g(v)`, stop. Otherwise exchange `v,w` and repeat.

The nearest-integer step uses only integer arithmetic. Ties may be
resolved consistently in either direction. Each exchange decreases
the norm of the first vector by a factor strictly smaller than `3/4`.
The first norm starts at `(VR)^2` and is a positive integer. Therefore
there are `O(log(VR))` exchanges. Size reduction does not increase
the second norm. The initial second norm is `V^2c^2+U^2`, with
`0<c<R`. These observations also bound every basis coordinate and
intermediate numerator by a polynomial number of input bits.

At termination, the basis satisfies

\[
 2|\langle v,w\rangle|\le g(v),
 \qquad 4g(w)\ge3g(v).
 \tag{9}
\]

Negating the second vector if needed makes `det(v,w)=R`. The weighted
Gram determinant is

\[
 g(v)g(w)-\langle v,w\rangle^2=(UVR)^2.
 \tag{10}
\]

Combining (9) and (10) yields

\[
 g(v)^2\le2(UVR)^2.
 \tag{11}
\]

This proof uses a fixed reduction margin. Merely stating that a basis
is reduced would not charge the number of iterations in the reference
algorithm. The explicit exchange potential supplies that charge.

## 3. The rectangle meets few parallel lattice lines

Every point in (7) has a unique representation

\[
 p=(t,0)+jw+kv,\qquad j,k\in\mathbb Z.
 \tag{12}
\]

For a point in the conditioned rectangle, its transverse index is

\[
 j=\frac{\det(v,p-(t,0))}{R}.
 \tag{13}
\]

The numerator is a linear function of `p`, so its minimum and maximum
over the rectangle occur at the four inclusive corners. Let these
extrema be `L,H`. It suffices to enumerate

\[
 \left\lceil L/R\right\rceil\le j\le\left\lfloor H/R\right\rfloor.
 \tag{14}
\]

The difference between the extrema is

\[
 \begin{aligned}
 H-L
 &\le |v_x|(V-1)+|v_y|(U-1)\\
 &\le \sqrt{2g(v)}\\
 &\le 2^{3/4}\sqrt{UVR}
 <2\sqrt{UVR}.
 \end{aligned}
 \tag{15}
\]

Consequently, the number of lattice lines is at most

\[
 1+2\sqrt{UV/R}.
 \tag{16}
\]

This is an explicit cost bound at arbitrary ratio precision. It does
not assert a polynomial cost for a small `R` and very large intervals.

Now impose `R^3>=P`. Let `ell_I,ell_J` denote the original interval
lengths. A residue class with at least two points requires
`ell_I>R` and `ell_J>R`. Hence

\[
 U\le\ell_I/R+1\le2\ell_I/R,
 \qquad V\le2\ell_J/R.
\]

Positivity of the interval endpoints gives `ell_I ell_J<=P`.
Therefore

\[
 UV\le\frac{4\ell_I\ell_J}{R^2}
       \le\frac{4P}{R^2}\le4R.
 \tag{17}
\]

Equations (16) and (17) give at most five lines per root, a convenient
integer upper bound. Thus a ratio class uses at most twenty lines in
total. Roots handled by singleton division only reduce this work.
The line bound is uniform in the interval aspect ratio and position.

## 4. Exact equality on a line is a quadratic

For each `j` in (14), intersect

\[
 (t,0)+jw+kv\in U_r\times V_s
\]

with the two coordinate intervals. Signed ceiling and floor divisions
give the complete permitted integer interval for `k`. A zero direction
coordinate is handled by an interval-membership test. If the
intersection is empty, discard the line.

The original row and column are affine functions

\[
 a=a_0+d_a k,\qquad b=b_0+d_b k.
\]

Thus equality is equivalent to

\[
 d_a d_b k^2+(d_a b_0+d_b a_0)k+(a_0b_0-N)=0.
 \tag{18}
\]

If the quadratic coefficient is nonzero, compute the integer square
root of the discriminant. A candidate is retained only when the
discriminant is a square, the quadratic-formula numerator is divisible
by the denominator, and `k` lies in the permitted interval. If the
quadratic coefficient is zero, solve the resulting linear equation by
exact division.

The polynomial in (18) cannot vanish identically on a nonconstant
line when `N>0`. If one coordinate is constant, an identically constant
positive product would force the other coordinate to be constant too.
If both vary, the quadratic coefficient is nonzero. Thus there are
at most two integer roots per line. Verifying their original product
also provides a direct final check.

The low residue roots give disjoint row classes. The representation
(12) gives disjoint lines and distinct coordinates on each line.
Therefore every factor pair is found exactly once, apart from the
harmless deduplication of a repeated quadratic root. The reference
implementation retains at most forty pairs during a terminal ratio
query, then streams their count and least row to the outer search.

All endpoints, reduced basis coordinates, line indices, quadratic
coefficients, and discriminants have `O(n)` bits when the original
endpoints have `O(n)` bits and `h=O(n)`, as in the stated factoring
regime. The reduction, divisions, square root, and final verification
therefore take polynomial bit time per ratio class. The algorithm does
not treat integer square root or lattice reduction as a unit-cost
operation.

## 5. Complete rectangle search

First check the positive minimum and maximum rectangle products. If
`N` lies outside that range, the count is zero. Empty rectangles are
also immediate. Otherwise choose

\[
 h=\max\!\left(1,\left\lceil\frac{\mathrm{bits}(P)}3\right\rceil\right).
 \tag{19}
\]

Then `R^3>P`, and the terminal bound applies. Enumerate all odd ratios
modulo `R`. For `h>=3`, only `c=N mod 8` can contribute: an odd square
is one modulo eight, so (3) requires `cN=1 mod 8`, and an odd integer
is its own inverse modulo eight. This reduces the enumeration to
`R/8` classes. For `h=1,2`, use the full list of odd residues.

Every factor pair has a unique ratio `a/b mod R`, so summing the local
counts gives the exact rectangle count. Maintaining the least local
row gives the least row of the entire rectangle. Neither operation
requires storing the complete factor list. The count itself has
`O(n)` bits under the endpoint assumptions.

Since corner pruning leaves `P>=N` and the theorem assumes `P<=C N`,
equation (19) gives `h=n/3+O(1)`. There are `O(2^(n/3))` ratio
queries. The polynomial cost per query and the constant local output
bound prove (2). All loops can stream their inputs; the largest active
objects have polynomial bit size.

## 6. Smallest prime factor, including unbalanced inputs

An ordinary proper factor need not be prime. To obtain the smallest
prime factor without a hidden recursion, use increasing powers of two
`A=2,4,8,...` with `A<=isqrt(N)`, and query

\[
 \begin{aligned}
 I_A&=[A,\min(2A,\lfloor\sqrt N\rfloor+1)),\\
 J_A&=[\lfloor N/(2A)\rfloor+1,\lfloor N/A\rfloor+1).
 \end{aligned}
 \tag{20}
\]

Every divisor `a` in `I_A` has its complementary divisor in `J_A`.
Every positive product in this rectangle is less than `2N`.
The exact rectangle routine therefore meets (2) with a fixed constant,
regardless of the factor sizes. At the first nonempty rectangle,
return its least row. All smaller possible nontrivial divisors were
covered earlier, so this row is the least nontrivial divisor of `N` and
is prime. If every rectangle is empty, return `N`. Even inputs return
two immediately.

There are `O(n)` rectangles. Each least-row answer already scans every
ratio class in its rectangle, so there is no uncharged prefix search.
The resulting smallest-prime-factor procedure obeys (2). If complete
factorization is desired, repeatedly remove the returned prime; at
most `n` extractions are needed, which preserves the stated bound.

## 7. Exact checks and evidence

Run:

```sh
python3 tools/check_circuit_box_ratio_lattice.py \
  --output evidence/circuit-box-ratio-lattice.json
```

The [production implementation](../tools/circuit_box_ratio_lattice.py),
[independent checker](../tools/check_circuit_box_ratio_lattice.py), and
[receipt](../evidence/circuit-box-ratio-lattice.json) contain:

- 1,000 weighted lattice reductions, with exact determinant, membership,
  Gram determinant, terminal inequalities, and exchange-potential checks;
- 1,200 local ratio factor lists and 1,200 complete rectangle counts
  compared against literal enumeration;
- 1,198 smallest-prime-factor results for inputs `2` through `1199`,
  compared against direct trial division;
- four single-class controls on inputs of 151, 301, 602, and 1,204 bits,
  with independently enumerated divisors from a supplied `3^u 5^v`
  factorization and exponentially long intervals;
- four empty rough-semiprime rectangles from the
  [inverse-window audit](circuit-box-inverse-windows.md), with complete
  ratio enumeration and exact trial-division prime certificates;
- three complete smallest-prime-factor examples at 21, 27, and 33 bits.

For the empty rectangles, the measured work is:

| Input bits | Ratio precision | Ratio queries | Lattice lines | Exact count |
|---:|---:|---:|---:|---:|
| 33 | 12 | 512 | 744 | 0 |
| 39 | 14 | 2,048 | 2,994 | 0 |
| 45 | 16 | 8,192 | 12,389 | 0 |
| 51 | 18 | 32,768 | 49,389 | 0 |

The 1,204-bit known-factorization control uses a 512-bit ratio modulus
and one lattice line. This is a local terminal control, not a complete
1,204-bit factorization benchmark. The receipt separately records all
work in the complete examples and hashes the exact checked sources.

## 8. What a further improvement must supply

The lattice terminal uses exact integer equality and the area bound
`R^3>=P`. It does not give a cheap query for arbitrarily many extra
product-congruence bits in a large rectangle. At a smaller ratio
modulus, (16) can require many lines, and that cost remains charged.

An improvement can come from a stronger local recovery theorem, a
representation that handles many ratios together, or a useful
arithmetic structure that narrows the candidate set. The preceding
positive lookahead count remains available for pruning, but it is not
needed for the bound proved here. No cancellation argument or heuristic
distribution of modular inverses is used.
