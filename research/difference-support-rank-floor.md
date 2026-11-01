# An all-code barrier for the elementary crossing coefficient span


## Weighted block theorem: arbitrary paired codes, every field

Let A and B be arbitrary nonempty subsets of {0,1}^m, of sizes a and b.
Use A on the X,V legs and B on the Y,U legs. Over ANY field, restrict

```
product_j ((1-z_j) I_j + z_j E_j)
```

to these physical alphabets. In characteristic two this is the earlier
I+z_j(E+I) convention. No downset, linearity, site symmetry, or equal-size
assumption is made. For D in A+B, define

```
S_D = {x in A : x xor D in B},
n_D = |S_D|,
h_D = |{x restricted outside D : x in S_D}|,
g_D = |{x restricted to D : x in S_D}|.
```

Every decomposition of a nonzero scalar multiple of the crossing endpoint
into full multivariate coefficient-span tensors satisfies

```
sum_i sqrt(r_i c_i) >= sum_(D in A+B) sqrt(h_D g_D)
                     >= a*b / sqrt(min(a,b))
                     >= (a*b)^(3/4).                 (W)
```

In particular, for a=b=q the lower bound is q^(3/2), and with s pieces
some r_i*c_i >= q^3/s^2. Polynomially many pieces therefore force
rho >= 3/4 in the maximum-rank-product interface. This uniformly excludes
ALL equal-size code restrictions of this elementary seed, including
arbitrary different codes on the two crossing pairs. The unequal-size
bound also defeats the ordinary paired-exchange rank-sum envelope.

This is a coefficient-span theorem, not an obstruction to new local seeds,
context-dependent auxiliary kernels, or arbitrary factoring algorithms.
The complete factoring score remains 1/2. No new Lean proof is asserted.

### Proof

For a coefficient functional L put beta_D=L(product_(j in D) z_j), now
including D=0. A swapped entry (u,v)=(y,x), where x xor y=D, is exactly
beta_D. Outside D all four bits coincide, so no extra factor is present.

In the row cut, nonzero entries preserve both the xor difference D and
the intersection x intersection y. For fixed D, these intersections are
exactly the h_D outside projections of S_D. Every one of these h_D
separate blocks contains a swapped entry beta_D. Thus

```
r(L) >= sum_(D: beta_D != 0) h_D.
```

In the column cut, nonzero entries preserve D=x xor u=y xor v. On D,
local multiset equality forces y=u and v=x. Each of the g_D orientations
x restricted to D can couple only to its complementary column orientation.
After permuting column blocks these are separate blocks, each containing
a swapped entry beta_D. Hence

```
c(L) >= sum_(D: beta_D != 0) g_D.
```

Cauchy--Schwarz yields sqrt(r(L)c(L)) at least the sum of sqrt(h_D*g_D)
over those active differences. The endpoint identity requires
sum_i beta_(i,D)=kappa != 0 for every D in A+B. Summing the bound over
pieces therefore counts every difference at least once, proving the
first inequality of (W). This uses the same beta_D in both cuts.

An element of S_D is uniquely specified by its inside and outside
projections, so h_D*g_D >= n_D. Moreover n_D <= min(a,b) and
sum_D n_D = a*b, since ordered pairs in A x B partition by their xor.
Therefore

```
sum_D sqrt(h_D*g_D) >= sum_D sqrt(n_D)
                    >= sum_D n_D/sqrt(min(a,b)).
```

Finally max(a,b) >= min(a,b) gives the last inequality of (W).
All block witnesses are literal nonzero entries. The argument uses no
order, positivity, algebraic closure, or division by the characteristic.

### Consequence for the current paired-code lane

The earlier unweighted difference count left staggered paired downsets
open. The weighted argument above closes their ENTIRE elementary
coefficient-span envelope too, without enumerating larger matching codes.
Their exact normal forms and checkers remain useful, but a positive pivot
must change the local kernel or the way auxiliary states couple to it.
Do not reinterpret this as ruling out all nonlinear-code constructions.

## Unweighted active-difference corollary

For a nonempty binary word code C, let
Delta(C) = {x xor y : x,y in C, x != y}. Restrict the elementary
characteristic-two crossing to C, permitting independent variables at
every position. Every decomposition of its endpoint into tensors from
the full multivariate coefficient span obeys

```
sum_i sqrt(r(A_i)c(A_i)) >= sqrt(2) |Delta(C)|.       (1)
```

For a decomposition with s pieces and nonzero endpoint scale, some piece
has r(A_i)c(A_i) >= 2|Delta(C)|^2/s^2. Coefficients may lie in any
characteristic-two field. There is no symmetry, downset, or bound on
coefficient-functional support in these statements.

This extends the [difference-unique code bound](sidon-code-rank-obstruction.md),
which is the case |Delta(C)| = binom(|C|,2). It excludes the elementary
crossing on the full hard-core language and on all multiplication graphs
below. It is not a lower bound for arbitrary code-coupled auxiliary
crossings or for factoring. The complete factoring score remains 1/2.

## A same-functional witness for both cuts

The local polynomial is I + z_j(E+I). A nonzero entry requires equality
of the two local multisets. For a coefficient functional L and a set D of
positions, put beta_D = L(product_(j in D) z_j). Let a(L) be the number
of D in Delta(C) with beta_D nonzero.

**Row cut.** Rows (x,y) and columns (u,v) couple only when x xor y =
u xor v. Thus the row matrix is a direct sum over differences D.
For each active D choose x,y in C with x xor y = D. The swapped entry
(u,v) = (y,x) is exactly beta_D, so that entire D block is nonzero.
Therefore r(L) >= a(L).

**Column cut.** Rows (x,u) and columns (y,v) couple only when
x xor u = y xor v = D. At every position in D, local multiset equality
forces y=u and v=x. Hence row orientation a=x restricted to D can couple
only to column orientation D minus a. Permuting the column orientation
blocks gives a direct sum. For each active D, an ordered pair (x,u) and
its reverse (u,x) lie in two distinct row orientation blocks, since D is
nonempty. Each block contains a swapped entry beta_D. Both are nonzero,
so c(L) >= 2a(L).

The beta_D here is the SAME functional on the SAME monomial in both cuts.
No independent optimization of blocks is allowed or needed. It follows
that sqrt(r(L)c(L)) >= sqrt(2)a(L).

If sum_i A_i = kappa E_C, then sum_i beta_(i,D) = kappa != 0 for each D
in Delta(C). Every difference is active in at least one piece. Sum the
preceding bound to obtain (1), and use pigeonhole for the maximum-product
version. Cancellation among pieces does not invalidate this coverage step.

## 1. All no-adjacent-ones words are excluded, uniformly

Let C_m consist of the length-m words with no adjacent ones. Every binary
word is the xor of two such words: alternately color the ones within each
run, assigning the two colors to the two summands. Thus

```
|Delta(C_m)| = 2^m-1,
|C_m| = F_(m+2) = Theta(phi^m),  phi=(1+sqrt(5))/2.
```

Consequently every full-span decomposition has rank sum at least
sqrt(2)(2^m-1). With polynomially many pieces, its leading normalized
maximum rank-product parameter satisfies

```
rho >= 1/(2 log_2(phi)) > 5/7 > 1/sqrt(2).          (2)
```

The first strict inequality has an elementary rational check:
phi^10 = 55phi+34 < 128, since phi < 17/10. Hence log_2(phi) < 7/10.
The second follows by squaring: 25/49 > 1/2.

The exact-size increasing-binary prefixes used in
[the downset construction](downset-code-crossings.md) are covered
asymptotically too. A prefix of size N in the least sufficient C_m
contains the complete C_(m-1) with an initial zero, so its difference
set has at least 2^(m-1)-1 elements. The constant factor does not change
(2), because log_2(N) = m log_2(phi)+O(1).

This closes that language's scalar elementary coefficient-span route.
The weighted theorem above also handles unequal physical-code pairs.
Different local polynomials and new auxiliary dynamics remain outside it.

## 2. Multiplication graph codes are excluded at every field size

Let C_k = {(a,b,ab): a,b in F_(2^k)}, written in any fixed binary basis.
Then q=|C_k|=2^(2k). For a nonzero input difference (h,g), the output
difference is

```
(a+h)(b+g)+ab = a*g+b*h+h*g.
```

The linear part is surjective onto F_(2^k), so there are exactly 2^k
possible output differences for each of the 2^(2k)-1 nonzero input
differences. Zero input difference has zero output difference. Therefore

```
|Delta(C_k)| = (q-1)sqrt(q),
sum_i sqrt(r_i c_i) >= sqrt(2)(q-1)sqrt(q) > q^(3/2)
```

for q >= 4. The strict inequality is equivalent to
2(q-1)^2 > q^2. Since 3/2 > sqrt(2), this defeats the repeated-block
target at every k >= 1. With polynomially many pieces, the maximum-product
version gives leading rho >= 3/4. A change of binary field basis does not
alter the difference count. This upgrades the earlier k=2 diagnostic to
a uniform family exclusion without a larger search.

## 3. Unequal paired alphabets: a weaker but useful bound

For A on the X,V pair and B on the Y,U pair, the same proof gives

```
sum_i sqrt(r_i c_i) >= |A+B|,
```

where addition is xor, and zero is included if present. The row and column
cuts each have a nonzero block for every active cross-difference. The
factor two need not survive because reversing a pair may leave A x B.
The weighted theorem above supplies the stronger paired-code conclusion.
The factor two in the unweighted bound still requires equal codes.

## Evidence and scope

```
python3 tools/check_difference_support_floor.py > /tmp/difference-support.json
cmp /tmp/difference-support.json evidence/difference-support-floor.json
```

The exact F2 checker exhausts all 273 nonempty codes at lengths one to
three and all 65,532 multivariate coefficient functionals. It checks both
the unweighted and weighted block bounds. All 225 ordered pairs of nonempty
two-bit codes add 3,600 paired-code functionals. A separate signed dense
constructor and modular elimination over F3 check 1,215 functionals on all
15 nonempty two-bit codes, using (1-z)I+zE. Additional sumset controls cover
hard-core words through length 14 and multiplication graphs through k=4,
plus the rational phi bound.
These controls support the written all-size proof; they are not independent
review or a new Lean proof. No improved factoring bound is asserted.

The arbitrary-code entry formula and source credit follow
[downset-code-crossings.md](downset-code-crossings.md). The rank-sum target
and block amplification are in
[auxiliary-coefficient-span.md](auxiliary-coefficient-span.md).
This work used substantial OpenAI assistant support.

## Relation to the earlier full-conservation theorem

The [full-conservation theorem](full-conservation-rank-sum-barrier.md), retained from [PR #27](https://github.com/1337isnot1337/integer-factoring/pull/27), proves the same `d^(3/2)` rank-sum floor on the full binary cube for every kernel `[x+y=u+v] K(x+y,x+u)` over every field of characteristic two. Its matched Cauchy--Binet and Boolean-degree argument allows a broader kernel space than the elementary pencil span studied here.

The present result instead allows arbitrary coordinate-code pairs and every coefficient field, but retains the elementary multivariate coefficient-span hypothesis. Restricting to code subsets need not preserve the convolution blocks used by the full-conservation proof. Thus neither result generally supersedes the other, and the common full-cube bound is not claimed as new relative to PR #27. Both are restrictions on construction classes, not lower bounds on factoring.
