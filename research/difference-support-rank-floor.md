# Active differences force two-cut ranks in restricted elementary crossings

## Theorem

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
It does not exclude different local polynomials, new auxiliary dynamics,
or the separately constructed *unequal* physical-code pairs.

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
Do not apply the equal-code factor or the hard-core conclusion to the
staggered paired-code lane without checking its alphabets.

## Evidence and scope

```
python3 tools/check_difference_support_floor.py > /tmp/difference-support.json
cmp /tmp/difference-support.json evidence/difference-support-floor.json
```

The exact F2 checker exhausts all 273 nonempty codes at lengths one to
three and all 65,532 multivariate coefficient functionals. It checks
r >= active differences and c >= twice active differences for every
case. Additional sumset controls cover hard-core words through length
14 and multiplication graphs through k=4, plus the rational phi bound.
These controls support the written all-size proof; they are not independent
review or a new Lean proof. No improved factoring bound is asserted.

The arbitrary-code entry formula and source credit follow
[downset-code-crossings.md](downset-code-crossings.md). The rank-sum target
and block amplification are in
[auxiliary-coefficient-span.md](auxiliary-coefficient-span.md).
This work used substantial OpenAI assistant support.
