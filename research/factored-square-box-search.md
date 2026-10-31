# Exact square-box search and finite relation factorizations

## Status

This note supplies a deterministic, exact algorithm for the variable-cap
square-representation problem isolated in
[the decision reduction](smooth-square-decision-reduction.md). Its cost is
exponential in the number and size of the exponent caps. It is **not** the
polynomial representation oracle required by
[the conditional fixed-strip deduction](fixed-strip-smooth-square.md),
and it does not improve the proved complete factoring exponent.

The implementation nevertheless completes the arithmetic route on two
explicit inputs of 51 and 60 bits. It constructs all dominant relations,
derives an annihilator from their determinant, and extracts proper factors.
The inputs to the factoring driver are only the integer, the number of
small split primes, and the cap. Neither the factors, a class number, an
annihilator, nor a complete class-group enumeration is supplied.

Meet-in-the-middle search in quadratic class groups is an established
technique. The primitive supplied-norm version and its prior-art comparison
are recorded in [the factored-norm note](factored-norm-search.md), including
Mamah's [Theorem 5.4](https://arxiv.org/abs/2502.11402v2). The work here makes
the nonprimitive, variable-cap square problem executable, balances a
single large cap, reuses the index across the dominant-row targets, and
records complete finite factoring certificates.

## 1. Exact problem and its exponent box

Let \(\Delta<0\) be an order discriminant. Let \(p_1,\ldots,p_m\) be
distinct split rational primes coprime to \(\Delta\), with supplied roots
\(r_i^2\equiv\Delta\pmod{4p_i}\). Write

\[
g_i=\left[p_i,r_i,\frac{r_i^2-\Delta}{4p_i}\right].
\]

The caps \(t_i\) are nonnegative integers, and the target is a primitive
positive definite form \(Q\) of discriminant \(\Delta\). Put

\[
M=\prod_i p_i^{t_i},\qquad
V=\prod_i(2t_i+1).
\]

The ideal factorization argument in
[the universal-square note](universal-smooth-square.md) gives the exact
equivalence

\[
Q(x,y)=M^2\text{ for some integers }x,y
\quad\Longleftrightarrow\quad
[Q]=\prod_i g_i^{2e_i}\text{ with }-t_i\le e_i\le t_i.
\tag{1}
\]

The representation need not be primitive. For a word on the right, the
primitive represented norm is
\(\prod_i p_i^{2|e_i|}\). Multiplying its two coordinates by
\(\prod_i p_i^{t_i-|e_i|}\) gives a representation of \(M^2\).
Conversely, divide any representation by the coordinate gcd and factor
its primitive norm over the supplied primes. Each prime power has a
single orientation; its exponent is even because the original norm is
a square. This recovers a word in the stated box.

Set \(d_i=2t_i+1\), \(f_i=e_i+t_i\), and \(a_i=g_i^2\). The exact
search problem becomes

\[
\prod_i a_i^{f_i}=[Q]\prod_i a_i^{t_i},
\qquad 0\le f_i<d_i.
\tag{2}
\]

This is a bounded product search in the proper class group, with \(V\)
possible words. Both a positive answer and a negative answer are exact.

## 2. Balancing arbitrary caps

Simply dividing the primes into two sets can leave one side much larger
when one cap dominates. The implementation partitions the mixed-radix box
into at most two rectangles, each with balanced independent halves.

Let \(k=\lfloor\sqrt V\rfloor\). If \(V=1\), there is only the identity
word. Otherwise choose the first coordinate \(j\) for which the prefix
product \(P d_j\) exceeds \(k\), where
\(P=\prod_{i<j}d_i\le k\). Put

\[
L=\lfloor k/P\rfloor,\qquad d_j=qL+r,\quad 0\le r<L.
\]

The full rectangle writes \(f_j=Lu+v\), with
\(0\le u<q\) and \(0\le v<L\). Its left half contains the prefix and
\(v\); its right half contains \(u\) and the suffix. If \(r>0\), the
second rectangle writes \(f_j=qL+v\), with \(0\le v<r\). The fixed
offset \(qL\) is removed from that rectangle's target.

These rectangles are disjoint and cover the box exactly. If
\(R=\prod_{i>j}d_i\), their half sizes are respectively

\[
(PL,qR)\quad\text{and, if }r>0,\quad(Pr,R).
\]

The left sizes are at most \(k\). Also
\(k+1\le P(L+1)\le2PL\), hence

\[
qR\le\frac V{PL}<2(k+1),
\qquad R\le\frac V{Pd_j}<k+1.
\]

Every half therefore has \(O(\sqrt V)\) assignments, including when
there is only one enormous cap. Sharing the split coordinate causes no
dependence: its two contributions add, and their stated ranges produce
exactly the relevant digit interval.

## 3. Deterministic index, witnesses, and cost

For each half, enumerate products of the reduced forms \(a_i\), keeping
the associated digit assignments. Sort the canonical reduced forms and
deduplicate equal class values after each layer. One witness suffices per
class because the remaining half choices are independent. Canonical
reduction, composition, inversion, and powering use the exact arithmetic
from [the factored-norm implementation](factored-norm-search.md).

For a target, enumerate right-half values and binary-search the sorted
left table for the required complementary class. This uses exact class
equality. It does not require the lexicographic order on reduced forms to
be translation invariant, and it does not use randomized hashing.

One index can answer any number of targets with the same discriminant,
primes, and caps. Exhausting both rectangles proves a negative answer.
For a positive answer, the retained assignments give the bounded word.
The supplied-word coordinate constructor then gives actual integers
\((x,y)\), and the implementation checks \(Q(x,y)=M^2\) directly.

Let \(b=\lceil\log_2|\Delta|\rceil\), and let \(s\) denote the total
bit length of the supplied primes and caps. Building the index and each
target search have deterministic bit cost

\[
\widetilde O(\sqrt V)\operatorname{poly}(b,m,s)
\tag{3}
\]

and space bounded by the same expression. The logarithm from sorting is
included in the polynomial factors because \(\log V\) is bounded by the
cap-input length. Integer coordinate construction additionally costs a
polynomial in the expanded norm length
\(\log(M^2)=2\sum_i t_i\log p_i\). The decision and word API can omit
coordinate construction, so a cap given in binary need not first be
expanded into the integer \(M^2\).

The evidence fields named `index_group_multiplications` and
`query_group_multiplications` count composition calls in the index and
lookup loops. They exclude polynomial binary-power setup and coordinate
construction. They are useful finite operation counts, not complete
runtime measurements.

For a uniform cap \(T\), the charged exponential factor is
\((2T+1)^{m/2}\). A polynomial-bit universal square does not make this
factor polynomial. A box covering a subgroup \(H^2\) must have
\(V\ge |H^2|\), so this exact enumeration alone gives no route past the
general square-root search scale for such a covering box.

## 4. Dominant relations and factor extraction

For the finite factoring driver, take \(\Delta=-4N\) with odd \(N>1\).
Choose the first \(m\) split unramified rational primes by deterministic
trial enumeration and the least nonnegative splitting root used by the
implementation. Fix the uniform cap \(T\), then set

\[
K=Tm+1.
\]

For each \(i\), search the target \(g_i^{2K}\). A returned half-word
\(e^{(i)}\) yields the relation row

\[
R_i=2K\,\mathbf e_i-2e^{(i)}.
\]

Its diagonal entry is positive and exceeds the sum of the absolute
off-diagonal entries by at least two. Thus the integer matrix \(R\) is
nonsingular. If \(E=|\det R|\), the adjugate identity implies that every
supplied generator is annihilated by \(E\). The driver checks all the
group identities again directly.

It then runs [the early-stop extraction](early-stop-class-factoring.md).
The small bases used in these finite runs have no certificate of generating
the entire class group; none is asserted. A returned proper divisor is
checked directly. A failed bounded search or an unhelpful subgroup is
reported as such and is never used as a primality certificate.

There is also a possible earlier factor: compare the known square root
\(g_i^K\) with the recovered root \(\prod_j g_j^{e_j^{(i)}}\). Their
quotient is ambiguous, and a helpful ambiguous class gives coprime factors
of \(N\). The evidence records this check but still completes the full
relation matrix and annihilator route.

## 5. Completed finite constructions

| Quantity | 51-bit input | 60-bit input |
| --- | ---: | ---: |
| \(N\) | 1621982935400141 | 1000000016000000063 |
| Supplied \((m,T)\) | \((9,4)\) | \((12,3)\) |
| \(K=Tm+1\) | 37 | 37 |
| \(V=(2T+1)^m\) | 387420489 | 13841287201 |
| Left and right class counts | 19683 each | 117649 each |
| Index loop compositions | 39362 | 235295 |
| Total target-query loop compositions | 4628 | 38543 |
| Bit length of \(M^2\) | 286 | 361 |
| Minimum row dominance slack | 30 | 14 |
| Derived \(E\) | 63391378092879872 | 22962539960812291571712 |
| Bit length of \(E\) | 56 | 75 |
| Returned factorization | \(32452843\cdot49979687\) | \(1000000007\cdot1000000009\) |

The 51-bit basis is
\(3,5,7,11,23,31,37,41,43\). Its early-stop routine examines two
generators. The first produces the unhelpful class
\((2,2,810991467700071)\). The second has helpful top
\((41216265,17526844,41216265)\), whose equal-outer-coefficient formula
returns the stated factors. Four of its nine direct root differences
already give the same factorization.

The 60-bit basis is
\(3,7,13,19,31,47,59,61,67,79,89,101\). The first processed generator
has helpful top \((1000000007,0,1000000009)\), which immediately exposes
the factors. Five of its twelve direct root differences are helpful.

These are finite existence-and-construction certificates for the
arithmetic pipeline. They make no claim to outperform practical factoring
algorithms on these inputs or to establish an asymptotic score increase.

## 6. Independent checks and reproduction

[The exact-decision evidence](../evidence/factored-square-box-search.json)
compares all reduced target classes for the selected small discriminants
against an independent coordinate enumeration. It includes 3277 exact
decisions: 1951 positive and 1326 negative. Of the returned witnesses,
1062 are nonprimitive; 69 successful searches use the remainder
rectangle. There are also 100 nonreduced target controls.

One large-cap control uses a 260-bit discriminant, one split prime, and
cap 1000000. Its two rectangles have half sizes \((1414,1414)\) and
\((605,1)\), and it recovers the supplied exponent 123456 exactly without
constructing the 2000001-bit full norm. A separate 131-bit discriminant
uses sixteen split primes with cap one and verifies integer coordinates.
Those two targets are generated from supplied words, so they are checks
of extraction and arithmetic, not generic coverage claims.

The factoring certificates are
[the 51-bit record](../evidence/square-relation-factoring-51bit.json) and
[the 60-bit record](../evidence/square-relation-factoring-60bit.json).
[The independent replay](../evidence/square-relation-factoring-replay.json)
does not rerun the meet-in-the-middle solver. It recovers signed words
from the integer coordinates, verifies every dominant relation, computes
the determinant by rational elimination rather than the construction's
Bareiss algorithm, verifies every annihilation identity, and replays the
ambiguous-class extraction. Trial division certifies the returned factors
as prime only after extraction; those factors are never solver inputs.
The construction runs were also repeated and produced identical arithmetic
records.

From the repository root:

```sh
python3 tools/check_factored_square_box_search.py
python3 tools/run_square_relation_factoring.py --n 1621982935400141 --primes 9 --cap 4
python3 tools/run_square_relation_factoring.py --n 1000000016000000063 --primes 12 --cap 3
python3 tools/check_square_relation_factoring.py
```

The exact search implementation is
[`factored_square_box_search.py`](../tools/factored_square_box_search.py).
The remaining asymptotic gap is an algorithm that replaces its explicitly
charged exponent-box enumeration while preserving exact witnesses or the
negative-instance decisions needed by the search-to-decision reduction.
