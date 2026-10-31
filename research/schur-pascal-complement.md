# Greedy Pascal complements give every one-species mode profile

## Result

For an arbitrary finite set of nonnegative mode indices, the normalized
one-species infinity matrix has a Smith profile determined by binary
linear algebra. No polynomial elimination is needed to find it.

Keep the original indices. Scan the indices outside that set in increasing
order, and retain a new index whenever its Pascal column is independent
of the columns retained so far. Stop after retaining as many columns as
there were original indices. The union of the original indices and the
retained indices is exactly the Smith order list.

The algorithm can skip whole binary intervals whose columns are already
dependent. Its traversal depends on the bit lengths of the indices, so
large gaps do not require a scan through every missing integer.

This gives an exact construction tool for the
[one-species crossings](schur-one-species-interval.md), including gapped
sets and sets that are not Boolean down-sets. It also proves the signed-odd
interval profile for every interval length. No improved factoring
exponent is established here. The derivation and checks used substantial
OpenAI assistant support.

## 1. The finite matrix

Let \(S\) be a nonempty finite set of nonnegative integers, of cardinality
\(n\). All matrix rows and columns below are indexed by \(S\), in increasing
order. Work over \(\mathbb F_2[[t]]\), and put

~~~math
P_{ab}=\binom ab\bmod2=\mathbf1_{b\subseteq a},\qquad
D=\mathrm{diag}(t^a),\qquad
Q=\mathrm{diag}((1+t)^a),\qquad
O_{ab}=t^{a\mathbin{|}b}.                               \tag{1}
~~~

The matrix of interest is

~~~math
F_S(t)=
\begin{pmatrix}
QP^{\mathsf T}D&DP\\
D&(1+t)O
\end{pmatrix}.                                        \tag{2}
~~~

It is the same normalized matrix obtained from the explicit physical
gauges in the interval calculation. That normalization uses only the
displayed entry formulas, so it also holds for an arbitrary mode set.
No identity such as \(P^2=I\) is assumed here. The constant matrix \(P\)
is invertible because it is triangular with diagonal entries one.

For an index \(c\ge0\), define its Pascal column

~~~math
u_c=\left(\binom ca\bmod2\right)_{a\in S}
    =\left(\mathbf1_{a\subseteq c}\right)_{a\in S}.
~~~

Let \(c_1<\cdots<c_n\) be the greedy independent columns among indices
\(c\notin S\). Such a list exists. If \(L\) is a power of two exceeding
\(\max S\), the columns repeat with period \(L\), and one complete
period has row rank \(n\). All indices in \([L,2L)\) are outside \(S\).

**Theorem.** The Smith orders of \(F_S(t)\) are

~~~math
\boxed{\quad
\mathrm{sort}\bigl(S\cup\{c_1,\ldots,c_n\}\bigr).
\quad}                                                 \tag{3}
~~~

The two sets in (3) are disjoint, so every order is distinct.

## 2. The complement factorization

The regular elimination from the interval calculation gives

~~~math
F_S\ \sim\ \mathrm{diag}(D,K_S),\qquad
K_S=DP+(1+t)QP^{\mathsf T}O.                            \tag{4}
~~~

To check regularity, first subtract \(QP^{\mathsf T}\) times the second
row block from the first and interchange the row blocks. Then clear
\((1+t)O\) above \(K_S\) using a column operation. This operation is
regular because every entry of row \(a\) in \(O\) is divisible by \(t^a\).

The Schur complement has the formal-series expression

~~~math
(K_S)_{ab}
=(1+t)^{a+1}
 \sum_{c\notin S}
 \binom ca\,t^c t^{b\mathbin{\&}\neg c}.               \tag{5}
~~~

Only the parity of the binomial coefficient is used. The sum is well
defined formally: its term at index \(c\) has valuation at least \(c\).

For a proof, choose a power of two \(L>\max S\). Splitting the sum over
all \(c\ge0\) into periods of length \(L\), the binary-digit identity
for the full interval gives

~~~math
\begin{aligned}
&(1+t)^{a+1}\sum_{c\ge0}\binom ca\,t^{c\mathbin{|}b}\\
&\qquad=
\frac{(1+t)^L}{1+t^L}\,t^a\binom ab
=t^a\binom ab.                                        \tag{6}
\end{aligned}
~~~

The terms with \(c\in S\) are exactly the second summand of (4).
Removing those terms proves (5), since addition and subtraction agree
in characteristic two.

Write (5) as an infinite matrix product

~~~math
K_S=U\mathrm{diag}(t^c:c\notin S)V,
\qquad
U_{a,c}=(1+t)^{a+1}\binom ca,
\quad
V_{c,b}=t^{b\mathbin{\&}\neg c}.                       \tag{7}
~~~

The constant column of \(U\) at \(c\) is \(u_c\), and the constant
row of \(V\) at \(c\) is \(u_c^{\mathsf T}\).

## 3. Determinantal minima are greedy weights

For each \(k\le n\), the greedy set \(\{c_1,\ldots,c_k\}\) is the
unique independent \(k\)-set of minimum total index weight. Indeed, if
\(b_1<\cdots<b_k\) is any independent set, then the first \(b_j+1\)
columns have rank at least \(j\), so \(c_j\le b_j\). Equality in the
sum of weights forces equality at every position.

Apply Cauchy--Binet to a \(k\)-minor of (7). A nonzero minor of \(U\)
requires the corresponding Pascal columns to be independent, since its
row factors \((1+t)^{a+1}\) are units. Every nonzero term therefore has
valuation at least

~~~math
c_1+\cdots+c_k.                                        \tag{8}
~~~

Choose \(k\) rows \(I\subseteq S\) on which the first \(k\) greedy
columns have a nonzero minor. Use the same set \(I\) for the columns of
\(K_S\). The term belonging to \(c_1,\ldots,c_k\) has leading
coefficient

~~~math
\det((u_{c_j})_{I,j})^2=1.
~~~

Every other independent index set has strictly larger total weight, so
this coefficient cannot cancel. The minimum valuation of a \(k\)-minor
is exactly (8). Taking successive differences proves that the orders of
\(K_S\) are \(c_1,\ldots,c_n\). Together with the diagonal block \(D\),
this proves (3).

One may apply Cauchy--Binet to finite truncations in this argument. Terms
outside any fixed precision have valuation at least that precision, so
no convergence assumption beyond formal power series is needed.

## 4. Every initial interval has the odd profile

Take \(S=\{0,\ldots,n-1\}\), with arbitrary \(n\ge1\). The first
\(n\) available columns are indexed by \(n,\ldots,2n-1\). They are
independent because

~~~math
\binom{n+j}{a}
=\sum_{h=0}^{a}\binom n{a-h}\binom jh                  \tag{9}
~~~

factors their square matrix into two triangular integer matrices with
diagonal entries one. Thus the determinant remains one modulo two.
The greedy complement is exactly \(\{n,\ldots,2n-1\}\), and (3) gives

~~~math
0,1,\ldots,2n-1.                                      \tag{10}
~~~

Consequently the explicit gauges in
[the interval note](schur-one-species-interval.md) have positive odd
zero orders and signed odd infinity orders for every interval length,
using any ambient power of two that contains the interval. The
interval note's analytic pressure comparison then applies to every
odd or even length, rather than only the two special size families used
in its original finite elimination proof.

## 5. A binary traversal for large gaps

Consider the aligned block

~~~math
\{p,p+1,\ldots,p+2^k-1\},\qquad 2^k\mid p.
~~~

For a row index \(a\), its Pascal entry on this block is identically zero
unless its high bits are included in the high bits of \(p\). Otherwise
the entry depends only on the low \(k\) bits of \(a\).

Group the active rows by their low-bit masks. The span of all columns in
this block is exactly the span of the indicator vectors of these groups.
The corresponding distinct low-bit inclusion functions are independent:
they are distinct rows of the full binary Pascal matrix on \(k\) bits.

This gives a depth-first greedy algorithm. Maintain the span of accepted
columns, visit binary blocks in increasing order, and discard a block
when all its group indicators lie in the current span. Otherwise split
the block in half. At a leaf, skip the excluded indices in \(S\); accept
the column if it is independent.

Every expanded branch leads to an accepted column or an excluded leaf
whose presence prevented an earlier pruning. There are at most \(n\)
accepted columns and \(n\) excluded indices. The number of expanded
nodes is therefore \(O(n\log L)\), where \(L\) is a containing power of
two. A node uses at most \(n\) indicator vectors and ordinary binary
basis reductions. This avoids work proportional to the numerical value
of the largest index.

The checker compares this traversal with a separate linear scan on small
sets. It also verifies explicit examples with indices containing up to
256 bits.

## 6. A direct pressure formula for arbitrary geometry

Let \(\nu_1,\ldots,\nu_{2n}\) be the order list (3). With the original
one-species gauges, the zero orders are two copies of
\(\{2a+1:a\in S\}\), the infinity orders are
\(\{2\nu_i+1-2n\}\), and the restored heights are
\(G_a=8a-4(n-1)\).

One may also add a height \(h\) to both remaining plus boundaries at
infinity. This adds \(h\) to every infinity order and \(2h\) to every
restored mode height. Hence the exact product-type criterion becomes

~~~math
\begin{aligned}
\mathcal B_{S,h}(s)
&=\prod_{a\in S}(1+s^{-(2a+1)})^2
  \prod_{i=1}^{2n}(1+s^{-(2\nu_i+1-2n+h)}),\\
\mathcal Z_{S,h,R}(s)
&=\prod_{a\in S}
  \left(1+s^{-[8a-4(n-1)+2h]/R}\right),\\
\mathcal B_{S,h}(s)&<\mathcal Z_{S,h,R}(s)^R.
                                                               \tag{11}
\end{aligned}
~~~

Integer heights give polynomial or Laurent monomial gauges directly.
Rational heights can be realized after a common variable substitution.
Floating optimization of \(h\) or \(s\) is only a search procedure;
a proposed strict inequality still requires exact replay.

Formula (11) makes a large geometry search inexpensive: determine one
greedy binary basis, then optimize a product of \(O(n)\) scalar factors.
It is a restricted gauge family. An unsuccessful search in it would not
exclude other physical heights or a nonempty second species.

## 7. Reproduction

```sh
python3 tools/check_schur_pascal_complement.py \
  --output evidence/schur-pascal-complement.json
```

The [checker](../tools/check_schur_pascal_complement.py) compares the
greedy formula with exact polynomial Smith elimination on every nonempty
subset of \(\{0,\ldots,7\}\) and four additional gapped sets. It also
checks the formal complement expansion, uses independent block-Toeplitz
ranks on selected cases, verifies all interval complements through size
64, and tests the binary traversal on large index gaps. The
[evidence](../evidence/schur-pascal-complement.json) records all profiles
and source hashes.
