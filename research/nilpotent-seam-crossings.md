# Nilpotence on a braid seam gives complementary asymmetric channels

A four-state nilpotent auxiliary placed on the central seam of the
[bounded-path braid](quantum-braid-path-quotient.md) gives whole coefficient
families with explicit asymmetric rank bounds. The opposite seam marks
the physical block boundary. Both seams have n gates, independent of
the braid's total n-squared size.

For an L-letter alphabet at its least physical length n, the highest
allowed central-degree coefficients have row rank

```math
L^{\beta+o(1)},\qquad
\beta=\frac{\log2}{2\log\varphi}=0.72021\ldots<\frac34,
\qquad \varphi=(1+\sqrt5)/2,                    \tag{1}
```

and column rank at most `L^2`. The zero boundary-degree coefficients
have polynomial column rank and row rank at most `L^2`. These are
complementary asymmetric families, rather than bounds requiring both
ranks to be small. At eight sites one highest-degree coefficient has
the exact pair **(8,169)** on the full13-letter alphabet.

The crucial remaining obligation is a spanning theorem: these cheap
edge families do **not** span the crossing in the checked small
alphabets. Thus (1) is a uniform bound for selected coefficients,
not a rank-sum bound for a complete crossing. The strongest complete
factoring score remains `1/2`.

## 1. Two disjoint seams and a strict crossing

In the n-by-n block-braid grid, cell `(j,k)` has local index
`n+j-1-k`. Order cells by `j+k`, then by j. This changes the ordinary
braid word only by commuting disjoint local gates. For even n the
two sets

```math
A=\{j+k=n-1\},\qquad B=\{j=k\}
```

are disjoint and each has n cells. A is the central anti-diagonal:
its local indices are `0,2,...,2n-2`. B consists of the n occurrences
of `e_(n-1)`, exactly the gates crossing the physical block boundary.

Over F4 define

```math
Q=\begin{pmatrix}
0&0&0&1\\0&0&1&0\\0&1&\alpha&0\\1&0&0&\alpha
\end{pmatrix},\qquad
N=\begin{pmatrix}
0&1&0&0\\0&0&1&0\\0&0&0&1\\0&0&0&0
\end{pmatrix},\qquad \alpha^2=\alpha+1.
```

Then `Q^2=I+alpha Q`, `Q^5=I`, `N^4=0`, and
`N^3 Q e_0=e_0`. Q is two copies of the regular F16 braid parameter.
Use the following gates in the stated anti-diagonal order:

```math
\begin{array}{ll}
A:& (I+x(Q-I))\otimes I+(xI+x(1+x)N)\otimes e_i,\\
B:& Q\otimes I+yI\otimes e_i,\\
\text{elsewhere}:&Q\otimes I+I\otimes e_i.
\end{array}                                                   \tag{2}
```

The noncommuting auxiliary is introduced **after** selecting the
anti-diagonal circuit. Its gates away from x=1 need not commute on
disjoint physical sites. No such commutation is used.

Start in auxiliary state `e_0`. For `k=-n^2/2 mod5`, project by
`ell=(e_0^T+tau e_3^T)Q^k`, with
`tau=(0,1,0,alpha,alpha+1)[k]`. At `x=y=1` every gate is the regular
representation of `qI+e_i`. The earlier all-length braid theorem
therefore gives `G(1,1)=SWAP` exactly. All path matrices and projected
coefficients `G_ab` are in F4. The negligible-wall comparison applies
termwise to the physical TL words, including auxiliary matrix weights.

The bivariate matrix need not equal I at `(0,0)`. A strict univariate
crossing is explicitly

```math
H(t)=(1+t)I+tG\bigl(t,t^{2n+1}\bigr).            \tag{3}
```

There are no coefficient collisions: the x-degree is at most2n.
The coefficient at t is `I+G_00`; all other nonconstant coefficients
are individual `G_ab`. The degree is `O(n^2)`, both strict endpoints
hold, and the extra I term has rank pair `(L^2,1)`.

Take any common exact binary alphabet at the least even n, using the
existing suffix-count unranking. Its size L satisfies `L=Theta(phi^n)`.
Restriction cannot increase either rank below. The explicit path
propagation, coefficient selection and dense factors cost `2^O(lambda)`
for binary width lambda; all polynomial degrees and auxiliary sizes
have the required bounds.

## 2. A degree bound from nilpotent runs

Put `m=floor(n/4)` and `a_*=2n-m`. A degree-two central factor is
`N tensor e_i`. Expand every degree-one auxiliary factor into Q, I
and N letters. A nonzero word has at most three N letters between
successive Q letters, since `N^4=0`.

If k central gates supply degree two and q supply genuine Q letters,
then `k<=(q+1)3` and `k+q<=n`. Thus `k<=n-m`, and total degree is
at most `n+k<=a_*`. A term at degree `a_*` must have exactly
`n-m` degree-two gates and m degree-one gates; no degree-zero gate
can occur.

There are at most `2(n+1)^3` possible sets of degree-two positions
that can contribute to this highest degree. To prove the count,
write `n=4m+s`, `0<=s<4`. If all m degree-one positions supply Q,
the m+1 runs of N have total slack `3-s<=3` from their maximum
length three. Stars and bars bounds their possible lengths by
`(n+1)^3`. A genuine Q count of m-1 is possible only when s=0;
then the N runs are all saturated and the single other degree-one
letter must be I. Its possible insertion positions contribute at
most n more patterns. Fewer Q letters cannot support `n-m` N's.
The alternative N letter in that exceptional position would give
four consecutive N's and vanish.

For every such position set, the mandatory `e_i` factors project
`n-m` disjoint central physical pairs onto one-dimensional cups.
The remaining m pairs have total dimension `4^m`. All additional
central physical factors commute with these projections, so summing
the remaining auxiliary choices does not enlarge their image.
Including the four-state auxiliary, this central operator has rank
at most `4*4^m` for each position set.

The circuits before and after A are polynomials in y. Extracting
degree b is a convolution of at most b+1 left/right factors.
Projection onto the boundary frames cannot increase matrix rank.
Consequently the actual row cut satisfies the explicit all-length
bound

```math
\boxed{\quad r(G_{a_*,b})
 \le8(n+1)^3(b+1)4^{\lfloor n/4\rfloor}.\quad}  \tag{4}
```

This bound is proved on physical matrices over F16, and gives the
same rank for their F4 path entries: matrix rank is unchanged by
extension of the coefficient field. No extension-degree penalty is
needed for (4). The highest allowed coefficient can vanish: at six
sites `a_*=11` is zero and the actual degree is10. This is recorded,
not replaced by an unsupported bound on the lower coefficient.

## 3. Polynomial opposite-cut factors on the other edge

To form `G_a0`, every B gate supplies Q and no physical crossing.
All gates outside A then act within the two physical blocks. Their
auxiliary coefficients are polynomials in Q. Over F16 diagonalize Q;
the fixed input and output use its first regular two-dimensional
copy. There are at most two incoming and two outgoing eigenchannels.

On A, the even-index gates in the first physical block precede
those in the second block. Cutting between them leaves at most four
auxiliary states. For each of the at most `2*2*4=16` channels, the
whole physical operator factors as a first-block polynomial in x
times a second-block polynomial in x. Extracting degree a gives at
most a+1 products. This is exactly the column partition, so

```math
\boxed{\quad c(G_{a,0})\le16(a+1).\quad}         \tag{5}
```

The coefficient `I+G_00` in (3) therefore has column rank at most17.
Again rank over F4 equals rank after the temporary F16 extension.

## 4. The remaining spanning problem

Since `L=Theta(phi^n)`, equation (4) gives (1), while (5) is
`L^o(1)`. The exact inequality `phi^3=2+sqrt(5)>4` implies
`beta<3/4`. Thus the union of these two O(n)-size edge families has
the favorable rank budget

```math
\sum_{\text{edges}}\sqrt{r(T)c(T)}
 \le L^{1+\beta/2+o(1)}
 < L^{11/8+o(1)},\qquad 11/8<\sqrt2.            \tag{6}
```

Equation (6) is a budget **for the listed pieces**. To use the
[strict amplification criterion](geometric-cut-rank-budget.md), the
crossing must belong to their span, possibly after a changed seam
template or an enlarged family with the same rank bounds. If it did
at every sufficiently large binary width, exact dense linear solving
would construct the coefficients in `2^O(lambda)` time. Expose a
solution for `SWAP+I` in distinct positive monomial degrees and use
constant coefficient I. This makes a strict crossing with precisely
these cheap pieces and polynomial degree.

The current small edge spans miss SWAP. There is no claim that (6)
already applies to a whole crossing. The new positive object is the
pair of complementary asymmetric edge families. Controlling the
interior through another auxiliary or a spanning identity is the
open task.

## 5. Reproduction and finite scope

[The implementation](../tools/quantum_braid_nilpotent_seam.py) constructs
the complete bivariate tensors. The
[replay](../tools/check_quantum_braid_nilpotent_seam.py) checks endpoints,
degree bounds, constant auxiliary identities, independent scalar
propagation over F16, exact cut profiles and strip-span membership.
The [source-pinned receipt](../evidence/nilpotent-seam-crossings.json)
records full4/6/8-site results. These have coefficient-span dimensions
10,69,118. At eight sites the highest central edge includes the pairs
`(25,169)`, `(16,169)` and `(8,169)`; the opposite edge includes
`(134,2)`, `(110,4)` and `(144,6)`.

```sh
make verify-braid-nilpotent-seam
```

The optional replay needs NumPy. This extends the earlier
[matrix-auxiliary construction](quantum-braid-matrix-auxiliary.md)
by placing noncommutativity on one contiguous seam and combining
nilpotent run counting with the central cup factorization. It uses
the earlier braid and path-quotient proofs with their retained source
credit. No new Lean theorem or improved complete factoring score is
claimed.
