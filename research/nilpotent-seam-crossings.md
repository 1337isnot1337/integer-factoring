# Nilpotence on a braid seam gives complementary asymmetric channels

A four-state nilpotent auxiliary placed on the central seam of the
[bounded-path braid](quantum-braid-path-quotient.md) gives whole coefficient
families with explicit asymmetric rank bounds. The opposite seam marks
the physical block boundary. Both seams have n gates, independent of
the braid's total n-squared size.

For an L-letter alphabet at its least physical length n, the highest
allowed central-degree coefficients have row rank

```math
L^{1/2+o(1)},\qquad \varphi=(1+\sqrt5)/2,        \tag{1}
```

and column rank at most `L^2`. The zero boundary-degree coefficients
have polynomial column rank and row rank at most `L^2`. These are
complementary asymmetric families, rather than bounds requiring both
ranks to be small. At eight sites one highest-degree coefficient has
the exact pair **(8,169)** on the full13-letter alphabet.

These cheap edges cannot span the crossing by themselves at all
large widths: the existing minimum-cut budget forces additional
balanced-rank pieces. The remaining task is a mixed decomposition
combining the asymmetric channels with such middle pieces. Thus (1)
is a uniform bound for selected coefficients, not a rank-sum bound
for a complete crossing. The strongest complete factoring score remains `1/2`.

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
`n-m` disjoint central pairs onto cups. The relevant image dimension
is smaller than the ordinary physical bound `4^m`: it is exactly
the dimension of the **2m-strand vacuum quotient**.

Write `v_(2m)` for the number of height-three return paths of length
2m, with `v_0=1`. Insert the `n-m` adjacent cups to obtain
`C:V_(2m)->V_(2n)`, and remove them with caps to obtain D in the
opposite direction. These maps are well defined on the same
negligible-wall quotient: its defining ideal is killed by every
closed vacuum pairing, including after inserting or removing cups.
The explicit cup/cap identities give

```math
DC=\alpha^{n-m}I,\qquad
CD=\prod_{\text{mandatory sites}}e_i.
```

The scalar is nonzero, so C is injective, D is surjective, and the
displayed central projector has exact rank `v_(2m)`. This argument
uses the quotient representation, without an ordinary ambient
invariance assumption. All additional central physical factors
commute with the mandatory projectors, so their common image is
preserved. Including the auxiliary, each position set therefore
contributes rank at most `4*v_(2m)`.

The circuits before and after A are polynomials in y. Extracting
degree b is a convolution of at most b+1 left/right factors.
Projection onto the boundary frames cannot increase matrix rank.
Consequently the actual row cut satisfies the explicit all-length
bound

```math
\boxed{\quad r(G_{a_*,b})
 \le8(n+1)^3(b+1)v_{2\lfloor n/4\rfloor}.\quad} \tag{4}
```

The quotient matrices and their constant rank factors are over F4;
the earlier quotient comparison identifies their boundary blocks
with the ambient encoder calculation. No extension-degree penalty
is needed for (4). The highest allowed coefficient can vanish: at six
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

## 4. A mixed decomposition is necessary

Since `v_(2m)=O(phi^(2m))` and `L=Theta(phi^n)`, equation (4) gives
(1), while (5) is `L^o(1)`. Thus the union of these two O(n)-size edge families has
the favorable rank budget

```math
\sum_{\text{edges}}\sqrt{r(T)c(T)}
 \le L^{5/4+o(1)},\qquad 5/4<\sqrt2.            \tag{6}
```

Equation (6) is a budget **for the listed pieces**. The earlier
[minimum-cut theorem](minimum-cut-rank-budget.md) gives the independent
necessary condition

```math
\sum_i\min(r(T_i),c(T_i))\ge L                 \tag{7}
```

for every decomposition of SWAP. One can also see (7) directly by
fixing `Q=R=q`: SWAP becomes the L-by-L identity on P,S, whereas
each piece's slice rank is at most either of its two cut ranks.

The sum of the chosen small ranks of the two edges, including I,
is bounded by

```math
B_n=1+4(n+1)^4(n+2)v_{2\lfloor n/4\rfloor}
       +8(a_*+1)(a_*+2)=o(L).                  \tag{8}
```

Thus they **cannot** span SWAP for all large widths. Any polynomial
number of changed auxiliary templates having these same edge bounds
still gives `poly(n)B_n=o(L)` and cannot fix this by itself. This
consequence supersedes the earlier edge-only spanning target; no
larger scan of that target is warranted.

A complete decomposition must add middle pieces whose minimum ranks
sum to at least `L-o(L)`. With polynomially many such pieces, at least
one must have **both** ranks at least `L^(1-o(1))`. That is compatible
with the strict rank-sum threshold: a piece with both ranks O(L)
costs O(L), rather than `L^sqrt(2)`. The task is to find compatible
balanced pieces and exact coefficients, while retaining the cheap
asymmetric corrections and respecting all other cut budgets.

One available source of balanced pieces, already present in the
[finite-state pairing family](finite-state-pairing-crossings.md), is
the group-sum tensor. For any cyclic label group put
`K_sigma(P,Q,R,S)=[P+sigma*Q=S+sigma*R]`. Both cut ranks equal L,
because each cut groups rows and columns by their displayed linear
sum; coefficient one on P and S makes every class nonempty.
For binary labels XOR gives the analogous F2-vector-space version.
Such known cheap cores may be appended to a strict crossing by
`(t^k+t^(k+1))K_sigma`, beyond its old degree. Both endpoints are
preserved. This supplies an effective hybrid design space, not an
assertion that any chosen core completes the decomposition. In
particular, neither the pairing-only barriers nor an edge-only
exclusion settles a hybrid with the new path channels.

Once a mixed favorable spanning identity is found, expose a solution
for `SWAP+I` in distinct positive monomial degrees and use constant
coefficient I. Exact dense solving and factors have the allowed
`2^O(lambda)` cost. Finding that identity, with middle-rank mass and
the [other necessary budgets](geometric-cut-rank-budget.md) satisfied,
is the open positive task. No claim that (6) already applies to a
whole crossing is made.

## 5. Reproduction and finite scope

[The implementation](../tools/quantum_braid_nilpotent_seam.py) constructs
the complete bivariate tensors. The
[replay](../tools/check_quantum_braid_nilpotent_seam.py) checks endpoints,
degree bounds, exact ranks of all mandatory-cup placements at4/6/8sites,
constant auxiliary identities, independent scalar
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
