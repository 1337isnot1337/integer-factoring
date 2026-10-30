# The highest channel retains its rank on the actual binary labels

## Theorem and scope

The [highest-channel rank theorem](top-channel-rank.md) survives the
source's restriction to binary labels.

Let `h=hOf(lambda)` be the least integer at least two with
`p(h)>=2^lambda`. Let `S_lambda` contain the first `2^lambda` partitions
of `h`, in the source's lexicographic order. Define

```math
U_\lambda:\mathbb Q^{S_\lambda\times S_\lambda}\longrightarrow V_{2h},
\qquad e_{(\xi,\eta)}\longmapsto (t^\xi)[h]t^\eta.
```

Then, with `kappa=pi sqrt(2/3)/ln 2`,

```math
\log_2\operatorname{rank}_{\mathbb Q}U_\lambda
=\kappa\sqrt{2h}+o(\sqrt h)
=\sqrt2\,\lambda+o(\lambda).
\tag{1}
```

The proof also gives the quantitative bounds

```math
\sqrt2\,\lambda
-O\!\left(\frac{\lambda\ln\ln\lambda}{\ln\lambda}\right)
\le\log_2\operatorname{rank}U_\lambda
\le\sqrt2\,\lambda+O(\ln\lambda).
\tag{1a}
```

These bounds apply for sufficiently large `lambda`.

Using the classical full Hardy--Ramanujan asymptotic, Section 6 sharpens
the bounds to

```math
\begin{aligned}
\sqrt2\,\lambda-(3-2\sqrt2)\log_2\lambda-O(1)
&\le\log_2\operatorname{rank}_{\mathbb Q}U_\lambda\\
&\le\sqrt2\,\lambda+(2\sqrt2-2)\log_2\lambda+O(1).
\end{aligned}
\tag{1b}
```

This is a theorem about one factor map of the central piece.
It does **not** establish the rank of that complete piece.
Cancellations between different channels remain possible.
It gives no lower bound for factoring runtime.

The source is `anthropics/formal-math` at
[`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The source uses descending lists for partitions and orders these lists
lexicographically. Its declarations are
[`AppHalfPieces.desc`, `lexNum`, and `IsMu`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
The implemented order is identified in
[`GeneralPieces/Half/Order.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Half/Order.lean).

Sections 1--5 use only the source's coarse partition estimates.
Section 6 uses the classical full asymptotic, cited there separately.
Neither argument differentiates an asymptotic equivalence or requires
a new Lean formalization.

## 1. Bounded partitions lie in the actual prefix

Write `P_m(n)` for the number of partitions of `n` whose parts are at
most `m`. Put `a=pi/sqrt 6` and `c=2a=pi sqrt(2/3)`.
These constants use natural logarithms. The source's partition bounds imply

```math
\ln p(n)=c\sqrt n+O(\ln(n+2)).
\tag{2}
```

Fix any real `C>=4`. Set `m=floor(C sqrt h)` and `t=a/sqrt h`.
For sufficiently large `h`, both `m>=1` and `m<=h` hold.
The generating function for bounded partitions is

```math
Z_m(t)=\prod_{j=1}^m(1-e^{-tj})^{-1}.
```

Because its coefficients are nonnegative,
`P_m(h)<=exp(th) Z_m(t)`.
For `f(u)=-ln(1-exp(-u))`, monotonicity gives

```math
\sum_{j\ge1}f(tj)
\le\frac1t\int_0^\infty f(u)\,du
=\frac{\pi^2}{6t}.
```

The integral follows by expanding `-ln(1-z)` and integrating its
nonnegative terms: the result is `sum_(k>=1)1/k^2=pi^2/6`.
The omitted terms with `m<j<=2m` contribute at least
`m exp(-2tm)`, since `-ln(1-z)>=z` for `0<z<1`.
Therefore

```math
\ln P_m(h)
\le th+\frac{\pi^2}{6t}-m e^{-2tm}
\le c\sqrt h-\delta_C\sqrt h
\tag{3}
```

for some fixed `delta_C>0` and all sufficiently large `h`.
For example, `delta_C=(C/2) exp(-2aC)` suffices once
`m>=C sqrt(h)/2`.

Equation (2) gives `ln p(h-1)>=c sqrt h-O(ln h)`.
Thus (3) implies `P_m(h)<p(h-1)` eventually.
Minimality of `h` also gives `p(h-1)<2^lambda` when `h>2`.

Every partition with largest part at most `m` comes before every
partition with a larger first part. These bounded partitions therefore
form an initial segment of the source's order. Their number is smaller
than `2^lambda`, so every one belongs to `S_lambda`.

The constant `C` remains fixed throughout this argument.
The threshold can depend on `C`. We send `C` to infinity only after
establishing the limiting bounds for each fixed `C`.

## 2. Bounded parts retain nearly all leading partition entropy

For integers `n,m>=1`, replace each part greater than `m` by that many
ones. This maps partitions of `n` to partitions counted by `P_m(n)`.
Let `k=floor(n/(m+1))`. At most `k` parts are removed.

Given an output partition, its preimage is determined by the multiset
`tau` of removed parts. Remove `weight(tau)` ones from the output and
restore `tau`. Thus every fiber has size at most `B(n,k)`, the number
of partitions of weight at most `n` with at most `k` parts.

Conjugating partition diagrams identifies this count with partitions of
weight at most `n` whose largest part is at most `k`.
Suppose `k>=1` and `k^2<=n`. With `u=k/n`, positivity gives

```math
B(n,k)\le e^{un}\prod_{j=1}^k(1-e^{-uj})^{-1}.
```

For `0<=v<=1`, `1-e^(-v)>=v/2`.
Here `uj<=k^2/n<=1`, so

```math
B(n,k)
\le e^k\frac{(2n/k)^k}{k!}
\le\left(\frac{2e^2n}{k^2}\right)^k.
\tag{4}
```

The last inequality uses `ln(k!)>=k ln k-k`, which follows by
integrating `ln x` from one to `k`.
It follows that

```math
\ln P_m(n)
\ge\ln p(n)-k\ln\left(\frac{2e^2n}{k^2}\right).
\tag{5}
```

When `k=0`, every partition already has bounded parts and the assertion
is immediate without (4).

Apply (5) with `m=floor(C sqrt h)` and `n=2h-m`.
For fixed `C>=4` and increasing `h`,

```math
\frac n h\longrightarrow2,\qquad
\frac k{\sqrt h}\longrightarrow\frac2C,\qquad
\frac n{k^2}\longrightarrow\frac{C^2}{2}.
```

In particular, `k^2<=n` eventually. Equations (2) and (5) give

```math
\liminf_{h\to\infty}
\frac{\ln P_m(2h-m)}{\sqrt h}
\ge c\sqrt2-\frac2C\ln(e^2C^2).
\tag{6}
```

The final term tends to zero as `C` tends to infinity.

## 3. Complete each bounded partition with reserved ones

Take a partition `theta` of `2h-m` with all parts at most `m`.
Append exactly `m` ones, giving a partition `pi` of `2h`.
This construction is injective: removing exactly `m` ones recovers `theta`.

The partition `pi` splits into two partitions of weight `h`, each with
largest part at most `m`. To construct the split, take an initial segment
of the parts of `theta` with maximal sum at most `h`.
If that sum is less than `h`, the gap is smaller than the next part,
which is at most `m`. The reserved ones fill the gap.
The complementary multiset has weight `h` as well.

There is a next part unless the sum already equals `h`.
Indeed, `weight(theta)=2h-m>=h` for the range under consideration.
This also covers the endpoint `m=h`.

By Section 1, both halves belong to `S_lambda` for every sufficiently
large `h=hOf(lambda)`. Hence at least `P_m(2h-m)` distinct partitions
of `2h` occur as unions of two actual binary labels.

## 4. Independent leading monomials give the rank

The [channel theorem](top-channel-rank.md) proves the triangular identity

```math
(t^\xi)[h]t^\eta
=t^{\xi\cup\eta}
+\text{terms with fewer parts}.
\tag{7}
```

The creation-only term has coefficient one. Any annihilation lowers
the final number of parts by two relative to that term.
Choose one split for each union supplied by Section 3.
The corresponding images have different leading monomials and are
linearly independent, by descending induction on the number of parts.

Thus `rank(U_lambda)>=P_m(2h-m)` for each fixed `C` and sufficiently
large `lambda`. Equation (6), followed by `C -> infinity`, gives

```math
\liminf_{\lambda\to\infty}
\frac{\ln\operatorname{rank}U_\lambda}{\sqrt h}
\ge c\sqrt2.
```

The upper bound `rank(U_lambda)<=dim V_(2h)=p(2h)` and (2) give equality.
Finally, minimality of `h` and (2) imply
`lambda=(c/ln 2)sqrt h+O(ln h)`.
This proves (1).

The highest row channel of the source's central piece differs from
`U_lambda` only by nonzero signs on its input rows.
The rank is unchanged. It therefore needs
`sqrt(2) lambda+o(lambda)` bits within this factor map.
This does not prove that a factorization of the complete piece needs
the same number of bits.

Projecting the total row factor map onto its `i=2h` coordinates gives
this highest channel. The total factor map therefore has at least
the same rank. Its product with the other factor can still have
smaller rank.

## 5. A quantitative error term

The same argument gives (1a) without successive limits.
For sufficiently large `h`, choose

```math
C(h)=\frac{\ln h}{8a},\qquad
m=\left\lfloor\frac{\sqrt h\ln h}{8a}\right\rfloor,
\qquad t=\frac a{\sqrt h}.
```

Here `2tm<=ln(h)/4`. The omitted Euler-product block in (3) is
therefore at least

```math
m e^{-2tm}
\ge\frac{h^{1/4}\ln h}{16a}.
```

This dominates the `O(ln h)` error in the source's lower bound for
`ln p(h-1)`. Thus the bounded partitions still lie in the actual prefix.

Put `n=2h-m` and `k=floor(n/(m+1))` as before. Then

```math
k=\Theta\!\left(\frac{\sqrt h}{\ln h}\right),\qquad
\ln(n/k^2)=O(\ln\ln h).
```

Equation (5) loses only `O(sqrt(h) ln ln(h)/ln(h))`.
Also `sqrt(2h)-sqrt(n)=O(ln h)`, so shifting the partition weight
changes its logarithmic count by only `O(ln h)`.
The reserved-one construction and triangular argument therefore give

```math
\ln\operatorname{rank}U_\lambda
\ge c\sqrt{2h}
-O\!\left(\frac{\sqrt h\ln\ln h}{\ln h}\right).
```

Use `lambda=(c/ln 2)sqrt h+O(ln h)` to obtain (1a).
The upper bound follows from `rank(U_lambda)<=p(2h)` and (2).

## 6. A logarithmic error using the full partition asymptotic

The ordinary Hardy--Ramanujan theorem gives

```math
p(n)\sim\frac{e^{2a\sqrt n}}{4\sqrt3\,n},
\qquad a=\frac{\pi}{\sqrt6}.
\tag{8}
```

This is equation (1.41) of G. H. Hardy and S. Ramanujan,
[`Asymptotic formulae in combinatory analysis`, 1918](https://doi.org/10.1112/plms/s2-17.1.75);
the [original text is also available online](https://ramanujan.sirinudi.org/Volumes/published/ram36.html).
This section invokes (8) as an established ordinary theorem.
The pinned Lean source's coarse estimates alone were sufficient for
Sections 1--5, but do not supply (8).

### A uniform upper bound for large parts

Let `T_m(n)=p(n)-P_m(n)`. Removing one largest part gives

```math
T_m(n)=\sum_{j=m+1}^nP_j(n-j)
\le\sum_{j=m+1}^np(n-j).
```

For `j<=n/2`, formula (8) gives
`p(n-j)/p(n)<=C exp(-a j/sqrt(n))` for an absolute constant `C`
and all sufficiently large `n`. This follows from
`n/(n-j)<=2` and
`sqrt(n)-sqrt(n-j)>=j/(2sqrt(n))`.
The asymptotic comparison is uniform because every argument is at least
`n/2`. For `j>n/2`, monotonicity bounds every term by
`p(floor(n/2))`. Summing the first range as a geometric progression gives

```math
\frac{T_m(n)}{p(n)}
\le C_1\sqrt n\,e^{-am/\sqrt n}
+     C_2 n e^{-\beta\sqrt n},
\qquad
\beta=2a(1-1/\sqrt2)>0,
\tag{9}
```

for `1<=m<=n/2` and all sufficiently large `n`.
All constants are independent of `m`.

### The binary prefix omits only a small fraction

Put `D(n)=p(n)-p(n-1)`. This counts partitions without ones.
For `n>=2` it is nondecreasing: increase one largest part by one.
The new largest part is unique, so decreasing it recovers the input.
This is an injection between successive sets of partitions without ones.

Take `d=floor(sqrt(n))`. Telescoping and monotonicity give

```math
(d+1)D(n)\le\sum_{j=0}^d D(n+j)
\le p(n+d).
```

Formula (8) implies `p(n+d)/p(n)=O(1)`, and hence

```math
D(n)=O\!\left(\frac{p(n)}{\sqrt n}\right).
\tag{10}
```

This argument does not differentiate (8) or assume an error rate for
its relative asymptotic.

### A cutoff suitable at both weights

Choose `gamma=3/4` and

```math
m=\left\lfloor\frac{\gamma}{a}\sqrt h\ln h\right\rfloor,
\qquad n=2h-m.
```

For `m<j<=m+floor(sqrt(h))`, the largest-part classes are disjoint.
Equation (9), applied to weight `h-j` and cutoff `j`, gives

```math
P_j(h-j)=(1-o(1))p(h-j)
```

uniformly in this interval. Indeed, its relative upper-tail bound is
`O(h^(1/2-gamma+o(1)))=o(1)`. Formula (8) also gives the uniform ratio

```math
\frac{p(h-j)}{p(h)}
=\exp\!\left(-\frac{aj}{\sqrt h}+o(1)\right)
=\Theta(h^{-\gamma}).
```

The square-root expansion has error
`O(j^2/h^(3/2))=O((ln h)^2/sqrt(h))=o(1)`.
Adding these disjoint largest-part classes yields

```math
T_m(h)\ge c_0 h^{1/2-\gamma}p(h)
=c_0 h^{-1/4}p(h)
```

for some `c_0>0` and all sufficiently large `h`.
This exceeds `D(h)` by (10), so
`P_m(h)<p(h-1)<2^lambda`.
The bounded partitions again belong to the actual binary prefix.

At weight `n=2h-m`, equation (9) instead gives

```math
\frac{T_m(n)}{p(n)}
=O\!\left(h^{1/2-\gamma/\sqrt2+o(1)}\right)+o(1)
=o(1),
```

because `gamma>1/sqrt(2)`. Consequently,
`P_m(2h-m)=(1-o(1))p(2h-m)`.
Sections 3 and 4 therefore imply

```math
\log_2\operatorname{rank}U_\lambda
\ge \log_2 p(2h-m)+o(1)
=\kappa\sqrt{2h}-O(\log h).
```

Here `sqrt(2h)-sqrt(2h-m)=O(log h)`.
The opposite bound follows from `rank(U_lambda)<=p(2h)`.
The label calibration used in Section 4 already gives
`log2 rank(U_lambda)=sqrt(2) lambda+O(log lambda)`.

### The precise logarithmic coefficients

To obtain (1b), use a slightly smaller cutoff with a fixed additive shift:

```math
m=\left\lfloor
\frac{\sqrt h\ln h}{a\sqrt2}+B\sqrt h
\right\rfloor,\qquad n=2h-m.
```

Here `B>0` is fixed and sufficiently large. At weight `h`, the preceding
disjoint-class proof now has `gamma=1/sqrt(2)`. The shift changes its
positive constant, and it still gives

```math
T_m(h)=\Omega\!\left(
h^{\,1/2-1/\sqrt2}p(h)
\right)\gg D(h).
```

Thus the bounded partitions still lie in the actual prefix.
At weight `n=2h-m`, the first term on the right of (9) tends to
`C_1 sqrt(2) exp(-aB/sqrt(2))`. Choose `B` to make this limit smaller
than `1/4`. The second term tends to zero. Hence
`P_m(n)>=p(n)/2` for all sufficiently large `h`.
A positive retained fraction suffices; it need not tend to one.

Expanding (8) at this cutoff gives

```math
\ln p(2h-m)
=c\sqrt{2h}-\frac32\ln h-O(1),
\qquad c=2a.
```

Indeed, the shift in the square-root term is
`a m/sqrt(2h)=(ln h)/2+aB/sqrt(2)+o(1)`,
and the prefactor contributes `-ln h-O(1)`.
The reserved-one construction therefore gives this lower bound for
`ln rank(U_lambda)`. The ambient dimension gives the upper bound

```math
\ln\operatorname{rank}U_\lambda
\le\ln p(2h)
=c\sqrt{2h}-\ln h+O(1).
```

Minimality gives `p(h-1)<2^lambda<=p(h)`. Formula (8) implies
`p(h-1)/p(h)->1`, so the more precise label calibration is

```math
\lambda\ln2
=c\sqrt h-\ln h-\ln(4\sqrt3)+o(1),
\qquad
\ln h=2\ln\lambda+O(1).
\tag{11}
```

Combining these expressions gives (1b). The coefficient
`3-2sqrt(2)` uses this shifted threshold cutoff. The simpler cutoff
`gamma=3/4` above establishes the logarithmic error but has a larger
lower-bound coefficient.

## 7. Lossless compression of the existing factor vectors

The source's two factorizations, displays (54) and (55), use the vector
families

```math
F_r(\xi,\eta)
=\bigl((a^\xi)^{\langle i\rangle}[i-h]t^\eta\bigr)_{i=0}^{2h}
\in\bigoplus_{i=0}^{2h}V_i,
```

and

```math
F_c(\nu,\eta)
=\bigl((b^\nu)[h-i]t^\eta\bigr)_{i=0}^{2h}
\in\bigoplus_{i=0}^{2h}V_{2h-i}
```

for the central piece `s=0`. View each family as a linear map from
the vector space with its label pairs as a basis. The alternating signs
in the factorizations do not affect the following rank statements.
These are families in two different factorizations, not the two factors
of one matrix product.

The projection of `F_r` onto `i=2h` is `U_lambda`, up to the nonzero
input sign `(-1)^(h-length(xi))`. Indeed, `angle h (2h)` kills all
homogeneous parts below weight `h`, and leaves the top part of
`a^xi` with coefficient one. For `F_c`, the projection onto `i=0` is

```math
(b^\nu)[h]t^\eta
=\frac{h!}{z_\nu}(t^\nu)[h]t^\eta.
```

Every factor `h!/z_nu` is nonzero. Thus both projections have rank
`rank(U_lambda)`, and

```math
\operatorname{rank}F_r\ge\operatorname{rank}U_\lambda,
\qquad
\operatorname{rank}F_c\ge\operatorname{rank}U_\lambda.
\tag{12}
```

Here lossless linear compression of `F_r` means linear maps `L_r` and
`D_r` such that `D_r L_r F_r=F_r`, with the compressed vectors in
`Q^(M_r)`. Then `M_r>=rank(F_r)`. Define `M_c` in the same way for
`F_c`. If `w_r=ceil(log2 M_r)` and `w_c=ceil(log2 M_c)` are the
corresponding channel-index widths, (1b) implies

```math
w_r+w_c\ge 2\sqrt2\,\lambda-O(\log\lambda).
\tag{13}
```

The original ambient dimensions are both `sum_(j=0)^(2h) p(j)`, whose
logarithm is `sqrt(2) lambda+O(ln lambda)`. Therefore independent,
lossless linear compression of these existing vector families cannot
improve their leading width sum.

This conclusion specifically requires preserving each factor-vector
family. A quotient that discards coordinates invisible to the opposite
factor need not be lossless in this sense. A joint refactorization can
also change both vector families. Neither operation is excluded by
(12) or (13), and the rank of the complete central piece remains open.

The coarser error in (1a) gives the same leading conclusion without
invoking the full asymptotic (8).

## 8. Reproduction and limits

Run the finite implementation checks with

```sh
python3 tools/check_binary_prefix_channel.py
```

The checker verifies the replacement-map fibers, reserved-one splits,
lexicographic inclusion, largest-part classes, the injection for partitions
without ones, triangular terms, and small restricted matrix ranks.
Its recorded output is [binary-prefix-channel.json](../evidence/binary-prefix-channel.json).
It does not estimate the asymptotic threshold in Section 1.
The proof supplies the all-length asymptotic statement.

The source partition estimates are in
[`ClauseBound/Part/Bounds.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Part/Bounds.lean).
The operator construction and its source attribution remain unchanged.
This extension and checker used substantial OpenAI assistant support.
