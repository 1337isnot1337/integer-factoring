# Coordinate cores with a filled low-mode background

## Status

This note constructs rational local tensors whose restriction to matching
physical word types on each crossing pair is polynomial and has the exact
crossing endpoint. The
construction keeps the low mode occupied in a distinguished core state.
It applies to arbitrary finite mode sets, including sets that are not
downward closed. It provides an admissible family for the two-place Smith
pressure method. No improved factoring exponent is claimed.

The two separately constructed search inputs with mode sets
\(\{0,1,2\}\) and \(\{0,1,2,3\}\), both with ambient parameter \(M=4\),
agree coefficient by coefficient with the checker accompanying this note.
The local tensors have genuine poles at the crossing point. Their validity
depends on the paired-type restriction proved below.

## 1. The rational tensor and its coordinate maps

Work in characteristic two. Let \(S\subset[0,M)\) be a finite set containing
\(0\), let \(T=\max S>0\), and let \(M\) be a power of two. Write
\(c=1+w\), \(z=w^2\), and \(q=c^2\). The one-particle matrix of the
one-species Gaussian, with rows \((v,y)\) and columns \((u,x)\), is

\[
D=\begin{pmatrix}A&z^M P\\ I&B\end{pmatrix},\qquad
P_{ab}=[b\subseteq a],\quad
A_{ab}=c q^a[a\subseteq b],\quad
B_{ab}=c z^{M-1-(a\mathbin\vert b)}.
\tag{1}
\]

All indices lie in \(S\). The notation \(b\subseteq a\) means binary
inclusion. The scalar tensor \(G\) consists of all minors of \(D\), with
physical occupancy labels \((x,y,u,v)\). No closure property of \(S\) is
needed for (1).

Introduce the exterior monomial

\[
Q=x_0y_Tu_0v_T
\]

and restrict \(G(1+Q/c^2)\) by the following coordinate maps on the two
extreme modes. All other modes remain free.

| Logical core letter | Input \(x\), output \(v\) | Input \(u\), output \(y\) | Index height |
|---|---|---|---:|
| \(0\) | low mode only | low mode only | \(0\) |
| \(1\) | vacuum | vacuum | \(0\) |
| \(2\) | both extreme modes | top mode only | \(T\) |

Thus the physical alphabet has size

\[
d=3\cdot 2^{|S|-2}.
\tag{2}
\]

For a full physical letter \(a\), define its height by

\[
h(a)=\sum_{j\in S\setminus\{0,T\}}j\,n_j(a)
       +T[\text{the core letter is }2].
\tag{3}
\]

This is the actual occupied-index sum on every one of the four legs,
despite the different coordinate maps in the table. Let \(K\) denote the
resulting rational local tensor. Its explicit polynomial numerator is

\[
N=\operatorname{restrict}(c^2G+QG),\qquad K=N/c^2.
\tag{4}
\]

The checker implements (4) by a cached minors recursion. It also constructs
the same numerator independently by multiplying the literal exterior
factors of (1).

## 2. A regular ramified tensor

Make the substitution

\[
w=1+t^T,\qquad c=t^T,
\]

and multiply each coefficient by
\(t^{\Delta h}\), where
\(\Delta h=h(y)+h(v)-h(x)-h(u)\).
The resulting tensor is polynomial in \(t\).

To see this before the coordinate restriction, attach the index gauge
\(t^{a-b}\) to an edge from input index \(b\) to output index \(a\).
The \(A\)-edge has exponent

\[
T(2a+1)+a-b\ge0
\]

on its support \(a\subseteq b\). Equality occurs only for \(a=0,b=T\).
The \(B\)-edge has exponent \(T+a-b\ge0\), again with equality only for
\(a=0,b=T\). Its remaining factor
\(w^{2(M-1-(a\mathbin\vert b))}\) is a polynomial with constant term one.
The \(x\)-to-\(v\) edge has exponent \(a-b\ge0\) on the support of
\(P_{ab}\); equality means \(a=b\). The \(u\)-to-\(y\) edges are the
identity. Finally, the index increase of \(Q\) is \(2T\), which exactly
cancels its denominator \(c^2=t^{2T}\).

Consequently the constant term of the full ramified tensor is the
crossing on all modes, multiplied in the exterior algebra by

\[
(1+u_Tv_0)(1+x_Ty_0)(1+x_0y_Tu_0v_T).
\tag{5}
\]

The bulk modes are untouched. Applying the coordinate maps to the
extreme modes gives a nine-dimensional transfer matrix

\[
F=I+N_0,\qquad
N_0=E_{0,5}+E_{3,7},\qquad N_0^2=0.
\tag{6}
\]

Here rows are indexed by \(v+3y\) and columns by \(x+3u\). In explicit
core pairs, the two extra transitions are

\[
(x,u)=(2,1)\longmapsto(v,y)=(0,0),\qquad
(x,u)=(1,2)\longmapsto(v,y)=(0,1).
\tag{7}
\]

Each strictly lowers the total height by \(T\). All diagonal transitions
preserve height. Equation (6) can be verified directly from the seven
exterior factors in (5) and the crossing.

## 3. Independent types on the two crossing pairs

Fix a word length \(m\) and arbitrary physical word types
\(\pi_X,\pi_U\). Restrict the \(x,v\) boundary words of
\(K^{\otimes m}\) to \(\pi_X\), and restrict its \(u,y\) boundary
words to \(\pi_U\). The two types may differ.
Every coefficient of this restriction is a polynomial in \(w\), and at
\(w=1\) the tensor is the exact rectangular crossing on the two type classes.

Indeed, the total height \(\sum h(a)\) depends only on the word type.
The \(x,v\) heights cancel as one pair, and the \(u,y\) heights cancel
as the other. Hence the product of the ramified coefficient gauges is exactly one on
every restricted entry. The rational restricted entry, after substitution
\(w=1+t^T\), is therefore regular at \(t=0\) by Section 2. Its only
possible finite pole in the original variable is \(w=1\), because (4)
has denominator \(c^{2m}\). Ramification multiplies its valuation at this
point by \(T\). Regularity after ramification proves that the original
entry has no pole, so it is a polynomial in \(w\).

Its value at \(w=1\) is computed by the constant tensor (5). A term using
any transition in (7) strictly lowers total height. Since all other
constant transitions preserve height, such a term cannot survive the
paired-type restriction. Only the diagonal part of \(F^{\otimes m}\)
remains, together with the bulk crossing. This proves the stated
endpoint.

The proof uses actual index heights and therefore works for every finite
set \(S\) in Section 1. Neither \(P^2=I\) nor a down-set identity is used.
The polynomial degree of a restricted entry is at most \(m\deg N\), so
the number of coefficient tensors grows only linearly with \(m\) for a
fixed seed.

If equal physical dimensions are required, take a second copy with the
crossing pairs exchanged: swap \(x,u\) and \(y,v\). Its type-class
dimensions are exchanged, so their tensor product has the same dimension
\(D_XD_U\) on all four boundaries. The swap preserves each of the two
flattening ranks, up to permutations and transposition. Keep the two
coefficient indices separate when taking tensor products of the pieces.
Their coefficient-cost sums then multiply, while the external dimensions
multiply to \(D_XD_U\). The special case
\(\pi_X=\pi_U\) recovers the original common-four-leg-type statement.

## 4. Local pole control and Smith normalization

Every entry of \(N\) is divisible by \(c\). In fact the only possible
local poles occur when the \(x\)-core letter is \(0\).

The Gaussian contribution \(c^2G\) is already divisible by \(c^2\).
Consider a surviving \(QG\) contribution. Its core letters must satisfy

\[
x\in\{0,2\},\qquad u=0,\qquad y=v=2.
\]

After removing the four variables of \(Q\), the remaining minor contains
row \(v_0\) and no column \(x_0\). Every remaining \(x\)-column has zero
in row \(v_0\), since \(P_{0j}=0\) for \(j>0\). Every remaining
\(u\)-column has a factor \(c\) in that row. Expanding along it proves
that the minor is divisible by \(c\).

If the \(x\)-core letter is \(2\), column \(x_T\) also remains. The
\(v_T\) row was removed, so this column can only use a \(y\)-row, where
its entry has another factor \(c\). Thus these contributions are
divisible by \(c^2\). With no bulk modes, the correction vanishes
identically: after removing \(Q\), the only potential nonempty minor is
the zero entry from \(x_T\) to \(v_0\).

One may therefore use either of two exactly equivalent local inputs:

\[
K=N/c^2=(N/c)/c.
\tag{8}
\]

For two-place Smith profiles, division by \(c^2\) leaves the zero orders
unchanged and adds \(2\) to every infinity order of \(N\). If the
common factor is first canceled, division by \(c\) instead adds \(1\)
to every infinity order of \(N/c\). These are identical rational
profiles. The same restoration must be included when optimizing
physical gauges or word types.

The local pole is real in the checked bulk examples: the \(d=6\) seeds
have two entries with a simple pole, and the full four-mode \(d=12\)
seed has sixteen. Treating their local numerator as a polynomial
crossing without its denominator would therefore be invalid. Section 3
is the reason the rational profiles can be transferred to polynomial
coefficient tensors after the physical restriction.

## 5. Reproducible checks

Run from the repository root:

```sh
python tools/check_coordinate_core_crossings.py
```

The recorded [evidence](../evidence/coordinate-core-crossings.json) checks
the mode sets \(01,012,013,023,0123\), and \(025\), with valid ambient
parameters. It verifies all of the following exactly over \(\mathbf F_2\):

- The residual matrix (6), including \(N_0^2=0\).
- Equality of the minors constructor and a literal exterior expansion.
- Every local ramified coefficient is polynomial, with the predicted
  constant term.
- Every entry of every pair of independent two-copy full-letter types is divisible by
  the complete denominator and has the crossing endpoint.
- The common \(c\)-factor, the stronger \(x\)-core-2 divisibility, and
  the vanishing no-bulk correction.
An additional development replay compared the independently generated
search numerators: all 1,296 entries for \(S=012\) and all 20,736 entries
for \(S=0123\). This optional comparison accepts an explicit input path
through `--compare-seed`; the default evidence omits external search-file
comparisons. The `--seed-directory` option writes the
complete reproducible numerator seeds, with denominator degree and
logical-height metadata. The finite checks support the ordinary proof;
they make no numerical pressure or factoring-exponent claim.
