# Polynomial crossings from an asymmetric involutive core

An asymmetric three-letter restriction of the residual quartic tensor has
endpoint \(I+N\), where \(N^2=0\). This gives three constructive routes:
a polynomial tensor with a typed crossing endpoint in the original parameter,
a local output correction, and a two-layer construction with a strict local
crossing. The local alphabet is \(3\cdot2^{n-2}\) for \(n\) source modes.

The constructions and their endpoint proofs are exact. Their coefficient
pressures remain an optimization problem; no improved factoring exponent is
claimed here.

## 1. The asymmetric three-letter core

Use the residual tensor and mode conventions from
[`residual-quartic-codes.md`](residual-quartic-codes.md). The extreme states are
ordered as vacuum, low, top, both. On the two crossing pairs use these maps:

| Leg | Logical 0 | Logical 1 | Logical 2 |
|---|---|---|---|
| \(x\), input | vacuum | low | both |
| \(u\), input | vacuum | top | low + both |
| \(y\), output | vacuum dual | top dual | both dual |
| \(v\), output | vacuum dual | low dual | top dual + both dual |

The untwisted transfer has rows \((y,v)\), columns \((u,x)\), and pair index
\(v+3y\) or \(x+3u\). Direct expansion gives

\[
 F=I_9+N,\qquad N=E_{1,3}+E_{7,5},\qquad N^2=0,
 \qquad F^2=I_9. \tag{1}
\]

In explicit output pairs, the two nonidentity transitions are

\[
 (y=0,v=1)\leftarrow(u=1,x=0),\qquad
 (y=2,v=1)\leftarrow(u=1,x=2). \tag{2}
\]

The maps differ on the two crossing pairs. They therefore avoid the
same-code singularity obstruction in
[`residual-quartic-three-letter-obstruction.md`](residual-quartic-three-letter-obstruction.md).
The transfer in (1) is not itself a three-letter crossing, so it does not
contradict the general constant-code obstruction either.

## 2. A polynomial tensor in the original variable

Let \(M\) be a power of two and
\(S=\{0=a_0<\cdots<a_{n-1}=T\}\subset[0,M)\). No down-set condition is needed.
Put \(c=1+w\), and use the Gaussian one-particle matrix

\[
 D_S(w)=\begin{pmatrix}cQP^T&w^{2M}P\\I&cH\end{pmatrix},\quad
 Q_{aa}=c^{2a},\quad P_{ab}=[b\subseteq a],\quad
 H_{ab}=w^{2(M-1-(a\mathbin{|}b))}.
\]

Apply the core maps to

\[
 \Gamma(D_S(w))\left(1+\frac{x_0y_Tu_0v_T}{c^2}\right), \tag{3}
\]

but replace the \(u\)-input logical-two vector by

\[
 c\,|\{0\}\rangle+|\{0,T\}\rangle. \tag{4}
\]

Keep every bulk mode \(a_1,\ldots,a_{n-2}\) independently on all four legs.
The resulting logical alphabet has size \(q=3\cdot2^{n-2}\). A label is encoded
as its bulk bit mask plus \(2^{n-2}\) times its core letter.

### Proposition 1: entrywise polynomiality

Every coded entry of (3) with (4) belongs to \(\mathbb F_2[w]\).

**Proof.** The Gaussian part is polynomial. A quartic term forces the \(u\)
and \(v\) core letters to be two, the \(x\) core letter to be one or two,
and the \(y\) core letter to be one or two. After the quartic removes the four
extreme variables, only the optional \(u_T\) and \(v_0\) variables need be
summed over the two superpositions.

Their Gaussian mutual edge has coefficient \(c\). Write their remaining
linear couplings as

\[
 L_U=c\sum_{\substack{a\in S\\0<a<T\\a\subseteq T}}c^{2a}v_a,
 \qquad L_V=c\sum_{\substack{a\in S\\0<a<T}}u_a. \tag{5}
\]

The remaining Gaussian factor is independent of \(u_T,v_0\). The four
possibilities from (4) and the \(v\) projection contribute

\[
 c+L_U+cL_V+(c+L_UL_V)=L_U+cL_V+L_UL_V. \tag{6}
\]

Here \(c^3\mid L_U\), since every summation index in (5) is positive, and
\(c\mid L_V\). Thus (6) is divisible by \(c^2\). If the \(x\) core letter is
two or the \(y\) core letter is two, the extra variables \(x_T,y_0\) remain in
the independent Gaussian factor. They do not couple to the eliminated pair:
the relevant Pascal entry \(P_{0T}\) and identity entry \(I_{0T}\) are zero.
The same cancellation therefore applies in all cases. \(\square\)

## 3. Its endpoint is a crossing on every common word type

Let \(e(b)\) denote the sum of the indices of the occupied bulk modes. Define
four logical height functions by

\[
 \begin{aligned}
 h_x(b,j)=h_v(b,j)&=e(b)+(0,0,T)_j,\\
 h_u(b,j)=h_y(b,j)&=e(b)+(0,T,T)_j.
 \end{aligned} \tag{7}
\]

### Proposition 2: typed endpoint

Every nonzero endpoint entry has

\[
 \Delta=h_y(y)+h_v(v)-h_x(x)-h_u(u)\ge0.
\]

When \(\Delta=0\), the entry is exactly \([x=v][y=u]\). Therefore every
tensor power restricted to one common physical word type has the crossing
endpoint.

**Proof.** As in the earlier residual construction, set \(w=1+t^T\) and
apply the physical index-energy gauge before applying the constant core maps.
The resulting tensor \(\widetilde K(t)\) is polynomial. Its endpoint is the
core transfer \(F\) from (1), tensored with the bulk crossing. The original
polynomial tensor \(K\) satisfies

\[
 \widetilde K(t)_{xyuv}=t^{\Delta}K(1+t^T)_{xyuv}. \tag{8}
\]

Both extra transitions in (2) have height difference \(-T\). Thus the
zero-height part of \(F\) is exactly the crossing. Polynomiality of both sides
of (8) excludes a nonzero original endpoint entry with \(\Delta<0\), and the
zero-height entries equal the crossing.

On a common word type the local height differences sum to zero, since
\(h_x=h_v\) and \(h_u=h_y\) pointwise. Every nonzero local endpoint factor
has nonnegative difference, so all of them have difference zero. Their
product is exactly the word crossing. \(\square\)

The equality of the paired height functions is sufficient; the four functions
need not all be the same. This permits arbitrary common full-letter types in
the pressure search.

## 4. Monomial code parameters and output repair

The construction permits three integer parameters \(A,D,K\). Replace (4) by
\(c w^A|\{0\}\rangle+|\{0,T\}\rangle\), replace the \(v\)-output logical-two
row by \(\langle\{T\}|+w^D\langle\{0,T\}|\), and replace the quartic
coefficient by \(w^K/c^2\).

The numerator replacing (6) is

\[
 c(w^A+w^D)+L_U+c w^{A+D}L_V+w^D L_UL_V. \tag{9}
\]

It is again divisible by \(c^2\) in \(\mathbb F_2[w,w^{-1}]\):
\(w^A+w^D\) vanishes at \(w=1\), and the other terms have the same divisibility
as before. A common power of \(w\) clears all remaining Laurent poles.
Every monomial equals one at the endpoint, so Proposition 2 remains valid.

A fourth optional parameter \(k\in\mathbb Z\) acts on the output core pair by

\[
 I+c w^kN. \tag{10}
\]

Explicitly, retain the bulk bits and add \(c w^k\) times output pair
\((y=1,v=0)\) to \((y=0,v=1)\), and output pair \((y=1,v=2)\) to
\((y=2,v=1)\). Since \(N^2=0\), (10) is its own inverse. It preserves the
polynomial construction after a possible common \(w\)-clearing shift and
does not change its original-variable endpoint.

At infinity, \(c w^k\) is regular for \(k\le-1\). Thus, when the physical
gauge is constant on the core letters, (10) is an invertible regular column
operation on the \((x,u)\mid(y,v)\) flattening and preserves that Smith
profile. Arbitrary physical gauges require checking the transformed gate's
regularity, or computing the profile directly. The other cut can change.

## 5. Strict local endpoints from output correction or temporal composition

For the polynomial energy-gauged tensor \(\widetilde K\), a constant output
correction by \(F^{-1}=F\) gives a strict local crossing. In the original
variable this correction is exactly \(I+cN\), because both transitions of
\(N\) lower the output height by \(T\).

Alternatively, compose two copies of \(\widetilde K\) only on their core
wires. Sum over the first copy's output core letters \((i,j)\), feeding
\(i\) into the second copy's \(u\) core and \(j\) into its \(x\) core. Keep all
bulk occupations of both copies as separate external physical labels. At the
endpoint the core transfers multiply to \(F^2=I\), while both bulk endpoints
are crossings. The grouped tensor therefore has a strict local crossing on

\[
 3\cdot 2^{2(n-2)}
\]

letters. Its coefficients are polynomial by construction. This finite
composition does not by itself establish a favorable coefficient-rank cost;
that cost must be measured on the resulting scalar tensor.

## Reproduction and finite controls

```sh
python tools/check_involutive_core_crossings.py \
  --output evidence/involutive-core-crossings.json \
  --seed-output evidence/involutive-core-six-letter.json
```

The six-letter original-variable candidate uses \(S=\{0,1,2\}\), \(M=4\),
and has degree 24. The one-copy twelve-letter candidate on the full interval
\(S=\{0,1,2,3\}\), \(M=4\), has degree 32. The checker compares every
coefficient of both candidates with an independent exterior-factor expansion,
checks every local endpoint entry, and checks every two-copy common full-letter
type. It also tests noninterval modes, signed monomial exponents, output
repairs, and the strict twelve-letter temporal construction.

The reusable constructor accepts `--input-power`, `--output-power`,
`--quartic-power`, and `--repair-power`. For example:

```sh
python tools/involutive_core_crossing.py \
  --indices 0,1,2,3 --ambient 4 --repair-power -1 \
  --output /tmp/involutive-core-twelve-letter.json
```

All exports use polynomial bit masks at index \(x+qy+q^2u+q^3v\). The direct
mode explicitly reports a typed endpoint; `ramified-repair` and `temporal`
report a strict local endpoint.
