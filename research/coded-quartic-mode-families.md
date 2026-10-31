# Coded quartics on arbitrary mode sets

The polynomial construction in
[`residual-quartic-codes.md`](residual-quartic-codes.md) extends to arbitrary
finite mode sets. Its physical code can also vary with the parameter. Both
extensions preserve the exact typed crossing endpoint. They supply additional
families for optimization; no pressure improvement is claimed here.

## 1. Arbitrary modes

Let \(M\) be a power of two and let

\[
 S=\{0=a_0<a_1<\cdots<a_{n-1}=T\}\subset[0,M),\qquad n\ge2.
\]

No down-set condition is imposed. Index every matrix by \(S\), and use

\[
 D_S(w)=\begin{pmatrix}cQP^T&w^{2M}P\\I&cH\end{pmatrix},\quad
 c=1+w,\quad Q_{aa}=c^{2a},\quad
 P_{ab}=[b\subseteq a],\quad H_{ab}=w^{2(M-1-(a\mathbin{|}b))}.
\]

Apply the extreme-mode code from the preceding note to

\[
 \Gamma(D_S(w))\left(1+\frac{x_0y_Tu_0v_T}{c^2}\right).
 \tag{1}
\]

The logical alphabet has size \(q=2^{n-1}\). Its bulk bits correspond to
\(a_1,\ldots,a_{n-2}\), and its core bit encodes the two-letter code on \(0,T\).
Define the logical height

\[
 h(b)=\sum_{i=1}^{n-2}a_i[b_i=1]+T[b_{\rm core}=1]. \tag{2}
\]

### Proposition

The coded tensor is polynomial. Its endpoint has only nonnegative height
increase \(h(y)+h(v)-h(x)-h(u)\), and its zero-increase entries are exactly the
crossing. Consequently, every common fixed-height word restriction, including
every common physical word type, has the exact crossing endpoint.

**Proof.** After the quartic consumes both top output rows, a remaining top
input column is divisible by \(c\). For the \(u_T\) column the identity entry is
removed. For the \(x_T\) column, a surviving Pascal entry would require
\(T\subseteq a\) with \(a\in S\). Maximality of \(T\) forces \(a=T\), whose row
has been removed. Every other entry in these columns has an explicit factor
\(c\). Thus the same two-input argument cancels \(c^2\) termwise.

For the endpoint, put \(w=1+t^T\) and apply the index-energy gauge. The internal
entries acquire exponents

\[
 T(2a+1)+a-b\quad\hbox{and}\quad T+a-b.
\]

Both are nonnegative because \(0\le a,b\le T\). In each internal block the
only exponent-zero entry is \(a=0,b=T\); its remaining power of \(w\)
evaluates to one. The Pascal cross block has exponent \(a-b\), which is
nonnegative on its support and zero only on its diagonal. Therefore the
ramified endpoint is exactly the same residual extreme-mode tensor as before,
tensored with the bulk crossing. The same code corrects it, and the exact
relation

\[
 \widetilde K_S(t)_{xyuv}
 =t^{h(y)+h(v)-h(x)-h(u)}K_S(1+t^T)_{xyuv}
\]

proves the assertion. Neither \(T=M-1\) nor odd \(T\) is required. \(\square\)

The generic rank-loss identity for the unprojected full-interval quartic does
not automatically extend to arbitrary \(S\). The proposition establishes the
polynomial endpoint construction independently of that identity.

## 2. Parameter-dependent physical codes

Let \(h_x,h_u,g_y,g_v,k\in\mathbb F_2[w,w^{-1}]\) satisfy

\[
 h_x(1)=h_u(1)=g_y(1)=g_v(1)=k(1)=1. \tag{3}
\]

Use the following input columns and output rows on the two extreme modes:

| Leg | Logical zero | Logical one |
|---|---|---|
| \(x\), input | \(|\varnothing\rangle\) | \(c h_x(w)|\{0\}\rangle+|\{0,T\}\rangle\) |
| \(u\), input | \(|\varnothing\rangle\) | \(c h_u(w)|\{0\}\rangle+|\{0,T\}\rangle\) |
| \(y\), output | \(\langle\varnothing|\) | \(\langle\{T\}|+g_y(w)\langle\{0,T\}|\) |
| \(v\), output | \(\langle\varnothing|\) | \(\langle\{T\}|+g_v(w)\langle\{0,T\}|\) |

Replace the quartic coefficient in (1) by \(k(w)/c^2\). The resulting entries
are Laurent polynomials with poles only at \(w=0\). A common power of \(w\)
makes the entire tensor polynomial and leaves the endpoint unchanged.

Indeed, the termwise \(c^2\) cancellation is unchanged. All five additional
coefficients are regular at \(w=1\) and have value one, so the residual endpoint
and the relation with logical height remain valid. This proves the same typed
endpoint for every choice satisfying (3).

The four code coefficients may differ. This is useful when trying to align
different zero and infinity flags. For example, the inherited full-interval
four-mode Gaussian flags suggest zero orders

\[
 (\operatorname{ord}_0 h_x,\operatorname{ord}_0 h_u,
   \operatorname{ord}_0 g_y,\operatorname{ord}_0 g_v)=(-6,-7,-6,-7)
\]

and degrees \((6,5,7,6)\), respectively. These requirements are compatible with
(3). One can choose

\[
 \begin{aligned}
 h_x&=w^{-6}(1+w^j+w^{12}),&h_u&=w^{-7}(1+w^\ell+w^{12}),\\
 g_y&=w^{-6}(1+w^r+w^{13}),&g_v&=w^{-7}(1+w^s+w^{13}),
 \end{aligned}
\]

with \(1\le j,\ell<12\) and \(1\le r,s<13\). These are concrete search
parameters, not a claim that matching the individual flags improves the full
tensor pressure.

## 3. Reusable constructor and exact controls

The constructor exports polynomial bit masks with array index
\(x+qy+q^2u+q^3v\). It checks every local endpoint entry against (2), including
all possible additional positive-height entries. It records the common power
of \(w\) used to clear Laurent poles.

For the five-mode, sixteen-letter candidate:

```sh
python tools/build_coded_quartic_crossing.py \
  --indices 0,1,2,3,4 --ambient 8 \
  --output /tmp/coded-quartic-five-mode.json
```

This candidate has polynomial degree 80. The construction checks all
\(16^4=65{,}536\) local entries. Four-mode examples can use noninterval sets,
such as \(S=\{0,1,3,7\}\), without increasing the eight-letter alphabet.

Each optional coefficient argument is `POLYNOMIAL:SHIFT`, meaning
\(w^{\mathrm{SHIFT}}\) times the binary polynomial mask. For example, `7:-2`
denotes \(w^{-2}(1+w+w^2)\). The command rejects a coefficient whose value at
one is not one. Options `--hx`, `--hu`, `--gy`, `--gv`, and `--quartic` control
the five coefficients in (3).

The finite checks in [`coded-quartic-mode-families.json`](../evidence/coded-quartic-mode-families.json)
include interval and noninterval mode sets, the sixteen-letter candidate,
Laurent code changes, and nonconstant quartic numerators. These checks verify
the construction and its endpoint; numerical search results are separate.
