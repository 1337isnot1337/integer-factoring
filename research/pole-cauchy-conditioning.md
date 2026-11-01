# Balanced padding bounds the inverse precision of every pole rectangle

## Result

The [global Cauchy target](pole-cauchy-content.md) has potentially enormous
determinant content. This does not require enormous precision in its
normalized inverse. A different padding rule gives, for every factoring
rectangle, a square matrix whose every Smith order lies between zero and
`k-1`, where `k` is the displacement exponent in

\[
 M_k(a,b)=\frac{2^kw}{N+2^kw-ab},\qquad w^2+w+1=0.
 \tag{1}
\]

The bound is independent of the matrix dimension. In particular,
`2^(k-1) M_k^(-1)` is integral. Input entry precision
`p+2(k-1)` suffices for the inverse/Hadamard-square trace modulo `2^p`.
The displacement and this input precision both have numerical size
`O(n+p)` for `n`-bit factoring input.

There is also an effective entry oracle. Given one row and one column
by binary integers, it computes the exact valuation of the inverse
entry and its normalized unit modulo `2^p`, in bit cost polynomial in
the binary input length and `p`. It does not enumerate either domain.

The remaining global task is to sum or contract these entries at the
required cost. An inexpensive individual inverse entry does not make
that exponentially large contraction inexpensive. The complete
factoring score remains `1/2`.

Implementation:
[pole_cauchy_conditioning.py](../tools/pole_cauchy_conditioning.py).

## 1. Balanced residue sets and a cofactor identity

Let `A,B` each contain `m` distinct positive odd nodes. Choose `h>=1`
with `|N-ab|<2^h` for every pair and take `k>=h`. Assume `N` is odd.
For an odd residue `c modulo 2^ell`, put

\[
 r_{\ell,c}=\#\{a\in A:a\equiv c\pmod{2^\ell}\},\qquad
 s_{\ell,c}=\#\{b\in B:Nb^{-1}\equiv c\pmod{2^\ell}\},
 \quad \Delta_{\ell,c}=r_{\ell,c}-s_{\ell,c}.
 \tag{2}
\]

The [minor-energy identity](pole-smith-residue-flow.md) gives the exact
valuation of every minor as half the sum of its squared residue
imbalances. Apply it to the full determinant and the cofactor obtained
by deleting row `a` and column `b`. Let `d` be the largest level through
`k` at which `a` and `N/b` share a residue. Equivalently,

\[
 d=\min\bigl(k,v_2(N-ab)\bigr),
\]

where `v_2(0)=+infinity`. Since all nodes and `N` are odd, `d>=1`.
The two deletions cancel on their shared ancestor path. Below their
split, deleting the row changes an energy by `1/2-Delta`, while
deleting the column changes it by `1/2+Delta`. Therefore

\[
 v_2\bigl((M_k^{-1})_{b,a}\bigr)
 =\sum_{\ell=d+1}^k
   \left(1-\Delta_{\ell,a\bmod2^\ell}
             +\Delta_{\ell,Nb^{-1}\bmod2^\ell}\right).
 \tag{3}
\]

This is an exact cofactor calculation; it does not replace the inverse
by entrywise reciprocals. In particular, if every imbalance has
absolute value at most `epsilon`, then

\[
 v_2\bigl((M_k^{-1})_{b,a}\bigr)
   \ge -(2\epsilon-1)(k-d)
 \qquad(\epsilon\ge1).
 \tag{4}
\]

Call an odd node set of size `m` *balanced* if, at every level `ell`,
each odd residue class contains either
`floor(m/2^(ell-1))` or `ceil(m/2^(ell-1))` nodes. An ordinary block of
`m` consecutive odd integers is balanced. Multiplication by `N` after
inversion permutes the odd residue classes, so it preserves this
property. If both `A,B` are balanced, every imbalance in (2) is at most
one in absolute value. Equations (3)--(4) then give

\[
 \boxed{\quad 2^{k-1}M_k^{-1}\text{ is integral}.\quad}
 \tag{5}
\]

The matrix (1) is itself integral. Multiplying its Smith normal form
by the integral inverse bound in (5) shows that every Smith order is
in `[0,k-1]`. The bound can be attained: the evidence includes `N=29`,
rows `{3,5,7}`, columns `{5,7,33}`, and `k=12`, with orders `(9,9,11)`.

## 2. Balanced padding with binary input cost

Start from positive integer intervals and clip them to `[1,N+1)`.
For odd `N`, retain only their odd nodes. Empty domains have count zero.
Each nonempty domain is a progression

\[
 \{a_0,a_0+2,\ldots,a_0+2(r-1)\}.
\]

Set `m` to the larger of the two cardinalities. The longer domain is
already balanced. To pad a shorter domain, let `Q=2^K` be the least
power of two at least `2m`. Consider its missing continuation
`a_0+2r,...,a_0+2(m-1)`. Shift this entire tail by `tQ`, with

\[
 t=\max\left(0,
       \left\lfloor\frac{N-(a_0+2r)}Q\right\rfloor+1\right).
 \tag{6}
\]

All added nodes are greater than `N`, so no new positive factor pair
is created. The padded domain is the union of at most two disjoint
odd progressions.

At levels through `K`, its residue counts equal those of the ordinary
length-`m` continuation, hence are balanced. All `m` continuation
residues are distinct modulo `Q`: the difference between their indices
has absolute value less than `m<=Q/2`. Consequently the padded nodes
remain distinct at every higher level, where balanced counts are
exactly zero or one. This proves balancedness at every level.

All added nodes are at most `N+Q+2m<4N+3`, because `m<=(N+1)/2`.
The padding parameters have `O(n)` bits and are constructed with a
constant number of integer operations. The matrix dimension is still
potentially exponential, but neither the padding construction nor the
precision bound visits its nodes.

Combining this padding with (5) gives the same inverse bound for every
original rectangle. This improves the earlier arbitrary-tail padding
by preserving the required residue balance.

## 3. Precision for the normalized scalar is independent of dimension

Write `L=k-1`. Let `A=M_k+E` be any integral entrywise approximation with
every entry of `E` divisible by `2^P`, where

\[
 P=p+2L,\qquad p\ge1.
 \tag{7}
\]

Since `M_k^(-1)` loses at most `L` bits, `M_k^(-1)E` is divisible by
`2^(P-L)`, hence by two. The matrix `I+M_k^(-1)E` is therefore invertible
over the valuation ring. Thus `A` is invertible and

\[
 v_2\bigl(A^{-1}-M_k^{-1}\bigr)\ge P-2L=p,
 \tag{8}
\]

where a matrix valuation is the minimum entry valuation. This follows
directly by subtracting the two inverses; the inverse of
`I+M_k^(-1)E` is integral. Moreover the entrywise squares of `A` and
`M_k` differ by a matrix divisible by `2^P`. Matrix multiplication and
trace cannot lower the minimum valuation merely because there are many
summands. Consequently

\[
 \mathrm{tr}\bigl(A^{-1}(A\mathbin{\circ}A)\bigr)
 \equiv
 \mathrm{tr}\bigl(M_k^{-1}(M_k\mathbin{\circ}M_k)\bigr)
       \pmod{2^p}.
 \tag{9}
\]

For `k>=h+p`, the right side is the rectangle's factor count modulo
`2^p`, by the existing normalized trace identity. Taking `p` large
enough to exceed the binary length of the smaller original cardinality
recovers the count exactly. Here `h,p,k,P=O(n)`.

Equation (9) is a precision guarantee, not a free inversion operation.
One may choose ordinary integer representatives of both coefficients
in the quadratic basis `1,w` for each approximate entry, then use
exact linear algebra over `Q(w)`, with its numerical dimension cost
charged. A faster implicit algorithm must still account
for every global operation. No precision proportional to the enormous
determinant content is required by the input-conditioning problem.

## 4. The individual inverse entries are effective

Let `z=N+2^k w`, `t=2^k w`, and
`P_A'(a)=prod_(u in A,u!=a)(a-u)`, with the analogous notation for `B`.
The exact Cauchy cofactor formula gives

\[
 (M_k^{-1})_{b,a}
 =\frac{
       \prod_{u\in A}(z-ub)\,
       \prod_{v\in B}(z-av)}
      {t\,z^{m-1}(z-ab)P_A'(a)P_B'(b)}.
 \tag{10}
\]

The cofactor and two Vandermonde signs cancel. Every displayed product
has a fixed number of odd-progression factors. After indexing one such
progression, it is a product of affine terms `c+di+2^kw`, or an integer
affine product with its one zero term omitted.

The algorithm separates every product into a binary valuation exponent
and a unit. A true zero in the constant affine part contributes exactly
`2^k w`; otherwise group indices by the exact valuation of `c+di`.
There are only `O(n)` groups, found by linear congruences. Dividing each
group by its known power of two leaves a unit affine progression with
even step.

For unit `U` and even `V`, its normalized product modulo `2^p` is

\[
 \prod_{i=L}^{U_0-1}(U+iV)
 =U^{U_0-L}
   \sum_{j=0}^{\min(p-1,U_0-L)}e_j(L,\ldots,U_0-1)(VU^{-1})^j
       \pmod{2^p}.
 \tag{11}
\]

The [unit-jet algorithm](interval-unit-jets.md) computes the needed
elementary symmetric functions by exact power sums and Newton
identities. Its cost is polynomial in the binary endpoints and `p`;
the number of factors is not enumerated. The same valuation separation
is used in the [normalized interval jets](normalized-interval-jets.md).
Here it is applied at the constant quadratic unit `w`.

Multiply and divide the resulting units in the fixed quadratic ring,
combine the binary valuation exponents, and obtain (10). The unit
power `z^(m-1)` uses binary powering. This proves polynomial bit cost
for each requested inverse entry. Formula (3) independently computes
its valuation using `O(k)` residue-capacity queries, each a floor count
and modular inverse.

Returning an inverse entry in valuation-plus-unit form avoids allocating
its possibly large power of two. By (5), the scaled entry
`2^(k-1)(M_k^(-1))_(b,a)` is integral and can also be returned modulo
`2^p` at polynomial bit cost. These are individual query guarantees;
global trace summation is still the open operation.

## 5. Exact controls

Run:

```
python3 tools/check_pole_cauchy_conditioning.py \
  --output /tmp/pole-cauchy-conditioning.json
cmp /tmp/pole-cauchy-conditioning.json evidence/pole-cauchy-conditioning.json
```

The checker builds 200 padded finite matrices over `Q(w)`. It verifies
factor-count preservation and residue balance, computes their inverses
by exact Gaussian elimination, and checks 2,185 inverse entries against
both the cofactor-capacity formula and the normalized product oracle.
Two hundred rounded-entry traces use exactly the precision in (7) and
recover the expected factor count modulo sixteen.

The long control has a 402-bit input and padded domains of size
`2^200+13`, represented by three total progressions. It evaluates three
inverse entries at twelve unit bits and displacement `k=615`, without
enumerating the domains. Its exact valuations are `-16,-9,-234`, all
within the proved bound `1-k`. The receipt records the actual unit
coordinates, parameters, and finite sharpness example. This is an
implicit inverse-entry computation, not an evaluation of the global
factoring trace.
