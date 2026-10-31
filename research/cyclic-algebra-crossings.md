# Regular cyclic-algebra crossing transforms

## Result and scope

The regular-nilpotent crossing constructor remains valid after every
polynomial change of its nilpotent generator with nonzero linear term.
This gives explicit, generally dense operators, beyond the particular
controlled-monomial hypothesis in the
[predecessor](controlled-permutation-smith-floor.md). It preserves the
ideal column orders, while its row orders can change.

The first complete search covers 1,024 transforms of one full 64-cycle
at eight physical letters. Every tested transform has bare Smith
parameter at least `31/32`, far above the required `1/sqrt(2)`. This is
an exact finite exclusion, not an obstruction for all cyclic-algebra
transforms or for other interpolation schemes. The factoring exponent
remains `1/2`. No Lean source is changed.

## 1. Construction and exact endpoint

Let `q` be a power of two, `n=q^2`, and let `U0` be any full cycle on
`n` coordinate vectors over a field of characteristic two. Then
`S0=I+U0` is regular nilpotent of order `n`. Choose

```math
g(X)=X+g_2X^2+\cdots+g_{n-1}X^{n-1},\qquad
S=g(S_0),\quad U=I+S.                                      \tag{1}
```

Binary coefficients suffice for the implementation. Since
`g(S0)=S0 h(S0)` with `h(0)=1`, the factor `h(S0)` is invertible and
commutes with `S0`. Thus `rank S^j=n-j` for `0<=j<=n`; in particular,
`S` is again one regular nilpotent block. Equivalently, if `U=f(U0)`,
then the conditions are `f(1)=1`, `f'(1)=1`.

Use the predecessor's reshuffling

```math
\Phi(A)(x,y,u,v)=A_{(x,u),(v,y)}.
```

Define the degree-at-most-`n-1` polynomial tensor

```math
K_g(z)=\Phi\!\left(\sum_{j=0}^{n-1}z^{n-1-j}S^jU\right).
                                                               \tag{2}
```

Because `(I+S) sum_j S^j=I`, evaluation at one gives
`K_g(1)=Phi(I)=[x=v][y=u]`. This is the complete physical crossing
endpoint, without a restriction to a selected type or an endpoint repair.

The column flattening is the displayed operator polynomial, up to a
fixed permutation of columns. It equals
`z^n (zI+S)^(-1) U`. The Smith form of `zI+S` for one regular block is
`diag(1,...,1,z^n)`. Multiplication by `z^n`, inversion, and the
constant invertible factor `U` therefore give column orders

```math
(0,\underbrace{n,\ldots,n}_{n-1\text{ times}}).                 \tag{3}
```

The complete normalized reversal at infinity is

```math
H_g(w)=w^{n-1}K_g(w^{-1})
      =\Phi\left(U(I+wS)^{-1}\right).                         \tag{4}
```

Its row Smith orders `mu_i` are computed, rather than inferred from
(3). Equation (4) retains the full degree cost even if individual entries
have smaller degrees. Both cuts are generically invertible by the endpoint.

All fields, generators, and dimensions in a seed are fixed. Explicit
binary matrices compute every entry of (2); tensor powers, coefficient
extraction, and both rank factorizations then have `exp(O(lambda))`
construction cost in binary label width `lambda`. Thus any future strict
Smith witness from this construction can use the existing
[polynomial interpolation interface](polynomial-crossing-interpolation.md)
and [finite-field transfer](scale-one-modular-transfer.md). The finite
search below does not supply such a witness.

## 2. An exact two-jet exclusion certificate

The bare two-place product for (3)--(4) is

```math
B_g(s)=\left(\sum_i s^{n-1-\mu_i}\right)
                 \left(1+(n-1)s^{-n}\right),\qquad s\ge1,
\quad \rho_g={\inf_{s\ge1}\log B_g(s)\over4\log q}.           \tag{5}
```

Let `r` count zero orders, and `k` count order-one pivots. The
predecessor's two elementary inequalities give

```math
B_g(s)\ge nr+(n-1)k.                                      \tag{6}
```

For completeness, the respective contributions are
`s^(n-1)+(n-1)/s >= n` and
`s^(n-2)+(n-1)/s^2 >= n-1`. For the second, minimize at
`s^n=2(n-1)/(n-2)`; Bernoulli's inequality with exponent `n/2`
proves the stated lower bound for every `n>=4`.

If `A0,A1` are the first two coefficients of `H_g`, literal binary
elimination recovers the two multiplicities:

```math
r=\operatorname{rank}A_0,\qquad
k=\operatorname{rank}\begin{pmatrix}A_0&0\\A_1&A_0\end{pmatrix}-2r.
                                                               \tag{7}
```

Consequently (6) is a certificate using only two finite matrices. It
does not require trusting a floating-point minimizer or all later orders.

## 3. Completed eight-letter schedule

The reference full cycle is the controlled-transposition example from
the predecessor. On eight letters let `J(i)=i+1 mod 8`, let `T`
interchange zero and one, and set `L_j=T J` for `j=1,2,3`, and `L_j=J`
otherwise. Then `U0(i,j)=(j,L_j(i))` is a single 64-cycle.

The exact schedule is every binary polynomial

```math
g(X)=X+\sum_{j=2}^{11}b_j X^j,
\qquad (b_2,\ldots,b_{11})\in\{0,1\}^{10}.                  \tag{8}
```

There are 85 different full row profiles. The smallest exact right-hand
side of (6) is `3248`, at binary polynomial mask `3314`, where `r=35`
and `k=16`. Since

```math
3248^8>8^{31},
```

every member of (8) satisfies `rho_g > 31/32`. In particular, none can
meet the strict post-half threshold through (5).

A numerical minimization of the full profiles finds the best parameter
about `0.9991774294`, at mask `14`. This transform is merely
`U=U0^3`: `g(X)=X+X^2+X^3`. The numerical value is descriptive;
(6)--(8) supply the exact exclusion. Other transforms, such as mask
`6`, give dense operators and different profiles. No claim is made for
other cycles, larger degree supports, nonbinary coefficients, physical
gauges, or collective coefficient recombinations.

## 4. Reproduction and independent controls

```sh
make verify-cyclic-algebra
python3 tools/check_cyclic_algebra_crossings.py \
  --search-bits 10 --output /tmp/cyclic-algebra-search.json
```

The first command uses only Python's standard library. It checks 13
seeds at two, four, and eight letters. Literal matrix multiplication
independently evaluates (1), verifies selected ranks of `S^j`, compares
every coefficient in (2), and checks every crossing endpoint entry.
Binary block-Toeplitz ranks recover all 13 row profiles and the complete
column profiles in the two- and four-letter controls. Invalid constant,
linear, and degree conditions are rejected.

The optional search needs a C++17 compiler for the existing native Smith
engine. Every native profile is checked against independent exact
two-jet ranks. A representative of each distinct complete profile is
also recovered independently by block-Toeplitz ranks through its largest
order. The engine promises every pivot from the exact endpoint, and
raises if its selected precision fails to resolve one; missing pivots
are never interpreted as zeros. Saved evidence is in
[`cyclic-algebra-crossings.json`](../evidence/cyclic-algebra-crossings.json).

The regular-nilpotent construction and elementary pressure inequalities
retain the predecessor's attribution. The polynomial-generator extension,
literal constructor, and bounded search used substantial OpenAI Codex
assistance. The next constructive question is whether another cyclic
coordinate ordering and higher-degree generator can jointly increase
the row orders without sacrificing (3).
