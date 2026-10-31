# Polynomial-bit-time jets of binary-described interval products

## Constructive result and scope

The [interval orbit-gcd interface](interval-orbit-gcd.md) needs an
effective representation of products whose degrees can be exponential
in the input bit length. This note supplies their **local Taylor jets**
without expanding either input product.

For odd `N`, intervals `I,J` with positive binary endpoints, precision
`p>=1`, and jet degree `d>=0`, one deterministic algorithm computes

```math
P_I(\omega+Z),\quad Q_{N,J}(\omega+Z)
\quad\text{in }\quad
A_{p,d}=\frac{(\mathbb Z/2^p)[\omega,Z]}
 {(\omega^2+\omega+1,\ Z^{d+1})}.                 \tag{1}
```

Its bit cost is polynomial in `n,p,d+1`, where `n` bounds the bit
length of `N` and all endpoints. It never visits the interval elements.
It works more generally after substituting any given jet whose constant
residue in `F_4` lies outside `F_2`.

This is a local evaluation algorithm, not the unresolved global gcd
algorithm. In particular, both products in (1) are units, so computing
an ordinary ideal gcd inside this truncated ring cannot recover the
common global factor. The complete factoring score remains `1/2`.

## 1. Parity splitting makes every base factor a unit

The polynomial `X^2+X+1` is irreducible over `F_2`. Thus `omega` and
`omega-1` are units modulo two in the unramified quadratic ring. An
element `a+b*omega` has norm `a^2-ab+b^2`; if this is odd, its inverse
is computed by

```math
(a+b\omega)^{-1}
=\bigl((a-b)-b\omega\bigr)(a^2-ab+b^2)^{-1}\pmod{2^p}.
```

A jet with unit constant term has its inverse from the usual triangular
power-series recurrence, using `O((d+1)^2)` quadratic-ring operations.

Let `X` be a jet with constant residue outside `F_2`. Splitting an
integer interval by `a=r+2t`, `r in {0,1}`, gives factors

| Product | Unit base `u_r` | Remaining factor |
| --- | --- | --- |
| `P_I(X)` | `X-r` | `u_r-2t` |
| `Q_(N,J)(X)` | `N-rX` | `u_r-2tX` |

Both bases are units: for the second row, `N` is odd, and the two
constant residues are one and `1-X`, both nonzero. The parity-restricted
values of `t` form another ordinary integer interval with binary
endpoints computable by ceiling division.

## 2. Only the first `p` symmetric functions are needed

For any unit `u`, any jet `V`, and consecutive integers `A<=t<B`, put
`ell=B-A`. The finite product identity gives

```math
\prod_{t=A}^{B-1}(u-2tV)
=u^\ell\sum_{k=0}^{\ell}
 e_k(A,\ldots,B-1)(-2Vu^{-1})^k.
```

The `k`th term is divisible coefficientwise by `2^k`. Therefore modulo
`2^p` the exact expression is

```math
u^\ell\sum_{k=0}^{\min(\ell,p-1)}
 e_k(A,\ldots,B-1)(-2Vu^{-1})^k.                 \tag{2}
```

This truncation is independent of the numerical length `ell` and of
the Taylor degree `d`. It introduces no analytic approximation and
uses no division by two in the coefficient ring.

Compute the necessary symmetric functions from integer power sums.
For `S_j(H)=sum_(0<=t<H)t^j`, telescoping powers yields

```math
S_j(H)=\frac{H^{j+1}-\sum_{i<j}\binom{j+1}{i}S_i(H)}{j+1}.
```

Take `s_j=S_j(B)-S_j(A)`. Newton's identities then give

```math
e_0=1,\qquad
e_k=\frac1k\sum_{j=1}^k(-1)^{j-1}e_{k-j}s_j.     \tag{3}
```

All divisions in these two recurrences are exact divisions of ordinary
integers. They occur before reduction modulo `2^p`, so nonunits among
`j+1` and `k` cause no loss of information. An implementation that
instead divides these scalars in the residue ring would be invalid.

For `k<p`, the power sums and symmetric functions have `O(np)` bits:
`|s_j|<=ell*max(1,B)^j` and
`e_k<=binom(ell,k)*max(1,B)^k`. Consequently the `O(p^2)` integer
recurrence operations have polynomial bit cost. The power `u^ell` in
(2) takes `O(n)` binary-powering steps in the finite jet ring; it is
not computed as a huge unreduced integer. Combining the two parities
proves the stated uniform bit bound.

## 3. What this changes in the elimination question

Local jets of the two long input polynomials are now effective data.
They need no dense coefficient list and no uncharged long product.
The construction also supports their boundary identities under
`X -> X+1` and `X -> NX/(N+X)`: both substitutions preserve the
condition that the constant residue is outside `F_2`.

It remains necessary to identify the **global polynomial factor**.
Every candidate polynomial with unit value at the chosen point divides
each input as a formal power series, because a unit is invertible.
Thus local divisibility alone imposes no global factor condition. A
method must use the finite polynomial degrees, the two boundary
recurrences, or another global certificate. The small-output circuit
result audited in the companion note does not remove this obligation.

The next constructive question is whether those boundary recurrences
allow global elimination using these inexpensive jets while keeping
all degree, precision, and scalar costs within the desired bound.

## 4. Exact replay

```sh
python3 tools/check_interval_unit_jets.py --output /tmp/interval-unit-jets.json
cmp /tmp/interval-unit-jets.json evidence/interval-unit-jets.json
make verify-arithmetic
```

The [implementation](../tools/interval_unit_jets.py) uses parity,
exact integer power sums, Newton identities, and finite quadratic jets.
The [checker](../tools/check_interval_unit_jets.py) passed 256 symmetric
function controls, 1,152 comparisons with literal interval products,
and 576 inverse controls over binary precisions 1 through 12 and
Taylor degrees 0 through 5.

The [saved long-range cases](../evidence/interval-unit-jets.json) include
an interval of length `2^200+12`, with 241-bit endpoints, at precision
48 and Taylor degree eight. It is checked by splitting the interval
and by both exact boundary recurrences; the smaller long-range case
uses length `2^110+12`. These eight long-range identities are not
literal product expansions. The receipt saves the resulting jets and
parameters so their exact finite claims can be replayed.
