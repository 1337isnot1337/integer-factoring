# First-jet separation of constant-field carry weights

## Result

The constant-field weighted carry interface has a simpler form that
requires only an **even** digit width. Its carry-dependent scalar is
determined by the first derivatives of the operand digit polynomials
at one, modulo three. It separates into two operand cases and into
eighteen digitwise phase terms, all over the existing `F_4`.

Thus the digit width need only be rounded upward by at most one bit.
The original weighted-carry representation retains its unit Fourier
factors. Its new operand-only representation uses the ordinary
unweighted carry constraints and adds no vertices or wires: only the
existing local input weights change.

This is an exact interface improvement, not a faster contraction
algorithm. The first-operand sum and cyclic coefficient extraction
remain the unresolved operations in
[weighted-dyadic-carry-contraction.md](weighted-dyadic-carry-contraction.md).
The established score remains `1/2`.

The [constant-field note](constant-field-carry-weights.md) remains useful:
its alternative phases all have order nine and give unit operand-box
factors. Some phases in the simpler `F_4` expansion below equal one,
so those alternate second-operand Fourier factors can vanish.

## 1. Parameters and the exact integer polynomial identity

Use the adaptive `D=2h^2`, so `t=log2 D` is odd. Let `w` be the
least even integer at least the preliminary width of either adaptive
constructor. Set

```math
b=2^w,\qquad
m=D\left\lceil\frac{2w+t}{D}\right\rceil,
\qquad M=2^m+1.
```

The root field `F_(2^(2m))` contains a constructively obtainable
nontrivial cube root `beta`. Since `m` and `w` are even,

```math
M=2\pmod3,\qquad b=1\pmod3,
\qquad\gamma=\beta^2,\quad\gamma^D=\beta.
\tag{1}
```

The root-field construction, tensor recovery algebras, and all unit
proofs of the constant-field note still apply. This even width is
no larger than that note's width `2 mod6`. All its explicit written
table bounds therefore hold by monotonicity in `w` and `m`, including
`99n,117n/P` for the first constructor and `49n,57n/P` for the
operand-length constructor. Input capacity increases by at most `D`
per `D`-digit block. The modulus ceiling can still add one `D`-bit
block, so its total residue-volume overhead is `O(D^2)=o(n)` under
the stated cutoff.

For an accepted factor pair `pq=N`, write its integer digit
polynomials as

```math
P(X)=\sum_i p_iX^i,\qquad Q(X)=\sum_j q_jX^j,
\qquad N_b(X)=\sum_j n_jX^j,
```

and put

```math
R(X)=\sum_{j=1}^{T-1}r_jX^{j-1},
\qquad\rho=R(1).
```

The exact carry equations give the identity over the integers

```math
P(X)Q(X)=N_b(X)+(b-X)R(X).                     \tag{2}
```

Differentiate before reducing any coefficient. At `X=1`, (2) gives

```math
P'(1)Q(1)+P(1)Q'(1)
=N_b'(1)+(b-1)R'(1)-\rho.
```

Reduce this integer identity modulo three and use `b=1 mod3`:

```math
\rho\equiv N_b'(1)-P'(1)Q(1)-P(1)Q'(1)\pmod3.
\tag{3}
```

These residue tags are ordinary integers modulo three. The scalar
weights below live in the different, characteristic-two field `F_4`.
No derivative in characteristic two is being substituted for (3).

## 2. An exact rank-two coefficient kernel

Let

```math
s(z)=\sum_i z_i\pmod3,
\qquad j(z)=\sum_i i z_i\pmod3
```

be the value and derivative tags of the base-`b` digit polynomial.
The letter `j` denotes a tag here, not a digit position. Because
`b=1 mod3`, `s(z)=z mod3`.

The roughness pretest handles the factor three. For the remaining
input put `c=N mod3 in {1,2}` and `h=j(N)`. If `p=u mod3` and
`q=v mod3`, then `u,v` lie in `{1,2}` and `uv=c mod3`.
Equation (3), together with `-2=1 mod3`, gives

```math
\beta^{q+2\rho}
=\beta^{v+2h}\,
 \beta^{v j(p)}\,
 \beta^{u j(q)}.                               \tag{4}
```

There are exactly two choices for `u`, and each determines `v`.
Thus (4) is a sum of two separated operand weights, with the two
residue selectors attached to their corresponding operands.

The associated `9` by `9` matrix has row tag `(s(p),j(p))` and
column tag `(s(q),j(q))`, both in `(Z/3Z)^2`. Its entries are

```math
K_{c,h}((u,r),(v,s))
=[uv=c\pmod3]\,
 \beta^{v+2h+v r+u s}.                         \tag{5}
```

It has two disjoint nonzero `3` by `3` blocks. Each block is an outer
product with nonzero entries. Its rank over `F_4` is therefore
**exactly two**, for every `c in {1,2}` and `h in Z/3Z`.

The two tags require only nine states per operand. They update on a
digit at position `i` by adding `(z_i,i z_i)` modulo three. The
position-dependent coefficient has period three.

## 3. Eighteen phase terms on the existing input vertices

The residue selector satisfies, in `F_4`,

```math
[s(p)=u]=\sum_{a=0}^2\beta^{a(s(p)-u)}.
\tag{6}
```

For an equal residue this is three, which equals one in
characteristic two; otherwise the cube-root sum is zero.
Combine (6) with the derivative factor in (4):

```math
[s(p)=u]\beta^{v j(p)}
=\sum_{a=0}^2\beta^{-a u}
 \prod_i\beta^{(a+v i)p_i}.                    \tag{7}
```

The corresponding formula for the second operand has index `d`:

```math
[s(q)=v]\beta^{u j(q)}
=\sum_{d=0}^2\beta^{-d v}
 \prod_j\beta^{(d+u j)q_j}.                    \tag{8}
```

Consequently, the precise weighted prefix query

```math
F_{N,k,x}(\alpha)
=\sum_{\substack{p q=N\\\lfloor p/2^k\rfloor=x}}
 \alpha^p\beta^{q+2\rho(p,q)}
```

equals the sum of the following eighteen ordinary factor
contractions:

```math
\sum_{\substack{u\in\{1,2\}\\v=c u^{-1}\bmod3}}
\ \sum_{a,d=0}^2
 \beta^{v+2h-a u-d v}
 \sum_{\substack{p q=N\\\lfloor p/2^k\rfloor=x}}
 \left[\prod_i
  \left(\alpha^{b^i}\beta^{a+v i}\right)^{p_i}\right]
 \left[\prod_j\left(\beta^{d+u j}\right)^{q_j}\right].
\tag{9}
```

All carry constraints in the inner contractions are unweighted. The
scalar factors on input digits are nonzero and can be inserted into
the existing input vertex weights. Each vertex knows the positions
of its own digits. The original prefix and leading-bit tests remain
unchanged. There is no new coefficient extension, side channel,
residue-check vertex, or uncharged summation in (9).

The first-operand factor is explicitly positional:
`alpha^p=product_i(alpha^(b^i))^(p_i)`. No requirement `ell|w` has
been reintroduced. Writing and evaluating eighteen instances is a
constant multiplier on the previously stated costs.

## 4. Unit factors and the remaining algorithmic task

The original representation of the left side of (9), with fixed
`beta` and carry weight `gamma=beta^2`, still has unit second-operand
and carry Fourier factors at every frequency. Its first-operand
factors have the positional recovery weights just described.

In the operand-only expansion, the second-operand phase
`beta^(d+u j)` equals one at one of every three digit positions.
Its full-box factor at zero frequency is then `S_b(1)=0`. The
unweighted carry factors also have their usual zeros. These are
properties of this alternate representation, not errors to remove.

The first-operand factors in (9) are still units at every frequency.
An odd recovery degree makes the order of `alpha^(b^i)` coprime to
three and to `M`; multiplying it by a cube-root phase cannot cancel
it. For `c=alpha^(b^i) beta^e`, the exact identity
`(1+c^M) S_3(c^M)=1+alpha^(3M b^i)` reduces its inverse to a
nonzero scalar in the recovery field. Combine it with
`(1+Y)S_M(Y)=1+c^M` for `Y=c zeta^s`.

The exact identity can transfer an algorithm between two concrete
interfaces: the common all-unit weighted-carry quotient, and eighteen
period-three input-phase instances of the ordinary multiplication
network. It does not assert a fast evaluation of either interface.
In particular, the remaining first-operand sum and cyclic coefficient
operation cannot be charged only by the two or eighteen terms here.

## 5. Reproduction

Run:

```sh
python3 tools/check_first_jet_carry_weights.py \
  --output evidence/first-jet-carry-weights.json
```

The checker differentiates literal digit convolution identities,
checks all six tag kernels, verifies the eighteen-term equality on
actual factor prefixes and individual digit vectors, exercises even
widths in all three residue classes modulo six, and rechecks the
uniform constructor bounds for the smaller width. It includes the
zero-factor controls of the operand-only interface.
