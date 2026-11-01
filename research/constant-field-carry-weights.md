# Constant-field carry weights and exact operand separation

## What improves

The aligned weighted-carry construction does not need its digit radix
to depend on the sparse recovery field. A fixed nontrivial cube root
`beta` in `F_4` supplies all auxiliary weights. Round the digit width
only to `w=2 mod 6`, adding at most five bits per digit.

This has three concrete consequences.

- The multiplication network and its root field are built once for an
  input length and group count, independently of the recovery degree.
- Input capacity increases by `O(D)`, replacing the earlier
  `O(D ell)` adjustment. Both adaptive constructors retain their stated
  uniform written-table bounds and their leading constants.
- After the roughness pretest removes the factor three, the combined
  second-operand and carry weight has an exact rank-two separation
  into functions of the two operand digit sums modulo nine.

An optional expansion writes that separation as eighteen product
weights from a fixed field of 64 elements. The actual multiplication
and carry constraints remain in place. None of these identities alone
evaluates the remaining coupled first-operand sum or cyclic coefficient
in [weighted-dyadic-carry-contraction.md](weighted-dyadic-carry-contraction.md).
The established score remains `1/2`.

## 1. Use a fixed cube root and a constant width adjustment

Use the adaptive dimensions `D=2h^2`, where `h` is a power of two,
and put `t=log2 D`. Thus `t` is odd. From the preliminary width `w_0`
of either adaptive constructor, choose the least integer

```math
w\ge w_0,\qquad w\equiv2\pmod6.
```

Then `w-w_0<=5`. Put

```math
b=2^w,\qquad
m=D\left\lceil\frac{2w+t}{D}\right\rceil,
\qquad M=2^m+1.                                \tag{1}
```

As before, `M>Db^2`, `D|m`, and the exact transform and padded carry
proofs apply. Since `D>=2`, `m` is even. Hence

```math
M\equiv2\pmod3,\qquad b\equiv4\pmod9.
\tag{2}
```

Construct `F=F_(2^(2m))` and its primitive `M`th root `zeta` as in
the carry Fourier note. This field contains `F_4`. Find a nontrivial
cube root `beta` by enumerating nonzero `a in F` and testing

```math
\beta=a^{(2^{2m}-1)/3}\ne1.
```

The multiplicative group of a finite field is cyclic, so this
projection reaches both nontrivial cube roots. The exhaustive worst
case is `2^O(m)` bit operations, already within the charged root-field
construction budget. There is no factoring oracle or free field
embedding. The result satisfies

```math
\beta^3=1,\qquad\beta^2+\beta+1=0.
```

Because `w` is even, `beta^b=beta`. Because `t` is odd, set

```math
\gamma=\beta^2,
\qquad\gamma^D=\beta.                          \tag{3}
```

Thus the same two constants work at every recovery degree. If the
general construction uses a power-of-two `D` with even `t`, the
corresponding choice is `gamma=beta`; only the adaptive odd-`t` case
is needed here.

For an odd-prime recovery degree `ell`, use

```math
\mathcal A_\ell=F\otimes_{\mathbb F_2}\mathbb F_{2^\ell}.
\tag{4}
```

Its dimension is still `2m ell` over `F_2`. There is no additional
coefficient-field factor. Only the recovery variable `alpha` ranges
over `F_(2^ell)\{0,1}`. The two constants `beta,gamma` remain in the
already constructed copy of `F_4` inside `F`.

## 2. The recovery polynomial and all unit statements survive

For each divisor let `rho(d)` be the sum of its genuine internal
base-`b` multiplication carries. The univariate prefix polynomial is

```math
F_{N,k,x}(X)=
\sum_{\substack{d\mid N\\\lfloor d/2^k\rfloor=x}}
 \beta^{N/d+2\rho(d)}X^d.                       \tag{5}
```

Every coefficient is one of the three nonzero elements of `F_4`.
The polynomial and its coefficients are now independent of `ell`.
After the same roughness pretest, its support is at most `K` and all
its exponents are below `2^n`. The univariate punctured-field theorem
therefore gives exactly the same deterministic recovery family and
cost as before. Remove the known term `beta^N X` for divisor one.

The first-operand weight must be retained as `alpha^p`. In digits it
is

```math
\alpha^p=\prod_i\left(\alpha^{b^i}\right)^{p_i}.
\tag{6}
```

It is generally **not** `alpha^(sum_i p_i)`, because `ell` need not
divide `w`. Each `alpha^(b^i)` is still nontrivial: raising to a power
of two permutes the recovery field. The second-operand weight does
simplify to `beta^(sum_j q_j)`, since `beta^b=beta`.

Every second-operand and weighted-carry Fourier box stays a unit.
The order-three constants have order coprime to `M`, and

```math
1+\beta^M=1+\beta^2=\beta\ne0.
\tag{7}
```

The same argument applies to `gamma` and its power-of-two images.
For a first-operand box, its scalar is `alpha^(b^i)` in the recovery
field. Odd `ell` gives `gcd(M,2^ell-1)=1`, so the tensor-algebra
inverse formula from the sparse recovery note applies unchanged.
It does not require that `F` and the recovery field have trivial
intersection, nor that their tensor product be a field.

The common Frobenius quotient of the weighted-carry note is unchanged:
`sigma:f->f^b` fixes `beta,gamma`, and every binomial of `H_DG_p` has
scalar `beta`. The first-operand factor outside that quotient is
`alpha^p`, as in (6). Its common scalar denominator is especially
simple:

```math
(1+\beta^M)^\nu=\beta^\nu,
\qquad\nu=T-1+d_q.                             \tag{8}
```

Clearing it has the same support-preserving interpretation already
proved for the disjoint bit-length instances.

## 3. Complete capacity and written-table costs

The only input-capacity increase is at most `5D` per `D`-digit block,
and `O(D)` across both operands. The modulus-rounding ceiling can add
one further `D`-bit block compared with the old modulus. Consequently
the worst increase in total residue volume is `O(D^2)`, not necessarily
`O(D)`. Under `P^40<=n`, this is `O(n^(9/10))=o(n)`.

The original explicit constants have enough margin as well. For the
first adaptive constructor, `C=Dw<n+6D`, `D^2<=n/2`, and `D>=2^19`.
Its proved aggregate bounds therefore give

```math
U\le36C+20Dt+16D^2
<36n+(216+20t)D+16D^2<45n,
```

```math
S\le20C+10Dt+8D^2
<20n+(120+10t)D+8D^2<25n.                       \tag{9}
```

For all `D>=2^19`, both `216+20log2 D` and
`120+10log2 D` are less than `D`. This proves the strict numerical
inequalities in (9). The earlier grouping proof therefore still gives
`P^2 lambda<=99n` and `B_tab<=117n/P`.

For the operand-length constructor, its preliminary capacity is
`C_0<0.501n`. The same inequalities imply

```math
C<0.50101n,\qquad Dt/n<1/2^{15}.
```

At `P=2`, use `D^2/n<=1/4`; at `P>=3`, use the proved
`D^2/n<19/50`. Substituting into the same expressions for `U,S`
retains respectively

```math
(U/n,S/n)<(22.1,12.1)\quad(P=2),
\qquad(U/n,S/n)<(25,14)\quad(P\ge3).
\tag{10}
```

Hence its earlier uniform bounds `P^2 lambda<=49n` and
`B_tab<=57n/P` also remain valid. The leading constants of both
constructors are unchanged.

All field construction, cube-root projection, coefficient arithmetic,
table writing, and field points are charged. The complete transfer
cost remains

```math
\exp O\!\left(
 \log(A(n)+2)+n/P+u+L+m+\log n
\right),                                      \tag{11}
```

with the notation of the sparse recovery theorem. The root field and
its constants are reusable across all degrees. There is no `D ell`
capacity condition in this version.

## 4. The carry scalar has an exact rank-two operand separation

The roughness pretest checks three; handle that factor immediately
if present. For the remaining input assume `3` does not divide `N`.
Write `s_b(z)` for its base-`b` digit sum. Summing the exact carry
equations without radix weights gives

```math
(b-1)\rho(p)=s_b(p)s_b(q)-s_b(N),\qquad pq=N.
\tag{12}
```

Let `c=s_b(N) mod3`, so `c` is one or two. For each of the two
values `u in {1,2}`, choose `v in {1,2}` with `uv=c mod3` and write

```math
s_b(p)=u+3i,\qquad s_b(q)=v+3j,
\qquad s_b(N)=c+3k.
```

Since `b=4 mod9`, `(b-1)/3=1 mod3`. Equation (12) yields

```math
\rho(p)\equiv
 \frac{uv-c}{3}+u j+v i-k\pmod3.
\tag{13}
```

Also `q=v mod3`, because `b=1 mod3`. Thus its coefficient in (5)
separates as

```math
\beta^{q+2\rho(p)}
=\beta^{\,v+2(uv-c)/3-2k}
 \beta^{2v i}\,\beta^{2u j}.                   \tag{14}
```

The first variable factor depends only on `s_b(p) mod9`, and the
second only on `s_b(q) mod9`. The prefactor depends on `N` and the
chosen residue case.

More explicitly, define a `9` by `9` matrix at the two operand digit
sums modulo nine. Put zero when their residues modulo three are
incompatible with `c`; put the value in (14) otherwise. It has two
disjoint nonzero `3` by `3` blocks. Each block is an outer product
with nonzero entries, so its rank over `F_4` is **exactly two**.

Therefore the weighted factor query is a sum of two contractions with
separate operand weights. No carry-dependent scalar is needed in
those two contractions. The original exact carry constraint is still
required. The extra operand information is only its digit sum modulo
nine, which has a nine-state addition automaton.

## 5. Eighteen product phases, if a digitwise interface is preferred

Adjoin a primitive ninth root `eta` with `eta^3=beta`. The polynomial
`X^3-beta` has no root in `F_4`, so it is irreducible and defines
`F_64`. Extending `F` by this polynomial costs only a factor three
in coefficient dimension; it is allowed to split over `F`.
Including the recovery field gives dimension `6m ell` over `F_2`.

For a residue case `(u,v)`, set `lambda=eta^(2v)`. The first operand
factor in (14), including its residue selector, is

```math
f_u(r)=[r=u\pmod3]\,
 \beta^{2v(r-u)/3}
=\sum_{a=0}^2(\lambda\beta^a)^{r-u},
\qquad r\in\mathbb Z/9\mathbb Z.              \tag{15}
```

The right side vanishes for incompatible residues by the cube-root
sum. For compatible residues, all three terms have the required
value; their sum is that value because three equals one in
characteristic two. The second operand has the analogous formula
with `mu=eta^(2u)` and exponent `s-v`.

Substitution into (14) gives `2*3*3=18` terms, each a constant times
a product of digitwise operand phases. Every phase has order nine:
its exponent on `eta` is `2v+3a` or `2u+3a`, never divisible by three.

These operand phases remain compatible with Fourier unit formulas.
The modulus `M` is coprime to nine. A second-operand phase has order
nine, so multiplying it by an `M`th-root character never gives one.
For a first-operand phase multiplied by `alpha^(b^i)`, the recovery
order `2^ell-1` is coprime to both nine and `M`; the corresponding
geometric inverse is therefore also a unit in the tensor algebra.
This optional interface has constant field and term-count overhead.

The inverses have explicit identities even if that algebra has zero
divisors. For `c=lambda alpha^(b^i)`,

```math
(1+c^M)S_9(c^M)=1+\alpha^{9M b^i},
```

and the right side is a nonzero recovery-field scalar. For a second
operand phase alone,

```math
(1+\lambda^M)S_3(\lambda^M)=1+\lambda^{3M},
```

whose right side is a nonzero `F_4` scalar. Combine these with
`(1+Y)S_M(Y)=1+c^M` for `Y=c zeta^s`. The resulting inverses cost
`O(m)` algebra operations plus ordinary scalar inverses; they require
no decomposition of the tensor algebra.

Removing the carry weights changes the carry Fourier factors back
to the unweighted ones, which can vanish. The rank-two or eighteen-term
identities are exact alternate representations, not a claim that a
carry constraint or the coupled frequency sum can be omitted.

## 6. Exact checks

Run:

```sh
python3 tools/check_constant_field_carry_weights.py \
  --output evidence/constant-field-carry-weights.json
```

The checker constructs cube roots inside the root fields, tests
recovery degrees that do not divide the digit width, checks the unit
formulas and positional first-operand weights, verifies the uniform
constructor bounds with exact integer arithmetic, and compares the
rank-two and eighteen-phase identities entry by entry. Literal
multiplications check (12)--(14), including nonzero carries.
