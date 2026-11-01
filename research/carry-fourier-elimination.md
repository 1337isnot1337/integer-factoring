# Exact Fourier evaluation of the whole bounded carry tensor

## Result and role in a coupled evaluator

The whole carry-check tensor of the adaptive multiplication network has
an efficiently evaluable additive Fourier transform. The exponentially
many carry assignments can be summed with a constant number of field
elements, geometric sums, and exact integer comparisons.

Theorem. Let `D>=2`, `b=2^w>D`, `K=D(b-1)`, and `M=2^m+1`, with
`D b^2<=2^m`. Fix `T` target digits `0<=n_j<b`. There is deterministic
preprocessing of bit cost `2^O(m)` that constructs a field of `2^(2m)`
elements and an element `zeta` of exact order `M`. After this
preprocessing, the additive Fourier transform of the complete bounded
carry acceptance tensor at **any specified** frequency vector
`s in (Z/MZ)^T` is computable with `O(T m)` field operations. With
schoolbook field arithmetic, its bit cost is polynomial in `T,m`.

For the [adaptive constructor](adaptive-transform-constructor.md),
`m=O(n/P^18)`, so the preprocessing is below the already charged
`2^O(n/P)` table budget. The same conclusion holds for the
[operand-length refinement](operand-length-transform-constructor.md).
No root of unity or extension-field operation is treated as unit cost.

The result is a constructive primitive for eliminating the inverse and
carry layers together in a Fourier representation. It does not give a
new factoring exponent. Section 7 states the exact remaining bilinear
sum over bounded input digits. Section 8 also explains why efficient
evaluation at one frequency vector does not imply a small uniform
tensor factorization.

## 1. The complete carry tensor

For a coefficient tuple `a in {0,...,M-1}^T`, let `W_N(a)` be one if
there are carries `r_0,...,r_T` satisfying

```math
\begin{gathered}
r_0=r_T=0,\qquad 0\le r_j\le K,\\
a_j+r_j=n_j+b r_{j+1}\qquad(0\le j<T),
\end{gathered}                                             \tag{1}
```

and zero otherwise. This is precisely the acceptance tensor from the
adaptive constructor, including the outgoing carry bound. Carries, when
they exist, are unique because (1) determines them successively from
`r_0=0`.

Let `zeta` have order `M` in a field, and write
`chi_s(x)=zeta^(sx)` with exponents reduced modulo `M`. Define

```math
\widehat W_N(s)
=\sum_{a\in(\mathbb Z/M\mathbb Z)^T}
 W_N(a)\,\chi_1\left(\sum_j s_ja_j\right).
\tag{2}
```

For any carries in `[0,K]`, the candidate coefficient
`n_j+b r'-r` is automatically below `M` whenever it is nonnegative:

```math
n_j+b r'-r\le b-1+bK
=Db^2-(D-1)b-1<Db^2\le2^m<M.                  \tag{3}
```

Thus only its lower bound remains. Equation (2) is exactly the carry
path sum with transition matrices

```math
T_{n,s}(r,r')
=\chi_s(n+b r'-r)\,[\,r\le n+b r'\,],
\qquad0\le r,r'\le K.                         \tag{4}
```

The initial vector is supported at carry zero, and the final answer is
its carry-zero coordinate after multiplying all `T` transitions.

## 2. The first compression: a small rank at a fixed frequency

For fixed `n,s`, the character in (4) is a nonzero row scaling times a
nonzero column scaling. Its rank is therefore the rank of

```math
H_n(r,r')=[\,r\le n+b r'\,].                   \tag{5}
```

Every column with `r'>=D` is all ones because `n+bD>K`.
The earlier columns are nested nonempty initial segments. If

```math
q_n=\left\lceil\frac{K-n}{b}\right\rceil,
```

the first `q_n` columns are distinct proper initial segments, and all
remaining columns are all ones. Their successive differences have
disjoint nonempty supports. Consequently, over every field,

```math
\mathrm{rank}(T_{n,s})=q_n+1\le D+1.           \tag{6}
```

This already permits a polynomial-dimensional matrix computation for
one specified frequency vector. The next sections remove even the
polynomial number of state coordinates.

## 3. Geometric sums without division

All formulas below work in a commutative ring containing the required
character values. Define

```math
S_z(L)=\sum_{r=0}^{L-1}z^r,\qquad S_z(0)=0,
\qquad
J_L(z,\beta)=\sum_{r=0}^{L-1}z^rS_\beta(r).
\tag{7}
```

They are computable in `O(log(L+1))` ring operations by binary powering.
There is no division by `z-1`, so `z=1`, positive characteristic, and
zero divisors cause no exceptional case. For example,

```math
\begin{aligned}
S_z(L+Q)&=S_z(L)+z^L S_z(Q),\\
J_{L+Q}(z,\beta)
&=J_L(z,\beta)
 +z^L\bigl(S_\beta(L)S_z(Q)+\beta^LJ_Q(z,\beta)\bigr).
\end{aligned}                                             \tag{8}
```

The five-tuple
`(z^L,beta^L,S_z(L),S_beta(L),J_L(z,beta))` therefore has a
constant-cost concatenation operation. Binary powering starts from the
length-zero and length-one tuples.

## 4. A constant-size invariant for the carry distribution

After any number of transitions, the function of the current carry has
the following form, for seven stored field elements
`a,e,u,v,g,alpha,beta`:

```math
F(r)=
\begin{cases}
a,&r=0,\\
\alpha^r\bigl(u+vS_\beta(r-1)\bigr),&1\le r\le D-2,\\
e,&r=D-1,\\
g\alpha^r,&D\le r\le K.
\end{cases}                                               \tag{9}
```

For `D=2`, the middle interval is empty. Initially use
`a=1`, `e=u=v=g=0`, and `alpha=beta=1`; this is exactly the initial
carry-zero vector.

The separate endpoint `D-1` is necessary. For some target digits,
`n+b(D-1)` exceeds `K`, while for others it does not. Omitting that
endpoint distinction would make the following recurrence false.

Given the stored state and a new frequency `s`, define the weighted
prefix

```math
P_s(L)=\sum_{r=0}^{\min(L,K)}F(r)\chi_{-s}(r).
\tag{10}
```

Put `z=alpha chi_(-s)(1)`. The value in (10) is the sum of these at
most four contributions:

* `a` from carry zero;
* `u z S_z(q)+v z J_q(z,beta)`, with
  `q=min(L,D-2)`, when `q>=1`;
* `e chi_(-s)(D-1)` when `L>=D-1`;
* `g z^D S_z(min(L,K)-D+1)` when `L>=D`.

Each prefix uses `O(log(K+2))` ring operations by (8). All lengths in
these sums are nonnegative integers in the branches where they occur.

## 5. The exact state update

Fix the next target digit `n` and frequency `s`. Let

```math
\begin{gathered}
z=\alpha\chi_{-s}(1),\quad
\alpha'=\chi_s(b),\quad \beta'=z^b,\\
A_0=P_s(D-1),\qquad L_1=n+b-D+1.
\end{gathered}                                             \tag{11}
```

Because `b>D`, `L_1` is positive. Update the five coefficient values by

```math
\begin{aligned}
a'&=\chi_s(n)P_s(n),\\
e'&=\chi_s(n)\,(\alpha')^{D-1}
       P_s\bigl(\min(K,n+b(D-1))\bigr),\\
g'&=\chi_s(n)P_s(K),\\
u'&=\chi_s(n)\bigl(A_0+g z^D S_z(L_1)\bigr),\\
v'&=\chi_s(n)g z^{D+L_1}S_z(b).
\end{aligned}                                             \tag{12}
```

These formulas preserve (9) exactly. To prove it, multiplication by
(4) gives for every output carry

```math
F'(r')=\chi_s(n)(\alpha')^{r'}
              P_s\bigl(\min(K,n+b r')\bigr).       \tag{13}
```

The cases `r'=0`, `r'=D-1`, and `r'>=D` give the first three lines
of (12) directly. For `1<=r'<=D-2`, the bound `b>D` gives

```math
D\le n+b r'<K.
```

Thus its prefix contains all of the small carry block and an initial
part of the geometric tail:

```math
P_s(n+b r')=A_0+g z^D S_z(L_1+b(r'-1)).
```

The identity

```math
S_z(L_1+b(r'-1))
=S_z(L_1)+z^{L_1}S_z(b)S_{z^b}(r'-1)
```

is just a division-free grouping of consecutive geometric terms.
Substitution gives exactly the interior case of (9), with the last two
lines of (12). This completes the invariant proof.

After processing all `T` target digits, return the final `a`. The final
carry is thereby fixed to zero. Equations (10)--(12) use a constant
number of prefixes, powers, and geometric sums per digit. Since
`log(K+2)=O(m)` and character exponents have `O(m)` bits, the total is
`O(Tm)` ring operations.

## 6. Construct and charge the coefficient field

The special modulus permits a small extension degree:

```math
M=2^m+1\mid2^{2m}-1.
\tag{14}
```

Work in `F=F_(2^(2m))`. A fully deterministic, sufficient construction
is as follows.

1. Enumerate all monic binary polynomials of degree `2m` until an
   irreducible one is found. Test a candidate `f` using
   `X^(2^(2m))=X mod f` and
   `gcd(f,X^(2^((2m)/p))-X)=1` for every prime `p` dividing `2m`.
   Repeated squaring and polynomial gcd use polynomially many bit
   operations in `m` for each candidate. At most `2^(2m)` candidates
   are tested.
2. Factor the known integer `M` by exhaustive trial division, costing
   `2^O(m)` bit operations.
3. Enumerate nonzero field elements `a`, set `zeta=a^(2^m-1)`, and
   test `zeta^(M/p)!=1` for each distinct prime `p` dividing `M`.
   The first success has exact order `M`. A success exists because
   the finite field's multiplicative group is cyclic. At most
   `2^(2m)-1` elements are tried.

The polynomial irreducibility criterion is the standard deterministic
test underlying M. O. Rabin, *Probabilistic Algorithms in Finite Fields*,
SIAM J. Comput. 9 (1980), 273--280
([primary paper](https://doi.org/10.1137/0209024),
[MIT report](https://hdl.handle.net/1721.1/149500)). The algorithm here
uses exhaustive search and is deterministic; no expected-time claim is
being substituted for a worst-case bound.

The total preprocessing is `2^O(m)`. Once the field is built, an element
uses `2m` bits and schoolbook multiplication has polynomial bit cost.
All character values in the recurrence are computed by repeated
squaring. The field degree is `2m`, not `M`, and the existence of a
primitive root is accompanied by the explicit search cost above.

The recurrence also works in any commutative unital `F`-algebra. This
allows the cyclic polynomial rings used by
[sparse-divisor-nonzeroness.md](sparse-divisor-nonzeroness.md): the
character values lie in the coefficient field and remain units after
extension. Integer counts are not recovered from characteristic two;
weighted divisor polynomials retain their distinct exponent support.

For the adaptive parameters, `D>=2P^18`, `D<8P^18`, `n>=P^40`, and
`m<=2ceil(n/D)+t+D`. Hence `m=O(n/P^18)`. The preprocessing is therefore
absorbed by the existing `2^O(n/P)` top construction. This estimate
charges both the field and its primitive root.

The hypothesis `b>D` also holds in these constructors. The first
adaptive choice has `w>=n/D>=2D`; the operand-length choice has
`w>(n+1)/(2D)>=D`. Either gives `2^w>D`.

## 7. The exact remaining Fourier-dual problem

Let `p_i,q_j` be the two bounded input digit sequences, with supports
whose sum is at most `T`. Their convolution coefficients are exact
integers below `M`. Because `M` is odd, its image in the characteristic
two field is one, and Fourier orthogonality gives

```math
W_N(a)=\frac1{M^T}\sum_{s\in(\mathbb Z/M\mathbb Z)^T}
 \widehat W_N(s)\,\chi_1(-s\mathbin\cdot a).
\tag{15}
```

The scalar `1/M^T` is written to display the usual normalization; it
equals one in this field. For the weighted divisor-prefix polynomial,
substituting the exact convolution gives

```math
\frac1{M^T}\sum_s\widehat W_N(s)
 \sum_{p\in\mathcal P}X^{\sum_i p_i b^i}
 \sum_{q\in\mathcal Q}
       \chi_1\left(-\sum_{i,j}s_{i+j}p_iq_j\right).
\tag{16}
```

Here `P` and `Q` are the actual products of digit intervals, including
prefix and top-digit restrictions. For a fixed `s,p`, the `q` sum
factors into explicit geometric sums:

```math
\prod_j\left(\sum_{q_j\text{ allowed}}
       \chi_1\left(-q_j\sum_i s_{i+j}p_i\right)\right).
\tag{17}
```

The remaining matrix of linear forms is Hankel, with entry `s_(i+j)`.
The new primitive makes `widehat W_N(s)` affordable, but it does not
sum (16) over all `s,p`. The bounded digit intervals in (17) are
incomplete character sums. Replacing them by unrestricted residue sums
would replace the original input domain and change the problem.

A further exponent improvement would require an exact contraction
algorithm for this coupled family, with a closure property under its
recursive steps and a better leading recurrence coefficient. Merely
evaluating each geometric sum quickly leaves exponentially many
frequency and input tuples in (16). No such contraction theorem is
claimed here.

## 8. Why fixed-frequency rank does not compress every interface

The rank bound (6) applies after fixing `s`. It does not hold for the
single tensor flattening whose rows are `r` and whose columns are
the pairs `(s,r')`.

Indeed, fix `r'=K` and let `s=0,...,K`. All lower-bound tests in (4)
are then one. Up to nonzero column scalings, the resulting square matrix
is

```math
\bigl((\zeta^{-r})^s\bigr)_{0\le r,s\le K}.
\tag{18}
```

Because `K<M` and `zeta` has exact order `M`, these evaluation points
are distinct. The Vandermonde determinant is nonzero in the constructed
field. This flattening therefore has full rank `K+1`.

Likewise, the constant number of stored field elements in (9) is an
arithmetic evaluation invariant, not a constant tensor bond dimension:
the state parameters depend on previous frequency labels. This scoped
distinction prevents converting the oracle theorem into an unsupported
low-rank claim for all labels simultaneously.

There is also an exact rank formula for the two-frequency interaction
created by contracting adjacent fixed-frequency factors. For any set
`I` of distinct residues, put

```math
K_I(s,t)=\sum_{r\in I}\zeta^{(bs-t)r},
\qquad s,t\in\mathbb Z/M\mathbb Z.
\tag{19}
```

The complete Fourier matrices on either side are invertible; multiplying
frequency indices by `b` only permutes their rows because `gcd(b,M)=1`.
Thus (19) is their product with the diagonal mask of `I`, and

```math
\boxed{\mathrm{rank}(K_I)=|I|.}                 \tag{20}
```

In particular, a geometric interval sum with `b` terms has rank `b`
as a two-frequency kernel. Its short evaluation formula does not
provide a channel with only `O(log D)` bits. If `I` is a consecutive
interval of length `ell`, restricting both frequency sets to any sets
of size at least `ell` still gives rank `ell`: the two rectangular
Vandermonde factors have full column and row rank. A genuine improvement
through these kernels must therefore use an additional restriction or
joint contraction property supplied by the multiplication layers.

There is nevertheless a useful exact composition law when a coefficient
is fixed. For intervals `I,J` inside `[0,K]`, a target digit `0<=n<b`,
and an actual convolution coefficient `0<=a<=D(b-1)^2`, consider

```math
\sum_{u=0}^{M-1}K_I(s,u)\chi_{n-a}(u)K_J(u,t).
\tag{21}
```

Orthogonality restricts its terms to `r=bq+n-a mod M`. In the stated
ranges,

```math
-M<bq+n-a-r<M,
```

so this is the exact integer equality `r=bq+n-a`. The admissible `q`
form the interval

```math
Q=J\cap\left[
 \left\lceil\frac{\min I+a-n}{b}\right\rceil,
 \left\lfloor\frac{\max I+a-n}{b}\right\rfloor\right].
\tag{22}
```

Equation (21) equals
`M sum_(q in Q) zeta^(bs(bq+n-a)-tq)`. Its rank in `(s,t)` is exactly
`|Q|<=D`, by the same independent character columns as in (20).
This gives a cheap character-twisted composition when `a` is known.

Keeping the coefficient as an unresolved physical variable changes
the rank again. With `I=J=[0,K]`, `n=0`, and `a=0,...,b-1`, flatten
the result with row `s` and column `(t,a)`. Its rank is exactly `K+1`.
For every `r in [0,K]`, choose `q=ceil(r/b)` and `a=bq-r`; these are
legal and satisfy (22). Varying `t` isolates each surviving `q` by
Fourier orthogonality, so the column space contains all `K+1` independent
characters `s -> zeta^(bsr)`. The reverse inequality follows because
every row function is a combination of those same characters.

Thus the small conditional rank in (21) cannot be substituted uniformly
when `a` is still computed from unresolved input digits. A coupled
evaluator would have to share that coefficient dependence across the
multiplication and carry layers, rather than discard it.

## 9. Reproducible checks

Run:

```sh
python3 tools/check_carry_fourier_elimination.py \
  --output evidence/carry-fourier-elimination.json
```

The checker constructs irreducible binary polynomials and roots of the
required orders, checks the division-free geometric-sum identities,
compares every propagated state with a direct carry-vector computation,
and checks the original sum over coefficients on tiny complete domains.
It includes target digits on both sides of the `D-1` clipping boundary,
zero frequencies, coincident character ratios, composite `M`, and the
full-rank flattening in (18), interval-kernel ranks in (20), and the
conditional composition (21), including restoration of full rank when
its coefficient remains unresolved. It also checks the complete
Fourier-dual factor count in a tiny model. These are finite controls for the stated
uniform proofs, not evidence of a score increase.
