# A fixed weight-one coefficient target for global divisor residues

The global sparse-residue interface has a modular-form realization with
**weight fixed at one**. Odd Dirichlet-character coefficients determine
the signed divisor residue polynomial by exact Fourier inversion. Only
logarithmic 2-adic coefficient precision is needed. Its support has at
most twice as many terms as the original divisor set, so
the same trial-budget transfer recovers the smallest prime factor.

This supplies a second coefficient target alongside the
[level-one growing-weight moments](divisor-moment-recovery.md). Here the
level varies over small primes and the weight stays one. A uniform
fast algorithm for the coefficient at the original composite index N
remains unproved. The established factoring score stays `1/2`.

## 1. The coefficient family

Let q be an odd prime, m=q-1, and g a primitive generator modulo q.
In Q(zeta_m), define the character chi_j by

```math
\chi_j(g^e)=\zeta_m^{je},\qquad
\chi_j(a)=0\text{ if }q\mid a.
```

For odd j, 1<=j<m, it is an odd nontrivial primitive character of
conductor q. Use the Eisenstein series with trivial first character:

```math
\mathcal E_{1,1,\chi_j}(z)
 =-\frac{B_{1,\chi_j}}2
  +\sum_{n\ge1}\left(\sum_{d\mid n}\chi_j(d)\right)z^n.
                                                               \tag{1}
```

This is a holomorphic modular form of weight one, level q and character
chi_j. The coefficient formula and parity condition are the specialization
of William Stein's [*Modular Forms, A Computational Approach*, Chapter 5,
formula (4) and Theorem 5.8](https://wstein.org/books/modform/modform/eisenstein.html)
with k=1 and trivial first character. The algorithm requests only

```math
M_j(N)=[z^N]\mathcal E_{1,1,\chi_j}(z)
      =\sum_{d\mid N}\chi_j(d).                         \tag{2}
```

The coefficient index N need not be coprime to the level. Character
values at multiples of q are zero, as already specified.

## 2. Odd characters determine the signed cyclic image

For a in Z/qZ let c_a count divisors of N congruent to a. The nonzero
residues are g^e, 0<=e<m, and (2) is their multiplicative Fourier transform:

```math
M_j(N)=\sum_{e=0}^{m-1}c_{g^e}\zeta_m^{je}.
```

The vector b_a=c_a-c_(-a) has b_0=0. Translation e -> e+m/2 multiplies
the j-th Fourier coefficient by (-1)^j. Therefore

```math
b_{g^e}=\frac2m\sum_{\substack{1\le j<m\\j\text{ odd}}}
                   M_j(N)\zeta_m^{-je}.                 \tag{3}
```

This is an exact cyclotomic calculation whose result is an ordinary
integer. Reducing these integers modulo two supplies the complete
cyclic image of the Laurent polynomial

```math
H_N(Z)=\sum_{d\mid N}(Z^d-Z^{-d})
\quad\text{in }\mathbb F_2[Z,Z^{-1}]/(Z^q-1).           \tag{4}
```

If q divides a divisor, its two contributions to position zero cancel;
that is why no coprimality restriction on N was needed. At q=2 all
opposite exponents coincide, so the image is simply zero and no character
query is needed.

As a formal Laurent polynomial, (4) has exactly 2*tau(N) distinct
nonzero terms, even in characteristic two: the positive and negative
integer exponents are distinct. Cancellation in individual cyclic
images is handled by sparse interpolation rather than treated as an
empty divisor set.

## 3. Signed exponent recovery and the complete transfer

After trial division through 2^t without a factor, use
D=2^ceil(b/t)>=tau(N) and S=2D. All signed exponents have absolute
value at most N, so a collision difference is less than 2^(b+1).
The [global sparse-residue proof](sparse-divisor-residues.md) applies with
collision bit bound b+1, sparsity cap S and

```math
L=(b+1)S(S-1)/2+S(b+2)+3                               \tag{5}
```

prime moduli. Maximal support identifies the primes where all signed
exponents are separated. Pair the nonzero positions a,-a and square one
representative from each pair. Distinct pairs have distinct squares modulo
q; the squared positions are exactly the residues of the positive divisor
squares. Their monic product reduces

```math
P_N(Y)=\prod_{d\mid N}(Y-d^2).                         \tag{6}
```

Its coefficient height is at most (N^2+1)^D<=2^(D(2b+1)), and the product
of retained primes exceeds twice that bound. Signed coefficientwise CRT
recovers the same degree-tau(N) polynomial as the moment interface.
Remove Y-1, isolate the least root in [2,N^2], and take its integer square
root to obtain the smallest prime factor. This covers prime powers,
squares and arbitrary nonsquarefree inputs without separate promises.
The paired-residue simplification was suggested in the
[coordination discussion](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19324705).

Suppose a uniform deterministic algorithm supplies each (2), for the
characters specified above, in F(b)poly(b,q) bit operations, with
F(b)>=1 and fixed polynomial degrees. Equations (3)--(6) give complete
factoring cost

```math
2^t\operatorname{poly}(b)+F(b)\operatorname{poly}(b,2^{\lceil b/t\rceil}).
                                                               \tag{7}
```

Indeed all q are polynomial in b,D, there are O(q) characters per
level, and the exact inverse transform has polynomial cost in b,q.
Choosing t for any fixed positive rational beta as in the trial-budget
note preserves F(b)=2^O(b/log^beta(b+2)). For beta<=1 all levels and
the number of requested weight-one coefficients are polynomial in b.
For larger fixed beta, their total overhead is still absorbed into
the target bound. No stronger factoring exponent is proved without
the coefficient algorithm in this hypothesis.

Only inputs that pass the chosen small-factor trial phase need this
engine. It may use the certified promise that every prime factor of N
exceeds 2^t; the divisor and level bounds above already exploit that promise.

The same conclusion follows from coefficient residues at the small
precision in Section 4.1; full exact coefficients are a stronger input
than the recovery needs.

## 4. Encoding and arithmetic are charged

The character input is (q,g,j), and output is an exact integer
coefficient vector in the power basis of Q(zeta_(q-1)), reduced modulo
the cyclotomic polynomial. The degree is at most q-1. Since each output
is a sum of at most tau(N)<=N roots of unity, a polynomial bound in q,b
on its encoded bit length suffices. For example the elementary
coefficient bound on cyclotomic polynomials and monic reduction give
polynomial-in-q heights for reduced root powers; summing at most N
terms adds at most b bits.

All auxiliary work has polynomial bit cost in q,b. Factor q-1 by small
trial division, find g by testing candidates against its distinct prime
factors, construct the cyclotomic polynomial by exact division of
X^(q-1)-1, and reduce every product in its power basis. These procedures
do not require factoring N, a large-field primitive element or numerical
approximations. The division by m in (3) is exact integer division after
the nonconstant coordinates have cancelled.

There is a baseline coefficient engine at the already established
half-rate bound: factor N recursively, then use

```math
M_j(N)=\prod_{r^e\parallel N}\left(\sum_{h=0}^e\chi_j(r)^h\right).
```

This has uniform polynomial overhead in b,q. The open target is a
faster composite-index coefficient engine. Computing coefficients only
at prime indices does not meet it: at a prime r, (2) is just 1+chi_j(r).
Algorithms requiring the prime factorization of their coefficient index
also do not supply the hypothesis. Constructing an initial q-expansion
through index N costs exponentially in b and is likewise insufficient.

One concrete software path was inspected at Sage commit
`36e49d91b0568f07ad4373362c15ce8970587e67`: the general Eisenstein
coefficient method explicitly obtains `divisors(m)` before summing the
character values. See the [pinned source, lines 3150--3192](https://github.com/sagemath/sage/blob/36e49d91b0568f07ad4373362c15ce8970587e67/src/sage/modular/modform/element.py#L3150-L3192)
(blob `6a3318347ea27971b5ffca25722396d7156a3c7b`). Its divisor
enumeration remains part of the cost here. This is a scoped source
inspection, not an audit of every available modular-form algorithm.

### 4.1 Only v2(q-1) coefficient bits are needed

Write m=q-1=2^v u with u odd and v>=1. For one residue position let

```math
T_e=\sum_{j\text{ odd}}M_j(N)\zeta_m^{-je}.
```

Equation (3) says exactly T_e=2^(v-1) u b_(g^e). In the standard
cyclotomic power basis this is an ordinary integer, so all its
nonconstant coordinates are zero. Therefore

```math
b_{g^e}\bmod2
 =\left(\frac{T_e\bmod2^v}{2^{v-1}}\right)\bmod2.       \tag{8}
```

The odd factor u reduces to one modulo two. Thus residues of the
M_j coefficients modulo 2^v suffice: polynomial multiplication and
monic reduction modulo Phi_m commute with coefficient reduction. The
quotient ring may have nilpotents; its free power basis still gives
canonical coordinates, and no nonconstant element is inverted. The
division by 2^(v-1) is performed on the constant coordinate, whose
required divisibility follows from (3).

Consequently the engine in (7) may supply only coefficient vectors in

```math
(\mathbb Z/2^v\mathbb Z)[T]/(\Phi_{q-1}(T)),
\qquad v=v_2(q-1)\le\log_2q.
```

Their encoded length is O(q log q), independent of the large coefficient
index except for the time to compute them. When q=3 mod 4, even precision
modulo two suffices. No restriction to primes in that progression is
needed by the transfer.

For a coefficient engine requiring an integral modular form at the
constant term too, use G_j=2q E_(1,1,chi_j). For the nontrivial character,
q B_(1,chi_j)=sum_(a=1)^(q-1) a chi_j(a) is integral, so every coefficient
of G_j lies in Z[zeta_m]. Its N-th coefficient modulo 2^(v+1) determines
M_j(N) modulo 2^v by division by two and inversion of the odd q. Thus
full-form normalization costs one additional precision bit. Both versions
retain the same conditional runtime; neither coefficient engine is supplied.

## 5. Exact controls

Run `make verify-odd-character-residues`, or

```sh
python3 tools/check_odd_character_residues.py > evidence/odd-character-residues.json
```

The [conversion](../tools/odd_character_residues.py) implements cyclotomic
arithmetic and the full odd-character inverse transform. The
[checker](../tools/check_odd_character_residues.py) supplies coefficients
by summing characters over explicit divisor lists and compares the
result with an independent signed integer residue count. Character
multiplication and odd parity are checked exactly, including coefficient
indices divisible by the prime level. Its large inputs use the same
independently trial-proved prime factors as the preceding controls.

The signed option in the
[sparse recovery core](../tools/sparse_divisor_residues.py) separately
recovers the full polynomial (6) from supplied cyclic residues, comparing
it with an independent product over every squared divisor.
This includes completely cancelled images and one 60-bit recovery.
Exact counts and source scope are in the
[receipt](../evidence/odd-character-residues.json). These are component
checks of the reduction; no fast weight-one coefficient provider or
new uniform factoring bound is implemented.

The receipt records 1,040 exact character conversions, 5,850 supplied exact
coefficients and 184 completely zero signed images. An additional 520
finite-precision conversions use 2,925 coefficient residues at precisions
one, two and four bits and agree with independent integer residue counts.
The separate signed
recovery controls check ten full support polynomials, including the
60-bit case through 2,207 supplied cyclic queries of degree at most 19,463.

The ordinary argument has not received independent mathematical review.
This work used substantial OpenAI Codex assistance and no Lean.
