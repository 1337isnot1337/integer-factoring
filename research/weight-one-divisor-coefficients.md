# A fixed weight-one coefficient target for global divisor residues

The global sparse-residue interface has a modular-form realization with
**weight fixed at one**. Odd Dirichlet-character coefficients determine
the signed divisor residue polynomial by exact Fourier inversion. Its
support has at most twice as many terms as the original divisor set, so
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

The receipt records 1,040 character conversions, 5,850 supplied
coefficients and 184 completely zero signed images. The separate signed
recovery controls check ten full support polynomials, including the
60-bit case through 2,207 supplied cyclic queries of degree at most 19,463.

The ordinary argument has not received independent mathematical review.
This work used substantial OpenAI Codex assistance and no Lean.
