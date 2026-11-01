# A p-typical Witt state space with perfect modular duality

## What is proved

The polynomial algebra on generators of degrees `1,2,4,...` can retain
its binary-partition state count **and** a nondegenerate homogeneous
pairing in characteristic two. The correct integral lattice is the
p-typical Witt subalgebra of symmetric functions, rather than the
independent additive-translation model.

This note proves two facts for every weight `h`:

1. Its Hall Gram determinant is an explicit positive odd integer.
2. The charged and neutral vertex-derivative fields defined below have
   invertible homogeneous vacuum Gram matrices modulo two.

These supply actual duals on a low-growth Fock space. They do **not**
yet supply a crossing decomposition or a factoring improvement. A
direct finite exchange test already fails for these fields at `h=2`
with a common scalar angle through degree 12. Length-dependent angles
and a different extraction remain separate possibilities.

The defining ghost and Witt formulas are documented in the primary
[SageMath Witt implementation](https://github.com/sagemath/sage/blob/36e49d91b0568f07ad4373362c15ce8970587e67/src/sage/combinat/sf/witt.py).
The Artin--Hasse definition and algorithms are documented in its
[p-adic implementation](https://github.com/sagemath/sage/blob/36e49d91b0568f07ad4373362c15ce8970587e67/src/sage/rings/padics/padic_generic_element.pyx).
The determinant and vacuum arguments below are derived here. This work
used substantial OpenAI assistant support and includes no Lean work.

## 1. The integral algebra and the Hall determinant

Let `Lambda_Z` be the ring of symmetric functions, with Hall pairing
for which Schur functions are orthonormal. Its power sums satisfy

```math
\langle p_\lambda,p_\mu\rangle=\delta_{\lambda\mu}z_\lambda,
\qquad z_\lambda=\prod_n n^{m_n}m_n!.
```

Let `w_n` be the integral Witt symmetric functions, characterized by

```math
\prod_{n\ge1}(1-w_n t^n)^{-1}=\sum_{r\ge0}h_r t^r,
\qquad p_n=\sum_{d\mid n}d w_d^{n/d}.
```

Work over `Z_(2)`, the rationals with odd denominator, and put

```math
A=\mathbf Z_{(2)}[w_1,w_2,w_4,\ldots].
```

For the rest of this note write `x_j=w_(2^j)` and `P_j=p_(2^j)`.
Then

```math
P_j=\sum_{i=0}^j2^i x_i^{2^{j-i}},\qquad \deg x_j=2^j. (1)
```

The degree-`h` basis consists of monomials `x^m` with
`sum_j 2^j m_j=h`; its size is the binary partition number `b(h)`.
The Witt addition polynomials show that `A` is a graded Hopf
subalgebra. Equivalently, solve (1) recursively after making every
`P_j` primitive; the resulting Witt coproducts have integral
coefficients.

Order monomials so that the triangular change from `x^m` to `P^m`
has diagonal entry `2^(-sum_j j m_j)`. The Gram matrix in the latter
basis is diagonal with entry
`prod_j (2^j)^(m_j)m_j!`. Therefore the Hall Gram determinant in the
Witt monomial basis is exactly

```math
\det G_h=\prod_{m:\,\sum 2^jm_j=h}
                \prod_j \frac{m_j!}{2^{j m_j}}.       (2)
```

This determinant is integral because each basis vector is an integral
symmetric function. In fact the powers of two cancel completely.
Let `F(q)=prod_j(1-q^(2^j))^(-1)`. Legendre's formula gives

```math
\begin{aligned}
\sum_h q^h\sum_{m:\,\sum2^jm_j=h}\sum_j v_2(m_j!)
 &=F(q)\sum_{j\ge0,k\ge1}\frac{q^{2^{j+k}}}{1-q^{2^{j+k}}}\\
 &=F(q)\sum_{r\ge1}\frac{r q^{2^r}}{1-q^{2^r}}\\
 &=\sum_h q^h\sum_{m:\,\sum2^jm_j=h}\sum_j j m_j.
\end{aligned}                                                        (3)
```

Combining (2)--(3),

```math
\boxed{\det G_h=
 \prod_{m:\,\sum2^jm_j=h}\prod_j\operatorname{oddpart}(m_j!)
 \quad\text{is odd}.}                                  (4)
```

Thus `A_h` is self-dual over `Z_(2)` and over every characteristic-two
field. For example, at `h=2`, in basis `x_0^2,x_1`,

```math
G_2=\begin{pmatrix}2&-1\\-1&1\end{pmatrix}
 \equiv\begin{pmatrix}0&1\\1&1\end{pmatrix}\pmod2.
```

The same determinant argument works for any prime `p`, replacing
`2^j` by `p^j` and `oddpart` by the prime-to-`p` part. We use only
the characteristic-two case below.

### Effective duality

Put `lambda=ceil(log2 b(h))`. The homogeneous state basis, Hall Gram,
and its inverse modulo two can be constructed in `2^O(lambda)` bit
operations. The vacuum bases from Sections 3--4 have the same bound.

Enumerate the binary partitions and expand the Witt generators by (1).
Each homogeneous ghost polynomial has at most `b(h)` monomials.
A coarse induction in (1) bounds the numerator and denominator bit
lengths of all coefficients at weight at most `h` by `O(h log(h+2))`:
raising a weight-`2^i` polynomial to `2^(j-i)` scales its coefficient
logarithms by that factor, and division by `2^j` adds only `j` bits.
An extra polynomial factor in this bit bound would also suffice.
Products of basis generators and the diagonal ghost norms have
polynomial bit length in `h`. Thus dense polynomial arithmetic and
Gram construction take `b(h)^O(1) poly(h)` bit operations.
Exact inversion over `F_2` has the same form.

The binary-partition estimate `log b(h)=Theta((log h)^2)` implies
`poly(h)=2^O(sqrt(lambda))`, absorbed by `2^O(lambda)`. Computing
the coefficients in (8) or (11) uses at most `h b(h)` monomials
across all lower weights and a polynomial number of truncated
multiplications. Unitriangular inversion then produces the actual
homogeneous operator duals within the same bound. This proves an
effective low-growth dual interface; it supplies no exchange identity.

## 2. An integral coherent series and exact scalar commutation

Write

```math
E_2(t)=\exp\left(\sum_{j\ge0}\frac{t^{2^j}}{2^j}\right),
\qquad
\mathcal H(t)=\exp\left(\sum_{j\ge0}\frac{P_j t^{2^j}}{2^j}\right)
             =\prod_{i\ge0}E_2(x_i t^{2^i})
             =\sum_{n\ge0}H_n t^n.                    (5)
```

Every coefficient belongs to `A`. One self-contained integrality
argument uses

```math
E_2(t)=\prod_{n\text{ odd}}(1-t^n)^{-\mu(n)/n}.
```

Taking logarithms proves the identity. For an exponent with odd
denominator, the generalized binomial coefficients are 2-integral:
in a product of consecutive terms with odd step, the number divisible
by `2^k` is at least the corresponding count in the factorial
denominator. Each coefficient of the displayed product uses finitely
many factors.

Let `Gamma_+(t)` be multiplication by `H(t)` and let `Gamma_-(s)`
be its Hall adjoint. The perfect pairing from Section 1 makes every
coefficient of both operators integral on `A`. Over the rational ghost
basis, `Gamma_-(s)` translates `P_j` to `P_j+s^(2^j)`. Hence

```math
\Gamma_-(s)\Gamma_+(t)
 =E_2(st)\Gamma_+(t)\Gamma_-(s).                       (6)
```

This equality descends to characteristic two. Unlike ordinary
independent additive translations on `x_j`, it does not collapse to
one vacuum functional per energy: the Witt coordinates carry
nonlinearly, and the homogeneous pairing is perfect.

The scalar series itself remains a difficulty. Kramer--Miller,
[arXiv:2404.06968](https://arxiv.org/abs/2404.06968), proves that the
Artin--Hasse exponential modulo `p` is transcendental. Schmidt,
[arXiv:2106.06084](https://arxiv.org/abs/2106.06084), proves that its
coefficients `u_n` have nonzero minors
`det[u_(pi-j)]_(1<=i,j<=ell)=prod_(k=1)^ell k!p^k/(pk)!`.
These facts prevent a bounded-degree algebraic scalar replacement;
they do not settle a finite-window contraction construction.

## 3. Precise fields and their homogeneous vacuum maps

In characteristic two, set

```math
V(z)=\Gamma_+(z)\Gamma_-(z^{-1}).
```

For `n` a power of two, define its divided derivative by expansion in
an auxiliary variable:

```math
V_n(z)=[a^n]V(z+a).
```

For a binary partition `lambda` of `h`, take the normally ordered
product of its `V_n`. Normal ordering places every multiplication
operator to the left of every adjoint multiplication operator and
introduces no contraction factor.

More explicitly, for one part the normal-ordered symbol is

```math
G_n=\sum_{r,s\ge0}{r-s\choose n}H_r\otimes H_s,          (7)
```

where `f tensor g` acts as `M_f g^perp`. A term of degrees `r,s`
has `z` exponent `r-s-n`. For a label `lambda`, multiply the symbols
`G_n` in the commutative algebra `A tensor A`; its total exponent is
`degree(output)-degree(input)-h`. Each fixed matrix entry is a finite
sum: the annihilation degree is at most the input energy and the
creation degree is at most the output energy.

The constant vacuum-creation state is `H_lambda=prod_(n in lambda)H_n`.
These form a basis: `H_(2^j)=x_j` plus monomials of larger ordinary
monomial length. Here length means the sum of exponents of the Witt
variables, not weighted degree.

Define

```math
B_n(t)=\sum_{k\ge1}{-k\choose n}H_k t^k,
\qquad U_\lambda=[t^h]\prod_{n\in\lambda}B_n(t).         (8)
```

The homogeneous vacuum pairing of the operator labelled `lambda`
with a state `v in A_h` is `Hall(U_lambda,v)`.

The `U_lambda` are a basis for **every** `h`. The length-one part of
`H_k` is `x_j` only when `k=2^j`, with coefficient one. For two
positive powers of two `m,n`, Lucas' rule gives

```math
{-m\choose n}\equiv[m\le n]\pmod2.
```

Thus the length-one part of `B_n(t)` is
`sum_(m<=n,m power2) x_(log2m)t^m`. Every monomial in the product
in (8) has length at least `ell(lambda)`. To have exactly that length,
each factor must contribute a length-one term of weight at most its
labelled part. The total energy is `h=sum_(n in lambda)n`, so equality
forces every selected weight to equal its labelled part. Therefore

```math
U_\lambda=x^\lambda+
       \{\text{monomials of strictly larger ordinary length}\}.       (9)
```

The change of basis is unitriangular by length. Combined with (4),
this proves invertibility of the actual homogeneous vacuum Gram and
constructs homogeneous dual operator parameters.

## 4. Neutral ratio fields preserve the duality

The charged fields are not ordinary Heisenberg currents. A separate
candidate is the neutral ratio

```math
W(z,a)=:V(z+a)V(z)^{-1}:,
\qquad W_n(z)=[a^n]W(z,a).                             (10)
```

The inverse in this normal-ordered notation means the product of the
inverse creation and annihilation series, followed by normal ordering.
Its symbol is `G_n` from (7) multiplied by
`H(1)^(-1) tensor H(1)^(-1)`, interpreted in the completed graded
algebra. For a label of length `ell`, the symbol is multiplied by
the `ell`th power of this factor.

The creation state remains `H_lambda`. Its vacuum parameter is

```math
U_\lambda^{\rm neu}
 =[t^h]\mathcal H(t)^{-\ell(\lambda)}
                  \prod_{n\in\lambda}B_n(t).          (11)
```

Since `H(t)^(-1)` has constant term one and all other terms have
positive Witt monomial length, the proof of (9) is unchanged. Neutral
fields therefore also have full homogeneous duals at every weight.

Their formal scalar exchange factor, before extracting jets, is

```math
\frac{E_2((w+b)/(z+a))E_2(w/z)}
     {E_2(w/(z+a))E_2((w+b)/z)}.                       (12)
```

Over `Q_2`, its logarithm has `j`th term

```math
\frac{((w+b)^{2^j}-w^{2^j})
       ((z+a)^{-2^j}-z^{-2^j})}{2^j}.                 (13)
```

For fixed positive jet orders `k,l`, coefficients in the **logarithm**
(13) have valuation at least `j-v2(k)-v2(l)` once `2^j>=k,l`.
The logarithm thus displays a potential frequency `h^2` when
`k=l=h` is a power of two. This is not an established necessary
energy cost or an exchange-kernel cutoff: exponentiation and its
factorial denominators must be treated separately. Moreover finite
Laurent terms are not automatically poles at
`z=w`; the source's delta/residue crossing argument must be supplied
anew. No locality or channel bound is inferred from (13).

## 5. Exact finite checks and the unresolved exchange step

`witt_artin_hasse.py` computes the Artin--Hasse coefficients over
the rationals, checks their odd denominators, expands Witt generators
in the ghost basis, computes integral Hall Gram matrices, and computes
the vacuum map (8) exactly over `F_2`.

Both ranks are full for every `1<=h<=16` and for `h=32`. At `h=32`
both equal `202`. The all-size duality statements are proved above;
these computations only verify their implementation.

`witt_exchange_jets.py` builds the full operator matrices at `h=2`
on all states through energy 16. It constructs the dual fields from
the actual vacuum Gram and sets up the source-shaped equations

```math
Q_s^{\rm row}=\sum_{i=0}^L a_i T_{s+i}^{\rm row},\qquad
Q_s^{\rm col}=\sum_{i=0}^L a_i T_{4-s-i}^{\rm col},
\quad -L\le s\le4,
```

requiring equality in both cuts and `sum_s Q_s=E_2`. Every equation
uses an energy at most `L+4`, so the test is exact for `L<=12`.
No common scalar angle solves these equations, for either the charged
or neutral fields. This is a finite rejection of that interface, not
an all-window impossibility theorem.

The neutral energy tensors happen to vanish at every tested energy
from 7 through 16. The charged tensors have further nonzero entries.
Neither observation is extrapolated. Partitioning by the two operator
lengths, choosing different regularizers, changing the fields, or using
a different crossing extraction requires a new exact identity.

### The necessary rank budget closes the small-window shortcut

A full-entropy crossing on `d=b(h)` labels has row rank `d^2`.
Suppose there are only `poly(h)` pieces and each row factor uses
at most `poly(h)` copies of the state basis through energy `H`. Rank
subadditivity then requires

```math
d^2\le\operatorname{poly}(h)b(H).                     (14)
```

For binary partitions,
`log b(h^c)=c^2 log b(h)+o(log b(h))` for fixed `c>0`.
Consequently `H=h^(c+o(1))` in this interface necessarily has
`c>=sqrt(2)`. In particular, a polynomial-piece crossing whose
single-state channels have energy at most `O(h polylog h)` is impossible, irrespective
of the exchange identity. The small state count cannot be substituted
for the missing rank factorization.

At such a necessary larger window, dividing the energy by a fixed
constant has no effect on its leading logarithmic state count.
The naive two-cut bound `b(H/2)` in each cut therefore has
`rho=c^2/2>=1`. This is a barrier to the dimension-only estimate,
not a lower bound on all actual ranks. A useful construction would
need stronger actual rank compression, a new asymmetric factor
profile, more efficient recombination, or a different interface.
The perfect duals proved here do not bypass (14).

```sh
python3 tools/witt_artin_hasse.py --h 32 --gram --output evidence/witt-gram-h32.json
python3 tools/witt_exchange_jets.py --window 16 --output evidence/witt-exchange-h2-window16.json
python3 tools/witt_exchange_jets.py --window 16 --neutral --output evidence/witt-neutral-exchange-h2-window16.json
```
