# One 2-adic completion determines the integral divisor coefficient

## Result and scope

The [faithful-character interface](sparse-faithful-character.md) uses the
whole algebra F2[T]/Phi_u(T). Its field factors cannot generally be
discarded modulo two. A different input contract suffices for **all**
integers: one chosen field factor, lifted to explicitly polynomial
2-adic precision, determines the whole bounded integral coefficient by
deterministic polynomial-time lattice decoding.

This is a complete **conditional** factoring reduction. It supplies no
fast coefficient algorithm at the original composite index. The complete
proven factoring exponent remains 1/2. No semiprime, squarefree,
one-inversion-orbit or favorable residue-degree promise is imposed.

The change is precision, rather than a claim of injectivity of one
mod-two field. It connects the global weight-one interface to a single
unramified coefficient ring, with all acquisition and reconstruction
costs charged. Those costs are polynomial in the numerical level q,
not in log(q).

The parent already gives a complete precision-one **single-field** query
contract by changing the primitive generator and requesting every Frobenius
orbit. Here one fixed character at the fixed canonical generator suffices,
with increased precision instead of those additional queries. Neither
change of encoding supplies an exponent gain. The potential use is an
evaluator that works in one connected unramified ring but cannot cheaply
vary its character or contract the whole product algebra.

## 1. The integral target and its coefficient bound

Let q>=7 be prime, q=3 mod4, u=(q-1)/2 and q not divide N. Acquire a
primitive generator g modulo q. Put zeta=zeta_u and chi(g)=-zeta. Because
u is odd, chi has order 2u and is an odd primitive character of conductor q.
The positive-index coefficient of the existing
[weight-one Eisenstein series](weight-one-divisor-coefficients.md) is

```math
a_N=\sum_{r\mid N}\chi(r)\in A=\mathbb Z[T]/\Phi_u(T).
```

Its constant Fourier term is not requested. In particular, a possible
denominator in that constant term supplies no division in this input.
Let n=phi(u). Use the n-dimensional integral power basis of A.
Compute all u monomial remainders T^j mod Phi_u and let C_u>=1 be their
largest absolute coefficient. If D>=tau(N), every coefficient of a_N
has absolute value at most B=D*C_u: it is a sum of at most D signed
monomial remainders.

The acquisition is polynomial in u. Since the roots of Phi_u have
absolute value one, its coefficients have absolute value at most 2^n.
Monic long division of the u monomials gives the coarse bound
log2(C_u)=O(u(n+1)). Thus computing C_u exactly does not hide exponentially
many bits. No divisor enumeration is used to obtain B.

## 2. Acquire the chosen completion

Factor Phi_u modulo two by the existing deterministic Berlekamp procedure.
Choose **any** irreducible factor f of degree d. Odd u makes Phi_u
squarefree modulo two. For any requested K, acquire its monic Hensel
factor F modulo 2^K, and the complementary monic factor G, so that

```math
\Phi_u=FG\pmod {2^K},\qquad F\bmod2=f.
```

Here is an explicit polynomial algorithm. If F,G are correct modulo 2^k,
form E=(Phi_u-FG)/2^k modulo two. With g0=(Phi_u mod2)/f, solve

```math
\delta F=E g_0^{-1}\bmod f,\qquad
\delta G=(E-\delta F g_0)/f
```

over F2, then add 2^k deltaF and 2^k deltaG. Their degrees stay below d
and n-d. The inverse exists by squarefreeness. Repeat for k=1,...,K-1.
All polynomial degrees are at most n, all coefficient words have O(K)
bits, and there are K-1 lifts. This deliberately elementary implementation
has polynomial bit cost in u,K.

The provider returns just a_N in

```math
R_K=(\mathbb Z/2^K\mathbb Z)[T]/F(T),
```

with chi(g)=-T. Its complete output has dK bits. It must return the
coefficient at N itself; factoring N before making this request is not
part of the contract.

## 3. Kernel separation by an integer norm

Let psi:A->R_K be reduction and Lambda its kernel in Z^n. It is
surjective, since the first d monomials give the power basis of R_K.
Therefore [Z^n:Lambda]=2^(Kd).

For nonzero W in Lambda, the principal ideal (W) lies inside Lambda.
Its finite index is the nonzero integer |Norm(W)|, since Phi_u is
irreducible over Q. Consequently

```math
2^{Kd}\mid\operatorname{Norm}(W),\qquad
2^{Kd}\le |\operatorname{Norm}(W)|
 \le (\sqrt n\|W\|_2)^n.
```

The last inequality evaluates W at the n complex primitive u-th roots
and uses Cauchy--Schwarz at each root. Hence

```math
\lambda_1(\Lambda)\ge 2^{Kd/n}/\sqrt n.                 \tag{1}
```

This index argument needs no unknown prime ideal, class number or
algebraic-number oracle. The acquired finite factor identity suffices.

## 4. Exact polynomial-time reconstruction

Compute r_i=T^i mod F modulo 2^K for 0<=i<n. For i<d these are the
standard coordinate vectors. An explicit row basis of Lambda consists of

```math
2^K e_i\ (i<d),\qquad e_i-\sum_{j<d}(r_i)_j e_j\ (d\le i<n).
```

For a supplied residue y, put a0=(y0,...,y_(d-1),0,...,0). Its preimage is
exactly a0+Lambda. Apply the existing exact delta=3/4 LLL reduction and
backward nearest-plane construction from
[the boundary decoder](lattice-decoded-boundary-constructor.md#4-exact-nearest-plane-recovery-of-the-bounded-input),
retaining its classical Lenstra--Lenstra--Lovasz attribution and polynomial
bit-cost theorem.

Write h=ceil(log2 B), and choose

```math
K=\left\lfloor\frac{n(h+2n+2)}d\right\rfloor+1.       \tag{2}
```

LLL's reduction inequalities and (1) give

```math
\min_i\|b_i^*\|_2
 \ge2^{-(n-1)/2}\lambda_1(\Lambda)
 >2B\sqrt n.
```

If a_N has its promised coefficient bound, its norm is at most B sqrt(n).
Run backward rounding on the target -a0. Compared with the true lattice
vector a_N-a0, the error is -a_N. At every step its projection has size
less than half the relevant Gram--Schmidt length; higher coefficients
have already been removed and lower basis vectors are orthogonal to that
Gram--Schmidt vector. Every rounding is therefore exact. Adding a0
recovers the entire integral coefficient. Check its coefficient bounds
and its supplied local residue before accepting.

The basis has dimension n and O(K)-bit entries. Formula (2), exact LLL,
rational Gram--Schmidt, rounding, and the final checks have polynomial
bit cost in u,log(D+1). This is not an exact-CVP oracle or a sparse-support
enumeration. The generous precision avoids any optimization of a decoder
constant.

## 5. Complete all-input transfer and costs

Reduce the recovered integral coefficient modulo two. The sign in
chi(g)=-zeta disappears, yielding the full faithful observation

```math
P_N(T)\bmod\Phi_u(T),\qquad
P_N=\sum_{r\mid N}T^{\log_g(r)\bmod u}.
```

The [existing sparse decoder and lower-sieve level acquisition](sparse-faithful-character.md)
apply unchanged. For a b-bit input, first complete trial division through
2^t. If it has not returned a factor, D=2^ceil(b/t) bounds tau(N), with
prime multiplicities included. Set S=2D and

```math
L=(b+1)S(S-1)/2+S(b+2)+3.
```

Acquire L+b+1 distinct q of
size poly(b,D), q=3 mod4 and leastprime(u)>2D; omit q3, and terminate on a
proper level gcd. The existing polynomial acquisition theorem supplies
the extra b+1 levels needed for exceptions. Its analytic hypotheses and
source attribution are retained; this note supplies no new prime theorem.

Each reconstructed faithful observation is decoded into the signed
divisor residues. The existing maximal-support selection, symmetric
polynomial CRT, height bound and least-root extraction handle arbitrary
inputs, squares and repeated prime factors. No prime factor count or
semiprime test is needed.

Suppose one uniform deterministic provider supplies the selected local
coefficient above in F(b)*poly(b,q,K) bit operations, including its complete
output and any internal precision losses. The bound is uniform for all
subinstances of the original bit size b; alternatively take a nondecreasing
F. All preprocessing, level/field
acquisition, Hensel lifting, lattice decoding, residue conversion and
integer recovery then cost

```math
(2^t+F(b)+1)\operatorname{poly}(b,2^{\lceil b/t\rceil}).
```

For t=ceil(sqrt b), this is (F(b)+1)exp(O(sqrt b)). Every fixed positive
logarithmic-saving exponent of that provider transfers to a complete
factoring machine. At most b successive least-factor extractions give the
entire prime factorization, within the absorbed polynomial factor.
A provider at factored indices, at precision one only,
or with uncharged exponential dependence on q or K does not meet this
contract.

## Controls, limits and credit

The q47 sparse and q503 reciprocal collisions in the parent both separate
already modulo four in the chosen lifted factor. They remain collisions
modulo two. The latter words are reciprocal sparse controls, not claimed
actual divisor data. This precision change does not withdraw either
earlier counterexample or the stronger precision-one semiprime result.

```sh
make verify-cyclotomic-completion
```

The [implementation](../tools/cyclotomic_completion.py),
[checker](../tools/check_cyclotomic_completion.py) and
[receipt](../evidence/one-cyclotomic-completion.json) record 681 actual
supplied divisor coefficients across every field at q7,11,19,31,43,47;
97 bounded-vector controls; nine independently constructed integral norm
and lattice-determinant controls; and eight invalid-input guards. Four
complete supplied-data CRT factorizations (8,9,11,15) use one completion
at each of six eligible levels, reaching global dimension 40/residue
degree 20. The fixture provider explicitly lists divisors and is not a
fast coefficient engine. The q503 control checks separation modulo four,
without running its 250-dimensional lattice decoder.

An optional private 48-dimensional exact Python LLL pilot at order105
was deliberately stopped after about seven minutes before output; its
exit143 is recorded as incomplete. The permanent control verifies the
nontrivial monomial coefficient bound at that order, without claiming
its lattice replay. No larger field or reconstruction sweep is proposed.
Focused checks pass; full default verification is pending at this edit.
Further decoder optimization is not the next factoring mechanism: the
open task is the stated original-index evaluator.

The construction reuses the existing cyclotomic/Berlekamp tools and exact
LLL implementation. Norm separation and bounded lattice decoding are
classical techniques; the new contribution is the fully charged use of
one chosen completion in this complete coefficient interface. Developed
by Sol worker 2 with substantial OpenAI Codex assistance. No new Lean or
independent mathematical review is claimed.
