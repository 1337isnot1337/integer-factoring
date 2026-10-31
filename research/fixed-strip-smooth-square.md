# A fixed zero-free strip gives a polynomial-bit universal smooth square

Under the explicit analytic inputs in
[zero-free class generators](zero-free-class-generators.md), the
arithmetic representation gate can be made entirely polynomial in the
discriminant bit length. More precisely, those inputs give a split-prime
factor base of polynomial size and a word of length O(b/log b) for every
class. One fixed smooth square of polynomial bit length then represents
every class in the square subgroup.

This is a conditional **existence and reduction theorem**. It does not
find the short words or solve the prescribed square-norm decision problem.
Combined with the existing decision reduction and early-stop extraction,
it shows that a polynomial-time exact decision algorithm for the specified
variable-cap squares would give polynomial-time integer factoring under
the analytic inputs. No such decision algorithm is supplied here. The
external zero-free assertions also retain the unverified status stated in
the linked note.

The spectral mechanism is classical. Jao, Miller, and Venkatesan prove
class-group Cayley expansion under GRH in
[*Expander graphs based on GRH with an application to elliptic curve
cryptography*, arXiv:0811.0647v2](https://arxiv.org/pdf/0811.0647v2),
Theorem 1.1, Remark 1.2(a), and Corollary 1.3, pages 3--4; their
character-sum proof is on pages 9--10. Their theorem covers narrow ray
class groups and their quotients. The argument below adapts the same
character-sum mechanism from the GRH width 1/2 to a fixed width eta,
uses smoothing to keep the error uniform, and spells out the consequence
for this repository's fixed-square representation gate. It makes no
novelty claim for the expander or mixing argument.

## 1. Analytic inputs and conductor dependence

Let Delta=f^2 D_K<0 be an arbitrary quadratic order discriminant, and let
G=Pic(O_f). Put b=bits(|Delta|), so b>=2. Fix a rational constant
0<eta<=1/2. Assume the zero-free input and effective uniform smoothed
explicit-formula estimates described in
[zero-free class generators](zero-free-class-generators.md).

Each character chi of G is a finite-order Hecke character of K whose
conductor divides f O_K. Its primitive analytic conductor satisfies

```math
q_\chi=|D_K|N\mathfrak f_\chi\le |D_K|f^2=|\Delta|.       \tag{1}
```

The full ring class field is one of the extensions allowed in the
zero-free input. Its Dedekind zeta function factors into the Hecke
L-functions of these characters. Away from the principal pole at one,
a zero of any factor cannot be cancelled by another factor's pole.
Standard nonvanishing at one handles that point. Thus the same fixed
strip applies to every character of G.

Fix a computable nonnegative nonzero smooth function phi supported in
(1,2), and put c_phi=integral_1^2 phi(t) dt. The Hecke version of the
smoothed explicit formula gives, uniformly in chi,

```math
\sum_{\mathfrak p}\sum_{j\ge1}
 \chi(\mathfrak p)^j\log N\mathfrak p\,
 \phi((N\mathfrak p)^j/x)
 =\mathbf1_{\chi=1}c_\phi x
   +O_{\phi,\eta}\!\left(bx^{1-\eta}\right).              \tag{2}
```

At conductor primes use the usual omitted Euler factors. The effective
implied constant must be independent of the order and the character.
This follows by the same Mellin-transform argument used for the
Dedekind formula: the unit-height zero count for a degree-two Hecke
L-function is O(log q_chi+log(|t|+2)), and (1) bounds its conductor by
the original order discriminant. The use of a potentially huge class
field to infer zero exclusion does not enter this conductor estimate.

As in the generator note, the fixed strip and the effective uniform
explicit-formula constants are analytic inputs. A field-dependent
unspecified constant would not justify the following uniform algorithmic
consequences.

## 2. A positive split-prime walk with small nontrivial eigenvalues

For each unramified split rational prime p, choose one prime-form class
g_p; its conjugate has class g_p^-1. Define

```math
A_\chi(x)=
 \sum_{\substack{x<p<2x\\p\nmid\Delta\\p\ \mathrm{split}}}
   (\log p)\phi(p/x)
   \left(\chi(g_p)+\chi(g_p)^{-1}\right).                 \tag{3}
```

Deleting the terms in (2) that do not occur in (3) costs

```math
O_\phi(\sqrt{x}\log x+\log|\Delta|).
```

Indeed, all higher prime-ideal powers have underlying norm at most
sqrt(2x). The degree-two prime ideals over inert rational primes also
lie over primes at most sqrt(2x). Their total absolute contribution has
the displayed elementary bound. At rational primes dividing Delta,
at most one norm power for each prime ideal lies in the support interval
(x,2x); summing their logarithms costs O(log|Delta|). Consequently

```math
A_\chi(x)=\mathbf1_{\chi=1}c_\phi x+
 O_{\phi,\eta}\!\left(bx^{1-\eta}+\sqrt{x}\log x+b\right). \tag{4}
```

Choose any fixed integer B>1/eta and set

```math
\sigma=B\eta-1>0,\qquad
x=(C_{\eta,B}b)^B,\qquad Y=2x,                            \tag{5}
```

where C_{eta,B} is a sufficiently large effectively computable integer.
The main relative error in (4) is O(b x^-eta)=O(b^-sigma).
The other terms are no larger on this scale: because eta<=1/2,

```math
B/2-\sigma=1+B(1/2-\eta)\ge1.
```

Increasing the fixed constant in (5) therefore makes A_1(x)>0 and

```math
\left|\frac{A_\chi(x)}{A_1(x)}\right|
 \le \frac{1}{4b^\sigma}
 \quad\text{for every nontrivial character }\chi,\ b\ge2. \tag{6}
```

This choice is effective from a uniform constant in (4), rather than a
separate threshold for each field. For example, after dividing the error
by x b^-sigma, the three terms are bounded by constants times
C^(-B eta), C^(-B/2)(1+log C), and C^-B. These tend effectively to zero
as the integer C increases. Choosing their sum below c_phi/16 proves (6),
with room for the normalization by A_1.

Normalize the positive weights in (3) to a probability measure mu on G,
giving the two orientations of each split prime equal weight. Its
character eigenvalues are exactly A_chi/A_1. Coincident classes and
self-inverse prime classes simply have their weights added. In
particular, the walk is symmetric and (6) is an absolute eigenvalue
bound, including possible negative eigenvalues.

The measure is used to prove support and word length. The eventual
algorithm does not need to evaluate phi or these real weights. It
enumerates every unramified split prime p<=Y, so its factor base contains
the support of mu.

## 3. An explicit logarithmic word-length bound

Let ell=floor(log_2 b), and choose the integer

```math
T=\left\lceil\frac{b+1}{2+\sigma\ell}\right\rceil.        \tag{7}
```

The elementary reduced-form count gives |G|<=|Delta|<2^b. Fourier
inversion and (6) imply, for every h in G,

```math
\left|\mu^{*T}(h)-\frac1{|G|}\right|
 \le \left(\frac1{4b^\sigma}\right)^T
 \le 2^{-(2+\sigma\ell)T}
 \le 2^{-(b+1)}
 <\frac1{2|G|}.                                         \tag{8}
```

Thus mu^{*T}(h)>0. Every class is a product of T supported prime-form
steps, with either orientation allowed. Combining repeated steps gives
a signed word e satisfying

```math
h=\prod_{p\le Y}g_p^{e_p},\qquad
\sum_p|e_p|\le T=O_{\eta,B}(b/\log b).                   \tag{9}
```

This simultaneously proves that the chosen prime forms generate the
full class group. It is stronger than generation alone, but it is still
an existence statement for a target word. Rapid mixing or a positive
transition probability does not find a path to a specified target in
polynomial time.

## 4. Polynomial-size fixed-square and relation data

Let P be all unramified split rational primes at most Y, let m=|P|,
and fix the cap T from (7) before choosing any relation target. Set

```math
M=\prod_{p\in P}p^T.
```

Every h in G^2 has a square root g in G. Applying (9) to g gives
h=product g_p^(2e_p) with |e_p|<=T. The algebraic construction in
[universal smooth squares](universal-smooth-square.md) therefore shows
that every such class represents this same M^2, with nonprimitive
representations allowed. This conclusion uses neither a smooth-leading-
coefficient density estimate nor a radius R.

The prime count satisfies m<=Y=O_{eta,B}(b^B). Since T=O(b/log b)
and log Y=O(log b), the bit length L of M^2 satisfies

```math
L\le 1+2Tm\log_2Y
  =O_{\eta,B}(b^{B+1}).                                 \tag{10}
```

The prime list and its splitting roots can be constructed deterministically
in polynomial time for fixed eta and B, for example by a sieve and finite
root search modulo 4p. Primes dividing Delta are detected by gcd and
omitted. For Delta=-4N, a proper gcd with N may be returned directly.
No factorization of Delta or class number is required.

Put K=Tm+1. Each target g_i^(2K) lies in G^2 and hence represents M^2.
Any returned representation normalizes to an even word v with
sum|v_i|<=2Tm<2K. The rows 2K e_i-v give the strictly diagonally dominant
relation matrix from the universal-square note. The determinant has a
polynomial bit bound, for example O(m(log K+log m)). Thus the factor
base, fixed norm, targets, relation entries, and annihilator all have
polynomial encoded size under the analytic inputs.

This statement does not compress an unknown witness into a usable
algorithm. The representation must actually be found. A prescribed long
word for g_i^(2K) cannot merely be padded into M^2, because its exponent
2K exceeds the fixed allowed exponent 2T.

## 5. The exact remaining decision problem

The [square decision reduction](smooth-square-decision-reduction.md)
uses at most O(m log(T+1)) exact queries to recover one representation.
Every queried norm has the form

```math
\left(\prod_{p\in P}p^{t_p}\right)^2,
\qquad 0\le t_p\le T,
```

with the same supplied split-prime list and roots. The target is always
in G^2, but the variable-cap instance may have no representation. These
negative answers are essential. An oracle restricted to the original
universal YES instances gives no search information.

If that exact decision problem admits a uniform algorithm polynomial in
b, L, and m, all m dominant rows can be found using
O(m^2 log(T+1)) such queries, plus polynomial work. The
[early-stop class-group routine](early-stop-class-factoring.md) then
finds a proper factor for Delta=-4N in polynomial time from the complete
generators and the relation determinant. Even inputs and perfect powers
can be handled separately; recursive splitting yields complete factoring.

**Conditional algorithmic gate.** A fixed effective zero-free strip with
the uniform smoothed estimates above, together with a polynomial-time
exact variable-cap smooth-square decision algorithm, implies a uniform
deterministic polynomial-time integer factoring algorithm.

Neither additional input is silently certified here. In particular,
supplied-norm primitive orientation search in
[factored norm search](factored-norm-search.md) does not implement the
variable-cap nonprimitive oracle: it searches only the two full
orientations at each active prime, whereas each cap here permits all
intermediate even exponents. Its dependence on the number of distinct
primes also cannot be hidden inside the discriminant bit length.

## 6. Dependence on the zero-free width

The polynomial statements treat eta as a fixed positive constant. Their
degrees deteriorate as eta decreases. The least integer choice B>1/eta
gives the following representative bounds for the **norm's bit length**,
not for the total algorithm runtime:

| Fixed width eta | Integer B | sigma=B eta-1 | Upper-bound scale for L |
| --- | ---: | ---: | --- |
| 1/2 | 3 | 1/2 | O(b^4) |
| 1/32 | 33 | 1/32 | O(b^34) |
| 10^-6 | 1,000,001 | 10^-6 | O(b^1,000,002) |

The constants also depend on eta, B, and the effective explicit-formula
constants. These are asymptotic conditional reductions, with no claim of
practical feasibility at the very small stated widths. The stronger
GRH case, the claimed recent fixed strips, and an independently verified
analytic theorem must remain distinct.
