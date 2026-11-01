# A partial compiler for original-index theta programs

The [phase-seven identity](q7-cubic-parameter.md) extends to a sound,
partial construction for arbitrary odd prime phases. A successful output
has an ordinary modular/Sturm certificate for every odd original index.
The implementation succeeds at phases 3,7,11,19,23,31. Its specified
polynomial module fails at 5,13,17, with short exact separating covectors.
There is no theorem that it succeeds at every prime congruent to three
modulo four. Its explicit coefficient evaluator remains exponential in the
input bit length, so the complete factoring exponent remains 1/2.

This construction uses standard Tate coordinates, coset theta modularity,
finite linear algebra and Sturm's congruence theorem. The compiler and
source checks were developed with OpenAI Codex assistance. These ordinary
proofs and certificates have not received independent mathematical review;
no new Lean result is claimed.

## 1. The candidate module and output contract

For an odd prime q let A_q be the intrinsic finite etale coefficient algebra
from the [Tate compiler](tate-coordinate-coefficients.md), of degree
h=(q-1)/2. It can be a product of fields. The construction does not replace
it by one field component. Its distinguished basis is 1,c,...,c^(h-1).
Let X(t)=X_q(Delta(t)) be the full Tate-coordinate series.
For 1<=a<=h set

```math
\theta_a(t)=\sum_{\substack{n\in\mathbb Z\\ n\equiv a\pmod q}}t^{n^2}
 =\sum_{\substack{n>0\\ n\equiv\pm a\pmod q}}t^{n^2},
\qquad V_q(t)=\sum_{n>0}(t^{n^2}+t^{qn^2}).
```

Write O for projection onto odd exponents. The compiler seeks

```math
O(X)=\beta O(V_q)+
 O\left(\sum_{|e|\le4}\alpha_e\prod_{a=1}^{h}\theta_a^{e_a}\right),
\qquad \alpha_e,\beta\in A_q.                       \tag{1}
```

There are 1+binomial(h+4,4) scalar columns, including the trace column
and constant monomial. Exact binary elimination treats every native
coordinate separately. The returned program combines their coefficient
vectors, so each monomial is evaluated once.

The three output statuses have distinct meanings:

- `certified_by_ordinary_modular_Sturm_proof`: the actual finite congruence
  covers the bound proved below, and (1) holds at every odd index.
- `bounded_module_excluded`: a finite covector annihilates every scalar
  column but detects the target. There is no full-coordinate identity in
  this specified module. Partial coordinate programs are retained.
- `validation_failed`: a training fit exists, but this particular program
  fails the larger check. This status excludes neither the module nor a
  different program.

Only the first status is accepted by the original-index evaluator. No
finite training fit alone is promoted to an all-index identity.

## 2. An odd-index modular lift of every native coordinate

The finite etale trace pairing on A_q is nondegenerate. Thus every native
binary coordinate is lambda(X)=Tr_Aq(eta X) for some eta in A_q. This
statement uses the full algebra trace, including all field components.

Let zeta be a complex primitive qth root. The h nodal values
zeta^a/(1-zeta^a)^2, indexed by nonzero residues up to sign, are distinct
modulo two: in characteristic two the equation z/(1+z)^2=c has exactly
the pair z,z^-1 as its possible nonzero solutions. They lift the roots of
the intrinsic nodal polynomial. Consequently the real cyclotomic order,
localized away from q, reduces to A_q, and its trace reduces to the algebra
trace. A polynomial lift of eta gives weights w_a whose weighted torsion
sum reduces coordinatewise to Tr_Aq(eta X).

Let L_a be the integral Tate x-coordinate at zeta^a. With the period lattice
Z+tau*Z, the normalized Weierstrass function satisfies

```math
\frac{\wp(a/q;\tau)}{(2\pi i)^2}=L_a(t)+\frac1{12}.
```

This follows directly from the Fourier and quasi-period expansions
[DLMF 23.8.1 and 23.8.5](https://dlmf.nist.gov/23.8).
Its positive coefficient is
sum over d dividing n of d*(zeta^(ad)+zeta^(-ad)-2).
The positive coefficients reduce to the signed odd-divisor data.

Each normalized torsion p-value is a holomorphic weight-two form on
Gamma_1(q). The only common shift relative to L_a is constant in tau.
Odd projection kills it exactly. Therefore

```math
O\left(\sum_a w_a\frac{\wp(a/q;\tau)}{(2\pi i)^2}\right)
```

is an integral-at-two modular lift of the selected odd native coordinate,
including coordinate zero. Its Fourier constant is zero. A fixed odd
power of q clears the remaining denominators. This argument does not
assume that the characteristic-zero weights sum to zero, or that Tr(1)
is nonzero. Holomorphy at cusps follows from the torsion Tate expansions,
as in the phase-seven proof.

The scalar O(V_q) has the lift G(t)+G(t^q), where
G=(-E_2(t)+3E_2(t^2)-2E_2(t^4))/24 is integral of weight two on Gamma_0(4).
It lies at level 4q.

## 3. The subgroup and a smaller uniform congruence bound

Let H(t)=sum over n in Z of t^(n^2), with reduction 1. Homogenize every
monomial of degree d in (1) by H^(4-d). Each homogeneous term is a finite
sum of theta series for the rank-four lattice Z^4, Q=sum of four squares,
with coordinate congruences modulo q. The lattice Gram matrix is I_4,
its level and discriminant are both 1 in the convention of
[Kane--Kim, Proposition 2.3](https://arxiv.org/html/2211.03987v2#S2.SS5).
The resulting forms have weight two, trivial character and group
Gamma_0(4q^2) intersect Gamma_1(q).

Odd projection increases the lower-left divisibility by at most four.
For a matrix gamma=(a,b;c,d), conjugating by tau to tau+1/2 gives

```math
\begin{pmatrix}a+c/2&b+(d-a)/2-c/4\\c&d-c/2\end{pmatrix}.
```

If c is divisible by 16q^2 and a,d are 1 modulo q, this is integral and
belongs to Gamma_0(4q^2) intersect Gamma_1(q). The projected theta forms,
the projected torsion forms and the scalar trace lift therefore share

```math
\mathcal H_q=\Gamma_0(16q^2)\cap\Gamma_1(q).
```

Reduction of the bottom-right entry modulo q maps Gamma_0(16q^2)
surjectively onto (Z/q)^*. Its kernel is this intersection. Thus

```math
[\mathrm{SL}_2(\mathbb Z):\mathcal H_q]
 =16q^2(1+1/2)(1+1/q)(q-1)=24q(q^2-1).
```

Sturm's congruence theorem for a weight-two form integral at a prime above
two consequently gives the sufficient inclusive bound

```math
B(q)=4q(q^2-1).                                    \tag{2}
```

See [Stein, Theorem 9.18](https://wstein.org/books/modform/stein-modform.pdf#page=185)
for the general prime-ideal statement. Clearing a fixed odd denominator
puts our projected forms under its algebraic-integer hypothesis. Their
constant and even coefficients vanish identically, so checking all odd
coefficients through B(q) covers the entire required congruence. We retain
the full SL_2 index; using the even weight to adjoin -I could reduce it,
but is unnecessary here.

This proves a uniform **certification rule**. It does not prove uniform
success of the chosen polynomial module.

## 4. Construction, controls and exact exclusions

[The compiler](../tools/theta_coordinate_program.py) inverts the intrinsic
embedding by h binary columns. It does not enumerate all 2^h coefficient
values. A literal odd-divisor sieve supplies the finite training and
verification prefixes, whose lengths depend only polynomially on q.
The initial training length is min(B(q),max(4096,twice the column count)).
A successful fit is checked through B(q)+128. An unsuccessful fit includes
an exact separating covector, supported on odd indices only.

The successful phases and their sufficient bounds are:

| q | h | B(q) |
| --- | --- | --- |
| 3 | 1 | 96 |
| 7 | 3 | 1,344 |
| 11 | 5 | 5,280 |
| 19 | 9 | 27,360 |
| 23 | 11 | 48,576 |
| 31 | 15 | 119,040 |

At phases 5,13,17 exactly native coordinate zero is missing from the
training span, and the returned covectors prove that failure for all-index
identities of form (1). For q=5 the covector is just coefficient 15:
every allowed monomial and V_5 has coefficient zero there, while the target
is 1. This does not contradict the different proved phase-five program,
which uses products involving the scaled square series at phase five.
The receipt lists the five-entry and thirteen-entry covectors at 13 and 17.
The congruence-class pattern of these nine phases is an observation, not a
uniform spanning or exclusion theorem.

A regression deliberately fits q=31 using only 4,096 coefficients. The
binary rank reaches all 2,048 odd rows, so arbitrary target training data
can fit. This program fails the larger prefix check and is refused as an
evaluator. The final q=31 compilation instead trains through 7,754 and
passes the full bound. This records a concrete failure of underdetermined
training rather than silently treating every fit as an identity.

Related prior work: Paul Monsky's
[2013 author post on level-p characteristic-two modular forms](https://mathoverflow.net/questions/121506/level-p-characteristic-2-modular-forms-and-thetas)
uses these same unary-square coset series, relates their fraction field to
an Igusa modular-function field, and gives explicit rescaled formulas for
small phases. His post distinguishes proved theta-ring assertions from
conjectural modular-form-ring generators and corrects one displayed
phase-19 equation in a comment. It is a lead for structural continuation;
the programs above are certified by their stated modular lifts and finite
congruences.

## 5. Complete cost and replay

Preprocessing is polynomial in q: the division polynomial, embedding
inversion, bounded-degree monomial family, prefixes through O(q^3), exact
elimination and validation all have polynomial sizes and operation counts.
This says nothing about the size of the original index N.

To evaluate a certified program at odd N, Frobenius turns exponent bits
into power-of-two theta scales. A degree-four monomial has at most four
positive quadratic variables. Enumerate all but its last variable and use
a square test and residue test for that variable. The resulting query
bound is O(poly(q)*N^(3/2)*poly(log(N+1))) time, with polynomial-in-q program
storage and O(poly(q)+log(N+1)) working bits once the program is supplied.
The checker returns the actual enumeration counts. It supplies no
unrestricted factorization or quadratic-representation oracle.

Run `make verify-theta-programs`. The checker regenerates every listed
program, checks the deterministic
[compressed program artifact](../evidence/theta-coordinate-programs.json.gz),
replays the native embedding, checks small enumerative original queries,
and performs 18 original phase-three controls through 1,032 bits using
only its unary program. The
[receipt](../evidence/theta-coordinate-programs.json) records the finite
congruences, covectors, source hashes, artifact hashes and underfit failure.
Use `python3 tools/check_theta_coordinate_program.py --write-programs` only
to regenerate the artifact deliberately after a source change.

A general proof of successful spanning, a usable alternative for the
missing coordinate, and a sufficiently fast high-index evaluator remain
open. More numerical fits do not resolve any of those obligations.
