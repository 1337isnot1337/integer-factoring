# Division-free trace-dual counts at sparsity-sized precision

## Result and changed representation

A trace-dual basis replaces the nonuniform scalar-trace word map in the
[unary inverse query](one-trace-inverse-words.md). One Teichmuller input
now represents each odd half-width candidate exactly once. The pointwise
circuit uses several coordinate traces and a bit-serial odd-division circuit
made of integer polynomial Boolean gates. It retains the certified no-wrap
test, with no exact-halving gates at all.

After completed trial division through `2^t`, its **unweighted** prefix
count needs only `min(ceil(b/2)-1,ceil(b/t))` output bits. In particular,
trial division below any fixed positive-score budget can reduce output
precision to a power of log b. Starting ring precision is the same K.
This is a constructive count interface, not a faster count algorithm or a
new complete factoring bound. The proved exponent remains **1/2**.

The [registered question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19361839)
changes the word map and its fibre weights. It does not claim to evade
the [pointwise trace-selection restrictions](trace-selection-weil-budget.md)
by a few low-degree equations. The operation requirement changes: the new
circuit is division-free, but has unrestricted degree. The original explicit-product polynomial carry
backend does not automatically contract its shared circuit DAG.

## 1. Acquire a bijective coordinate map

For odd N>=5 write `b=bits(N)`, `a=ceil(b/2)` and `d=a-1>=1`.
Acquire a model of `F_(2^d)` and a basis `beta_0,...,beta_(d-1)` over F2.
Its trace pairing is nondegenerate at **both even and odd degrees**.
Indeed the trace polynomial `X+X^2+...+X^(2^(d-1))` is nonzero and has
degree below `2^d`, so is not the zero function. Multiplication by any
nonzero field element is bijective; thus `c -> Tr(c v)` is a nonzero
functional for every nonzero v. This proves nondegeneracy.

Binary elimination of the d-by-d Gram matrix acquires the trace-dual
basis gamma_i, satisfying `Tr(gamma_i beta_j)=[i=j]`. Polynomially many
field operations suffice. The implementation uses the monomial basis
of a supplied irreducible model; it does not enumerate the field to
acquire its dual. General deterministic model acquisition retains the
finite-field construction input and credit of the
[preceding word construction](witt-trace-native-word-circuits.md#4-the-exact-missing-evaluator-and-complete-transfer).

At precision P use the same unramified ring `W_P(F_(2^d))`. Lift gamma_i
to Teichmuller elements, and for its one Teichmuller input U put

```math
e_i(U)=H_P\bigl(\operatorname{Tr}_W([\gamma_i]U)\bigr),
\qquad x(U)=1+2\sum_{i=0}^{d-1}2^i e_i(U).
```

Here the acquired integer parity polynomial H_P returns exactly 0 or 1
modulo `2^P` on every scalar. Reduction of the lifted trace modulo two
is the ordinary field trace. For `u=sum_i u_i beta_i`, therefore
`e_i([u])=u_i`. Thus U maps **bijectively** to the odd integers
`1<=x<2^a`. There is no Hasse--Weil fibre argument, unknown weight,
odd-degree requirement or assumption that a single lifted trace is uniform.
There are d coordinate traces; this is not the preceding single-trace map.
On the stipulated Teichmuller inputs each lifted trace is simply
`sum_(j=0)^(d-1) ([gamma_i]U)^(2^j)`, implemented by repeated squaring.
Thus its pointwise implementation uses ordinary ring arithmetic as well;
it does not require a separate variable-coordinate or Frobenius oracle.

## 2. Compute the cofactor bits without exact divisions

Suppose trial division has found no factor through `2^t`, t>=1. Set

```math
K=\min\{d,\lceil b/t\rceil\},\qquad k=a+b,\qquad P=K.
```

Integer polynomial gates preserve Boolean values exactly. For bits A,B,C,
use

```math
\operatorname{XOR}(A,B)=A+B-2AB,
```

and the borrow of `A-B-C` is

```math
(1-A)(B+C-BC)+ABC.
```

Both expressions are exactly zero or one on bit inputs, over the integers
and hence modulo `2^K`, including K=1. The difference bit is the XOR of
all three bits. AND is multiplication, and NOT is `1-A`.

Start with the public constant bits of `R=N mod 2^k`. At step i, let y_i
be the current bit R_i, and subtract `2^i y_i x` from the residual word
modulo `2^k`, using those polynomial gates. Candidate bits are `(1,e_0,...,
e_(d-1))`, padded by zero. There are `k(k+1)/2` bit-subtraction positions.
Every register remains a Boolean scalar at precision K; no variable bit
is extracted by a native integer operation or division.

Inductively, before step i the residual is

```math
R=N-x\sum_{j<i}2^j y_j\pmod{2^k},\qquad R=0\pmod{2^i}.
```

Since x is odd, the selected y_i clears bit i. After k steps the residual
is zero, so these bits specify exactly `y=N x^(-1) mod 2^k`, in its
canonical range. This is an ordinary bit-serial odd-division algorithm,
not the assumption of an inverse-bit oracle.

Use the same product mask for `x>1`, `2<=y<2^b`, and the requested prefix
of the d-bit word `w=(x-1)/2`, most significant bit first. As in the
predecessor,

```math
x>1,\quad 2\le y<2^b\quad\Longleftrightarrow\quad
x>1,\quad y>1,\quad xy=N.
```

For the forward direction the congruence is modulo `2^(a+b)`, while both
xy and N lie below that modulus. The converse is immediate. The modulus
is unchanged; the wrapped prime-N=7 false positive at modulus eight
remains excluded. Let the resulting 0-or-1 scalar mask be A(U).

All basis acquisition, lifted traces, parity Horner loops, Boolean gates
and prefix masks have polynomial construction and pointwise bit cost in b.
There are d traces of d Frobenius conjugates, O(b^2) bit-subtraction
positions, ring degree d and **only K bits of ring precision**. An element
uses O(bK) bits. The circuit contains no variable or nonunit divisions.
Its polynomial coefficients are integers, acquired uniformly as in the
preceding parity construction.

This does not bound the cost of summing its domain. Shared polynomial
circuits can have exponentially large expanded degree and support. The
existing product-of-explicit-Frobenius-monomials carry API does not supply
a fast evaluator for this unrestricted-degree circuit DAG.

## 3. The unweighted count fits the smaller precision

Every accepted U corresponds to one distinct proper divisor x. Since the
excluded candidate x=1 has one domain point, every prefix count is
strictly below `2^d`. Also let Omega(N) count prime factors with multiplicity.
All prime divisors exceed `2^t`, and `N<2^b`, so

```math
t\Omega(N)<b,\qquad
\tau(N)\le2^{\Omega(N)}<2^{b/t}\le2^{\lceil b/t\rceil}.
```

Consequently

```math
C_{N,prefix}=\sum_{U\in\mathcal T_d}A_{N,prefix}(U)<2^K.
```

Its residue modulo `2^K` is the exact integer count. This conclusion
uses the actual trial pretest. Merely supplying a large t does not
certify the smaller output precision.

The [implemented conditional consumer](../tools/coordinate_inverse_words.py)
performs the trial pretest first. Even inputs return two, and N=3 is a
separate constant branch. If no factor is found, query the empty prefix.
A zero result proves primality, since every odd composite has a prime
divisor below `2^a`. Otherwise at most d more prefix queries recover
the smallest valid x. There are at most `d+1=a` queries in total.
The pretest costs `2^t poly(b,t)` bit operations. It returns the smallest
factor when it terminates early. A uniform count provider with cost
`F(b) poly(b,t)` gives total cost `(2^t+F(b))poly(b,t)`.

For any fixed target beta>0 choose a fixed integer L>beta and, for large b,
`t=max(1,floor(b/(ceil(log_2(b+1)))^L))`. Then `K=O((log b)^L)` and
the pretest costs `exp(O(b/(log b)^L))`, below the intended
`exp(O(b/(log b)^beta))` scale. Finite small b is handled directly.
This removes linear **output** precision from this conditional interface;
it also reduces starting ring precision to K and removes exact-halving
gates. No count provider at the intended cost is supplied.

## 4. Exact controls, limits and credit

```sh
make verify-coordinate-inverse-words
```

The [checker](../tools/check_coordinate_inverse_words.py) uses the actual
lifted trace, acquired parity polynomial and division-free bit circuit.
It compares against independent literal candidate divisibility and complete
finite field enumeration. The [receipt](../evidence/coordinate-inverse-word-controls.json)
records all 126 points of six field models, including even-degree trace
pairings; every word has exactly one preimage. It records 1764 pointwise
controls and 67 complete minimum-factor recoveries, with full prefix panels
on fixed targets and five rough-input controls. Separate complete small
1960 bit-circuit controls compare every returned cofactor bit against modular
inversion, without using this inverse inside the implemented circuit.

The roughness guard is N=105 with an **unjustified** t=4: four candidates
would wrap to zero at the incorrectly shortened two-bit count precision.
The actual consumer instead finds factor three before any count query.
For N=323 and valid t=4, the new plan has degree four, K=P=3;
the older weighted-trace plan has degree 21, K=22 and P=41. These are
resource invoices for different represented sums, not measured speedups
of a fast contractor. The exhaustive checker oracle is not such a provider.

Trace duality, bit-serial odd division and polynomial lifts of Boolean
gates are standard algebraic ingredients. The ring/Frobenius, Hermite parity
and no-wrap proof retain their preceding source and contributor credit.
The work here is their division-free bijective composition, charged smaller
precision and implemented consumer, developed with OpenAI Codex assistance. No new Lean,
independent mathematical review or improved complete exponent is claimed.


The [full default receipt](../evidence/coordinate-word-full-verification.json)
pins clean source `e764b055dff1d6bf6f7b52e995b2694314618af7` and its
exit-zero 1214.86-second detached replay, fresh private-mount control and
complete compressed log. The [initial interrupted run](../evidence/coordinate-word-full-interrupted.json)
is retained separately and is not a pass. Subsequent archive/integration
prose has focused checks, without another full replay.
Sol worker 1's [inspection](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19362216)
and Sol worker 2's [separate inspection](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19362256)
at that frozen source found no missing condition in the coordinate,
bit-division, no-wrap and precision arguments. They did not rederive all
imported field/ring/Hermite primitives or supply the global sum. Their
scope is narrower than independent review of the entire reduction.
