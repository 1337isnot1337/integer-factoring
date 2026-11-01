# A cyclic Pauli crossing and unrestricted physical types

**Structural question, 2026-11-01.** Can noncommuting qutrit words with
rank-one leading coefficients support a stronger complete crossing bound
after coupling modes and assigning unrestricted whole-symbol potentials?
The desired support is `R<1`, with all physical costs and a uniform transfer;
the complete factoring score remains **1/2**. No new Lean work is attempted.

The [new research priority](../AGENTS.md#research-priority-cross-the-barrier-not-polish-a-blocked-family)
requires a mechanism and a stopping test. This construction changes the
representation from the reflection/Fock source to a finite Weyl algebra.
Its charge blocks are literal polynomial restrictions; the physical
potentials need not preserve mode occupations or split over digits.
The 96-letter protected-core hypotheses do not hold. The binary
coordinate-support theorem is not asserted for this ternary charge model.
The fixed-gauge axis theorem still applies; two trials fix the balanced
angle and one uses the registered bounded interior interval.

Bare tensor products are proved blocked below, and a coupled control has
an independent bare Jensen obstruction. **Those bare searches are skipped.**
The diagnostic changes that assumption by freeing both physical types
and all 208 nonreference heights. A useful result is an independently
verifiable negative support margin or a structural profile mechanism
that can cross it. If neither appears, do not enlarge the loop count,
field or cutoff merely to continue optimizing this family.

## Polynomial factors and source attribution

[Rowell, arXiv:1410.8863, Section 2](https://arxiv.org/abs/1410.8863)
gives rational coefficients for powers of a finite-order quantum-torus
generator. We use its unnormalized coefficient formula as an algebraic
source. At order three in F4 it specializes to

```
X0=1,       X1=X2=(a+1)/(omega+a*omega^2).
```

Our calculation clears that denominator and changes the parameter by
`a=omega+omega^2*z`. For `U^3=I` and `P=I+U+U^2`, this gives
`F_U,c=z I+c(1+z)P`, with `c=omega^2`; conjugating field coefficients
gives `c=omega`. The arbitrary words below require only these polynomial
identities. No Yang-Baxter relation for an arbitrary chosen word is assumed.

Use symmetric Weyl operators on `2k` ternary coordinates:

```
U(a,b)|s> = omega^(2*a.b+b.s)|s+a>,
```

with arithmetic modulo three in the exponents and labels. Each nonzero
Pauli vector has three equal eigenspaces, so `P` is an idempotent of rank
`n/3`, where `n=3^(2k)`. The polynomial factor has

```
F_U,c(1)=I,   F_U,c(0)=c P,
leading(F_U,c)=I+cP,
ord_0 det(F_U,c)=2n/3.
```

Its leading operator is a qutrit Clifford phase. A second allowed factor,

```
H_U,d(z)=z^2 I+d(1+z)P,       d in F4 nonzero,
```

has endpoint and leading operator `I`, zero coefficient `dP`, and zero
determinant order `4n/3`. This degree-two factor can couple modes while
leaving the leading channel unchanged. Its degree and determinant cost
are retained; no denominator or scalar is discarded.

## A three-gate cell

For one qutrit on each physical party take

```
U1=D((1,1),(0,0)),
U2=D((0,0),(1,2)),
U3=D((1,1),(1,1)),
T=F_U1,omega^2 * F_U2,omega * F_U3,omega^2.
```

`U1,U2` commute and their two projectors give the Bell line. Thus the
zero coefficient of `T` has rank one. Its row-vector matrices can be
chosen as `omega^2 I3` and `diag(1,omega^2,omega^2)`, both invertible.
The leading coefficient is the full transpose channel `M -> A M^T B`:

```
A = [[omega,0,0], [0,0,1], [0,1,0]],
B = [[1,0,0], [0,0,omega^2], [0,omega^2,0]].
```

This explicit identity is checked coefficientwise. Consequently the
opposite leading matrix also has rank one. Tensoring `k` cells gives
physical alphabet `q=3^k`, degree `3k` and two nondegenerate leading
rank-one coefficients.

In the implementation, all `2k` Bell-prefix factors precede the `k`
completion factors. Insert degree-two `H` factors after the prefix.
The zero row rank remains at most one and the leading opposite matrix
remains rank one. Vanishing or degenerate zero vector matrices are
reported; seeds are not filtered out for that outcome.

## Exact blocks with arbitrary physical potentials

Choose loop shifts `a_X=a_Y`. Every factor then commutes with all
`Z_Xi Z_Yi^(-1)`. The exact tensor support satisfies

```
x-y=v-u,     componentwise modulo three.
```

The zero R cut splits by `x-y`; the infinity C cut splits by `x+u`.
Each has `q` square blocks of order `q`. The endpoint is the identity
in both arrangements, so every block has full generic rank. Arbitrary
monomial heights on the four physical legs preserve these polynomial
zeros. Both full weighted profiles are the union of all block profiles.

The [backend](../tools/cyclic_pauli_pressure.py) merges nested block
minors and their intrinsic prefix orders, then reuses the existing
asymmetric/Gibbs derivative. It charges the complete two physical energy
lists. No Gaussian substitution, selected sector, guessed zero tail or
height-homogeneity assumption is used.

## Bare obstructions and the changed assumption

The three-gate cell has normalized infinity profile
`r=(0,1,1,1,1,1,1,2,2)` and zero profile
`s=(0,1,1,2,2,3,3,3,3)`, with degree three. Its raw infinity orders
are `r-3`. Put `L=log 9`, and let `C(t,a,b)` be its bare support numerator.
Jensen and the minimum terms give, for every positive scale,

```
C >= (a+b)L-t/9,
C >= aL+17t/9,
C >= bL+t.
```

For norm-one positive `a,b`, `a+b-1 >= min(1-a,1-b)`.
If `a>=b`, use the second bound for `t>=9(1-a)L/17` and the first
below that point. If `b>=a`, use the third above `(1-b)L` and the first
below. These bounds give **`C>=L`**, hence `R>=1`, at every angle and
scale. Tensor powers multiply both numerator and denominator by `k`.
The entire uncoupled bare family is therefore excluded, not just a finite
screen or a selected scale.

More generally, if complete raw orders have
`sum(infinity)+sum(zero)<=0`, the two Jensen bounds give
`C>=(a+b)log(n)`, also excluding bare improvement at every angle and scale.
The independently checked nine-letter and 27-letter coupled fixtures have
raw sums **-36** and **-335**. The pilot records and obeys that condition.
These statements do not cover unrestricted physical potentials.

A single separable spectral generator is also blocked: if constant local
bases diagonalize a polynomial path, each matched scalar has
`val_0(f)-deg(f)<=0`. Holder's inequality gives bare support `R>=a+b`.
Noncommuting local Weyl generators and nonadditive physical potentials
are precisely the assumptions changed here.

## Independent controls and registered diagnostic

The [checker](../tools/check_cyclic_pauli_crossings.py) compares 1,782
complete polynomial entries using independent packed Weyl actions, all
6,561 entries of the two-cell tensor identity, and **1,638 complete
untruncated weighted orders**. That includes all 1,458 orders of the
27-letter coupled fixture, computed block by block with fraction-free
Python arithmetic at new integer heights. Exact off-block zeros and
endpoints justify assembling those complete profiles. It also compares
the full 81-channel native engine with the block backend and checks
84 real derivatives; maximum error is below `4e-10`.

The [receipt](../evidence/cyclic-pauli-crossing-controls.json) records all
source words, gauges, independent orders and source hashes. These are
finite implementation controls; they are not a negative support certificate.

The [manifest](../evidence/cyclic-pauli-typed-diagnostic-manifest.json)
freezes three 27-letter sources: the tensor-product control, the checked
one-loop fixture, and one OS-seeded four-loop source. The latter uses
shared shifts and independently chosen ternary phases, with every outcome
retained. For each source attempt fixed-balanced near-zero and graded
starts, then a free-angle normalized-gap start. All 208 heights are free;
the graded point is only an initialization. Each trial has at most 200
SLSQP iterations and 30 seconds; the full diagnostic has 420 seconds.
No bare optimization is run when its recorded obstruction applies.

Bounds are heights `[-128,128]`, log scale `[-9,log(4)]`, log tangent
`[-2,2]`; balanced trials fix log tangent zero. Starts have recorded jitter
and scale at least `.05`. Every rejected query retains coordinates and
reason. Partial statuses, complete best/last profiles and source words
remain explicit. An improving floating point requires fresh independent
exact profiles and outward pressure before any complete-bound claim.

```
make verify-cyclic-pauli-crossings
python3 tools/screen_cyclic_pauli_crossings.py --manifest evidence/cyclic-pauli-typed-diagnostic-manifest.json --output /path/typed-diagnostic.json
```

Directed by Cinna Davis / 1337isnot1337 with substantial OpenAI assistant
support. The credited source formula, existing Weyl arithmetic, complete
Smith engines and Gibbs/support machinery are reused. This is a literal
source and structural diagnostic, without a claim of external novelty or
independent mathematical acceptance.
