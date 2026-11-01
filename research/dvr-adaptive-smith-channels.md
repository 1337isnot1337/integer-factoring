# Complete regular channel-basis choices at fixed lattice precision

## Result and scope

At a fixed full-tail precision, every regular unimodular change of Smith
channels can be represented by finitely many polynomial jets. A channel's
valuation, target lattice, and residue vector determine its cost. The
best basis is a minimum-weight vector-matroid basis on the residue
vectors. This is a complete oracle for the declared class of local basis
changes, rather than a search through selected polynomial shears.

The oracle is independent of the chosen Hermite basis of a state and
the chosen minimum-valuation Schur factorization. It extends the
[constant-basis construction](matroid-adaptive-smith-channels.md).

For the same binary companion crossing, ten states at precision two
and the rational point `s=5` give

```math
\beta_\infty\le\frac{266}{25},\qquad
\beta_0\le\frac{241}{200},\qquad
B\le\frac{32053}{2500}=12.8212.                         \tag{1}
```

There is a strict comparison within this state model. Every constant
channel-basis policy, including changes across unequal Smith orders,
has zero-side word rate at least `603/500=1.206`. The new zero-side
upper bound is `241/200=1.205`. Thus parameter-dependent basis changes
add genuine certificate strength.

The previous note separately proves that the original fixed-factor
word-pressure product exceeds `12.9` at every scale `s>=1`, so (1)
also improves that global comparison. These numbers remain above the
elementary scalar pressure `256/27`. No stronger factoring exponent is
claimed.

## 1. Every channel basis is a basis of one column lattice

At a state `M`, flatten the regularized physical action `B_M K` into a
matrix `A`. Its rows contain the input lattice coordinate, the selected
physical half, and the outgoing auxiliary index. Its columns contain
the opposite physical half.

Minimum-valuation Schur elimination gives

```math
A=UR,\qquad U_i=\tau^{\nu_i}\ell_i,\quad 1\le i\le r.    \tag{2}
```

All entries of `ell_i` and `R` are regular. The columns `ell_i` are
jointly independent modulo `tau`, even when their orders differ:
restriction to the ordered pivot rows gives a triangular matrix with
unit diagonal. The rows of `R` have a triangular unit minor on the
ordered pivot columns.

Consequently `R` has a regular right inverse. The columns `U_i` are an
`O=F_2[[tau]]` basis of the column lattice of `A`, since

```math
\mathrm{col}_O(A)=\mathrm{col}_O(U).
```

Every other basis of this column lattice has the form `U S^T` for
`S in GL_r(O)`, and

```math
A=(US^T)((S^{-1})^T R).                                  \tag{3}
```

The new right factors stay regular. After stripping each new left
column's common valuation and taking its auxiliary row lattice, (3)
has exactly the form required by the
[state-adaptive coefficient-rank theorem](state-adaptive-smith-channels.md).

This class allows all bases of the existing column lattice. It does
not assert optimality among decompositions through a larger column
lattice, different auxiliary representations, or arbitrary coefficient
recombinations.

## 2. A finite jet captures every primitive column

A vector `v in O^r` is primitive when its residue `b=v mod tau` is
nonzero. Such a vector can occur as a row of a unimodular basis change.
For its left column, define

```math
a(v)=\mathrm{ord}_\tau(Uv)
=\min_i\bigl(\nu_i+\mathrm{ord}_\tau(v_i)\bigr).          \tag{4}
```

The equality follows from residue independence of the `ell_i`: the
terms of least order cannot cancel. Primitiveness gives

```math
\nu_{\min}\le a(v)\le\nu_{\max}.                         \tag{5}
```

Reshape `tau^(-a(v)) Uv` into auxiliary rows. Let `M(v)` be their
`O`-span enlarged by `tau^P O^D`. The normalized column modulo `tau^P`
determines this full-tail lattice exactly.

Put

```math
J=P+\nu_{\max}-\nu_{\min}.                              \tag{6}
```

Replacing every coordinate of `v` by its polynomial residue modulo
`tau^J` preserves its residue vector, the value (4), and the normalized
column modulo `tau^P`. Indeed the change in `Uv` has order at least
`nu_min+J=nu_max+P>=a(v)+P`. Since `P>=1`, it cannot change (4).

Therefore the triples

```math
\bigl(b(v),a(v),M(v)\bigr)                               \tag{7}
```

are a finite, exactly computable set. Every triple has a polynomial
representative of degree less than `J`. This proves completeness over
all regular power-series vectors, including those that are not rational.
The finite representatives used by the certificate are polynomial, so
their inverse basis matrices and all factor identities are rational.

### Efficient normalized enumeration

The engine avoids enumerating irrelevant high coefficients. For each
possible `a` in (5), write `w_i=tau^(nu_i-a)v_i` and enumerate the
`P`-jets of the regular `w_i`. At least one `w_i` has nonzero constant
coefficient. If `nu_i>a`, impose divisibility by `tau^(nu_i-a)`.
If `nu_i<a`, the residue of `v_i` must be zero. If `nu_i-a>=P`, the
residue of `v_i` is invisible in the normalized jet and is enumerated
separately. These are precisely the conditions obtained from (4).

Each resulting triple is stored once with a concrete polynomial
representative. An independent checker also enumerates the original
`v`-jets directly at precision one and obtains the identical action set.

## 3. Greedy optimization on residues

Fix a positive state potential `h` and a rational `s>1`. The cost of a
candidate (7) is

```math
c(v)=s^{\epsilon L-a(v)}h_{M(v)},                         \tag{8}
```

where `epsilon=0` at zero and `epsilon=1` at infinity.

Vectors `v_1,...,v_r` form a unimodular matrix exactly when their
residues form a basis of `F_2^r`. For each nonzero residue, first keep
the least-cost lift in the finite catalogue. Then apply the ordinary
greedy minimum-basis algorithm to these `2^r-1` residue vectors.
Equivalently, all lifts are parallel elements in the same vector
matroid. There is no compatibility constraint between the individually
chosen lifts beyond residue independence.

This computes the exact smallest potential image obtainable from any
`GL_r(O)` change at that state. Independently selecting the minimum
basis at every state gives a finite stationary matrix. The positive
potential certificate and finite-policy Collatz--Wielandt formula from
the preceding note apply without change.

### Independence of implementation bases

Changing the Schur factorization in (2) changes one `O`-basis of the
column lattice to another. The complete set of bases in (3) is therefore
unchanged.

Changing the state basis to `B'_M=W B_M`, for `W in GL(O)`, multiplies
the flattening on the left by the regular unimodular operator that acts
by `W` on input lattice coordinates and identically on the other row
indices. It preserves every column valuation. After reshaping, it takes
invertible combinations among the input-coordinate rows at each fixed
physical half, so it preserves their auxiliary row span. Thus it
preserves the target lattices in (7) as well.

The complete local action set depends on the actual lattice state and
physical cut, not on which Hermite basis or Schur pivots the implementation
uses. This invariance concerns the complete oracle; an arbitrary fixed
subset of basis changes need not have it.

## 4. Exact companion certificate and strict comparison

Use the ten precision-two states, initial state, terminal orders, and
companion crossing from the preceding note. The new infinity potential
is unchanged:

```math
h_\infty=(10000,5320,10000,10000,5320,
               10000,5320,10000,10000,10000).
```

At zero use

```math
h_0=(100000,94634,63719,62303,94634,
           62828,94634,63719,54764,61230).
```

The complete greedy oracle gives

```math
T_\infty h_\infty\le\frac{266}{25}h_\infty,
\qquad T_0h_0\le\frac{241}{200}h_0.                      \tag{9}
```

The largest zero-side row ratio is exactly `1876877/1557575`, below
`241/200`. The boundary prefactors are `1` and `47317/31414`. The
state-adaptive theorem now proves (1) for every word length.

For example, one zero-side state's selected basis change has rows

```math
S=\begin{pmatrix}
\tau&1+\tau&0&0\\
0&\tau&1&0\\
0&\tau&0&1\\
1&0&0&0
\end{pmatrix}.
```

Its determinant is `1+tau`, a unit. The inverse has rational regular
entries. The poles it may have away from zero are harmless: the local
identity is exact over `F_2(tau)`, and the coefficient-rank proof uses
only the regular expansion at the selected endpoint.

To compare all constant bases, retain the earlier positive vector

```math
g=(10000,9365,6316,6168,9365,6211,9365,6316,5423,6062).
```

Let `F_const(g)` be the rowwise minimum over **all** constant invertible
basis changes, including combinations across different orders. Exact
enumeration of their fifteen nonzero vectors at each state gives

```math
F_{\mathrm{const}}(g)\ge\frac{603}{500}g.                 \tag{10}
```

The smallest row ratio is `163519/135575`, which is strictly larger
than `603/500`. Every allowed constant-basis policy `T` satisfies
`Tg>=F_const(g)`, so `rho(T)>=603/500`. Since all terminal weights
are positive in these full-tail states, the same lower rate applies
to the counted words from the declared initial state, not only to a
possibly inaccessible spectral component.

By (9), the regular-basis policy has zero rate at most `1.205`, whereas
(10) forces constant-basis rates at least `1.206`. This is a lower bound
on a class of sufficient word certificates, not on actual coefficient
ranks.

## 5. Replay and limits

```sh
python3 tools/check_dvr_smith_channels.py \
  --output evidence/dvr-smith-channels.json
```

The standard-library checker verifies the full regular factor identities,
all selected lattice memberships, both upper potentials, and the lower
potential for every constant basis choice. It compares 4,335 directly
enumerated primitive jets with the compressed oracle at precision one.
Literal polynomial products separately check 4,368 endpoint entries
and eighteen rank inequalities through length three.

The implementation reports state or jet-budget exhaustion as inconclusive.
Completeness concerns all local unimodular channel-basis choices at the
declared finite precision. It does not claim convergence to the true
coefficient ranks, a finite certificate attaining an unrestricted optimum,
or monotonicity under unrelated changes of precision and representation.

Files:

- [Complete jet engine](../tools/dvr_smith_channels.py).
- [Independent exact checker](../tools/check_dvr_smith_channels.py).
- [Exact evidence](../evidence/dvr-smith-channels.json).

The constructive objective is to apply this stronger oracle to a local
crossing family whose pressure is already competitive. The companion
example establishes strict improvements between certificate classes;
its absolute pressure remains above the threshold needed for factoring.
