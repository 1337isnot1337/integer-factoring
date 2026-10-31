# Exact endpoint lifts after a fixed field extension

## Result and scope

The coefficient constraints used in polynomial filter synthesis can be
solved over `F4`, then lifted to `F16` to obtain nonzero values in every
required endpoint coordinate. For the eight diagonal endpoint coordinates
of a three-mode kernel, the lift exists exactly when no coordinate is
forced to vanish by the homogeneous equations. The output is an explicit
finite polynomial coefficient assignment.

This enlarges the available coefficient field and removes one finite-field
feasibility obstruction. It does not, by itself, improve either cut profile
or the factoring exponent. Every degree cap, weighted profile, and endpoint
repair still has to be charged in the chosen construction.

## 1. Endpoint images and the extension criterion

Let `A` be a matrix over a finite field `K`, let `V=ker A`, and let
`e_1,...,e_s` be linear endpoint functionals on `V`. Extend scalars to a
finite field `L` containing `K`, with `q=|L|`. If `s<=q`, there is a vector
`v` in `V tensor_K L` satisfying

\[
e_j(v)\ne0\qquad(1\le j\le s)
\]

if and only if every functional is nonzero on `V`.

Necessity is immediate. For sufficiency, the bad vectors belong to at
most `s` proper linear hyperplanes in the extended solution space. If
its dimension is `d>=1`, each hyperplane contains `q^(d-1)` vectors.
For `s<q`, their union has size at most `s q^(d-1)<q^d`. At `s=q`,
all hyperplanes contain zero, so the union has size at most
`q q^(d-1)-(q-1)<q^d`. A proper kernel of smaller dimension only
decreases this bound. Thus some vector avoids every bad hyperplane.
Extension of scalars preserves the ranks and nonzero functional tests.

When `s>q`, separate coordinate tests are insufficient. The seven nonzero
binary linear forms on three variables, with one coordinate repeated,
give eight endpoint coordinates that are each nonzero somewhere, while
their `F4` extension has no fully nonzero vector. The exact two-vector
cover argument appears in
[the reduced Smith-lattice note](weyl-lattice-filter-synthesis.md).
The same endpoint map admits a fully nonzero vector over `F16`.

## 2. Construct the coefficient lift

Polynomial degree and valuation requirements give homogeneous linear
equations in a finite coefficient vector. Packed elimination over `F4`
stores each equation in two binary planes, for the coefficients of
`1` and `omega`, with `omega^2=omega+1`. Reducing the endpoint
functionals modulo the equation row space computes their restrictions
to the free variables exactly.

The implementation computes a basis of the small endpoint image and
retains a free-variable lift of each basis vector. When `s<q`, it grows
the set of nonzero endpoint coordinates. At each step, at most `s`
scalars would cancel an already nonzero coordinate or leave a newly
available coordinate zero. A remaining scalar adds the next basis
vector while preserving every coordinate already obtained.

The strict inequality in this greedy argument matters. At `s=q`, a
particular greedy sequence can get stuck even though a solution exists.
The implementation uses a complete finite endpoint-image search at
equality or above. It normalizes the first nonzero image coefficient to
one, removing scalar duplicates without changing feasibility. A separate
four-coordinate `F4` regression exercises this boundary case.

After finding a fully nonzero image vector, the retained basis lifts and
back substitution give its complete coefficient assignment. Every
homogeneous equation is checked again against that assignment. An affine
system can be handled by adding a homogenizing coordinate, requiring it
to be nonzero, and dividing the final assignment by that coordinate.

## 3. The concrete extension and factoring transfer

Use

\[
\mathbf F_{16}=\mathbf F_4[u]/(u^2+u+\omega).
\]

The polynomial is irreducible: in `F4`, the image of `x^2+x` is
`{0,1}`, which does not contain `omega`. An element is represented as
`a+u b`, with `a,b` in `F4`. Multiplication is

\[
(a+ub)(c+ud)=(ac+\omega bd)+u(ad+bc+bd).
\]

For eight endpoint coordinates, `8<16`, so the constructive greedy
criterion is sufficient. The coefficient equations may already contain
nonbinary `F4` row flags; they need not be defined over `F2`.

For diagonal filters, nonzero endpoint coordinates give the constant
physical diagonal repair. For occupation-transfer kernels, the required
zero endpoint entries are imposed as linear equations and the nonzero
diagonal coordinates complete the energy-filtered triangular repair.
The respective construction notes establish why those repairs preserve
both weighted profiles.

The field degree is fixed independently of input length and tensor-power
length. The [fixed-field descent theorem](fixed-field-descent.md) accounts
for descent after the full typed tensor power. It does not charge a new
extension factor independently at every tensor position. A larger
coefficient field therefore remains compatible with that transfer, once
a strict complete pressure certificate is obtained.

## 4. Implementation and independent checks

The [linear helper](../tools/f4_lattice_linear.py) exposes packed field
operations, elimination, endpoint restrictions, explicit `F4` or `F16`
lifts, and assignment verification. Its
[checker](../tools/check_f4_lattice_linear.py) compares 360 independent
systems with 98,208 complete `F4` coefficient assignments. It also checks
all 4,096 triples for associativity and distributivity in `F16`, every
nonzero inverse, and the embedding of `F4`.

The [saved evidence](../evidence/f4-lattice-linear.json) includes the
eight-coordinate obstruction over `F4`, its explicit successful `F16`
lift, and the four-coordinate greedy-boundary regression. Reproduce it
with:

```sh
python3 tools/check_f4_lattice_linear.py --output /tmp/f4-lattice-linear.json
```
