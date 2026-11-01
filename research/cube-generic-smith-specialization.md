# Exact generic Smith profiles from a homogeneous finite-field specialization

Specializing the other two cube variables to one can increase Smith
orders. Such a computation is useful for rejection screens but cannot
certify a positive bound. Homogeneity permits a finite specialization
that preserves every generic order exactly, without a large Kronecker
exponent.

The [implementation](../tools/cube_generic_smith128.py) covers binary
polynomial matrices of order at most 64 whose complete minor degree
bound is below 128. Its [native engine](../tools/cube_generic_smith128.cpp)
uses exact field and truncated-series arithmetic. These are explicit
implementation limits; the following lemma holds at arbitrary size.

## 1. Homogeneous specialization lemma

Let `A(x,y,z)` be a polynomial matrix over `F_2`. Suppose row weights
`r_i >= 0` and column weights `c_j >= 0` satisfy

\[
\deg A_{ij}=r_i+c_j
\]

for every nonzero entry. Half-integral weights are allowed when all
entry degrees are integers. For a square minor on row set `I` and
column set `J`, every determinant term has degree

\[
d_{I,J}=\sum_{i\in I}r_i+\sum_{j\in J}c_j.
\]

Let `Gamma` bound the total degree of every minor. Choose an element
`theta` of degree `e > Gamma` over `F_2`. Then

\[
A(x,1,\theta)
\]

has exactly the same determinantal valuations at `x` as the generic
matrix over `F_2(y,z)[[x]]`.

To prove this, write any nonzero homogeneous minor as

\[
\Delta(x,y,z)=\sum_a x^a f_a(y,z).
\]

The coefficient `f_a` is homogeneous of degree `d_(I,J)-a`. Substitution
`y=1` is injective on polynomials of this fixed homogeneous degree:
distinct terms have distinct powers of `z`. Its resulting degree is
at most `Gamma < e`, so substituting `z=theta` is also injective.
Consequently the first nonzero coefficient in `x` survives. This holds
for every minor, and therefore for every Smith order. The proof applies
to either other axis after permuting variables.

For a square matrix with nonnegative weights, one may take

\[
\Gamma=\left\lfloor\sum_i r_i+\sum_j c_j\right\rfloor.
\tag{1}
\]

Every subminor uses subsets of these weights. For an invertible matrix,
the complete determinant has integral degree and the floor is exact.

The input coefficient field matters. The retained engine assumes
coefficients in `F_2`; merely putting coefficients from a larger field
into the same degree-128 extension would not prove the same statement.
For a larger base field, the extension degree must exceed `Gamma` over
that base field.

## 2. Complete absolute precision

Assume the matrix is generically invertible. Its nonnegative Smith
orders sum to the valuation of its full determinant, which is at most
`Gamma`. Every individual order is therefore at most `Gamma`.
Computation modulo `x^(Gamma+1)` determines the complete profile.

The native engine repeatedly selects an entry of minimum valuation,
normalizes its unit factor, and eliminates its row and column over the
truncated power-series ring. These are invertible operations over the
valuation ring. It records all pivots and rejects an incomplete profile.
The wrapper also checks their sorting, nonnegativity, and total budget.

The field is

\[
\mathbb F_2[z]/(z^{128}+z^7+z^2+z+1).
\]

The checker verifies irreducibility by the exact Rabin criterion and
compares the native field multiplication with an independent polynomial
implementation. A portable multiplication path is included; processors
with carryless multiplication use an equivalent faster path.

## 3. Cost and intended use

With an explicit `n` by `n` matrix, degree bound `Gamma`, and field
extension degree `e`, direct series elimination uses polynomially many
bit operations in `n`, `Gamma`, and `e`. Naive series products give
`O(n^3 Gamma^2)` field operations, apart from polynomial preprocessing
and field inversion. This cost does not hide an oracle for generic
rational-function rank.

For a Gaussian cube tensor, the one-particle matrix has order `4M`.
Its full Fock cut orders then follow by taking subset sums. For a proper
physical restriction, one must construct and analyze the compressed
matrix itself, or prove an additional structural reduction. The lemma
does not remove that dimension cost.

The [checker](../tools/check_cube_generic_smith128.py) includes known
Smith profiles after polynomial unit transformations, all three cuts
of the scalar cube family, an independent complete-minor comparison
for a coupled weighted kernel, and a control where setting both other
variables to one gives an incorrect higher valuation. The
[receipt](../evidence/cube-generic-smith128.json) records the exact checks.

Run `python3 tools/check_cube_generic_smith128.py` to reproduce them.
