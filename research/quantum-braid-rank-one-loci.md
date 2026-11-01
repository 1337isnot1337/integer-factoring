# All-extension rank-one loci of two small braid spaces

This note determines both rank-one loci of the complete six-site
coefficient spaces of the scalar strict F4 braid and the
[noncommuting matrix auxiliary](quantum-braid-matrix-auxiliary.md).
The result holds over **every extension of F4**, including an algebraic
closure. It uses polynomial minor equations and Frobenius, rather than
enumerating F4-valued coefficient combinations.

| Six-site space | Dimension | Row-cut rank-one locus | Column-cut rank-one locus |
| --- | ---: | --- | --- |
| Strict scalar braid | 29 | Top-coefficient line | Constant-coefficient line |
| Matrix auxiliary | 55 | Only zero | Constant-coefficient line |

Here a locus includes zero, and a line means its extension-linear span.
The matrix auxiliary still has its uniform leading coefficient of
rank `(L,L)`; the absence of rank one is compatible with that positive
construction. These finite-space conclusions do not exclude higher
low-rank pieces, other alphabets, longer paths or a better factoring
algorithm. The complete score remains `1/2`.

## 1. Pure-square consequences of minor equations

Let `M_1,...,M_k` be an F4 basis of one cut of a coefficient space.
Write

```math
M(c)=\sum_i c_i M_i.
```

For every pair of rows r,s and columns u,v, its two-by-two minor is

```math
\sum_i a_i c_i^2+\sum_{i<j}b_{ij}c_ic_j.        \tag{1}
```

Assemble those coefficient vectors over F4, putting the k square
monomials first. Exact highest-pivot row elimination identifies the
intersection of their row space with the square-coordinate subspace:
the echelon rows whose pivot is below k form a basis of the intersection.
No relation such as `c_i^2=c_i` or `c_i^4=c_i` is imposed.

Every pure-square consequence `sum a_i c_i^2=0` implies

```math
\sum_i\sqrt{a_i}\,c_i=0,                       \tag{2}
```

over every extension field. The square roots of the coefficients are
already in F4: `sqrt(a+b alpha)=(a+b)+b alpha` for binary a,b.
The square of (2) is the original consequence, and Frobenius is injective
in a field. Thus (2) is necessary without any finite-field-point
assumption or surjectivity requirement on the extension's Frobenius.

Restrict M to the kernel of these linear conditions and reassemble
its minors if needed. This gives successive necessary linear envelopes
for the rank-one locus. Early termination is valid: any subset of the
minor equations gives necessary consequences. A zero terminal envelope
proves there is no nonzero rank-one matrix. A one-dimensional terminal
envelope, whose generator is explicitly rank one, is exactly the locus.
All witness and equation coefficients lie in F4, so the same conclusion
holds after every field extension.

This is the F4 analogue of the pure-square row-space method already
used for the binary [energy Riccati reduction](joint-energy-polynomial-dressings.md)
and implemented in `check_joint_energy_riccati.py`. Here the polynomial
equations are two-by-two minors of the actual crossing cuts; the
earlier energy-gauge scope and its conclusions are not imported.

## 2. Exact certificates for the four loci

The scalar strict crossing has coefficient basis dimension29; the
auxiliary has dimension55. Each cut is a 25-by-25 matrix. The
[replay](../tools/check_quantum_braid_rank_one.py) constructs complete
coefficient bases and derives the envelopes above. It checks the
nonzero terminal witnesses against the constant or top coefficient,
including a separate elementary F4 matrix-rank computation.

Column cuts first reduce to an envelope of dimension two; reassembling
the minor equations in that envelope gives the constant line. The
scalar row cut reduces directly to its top line. The auxiliary row
cut reduces directly to zero. The
[source-pinned receipt](../evidence/quantum-braid-rank-one-loci.json)
contains each chosen coefficient basis, the exact square-root
constraints at every stage, the number of minor equations used and
the final coordinates in the original basis.

The packed eliminator is checked against elementary F4 elimination
and the independent dimension identity

```math
\dim(W\cap U)=\dim W+\dim U-\dim(W+U),
```

where U is the square-coordinate space. Literal F16 determinants
at extension-valued coefficient vectors independently check the minor
assembly, including mixed terms. At four sites the auxiliary spans
the entire 16-dimensional tensor space; both cut rank-one loci span
that whole space. This supplies a negative control against a false
general obstruction inferred merely from span enlargement.

```sh
make verify-braid-rank-one
```

NumPy is needed for these optional checks. The receipt supports exact
finite-dimensional algebra plus the all-extension argument above;
it does not provide an all-length rank-pressure theorem. In particular,
the next positive task is to organize higher-rank asymmetric channels
in the changed auxiliary family, whose leading identity channel is
already explicit. No new Lean theorem is claimed.
