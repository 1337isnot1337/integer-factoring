# Finite fused six-vertex crossings over F16

## Result and scope

Literal finite-field fusion supplies additional strict polynomial crossings
on three and four symbols. Two complete 72-start physical-gauge searches
were run, and every saved full profile and cost was replayed. Their best
ratios are approximately 0.98994 and 1.00038, well above the critical value
1/sqrt(2). These tests do not improve the factoring bound.

This note records the construction and the unsuccessful finite search so
that the same family need not be rediscovered. No all-spin fusion identity,
generic endpoint claim, or global gauge optimum is assumed. Each requested
finite instance must pass the literal invariant-subspace and endpoint checks
before it enters a search.

## 1. A completely evaluated finite construction

Work in F16=F4[u]/(u^2+u+omega), with omega^2=omega+1. For a nonzero
theta other than one, put c=theta+theta^(-1). On two bits, the checked
polynomial matrix has diagonal entries theta*z+theta^(-1) on 00 and 11,
and the following block on 01 and 10:

\[
 \begin{pmatrix}cz&z+1\\z+1&c\end{pmatrix}.
 \tag{1}
\]

For m bits, the weight-k vector is the sum of all weight-k words, with
coefficient theta raised to minus the inversion count. Its canonical word,
with zeros before ones, has coefficient one.

The constructor starts with each tensor product of two such weight vectors
in the full space on 2m bits. For i=m-1 down to zero and j=0 through m-1,
it applies (1) at adjacent positions i+j and i+j+1, with
z=w*theta^(2(i-j)). No intermediate projection is made. The output is
extracted at the canonical words and then reconstructed on the full bit
space. Exact equality proves that the selected subspace is invariant for
this specific instance.

The exact common polynomial gcd is removed from all entries. The resulting
matrix must specialize to a nonzero scalar identity at w=1. Dividing by
that scalar gives F(1)=I. Relabeling

\[
 T_{xyuv}=F_{yx,uv}
 \tag{2}
\]

gives the strict crossing T(1)=[x=v][y=u]. The complete polynomial tensor,
with its actual degree after common-factor cancellation, is the search
input.

Two checked examples are:

| m | Physical alphabet | theta encoding | Cancelled common degree | Final degree |
|---:|---:|---:|---:|---:|
| 2 | 3 | 2 | 2 | 2 |
| 3 | 4 | 4 | 6 | 3 |

The same attempt at theta=2 with m=3 or m=4 has an invariant subspace but
fails the required nonzero identity endpoint after common-factor removal.
Those instances are rejected. An invariant subspace alone is insufficient.

## 2. Projective coordinates and physical-symbol costs

The search uses six projective transformations fixing w=1, specified by
their images of zero and infinity:

\[
 (0,\infty),(2,3),(9,15),(14,10),(8,10),(\infty,0).
 \tag{3}
\]

After substitution, multiply by the full denominator to the actual tensor
degree and normalize its value at one. This preserves the strict endpoint
and retains the complete degree cost. Each map uses twelve starting points
for optimization over every independent physical-letter height at both
places and a positive scale.

The native F4 engine receives restriction of scalars from F16. Multiplication
by a+u*b is represented by

\[
 \begin{pmatrix}a&\omega b\\b&a+b\end{pmatrix}.
 \tag{4}
\]

Every F16 Smith order appears twice in that matrix. The two log partition
sums therefore subtract a total of 2*log(2). The physical alphabet and its
entropy remain q. Introducing extra physical letters to compensate for
scalar restriction would be incorrect.

## 3. Completed finite results

| Family | Completed starts | Optimizer-reported successes | Best replayed rho |
|---|---:|---:|---:|
| m=2, theta=2 | 72 | 68 | 0.98993968592235 |
| m=3, theta=4 | 72 | 67 | 1.00038395766880 |

Every start retains its best successfully evaluated point. The starts
without a success flag are included in the records with their explicit
optimizer message; they are not represented as converged optima. All 144
retained complete profile and root calculations were reconstructed from
the literal polynomial tensors.

The independent exact controls check 2,022 endpoint entries, 24 complete
Smith profiles after further restriction to binary matrices, and 34
block-Toeplitz ranks. These checks are separate from the floating gauge
optimizer. The two invalid endpoint examples are explicit negative controls.

```sh
make verify-fused-six-vertex
```

The constructor is
[`fused_six_vertex_crossings.py`](../tools/fused_six_vertex_crossings.py),
and the search driver is
[`search_fused_six_vertex_gauges.py`](../tools/search_fused_six_vertex_gauges.py).
The complete records are
[`fused-six-vertex-gauges-m2-theta2.json`](../evidence/fused-six-vertex-gauges-m2-theta2.json),
[`fused-six-vertex-gauges-m3-theta4.json`](../evidence/fused-six-vertex-gauges-m3-theta4.json),
and [`fused-six-vertex-controls.json`](../evidence/fused-six-vertex-controls.json).

These finite results leave other fields, larger valid fused instances,
additional polynomial operations, and different constructions open. They
provide no obstruction to those wider families.
