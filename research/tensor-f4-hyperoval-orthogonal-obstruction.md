# The six-column F4 hyperoval has no all-subset orthogonal scaling

## Result and scope

Let C be a three-by-six matrix over F4 such that every three columns are
independent. No nonzero scaling of its columns over F4 makes every
column-subset Gram matrix have its full expected rank. Some four-column
subset always has column rank three and Gram rank two.

This includes the six-column hyperoval matrix used by the
[F4 crossing construction](f4-hyperoval-crossing-bound.md). The obstruction
concerns the stronger all-subset condition needed for an orthogonal
paired-differential kernel. It does not invalidate that boundary crossing
or rule out its other physical retractions.

The field restriction is essential. The
[extension-field scaling construction](tensor-orthogonal-scaling.md)
still applies after enlarging F4. No improved factoring bound is asserted.

## 1. A missing weight in the kernel

Let H be the kernel of C. It has dimension three. A nonzero vector of H
has at least four nonzero coordinates, since every three columns are
independent. In fact it has either four or six nonzero coordinates.

To prove this, fix a coordinate j and let H_j consist of the kernel
vectors with coordinate j zero. The other five columns have rank three,
so H_j has dimension two. For each k different from j, the subspace H_jk
with both coordinates zero has dimension one, because the remaining four
columns have rank three.

These five one-dimensional subspaces H_jk are pairwise distinct. Their
pairwise intersections would otherwise contain a nonzero kernel vector
on at most three columns. A two-dimensional F4 space has exactly five
one-dimensional subspaces. Therefore the H_jk exhaust all nonzero vectors
of H_j. Every nonzero kernel vector with one zero coordinate has another
zero coordinate. Since three zeros are impossible, its weight is four.

This proves the absence of weight-five words. Nonzero column scaling
preserves every zero pattern, so the same conclusion holds for the kernel
of every scaled matrix.

## 2. Every scaling leaves an isotropic four-column witness

Take any nonzero column scaling and keep the name H for its kernel.
Intersect H with the ordinary coordinate-sum hyperplane:

```math
K=\{h\in H:\sum_{i=0}^5 h_i=0\}.
```

The dimension of K is at least two. Restricting any one coordinate to zero
therefore leaves a nonzero vector h. By the preceding argument, h has
weight four. Let S be its support.

In characteristic two,

```math
h^{\mathsf T}h=\sum_i h_i^2=\left(\sum_i h_i\right)^2=0.
```

The four-column submatrix C_S has rank three, so its kernel is exactly
the isotropic line spanned by h restricted to S. Its row space is the
orthogonal complement of that line. Because the line is isotropic, it
lies inside its own orthogonal complement. The row-space pairing is
therefore degenerate, with one-dimensional radical. Thus

```math
\mathrm{rank}(C_SC_S^{\mathsf T})=2<3
=\mathrm{rank}(C_S).
```

This proves the obstruction for every matrix under the stated three-column
independence hypothesis, not only for one enumeration of scalings.

Over a proper field extension E/F4, a two-dimensional E space has more
than five one-dimensional subspaces. The five further-zero lines no longer
cover it, and weight-five kernel words can occur. The proof consequently
does not contradict the deterministic extension-field construction.

## 3. Exact finite controls

The [checker](../tools/check_tensor_f4_hyperoval_obstruction.py) uses columns

```math
(1,t,t^2)^{\mathsf T}\quad(t\in\mathbb F_4),\qquad
(0,1,0)^{\mathsf T},\quad(0,0,1)^{\mathsf T}.
```

It checks all twenty independent triples and all 64 kernel words. Their
weights are zero once, four 45 times and six eighteen times. For each zero
coordinate, the checker reconstructs the five disjoint further-zero lines.

All 243 scalings normalized to have first coordinate one are checked.
They cover all 729 nonzero scalings before common-scalar normalization.
The [receipt](../evidence/tensor-f4-hyperoval-obstruction.json) supplies an
explicit isotropic four-column witness for every normalized scaling.
Of these scalings, 165 have three degenerate four-column subsets, 72 have
five, and six have fifteen. Every retained witness is checked by exact
matrix rank and its zero coordinate sum.

Run

```sh
python3 tools/check_tensor_f4_hyperoval_obstruction.py
```

This result extends the orthogonal differential-code interface at revision
`941cc1d8f56e6e80ec45938f1d660532de2f5046`. Substantial OpenAI assistant
assistance was used for the proof, implementation and text. It is an
ordinary proof with exact finite controls, without a new Lean claim.
