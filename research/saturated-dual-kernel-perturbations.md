# Low-degree crossing perturbations from saturated dual kernels

## Result and scope

Common polynomial factors in bounded-source kernel slabs can be removed
before the slabs are used as perturbations of another strict crossing.
Multiplying the resulting polynomial row by \(1+w\) restores its zero
endpoint. This produces seven degree-at-most-18 directions from the
weight-16 source, including one direction of degree one.

After independent logical placements and all allowed monomial shifts,
these directions form a 150-dimensional affine family over \(\mathbf F_4\)
around the five-symbol degree-18 reflection restriction. This is a new
exact construction space. Its completed one-coordinate scan and finite
two-coordinate prefix did not improve either displayed physical gauge.

The construction deliberately permits tensors outside the original
bounded-source operator image. Every resulting tensor is still polynomial
and has the literal strict crossing endpoint. No rational endpoint
specialization is used to justify it.

The source is the existing
[`bounded_dual_completion.py`](../tools/bounded_dual_completion.py) and its
[weight-16 record](../evidence/bounded-dual-completion-h16.json). The new
constructor is
[`saturated_dual_kernel_directions.py`](../tools/saturated_dual_kernel_directions.py).
The separate [full polynomial-module saturation](source-slab-module-saturation.md)
also uses cancellations between source generators and gives a larger space.
The present note concerns normalization of individual generators and the
reusable prescribed-direction search interface.

## 1. A polynomial normalization preserving the zero endpoint

Let \(K(w)\ne0\) be a binary polynomial slab with entries indexed by
\((y,u,v)\). Let \(g(w)\) be the monic greatest common divisor of its
entries. Define

\[
H(w)=(1+w)\frac{K(w)}{g(w)}.
\tag{1}
\]

Every entry of \(K/g\) is a polynomial. Thus \(H\) is polynomial and
\(H(1)=0\) entry by entry. This remains true when \(g(1)=0\): the
definition divides polynomials exactly before evaluating the new expression.
Its justification does not cancel a possibly zero value at the endpoint.

If the largest entry degree of \(K\) is \(d\), then

\[
\deg H=d-\deg g+1.
\tag{2}
\]

For any logical label \(a\), any nonnegative shift \(s\), and any strict
crossing \(B\), the tensor

\[
B(w;x,y,u,v)+\lambda [x=a]w^sH(w;y,u,v),
\qquad \lambda\in\mathbf F_4,
\tag{3}
\]

has the same strict endpoint as \(B\). If
\(s+\deg H\le D\) and \(\deg B\le D\), its degree is also at most
\(D\). Arbitrary linear combinations of these directions preserve both
properties.

## 2. The weight-16 kernel and the degree-18 family

The ten nonzero source slabs have the following exact normalization data.
Powers of \(1+w\) are taken in characteristic two.

| Source slab indices | Original degree | Common entry factor | Degree after (1) |
| --- | ---: | --- | ---: |
| 0 | 24 | \(w^{16}(1+w)^8\) | 1 |
| 3 | 30 | \((1+w)^{10}\) | 21 |
| 5 | 31 | \((1+w)^{11}\) | 21 |
| 7, 8 | 28 | \((1+w)^{12}\) | 17 |
| 11 | 31 | \((1+w)^{13}\) | 19 |
| 15, 16 | 30 | \((1+w)^{14}\) | 17 |
| 21, 23 | 31 | \((1+w)^{15}\) | 17 |

The first slab has one common polynomial on its twenty nonzero entries.
Its quotient is therefore a constant support tensor. The factor is
\(w^{16}(1+w^8)=w^{16}(1+w)^8\). This also agrees with the source's
identity-operator slab formula. The construction does not infer a new
degree-one tensor from sampled valuations; it divides every literal
polynomial entry exactly.

At degree bound \(D=18\), the seven admissible normalized slabs have
degrees

\[
1,17,17,17,17,17,17.
\]

The first contributes \(5\cdot18=90\) logical-placement and monomial-shift
directions. The others contribute \(6\cdot5\cdot2=60\). Exact binary
elimination on the complete coefficient vectors gives rank 150. The rank
is unchanged after extending scalars to \(\mathbf F_4\), since the matrix
entries lie in \(\mathbf F_2\) and its nonzero minors remain nonzero.

The [direction record](../evidence/saturated-dual-kernel-directions-h16.json)
retains every common factor, every normalized slab, the excluded degrees,
and every placement and shift. The baseline is the degree-18 strict crossing
on letters \(\{0,2,4,8,12\}\) described in the
[nonbinary reference-frame note](f4-recentered-crossing-lattices.md).

## 3. Projecting full row equations onto prescribed directions

The reusable adapter is
[`restricted_f4_crossing_fiber.py`](../tools/restricted_f4_crossing_fiber.py).
It applies to any literal endpoint-zero directions, not only the source
slabs in this note.

The complete coefficient constructor writes

\[
T=B+(1+w)Q,\qquad\deg Q_{xyuv}<D.
\]

For prescribed directions \(\Delta_j\), exact division gives
\(Q_j=\Delta_j/(1+w)\). Thus restricting to
\(T=B+\sum_j t_j\Delta_j\) replaces every full affine equation

\[
\ell(Q)=c
\]

by

\[
\sum_j\ell(Q_j)t_j=c.
\tag{4}
\]

Equation (4) is exact over \(\mathbf F_4\). The adapter first retains an
independent subset of the supplied literal directions, so its parameter
dimension is the dimension of the tensor image. It then projects every
coefficient equation in both complete row windows. Optional common floors
are eliminated in this smaller parameter space. The large complete-entry
floor system need not be solved first.

Sampled tensors are reconstructed in the original coefficients. The
checker verifies membership in the prescribed image, the strict endpoint,
the degree bound, the requested literal frame products, and both complete
profiles. A row-target bound alone is not treated as a pressure witness.

## 4. Sparse coefficient scans and Frobenius orbits

The source and all direction vectors are binary. Applying the field
automorphism \(\omega\mapsto\omega^2\) therefore conjugates only the
chosen scalar parameters. It preserves every minor's vanishing order at
zero and infinity. The complete profiles and both physical costs are
unchanged.

For one nonzero coordinate, the coefficient representatives are
\(1,\omega\). For two distinct nonzero coordinates they are

\[
(1,1),\ (1,\omega),\ (\omega,1),\
(\omega,\omega),\ (\omega,\omega^2).
\]

Thus \(2N\) one-coordinate queries cover all \(3N\) nonzero
one-coordinate assignments, up to this exact symmetry. A complete
two-coordinate scan has \(5\binom N2\) representatives. This symmetry
does not make a sparse scan exhaustive over higher-support assignments.

The driver
[`search_saturated_dual_directions.py`](../tools/search_saturated_dual_directions.py)
uses the complete native matrix profiles at two explicit paired-energy
choices: the natural \((0,10,6,2,8)\) and the two-zero
\((0,12,6,0,8)\) choices. These integer weights are supplied directly to
the profile engine. Pressure roots are exploratory numerical optimizations.

The [saved per-generator prefix](../evidence/saturated-dual-kernel-five-symbol-prefix.json)
contains all 300 one-coordinate representatives and 300 two-coordinate
representatives. It records 251 natural-gauge and 255 two-zero-gauge profile
classes, including the reference classes. Neither baseline ratio changed:
approximately \(0.7263479477\) and \(0.8335422390\), respectively.
The requested larger schedule was stopped when the full module saturation
provided a stronger family. Only saved completed counters are retained.

## 5. Independent controls

The [normalization controls](../evidence/saturated-dual-kernel-direction-controls.json)
verify 1,250 exact polynomial identities
\(gH=(1+w)K\), all 93,750 entries of the 150 endpoint-zero directions,
their exact rank, and twelve complete profiles through polynomial regular
frames. Twelve further profile comparisons check Frobenius conjugation on
the recorded examples.

The [adapter controls](../evidence/restricted-f4-crossing-fiber-controls.json)
verify three image ranks through binary restriction of scalars, sixty
literal span reconstructions, two full-parameter equation identities, and
three invalid-endpoint rejections.

```sh
python3 tools/check_saturated_dual_directions.py \
  evidence/saturated-dual-kernel-directions-h16.json \
  --output /tmp/saturated-dual-kernel-direction-controls.json
python3 tools/check_restricted_f4_crossing_fiber.py \
  --output /tmp/restricted-f4-crossing-fiber-controls.json
```

Saved search witnesses are replayed by
[`check_saturated_dual_search.py`](../tools/check_saturated_dual_search.py).
It reconstructs every sparse perturbation literally and certifies complete
profiles with polynomial regular frames, without calling the native Smith
oracle. Historical implementation hashes are resolved separately against
the worktree or exact local git blobs; recorded hashes are not rewritten.
