# Bounded Frobenius modules at ten sites

## Outcome

The explicit natural-module decomposition improves the ten-site column
certificate from 5,052 to **5,738** without exceeding the existing matrix budget.
Together with the [small-core row certificate](two-even-site-row-corner-barrier.md),
it gives

\[
R=44638912,\qquad C=5738,\qquad RC=256138077056
\]

for the 11,520-label two-even-site seed. These are valid common-functional rank
certificates over every characteristic-two field and at every tensor power.
Their value does **not** reach the target:

\[
256138077056^{16}<11520^{45},\qquad
(45/32)^2<2.
\]

Thus the cost exponent supplied by this certificate is below \(45/32<\sqrt2\)
(approximately 1.40448). This statement does not exhibit a cheap decomposition.
The largest audited matrix has dimension 5,120; no full ten-site matrix was
constructed.

A leading-coefficient calculation gives a column ceiling of 10,638 for the
chosen principal restrictions. The common coefficient \(A_8+B_6\) improves
that ceiling to 9,810. These ceilings remain too large to establish that every
possible common-functional compression on these restrictions must fail. The
bounded search stops with this precise unresolved gap.

## The free factor has an actual direct-sum decomposition

The normalized local column factor has the constant decomposition

\[
I_4+tN\cong\mathbf1\oplus\mathbf1\oplus U(t),
\qquad U(t)=\begin{pmatrix}1&t\\0&1\end{pmatrix}.
\]

The natural-cube identity from
[the folded-Gram construction](folded-gram-channels.md) is a constant module
isomorphism:

\[
V^{\otimes3}\cong V\otimes V^{[2]}
 \oplus2(\det V\otimes V).
\]

For \(U(t)\), the determinant is one and the Frobenius factor is \(U(t^2)\).
The new checker replays the generic 64-entry polynomial identity and its
invertible binary change of basis. It then composes these explicit changes of
basis to split every \(U(t)^{\otimes j}\) for \(0\le j\le6\). It checks the
complete polynomial intertwining identities, not only character multiplicities.

Write \(B_k(t)\) for the tensor product of \(U(t^w)\) over the recursive
Frobenius weights \(\mathcal W_k\). The relevant cases are

| Block | Frobenius weights | Dimension |
| --- | --- | ---: |
| \(B_4\) | \((1,1,2)\) | 8 |
| \(B_5\) | \((1,2,2)\) | 8 |
| \(B_6\) | \((1,1,2,2)\) | 16 |

For example,

\[
U(t)^{\otimes5}\cong B_5(t)\oplus4B_3(t)\oplus4B_1(t),
\]

\[
U(t)^{\otimes6}\cong B_6(t)\oplus4B_4(t)\oplus4B_2(t).
\]

All maps are constant and defined over \(\mathbb F_2\), so attaching a core
preserves the direct sums and their coefficient identities over every extension
field.

## Degree removes most summands

The selected functional has scalar polynomial

\[
\chi_{10}(G_{10})=(1+t)^8=1+t^8,
\]

or, in the original homogeneous variables,
\((q+r)^2(p+r)^8\). In particular, \(\chi_{10}(A_8)=1\).

The four-site weight-four column core has \(t\)-degree at most two. Among its
six free factors, only \(B_6\) can reach coefficient \(A_8\). It appears once,
and the attached matrix has dimension \(96\cdot16=1536\).

The five-site column cores have \(t\)-degree at most four. Among their five
free factors, only \(B_4\) and \(B_5\) can reach \(A_8\). Expanding
\((2\mathbf1\oplus U(t))^{\otimes5}\) gives multiplicities ten and one,
respectively. The weight-two matrices have dimensions
\(384\cdot8=3072\); the zero-syndrome matrices have dimensions
\(640\cdot8=5120\).

A discarded module has zero \(A_8\) matrix. It cannot have a positive-dimensional
scalar compression with this functional. This justifies omitting it from the
blockwise witness search, while preserving the actual direct-sum identity.

## Exact bounded search

The same tagged discovery routine proposes coordinate subspaces in the smaller
Frobenius bases. The independent low-pivot audit reconstructs every left map and
checks all entries of \(WM_iV=\chi_iI\).

The weight-four core contributes 16 column coordinates per physical block.
For five-site cores, the fresh weighted dimensions are recorded below. A
\(B_4\) witness is counted ten times; a \(B_5\) witness is counted once.

| Relative syndrome | Weighted column compression |
| ---: | ---: |
| 0 | 6 |
| 3 | 8 |
| 5 | 62 |
| 6 | 8 |
| 9 | 64 |
| 10 | 66 |
| 12 | 8 |
| 17 | 8 |
| 18 | 66 |
| 20 | 64 |
| 24 | 8 |

The physical anchor list chooses among these exact witnesses for each
weight-two syndrome. Summing the retained blocks gives 5,738. The report contains
the exact coordinates, module multiplicities, direct-sum maps, and the
reconstructed left-map hashes.

The checker also evaluates 1,536 complete ten-site polynomial embedding controls:
1,260 in entire weight-four blocks, 270 in the selected weight-two principal
restrictions, and six in the zero-syndrome restriction. The universal embedding
claim still rests on the all-size core factorization and explicit anchor maps.

## What the rank ceilings establish

Any scalar compression with \(\chi_{10}(A_8)=1\) has dimension at most the
rank of its \(A_8\) matrix. Applying this to the full chosen core restrictions,
then allowing the best physical anchors, gives the exact ceiling 10,638.

The same argument applies to

\[
X=A_8+B_6,\qquad\chi_{10}(X)=1.
\]

The checker computes its rank in every relevant Frobenius module. Both \(A_8\)
and \(B_6\) vanish in every discarded module: the largest core degrees are
\((2,0)\) in the four-site block and \((4,2)\) in the five-site blocks for the
\(A\)- and \(B\)-directions. Therefore the coefficient-rank ceiling applies
before discarding those summands, even to compressions that mix them. The
resulting uniform ceiling is 9,810.

These are ceilings on the selected core restrictions and the selected functional.
They do not apply to larger cores, other output sections, or different
functionals. They also leave room above the present 5,738-dimensional witness.
No optimality claim is made.

## Reproduce

```bash
python tools/check_two_even_site_frobenius_probe.py \
  --output /tmp/two-even-site-P10-frobenius.json
```

Add `--discover` to reproduce the coordinate choices. A certified five-site basis
can be supplied with `--basis5`; otherwise the checker regenerates it literally.
The [certificate](../evidence/two-even-site-P10-frobenius-certificate.json)
and [audit report](../evidence/two-even-site-P10-frobenius.json) record the
complete bounded witnesses and their verification.
