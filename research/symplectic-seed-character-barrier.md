# A common-character barrier for every power of the 24- and 26-label seeds

The full six-direction coefficient spaces of the 24- and 26-label symplectic
seeds cannot pass the finite-seed criterion over any field of characteristic
two. This holds for every tensor power, mixed products of the two seeds,
independent parameter curves in different factors, and arbitrary linear
recombinations in the resulting tensor-product space.

The obstruction uses one linear functional shared by the two flattenings.
The row cut has a large triangular block with that functional on its
diagonal. The column cut has a smaller constant matrix compression giving
the same functional times the identity. These two rank certificates
tensorize, so one nonzero scalar forces both ranks to be large.

This is a theorem about these two coefficient spaces. It is not a lower
bound for integer factoring, other label codes, or all crossing tensors.

## Statement

Let \(V_d\), for \(d\in\{24,26\}\), be the six-dimensional tensor space
spanned by the normalized seed generators

\[
E_d,A_d,B_d,C_d,D_d,F_d.
\]

The normalization is the one in the
[weighted seed construction](symplectic-seed-weighted-powers.md):

\[
K_d(u,v,s)=u^3E_d+u^2vA_d+uv^2B_d+v^3C_d+usD_d+vsF_d.
\]

Here \(E_d\) is the crossing tensor. The 24-label code is the exact logical
restriction of the 26-label code to label indices \(2,\ldots,25\).
In the exterior-algebra description, it keeps precisely the labels whose
syndrome has two single sites. For each pair of even sites \(a<b\), its
input is the product of the two single-site modes times \(D_a+D_b\), and
its output selects the term containing \(D_b\). There are six choices of
the even-site pair and four single-site orientations, hence 24 labels.

Let \(k\) be any characteristic-two field. Fix integers \(n,m\geq0\) with
\(n+m\geq1\), set

\[
\mathcal V=V_{24}^{\otimes n}\otimes V_{26}^{\otimes m},
\qquad d=24^n26^m,
\]

and extend scalars to \(k\). For every finite decomposition

\[
E_{24}^{\otimes n}\otimes E_{26}^{\otimes m}=\sum_i T_i,
\qquad T_i\in\mathcal V,
\]

the two-cut rank cost satisfies

\[
S:=\sum_i\sqrt{r(T_i)c(T_i)}
\ \geq\ 16512^{n/2}18748^{m/2}
\ >\ d^{3/2}
\ >\ d^{\sqrt2}.
\tag{1}
\]

The pieces \(T_i\) need not themselves be pure tensors. In particular,
the result allows the entire tensor-product coefficient span, which can
be larger than the coefficient span of a common-parameter tensor power.

## The shared functional

Define a linear functional on either seed space by

\[
\chi(E,A,B,C,D,F)=(1,0,1,0,0,0).
\tag{2}
\]

On the row coefficient algebra, this functional is a character. Its
multiplicity in a verified simultaneous triangular representation is

\[
a_{24}=384,\qquad a_{26}=436.
\tag{3}
\]

More explicitly, there are fixed changes of basis and fixed principal
index sets of these sizes such that the restricted matrix of each
generator is upper triangular, with the corresponding entry of (2) at
every diagonal position. Therefore, for every \(T\in V_d\), this
principal submatrix is upper triangular with diagonal \(\chi(T)\).
If \(\chi(T)\ne0\), it is nonsingular and

\[
r(T)\geq a_d.
\tag{4}
\]

The column coefficient algebra does not need to have \(\chi\) as a
character. The required certificate is only a constant matrix compression:
there are matrices \(W_d\) and \(V_d'\) with 43 rows and 43 columns,
respectively, such that

\[
W_d M_c(X)V_d'=\chi(X)I_{43}
\qquad(X=E_d,A_d,B_d,C_d,D_d,F_d).
\tag{5}
\]

Here \(M_c\) denotes the column flattening; the prime distinguishes the
right compression matrix from the tensor space \(V_d\).
By linearity, (5) holds for every tensor in the seed span. When its scalar
is nonzero, rank monotonicity gives

\[
c(T)\geq43.
\tag{6}
\]

All changes of basis and compression matrices are defined over
\(\mathbb F_2\). Their verified identities remain valid after extension
to any characteristic-two field.

## Explicit column witnesses

The 24-label column representation has a verified constant direct-sum
decomposition. The following module types supply one scalar compression
per copy. A decimal integer encodes a binary vector by its set bits, with
bit zero denoting the first coordinate. The type identifiers refer to
`families.full.types.column` in the
[module certificate](../evidence/symplectic-radical-sigma2-modules-certificate.json).

| Type | Module dimension | Copies | Right vector \(v\) | Left covector \(w\) |
|---:|---:|---:|---:|---:|
| 0 | 8 | 2 | 2 | 3 |
| 1 | 8 | 1 | 16 | 20 |
| 5 | 8 | 12 | 32 | 36 |
| 6 | 8 | 8 | 64 | 80 |
| 8 | 8 | 4 | 32 | 36 |
| 9 | 6 | 16 | 8 | 10 |

For each listed type, direct multiplication gives

\[
\bigl(w^{\mathsf T}Ev,w^{\mathsf T}Av,w^{\mathsf T}Bv,
      w^{\mathsf T}Cv,w^{\mathsf T}Dv,w^{\mathsf T}Fv\bigr)
=(1,0,1,0,0,0).
\]

There are \(2+1+12+8+4+16=43\) disjoint summands. Placing one witness in
each summand makes the off-diagonal compression entries zero, yielding
(5) for \(d=24\). The independent checker also assembles the actual
43 covectors and vectors in the adapted direct sum and verifies all
\(6\cdot43^2=11094\) compression entries explicitly.

The 26-label column certificate follows by restriction. Both flattenings
of every one of the six 24-label generators are exact principal
restrictions of the corresponding 26-label generator. Extending the
24-label left and right compression maps by zero therefore proves (5)
for \(d=26\) too. This restriction is checked in all 6,912 supplied bit
rows, and the logical inputs and selected outputs are checked separately.

## Tensorization and proof of the cost bound

Write \(\Xi=\chi^{\otimes(n+m)}\) on \(\mathcal V\), using the
appropriate seed space in each factor. Taking tensor products of the
row principal restrictions gives an upper triangular matrix of size

\[
a=384^n436^m
\]

whose every diagonal entry is \(\Xi(T)\), for arbitrary
\(T\in\mathcal V\). This remains true for arbitrary linear combinations
because the restrictions and diagonal extraction are linear. Hence
\(\Xi(T)\ne0\) implies \(r(T)\geq a\).

Taking tensor products of (5) gives a constant compression

\[
\left(W_{24}^{\otimes n}\otimes W_{26}^{\otimes m}\right)
M_c(T)
\left((V_{24}')^{\otimes n}\otimes(V_{26}')^{\otimes m}\right)
=\Xi(T)I_{43^{n+m}}.
\]

Thus the same nonzero scalar also implies \(c(T)\geq43^{n+m}\).
Since \(\chi(E_d)=1\), a decomposition of the crossing satisfies
\(\sum_i\Xi(T_i)=1\). At least one piece has \(\Xi(T_i)\ne0\), and
its contribution alone is at least

\[
\sqrt{384^n436^m43^{n+m}}
=16512^{n/2}18748^{m/2}.
\]

Finally, the exact integer inequalities

\[
16512>24^3=13824,\qquad18748>26^3=17576
\]

give (1). The strict comparison \(3/2>\sqrt2\) is equivalent to
\(9>8\); no floating-point threshold is used.

## What this closes and what remains open

Every coefficient obtained from scalar polynomial or rational parameter
curves in the two normalized seed families lies in the corresponding
tensor-product span after clearing denominators. Coefficient basis changes,
different local curves, and extension fields therefore cannot evade (1).
The result subsumes the earlier binary square and cube exclusions as well
as the checked mixed-curve exclusions. It also covers tensor powers that
would be impractical to enumerate.

The argument does not apply automatically to larger symplectic codes,
changed input or output maps, or pieces outside the specified six-generator
tensor-product spaces. For example, changing the number of sites in the
two-even-site construction changes the seed, so the four-site certificate
does not settle that family at all sizes.

## Verification and reproduction

The proof uses exact finite witnesses with a separate reconstruction audit,
not a search-based extrapolation. The checker reconstructs all generator
matrices from the literal seed blocks and verifies:

- the 24-label row and column algebras and their adapted representations;
- the 26-label row algebra and its adapted representation;
- all 2,590 algebra-basis products for those three representations;
- all 10,968 actual change-of-basis column identities;
- 4,920 row triangular diagonal checks;
- 3,456 full-generator column-module intertwining identities;
- all 11,094 entries of the assembled column compression;
- all 6,912 generator restriction bit rows and the exact integer margins.

The source hashes bind the input evidence. The checker treats the two word
encodings used by the earlier algebra generators equivalently; it does not
alter their frozen certificates.

Run from the repository root:

```sh
python tools/check_noncommutative_crossing_character_compression.py \
  --output /tmp/noncommutative-crossing-character-compression.json
```

The committed report is
[noncommutative-crossing-character-compression.json](../evidence/noncommutative-crossing-character-compression.json).
The 24-label code, coefficient algebra, and independent module replay are
provided by
[check_symplectic_radical_sigma2.py](../tools/check_symplectic_radical_sigma2.py)
and
[check_symplectic_radical_sigma2_modules.py](../tools/check_symplectic_radical_sigma2_modules.py).
