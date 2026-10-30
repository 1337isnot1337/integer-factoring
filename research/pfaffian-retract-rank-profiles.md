# Exact ranks of Pfaffian retract families

## Result and scope

The three-variable Pfaffian family admits explicit linear retracts that
cancel one or both noncrossing vertices. The complete coefficient spaces
of the finite retracts checked here still fail the factoring seed
criterion
\[
S=\sum_i\sqrt{r_i c_i}<d^{\sqrt2}.
\]

This note records the two cut ranks of every binary coefficient
combination. It separates the algebraic endpoint construction from the
remaining rank requirement. The optimal common-retract dimension theorem
is proved in
[noncommutative-crossing-optimal-retracts.md](noncommutative-crossing-optimal-retracts.md).
The universal family and its coefficient-space interpretation are defined
in [pfaffian-ternary-edge-family.md](pfaffian-ternary-edge-family.md).

The tables below concern recombinations over \(\mathbb F_2\); Section 6
also excludes three complete coefficient spaces over \(\mathbb F_4\).
These results do not exclude other label dimensions, other linear maps,
or untested coefficient fields. The retained finite constructions remain
valid over every characteristic-two field.

## 1. Conventions and exact endpoint checks

Use degree-\(L\) monomials on \(P\) complementary pairs, with
\[
D=\binom{2P}{L}.
\]
The involution \(J\) swaps the two modes within each pair. In local
occupation notation it acts by \(0\mapsto0,\ 1\leftrightarrow2,\ 3\mapsto3\).
It is not whole-set complementation.

For the universal tensor \(G_L(p,q,r)\), the pure coefficients are
\[
[p^L]G=[x=u][y=v],\qquad
[q^L]G=J_{x,y}J_{u,v},\qquad
[r^L]G=[x=v][y=u].
\]
A common input map \(I\) on \(y,v\) and output map \(P\) on \(x,u\)
therefore give
\[
[p^L]G_{P,I}=(PP^{\mathsf T})_{x,u}(I^{\mathsf T}I)_{y,v},
\quad
[q^L]G_{P,I}=(PJI)_{x,y}(PJI)_{u,v}.
\]
The crossing vertex remains exact when \(PI=I_d\).

The checker verifies the map identities and also checks the three
endpoints after expanding the complete polynomial tensor. Its two cuts
are
\[
(x,y)\mid(v,u),\qquad (x,u)\mid(v,y).
\]
The crossing is the identity matrix in both conventions.

### Sparse double cancellation

Choose one representative from each nonfixed \(J\)-orbit. Pair distinct
representatives \(a_i,b_i\), and use
\[
I_i=e_{a_i}+e_{b_i},\qquad P_i=e_{a_i}^{\mathsf T}.
\]
These maps satisfy
\[
PI=I_d,\qquad I^{\mathsf T}I=0,\qquad PJI=0.
\]
The canonical sparse examples pair representatives with the same
single-occupation syndrome. This retains the exact XOR block
decomposition of both cuts. A syndrome with an odd number of available
representatives leaves one unused.

### Triangular double cancellation

The triangular constructions from the optimal-retract theorem retain
\(m-1\) labels when there are \(m\) nonfixed orbits and no fixed monomial,
and \(m\) labels when a fixed monomial is available. They use more terms
per input and output map. The exporter also implements groups of
triangular blocks and distinct fixed anchors, so their support size can
be reduced without changing the endpoint identities.

### Keeping the \(q^L\) endpoint

Another sparse construction uses the \(J\)-orbit itself:
\[
I_a=e_a+e_{Ja},\qquad P_a=e_a^{\mathsf T}.
\]
Pairing fixed monomials in the same way gives a total of
\(\lfloor D/2\rfloor\) labels. Now
\[
PI=I_d,\qquad I^{\mathsf T}I=0,\qquad PJI=I_d.
\]
Thus the \(p^L\) vertex vanishes, the \(r^L\) vertex is the crossing, and
the \(q^L\) vertex is
\[
C(x,y,u,v)=[x=y][u=v].
\]
Its cut ranks are \((1,d^2)\). Keeping this inexpensive term is a different
finite family from double cancellation.

## 2. Complete binary coefficient-space certificates

A homogeneous degree-\(L\) tensor has at most \(\binom{L+2}{2}\)
coefficient tensors. The checker substitutes
\[
p=z,\qquad q=z^{L+1},\qquad r=1.
\]
Every monomial \(p^a q^b r^{L-a-b}\) receives the distinct exponent
\(a+(L+1)b\). This retains all coefficients, including every contact
term in the local Pfaffian formula.

After selecting an independent binary basis \(B_0,\ldots,B_{k-1}\),
with \(B_0\) the crossing, the checker forms all \(2^k\) combinations.
For mask \(v\), let
\[
w(v)=\sqrt{\operatorname{rank}_R B_v\,
                    \operatorname{rank}_C B_v}.
\]
A decomposition of the crossing is exactly a collection of masks
whose XOR is \(1\). Its minimum cost is the shortest path from \(0\)
to \(1\) in the additive group \(\mathbb F_2^k\), with edge weights
\(w(v)\).

The reported decimal cost is a floating-point shortest-path result
together with an exact-rank decomposition witness. The exclusions use a
separate integer computation with weights
\[
w_{\rm int}(v)=
\left\lfloor\sqrt{\operatorname{rank}_R B_v\,
                         \operatorname{rank}_C B_v}\right\rfloor.
\]
Its exact shortest path is a lower bound for every real-cost
decomposition. If the resulting integer \(B\) satisfies
\[
B^{70}>d^{99},
\]
then every decomposition fails the seed threshold, because
\(\sqrt2<99/70\). The conclusion does not depend on numerical comparison
near the threshold.

## 3. Finite results

### Sparse common retracts

| Mode pairs \(P\) | Length \(L\) | Labels \(d\) | Coefficient-span dimension | Found cost | Integer lower bound | Seed threshold |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 3 | 5 | 3 | 25.000000 | 25 | 9.738518 |
| 4 | 2 | 6 | 4 | 35.678945 | 34 | 12.602945 |
| 4 | 3 | 12 | 8 | 144.000000 | 144 | 33.588666 |
| 4 | 4 | 16 | 11 | 242.947428 | 242 | 50.452514 |

The full ranks, coefficient monomials, sparse maps, and decomposition
masks are recorded in
[pfaffian-sparse-double-endpoint-retracts.json](../evidence/pfaffian-sparse-double-endpoint-retracts.json).

### Triangular common retracts

| Mode pairs \(P\) | Length \(L\) | Labels \(d\) | Coefficient-span dimension | Found cost | Integer lower bound | Seed threshold |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 3 | 9 | 8 | 80.835348 | 80 | 22.361591 |
| 4 | 2 | 12 | 4 | 144.000000 | 144 | 33.588666 |
| 4 | 3 | 27 | 8 | 729.000000 | 729 | 105.743589 |

The larger dimension is algebraically useful, but these particular
triangular orderings retain high intermediate ranks. The complete finite
records are in
[pfaffian-triangular-retracts.json](../evidence/pfaffian-triangular-retracts.json).

For the middle grade on eight modes, the maximal triangular retract has
\(d=32\) and a thirteen-dimensional coefficient span. A complete native
screen checks all \(8191\) nonzero binary combinations. Exactly three have
cost below the integer gate \(135\):

| Binary basis mask | Row rank | Column rank |
|---:|---:|---:|
| 15 | 1 | 196 |
| 2321 | 225 | 1 |
| 8191 | 108 | 95 |

Their three-dimensional span does not contain the crossing, whose basis
mask is \(1\). Since \(32^{\sqrt2}<135\), any successful decomposition
would have to use only these three combinations. This gives a complete
finite binary exclusion without computing an optional unrestricted
minimum cost. Each eligible rank pair was independently checked by the
Python binary elimination routine. The certificate is in
[pfaffian-triangular-P4-F2-exclusion.json](../evidence/pfaffian-triangular-P4-F2-exclusion.json),
with the complete capped-screen events in
[pfaffian-triangular-P4-F2-screen.jsonl](../evidence/pfaffian-triangular-P4-F2-screen.jsonl).

### Invariant-orbit retracts retaining \(C\)

| Mode pairs \(P\) | Length \(L\) | Labels \(d\) | Coefficient-span dimension | Found cost | Integer lower bound | Seed threshold |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 2 | 3 | 2 | 7.898979 | 7 | 4.728804 |
| 3 | 2 | 7 | 5 | 38.749016 | 38 | 15.672891 |
| 3 | 3 | 10 | 8 | 77.082039 | 77 | 25.954554 |
| 4 | 2 | 14 | 5 | 145.757910 | 145 | 41.770513 |
| 4 | 3 | 28 | 9 | 550.568990 | 550 | 111.324439 |
| 4 | 4 | 35 | 13 | 804.170969 | 803 | 152.630726 |

These records are in
[pfaffian-invariant-orbit-retracts.json](../evidence/pfaffian-invariant-orbit-retracts.json).
All three tables fail the integer certificate above.

## 4. Exhausting output orientations in the ten-label example

The \(P=L=3\) ambient space has ten nonfixed \(J\)-orbits. Each input
is the orbit sum \(e_a+e_{Ja}\), but the output coordinate can be either
member of its orbit.

Simultaneously reversing every output representative leaves the tensor
unchanged. To see this, apply global \(J\) to all four legs of the
original tensor; the local pairing formula is invariant, and both input
legs are already \(J\)-invariant. It therefore suffices to check
\(2^9=512\) output choices.

The complete sweep gives a minimum integer shortest-path lower bound
of \(73\), above \(10^{\sqrt2}=25.954553\ldots\). The best found
real-cost decomposition is
\[
\sqrt{35\cdot38}+\sqrt{84\cdot9}+\sqrt{1\cdot100}
=73.964619227\ldots.
\]
Thus changing only these output representatives cannot produce a
successful binary seed at this size.

The orientation masks, individual lower bounds, rank-table digest, and
best decomposition are in
[pfaffian-invariant-orbit-orientations.json](../evidence/pfaffian-invariant-orbit-orientations.json).

## 5. Reproduction and independent contraction checks

The source is
[check_double_endpoint_retracts.py](../tools/check_double_endpoint_retracts.py).

The sparse evaluator expands four original monomial entries per
contracted tensor entry. The dense evaluator instead builds the original
coefficient matrices, applies the four binary linear maps with bitset
matrix operations, and reshapes the result into the second cut.
A built-in audit compares both methods on every matrix entry of both
cuts for the five-label sparse example.

The following commands regenerate the finite evidence:

~~~sh
python3 tools/check_double_endpoint_retracts.py \
  --output evidence/pfaffian-sparse-double-endpoint-retracts.json

python3 tools/check_double_endpoint_retracts.py \
  --cases full:3:3,full:4:2,full:4:3 --construction triangular \
  --output evidence/pfaffian-triangular-retracts.json

python3 tools/check_double_endpoint_retracts.py \
  --cases full:2:2,full:3:2,full:3:3,full:4:2,full:4:3,full:4:4 \
  --construction invariant \
  --output evidence/pfaffian-invariant-orbit-retracts.json

python3 tools/check_double_endpoint_retracts.py \
  --cases full:3:3 --construction invariant --orientation-sweep \
  --output evidence/pfaffian-invariant-orbit-orientations.json
~~~

Add \(\texttt{--export-dir DIRECTORY}\) to export independent coefficient
bases in both JSON and native bit-matrix formats. The durable small
bases are:

- [P3 L3 sparse basis](../evidence/pfaffian-retract-basis-P3-L3-sparse.json);
- [P3 L3 triangular basis](../evidence/pfaffian-retract-basis-P3-L3-triangular.json);
- [P4 L4 sparse basis](../evidence/pfaffian-retract-basis-P4-L4-sparse.json).

For a larger case, \(\texttt{--skip-optimization}\) exports its complete
basis without enumerating every binary combination. A native
coefficient-span screen can then test whether all pieces below the
required cost gate span the crossing.

No new asymptotic factoring exponent follows from these finite results.
The constructive advance is the availability of larger, explicitly
computable endpoint-cancelling label spaces; a successful family still
needs a substantially smaller complete two-cut rank sum.

## 6. Complete four-element-field coefficient searches

The native [four-element-field screen](../tools/screen_pfaffian_span_f4.cpp)
enumerates projective vectors of the complete coefficient space over
`F2[t]/(t^2+t+1)`. It stops a rank calculation only when the established
partial rank product already exceeds the allowed cost gate. The
[independent scalar checker](../tools/check_pfaffian_span_f4.py) recomputes
every eligible vector and additional controls using separate elimination
code and the opposite pivot direction.

| Retract | Labels | Span dimension | Projective vectors | Below integer gate | Dimension of their span |
|---|---:|---:|---:|---:|---:|
| P3 L3 sparse | 5 | 3 | 21 | 1 | 1 |
| P3 L3 triangular | 9 | 8 | 21,845 | 2 | 2 |
| P4 L4 sparse | 16 | 11 | 1,398,101 | 41 | 5 |

The integer gates are 10, 23, and 51, respectively; each is strictly
larger than the corresponding seed threshold. In all three cases the
span of every vector below that gate omits crossing. Any decomposition
with total cost below the seed threshold would use only such vectors,
which is impossible. This excludes arbitrarily many pieces and all
scalar recombinations over the four-element field in these spaces.

The [combined scalar audit](../evidence/pfaffian-retract-F4-audit.json)
records exact ranks, span bases, endpoint remainders, and source hashes.
The complete native logs are
[P3 sparse](../evidence/pfaffian-retract-P3-sparse-F4.jsonl),
[P3 triangular](../evidence/pfaffian-retract-P3-triangular-F4.jsonl), and
[P4 sparse](../evidence/pfaffian-retract-P4-sparse-F4.jsonl).
The scalar audit checks all 21 projective vectors in the smallest case,
and 42 and 81 vectors in the other two cases.
