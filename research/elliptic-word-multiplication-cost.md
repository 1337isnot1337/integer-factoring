# Costs of multiplication in cyclic elliptic word encodings

Nonlinear word encodings can change the native-field arithmetic question.
For a cyclic elliptic subgroup, successor is translation, an algebraic map
of fixed complexity. This avoids the assumption of a binary-linear word
encoding. It does not make the complete ring operations cheap.

This note gives two quantitative ordinary bounds for the encoded
**multiplication modulo `M=2^w`** interface:

* Coordinate charts need
  `sum_j max(9(d_(j,X)+d_(j,Y)),1+w/2) >= M`.
* A nonzero direct equation of separate degree at most `D` on the three
  curve words, vanishing on every multiplication witness, needs
  `M <= 246 D^3`.

Thus polynomially many charts of polynomial word degree, or constant-degree
direct gate equations, do not provide the intended packed multiplication
interface. No curve-parameter or field-size search follows this audit.
These are representation bounds, not factoring lower bounds. The complete
score remains **1/2**; no new Lean work is used.

Auxiliary equations whose **finite-field** solvability selects the graph
remain outside the direct-equation theorem. Their projection can be dense
over the algebraic closure. This distinction identifies the next possible
arithmetic mechanism; it is not a completed low-cost representation.

## 1. Exact encoding and standard geometric inputs

Work in `k=algebraic_closure(F2)`. Let `E` be an ordinary smooth plane
elliptic cubic, `P` have exact order `M=2^w`, and `G=<P>`. Fix any `Q0`
and encode the word `a mod M` by

```
iota(a)=Q0+[a]P.
```

Successor is translation by `P`. The proposed product output is
`iota(ab)`. All points and coefficients lie in some finite extension,
but the cost of finding a curve, generator, encoding and prefix interface
is not assumed free. The bounds below hold even if those are supplied.

Use three standard facts from J. S. Milne,
[*Abelian Varieties*, version 2.00, March 2008](https://www.jmilne.org/math/CourseNotes/AVc.pdf):
I Theorem 3.2 extends rational maps from smooth varieties to abelian
varieties; I Corollary 3.7 decomposes a rational map from a group variety
into a homomorphism and a translation; I Theorem 7.2 gives degree `a^2`
and separability for odd multiplication maps on `E` in characteristic two.
Printed pages 17, 19 and 33 contain these results. The quantitative finite-grid
arguments below are additional arguments, not statements quoted from Milne.

A nonzero section of a line bundle of degree `d` on a smooth projective
curve has at most `d` distinct zeros. For `L=O_E(1)`, `deg L=3`.
Ordinarity gives exactly `2^j` geometric `2^j`-torsion points. In particular
the geometric `M`-torsion is exactly `G`. For a standard ordinary Weierstrass model in characteristic two,
there is one nonidentity geometric two-torsion point; each doubling fibre
has at most two points. Induction and the given cyclic subgroup give the
same assertion at level `M`.

## 2. A rational map agrees with ring multiplication only sparsely

**Lemma.** A rational map `f:E*E -> E` agrees with `iota(ab)` on at most

```
M(1+w/2)
```

of the `M^2` pairs `(iota(a),iota(b))`.

The standard decomposition above writes its values on this grid as

```
[a]A+[b]B+C
```

for fixed points `A,B,C`. Absorb the encoding translation into `C`.
The finitely many points involved belong to one finite group of curve
points. Project that group to its two-primary subgroup; this fixes `P`.
For a fixed `a`, agreement requires

```
[b](aP-B_2) = aA_2+C_2.                                 (1)
```

If `B_2` is outside the geometric `M`-torsion subgroup `G`, its order is greater than `M`.
Subtracting `aP` preserves that order, so the `M` possible values of `b`
in (1) are distinct. There is at most one agreement per row.

Otherwise write `B_2=beta P`. Equation (1) has at most
`gcd(a-beta,M)` solutions in `b`, irrespective of its right side. Therefore

```
agreements <= sum_(a mod M) gcd(a-beta,M)
           = M + w M/2.                                  (2)
```

For each valuation `j=0,...,w-1`, there are `M/2^(j+1)` residues with
gcd `2^j`, contributing `M/2`. The remaining residue contributes `M`.
The bound is sharp: the constant map with value `iota(0)` agrees when
`ab=0 mod M`.

## 3. Piecewise coordinate formulas have an exponential tradeoff

A chart is a triple of bihomogeneous polynomials on `P2*P2`, all of the
same bidegree `(d_X,d_Y)`. At a point where the triple is nonzero it
specifies a projective output. Branch conditions may choose any chart;
their computation cost is ignored, which only strengthens this bound.

Compose the triple with the cubic equation of `E`. If the resulting
section is identically zero on `E*E`, the chart is a rational map into
`E` and the preceding lemma bounds its correct pairs by `M(1+w/2)`.
A triple identically zero on `E*E` has no defined pairs and can be omitted.

Otherwise that section has degrees `9d_X,9d_Y` on the two curve factors.
At most `9d_X` rows can be entire zero fibres: evaluate at one fixed
second-coordinate point giving a nonzero horizontal section. Each other
row has at most `9d_Y` zeros. Thus at most

```
9(d_X+d_Y) M
```

grid points can even have output on the curve, and hence can be correct.
This counts basepoints as potential zeros, so undefined outputs do not
invalidate the upper bound.

If `K` charts cover every product pair, summing these agreement bounds gives

```
sum_j max(9(d_(j,X)+d_(j,Y)),1+w/2) >= M.                 (3)
```

In particular, if both input degrees are at most `d`, then
`K max(18d,1+w/2)>=2^w`. This statement allows arbitrary branch selection;
it is not restricted to formulas defined at every grid point.

## 4. Direct polynomial graph equations need large degree too

**Theorem.** Let `s` be a section on `E*E*E` coming from a
trihomogeneous polynomial with each separate degree at most `D>=1`.
Suppose `s` vanishes on every

```
(iota(a),iota(b),iota(ab)),       a,b mod M.
```

If its restriction to `E*E*E` is nonzero, then

```
M <= 246 D^3.                                            (4)
```

Assume instead `M>246D^3`. Choose `3D+1` odd integers

```
a=2D+1,2D+3,...,8D+1.
```

They satisfy `a<=9D<M` and `a^2>3D`; their encoded first inputs are
distinct. For each, let `s_a` be the section on the last two curve factors
after fixing the first input to `iota(a)`. Define

```
h_a(Q)=[a]Q+(1-a)Q0,
```

so `h_a(iota(b))=iota(ab)`. Pulling `s_a` back along
`Q -> (Q,h_a(Q))` gives a curve section of degree at most

```
3D(1+a^2) <= 3D(1+81D^2) <= 246D^3 < M.
```

It has the `M` distinct zeros `iota(b)`, so it is identically zero.
For any fixed output point `R`, the odd isogeny `h_a` has `a^2` distinct
preimages of `R`. The section `Q -> s_a(Q,R)`, of degree at most `3D`,
vanishes at all of them. Since `a^2>3D`, it vanishes identically.
This holds for every `R`, hence `s_a=0` on the whole two-factor product.

Finally fix arbitrary second and third points. The remaining section in
the first factor, of degree at most `3D`, vanishes at the `3D+1` distinct
points `iota(a)` above. It is identically zero. Thus `s=0` everywhere,
contradicting the hypothesis and proving (4).

A conjunction of equations defining the multiplication graph on the encoded
grid must include one whose restriction to `E^3` is nonzero. Otherwise
every triple of encoded words satisfies all equations. Consequently at
least one separate degree is `Omega(2^(w/3))`. Curve equations themselves
are identically zero on `E^3` and do not select the multiplication graph.
This claim concerns direct equations in the three curve words, without
auxiliary existential variables or disjunctive graph clauses.

**Affine-chart variant.** Suppose ordinary affine coordinate polynomials,
of degree at most `D` in each curve word, vanish on every product triple
whose three points avoid infinity. A polynomial whose restriction to the
three affine curves is nonzero must satisfy

```
M <= 366 D^3 + 2.                                        (5)
```

Homogenize in each word. Choose `3D+2` odd slopes from `2D+1` through
`8D+3`, all at most `11D`. At most one has first point at infinity.
For every remaining slope, at most two of the `M` values of `b` have an
input or output at infinity: multiplication by odd `a` permutes residues.
If `M>366D^3+2`, the at least `M-2` zeros still exceed the pullback
degree `3D(1+a^2)<=366D^3`. The previous argument gives at least `3D+1`
entire zero slices and the same contradiction. Thus merely separating
the exceptional infinity cases does not remove the exponential degree cost.
The auxiliary/disjunctive restrictions remain as stated above.

## 5. The remaining arithmetic interface is different

These results stop the direct functional/chart and direct-equation proposals
under their stated assumptions. A uniform factoring use still needs the
complete word, carry, factor-prefix and evaluation costs. The small degree
of the successor map supplies none of those missing costs by itself.

For auxiliary equations, algebraic elimination may produce a very large
degree; (4) then constrains that degree rather than the written gate degree.
More significantly, a cover can project densely onto `E^3` over the
algebraic closure while its **finite-field rational fibres** select a small
set of triples. Neither (3) nor (4) excludes this arithmetic selection.
An Artin--Schreier-type fibre predicate is an example of the distinction,
not a supplied encoding of multiplication.

Further work on this representation must provide such an explicit alternative
and its complete variable/degree/precision invoice. No larger curve search,
more chart fitting or direct low-degree equation sweep is warranted.

## 6. Reproduction and credit

Run `make verify-elliptic-word-costs`. The
[checker](../tools/check_elliptic_word_costs.py) exhausts eight cyclic
product grids and all residues, checks 30,368 affine cases with slopes
outside the encoded subgroup, and audits six integer slice constants.
The [receipt](../evidence/elliptic-word-cost-controls.json) records these
finite controls and the source/version pointers. It does not verify the
geometric theorems by finite extrapolation; those use the ordinary arguments
and standard results stated above. Regenerate with:

```
python3 tools/check_elliptic_word_costs.py --generate
```

This is an additional interface audit alongside the
[native arithmetic framework](arithmetic-framework.md) and
[Sanftenberg source intake](sanftenberg-semiprime-intake.md). It does not
invalidate ordinary bit arithmetic in the packed coefficient ring.
Standard geometric results retain their attribution to Milne and their
original sources. The finite-grid bounds and cost separation were developed
with substantial OpenAI assistant support; no independent mathematical
review or new Lean formalization is claimed.
