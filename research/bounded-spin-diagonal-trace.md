# Diagonal trace forms for the bounded-spin primary decomposition

## Result and scope

For a commutative coefficient algebra acting on a module, an idempotent
that is a projection has rank equal to its trace. In the normalized
path-pair basis used by the bounded-spin cut operators, the diagonal
entry at the state `(a,b)` is the original-pencil evaluation on the
four paths `(a,b,b,a)`. This note records those **trace forms**, checks
the published ten-site primary dimensions modulo five, and supplies the
twelve-site forms. It does not prove a twelve-site rank bound and does
not change the complete score.

## 1. Trace identity

Let `e` be an idempotent matrix over `F_5` with `e^2=e`. Then
`trace(e)=rank(e)` over `F_5`. In the cut matrix from
[the six-vertex constructor](../tools/bounded_spin_six_vertex.py), the
row state is `(x,y)` and the column state is `(v,u)`; the diagonal
entry is the entry at `(x,y,y,x)`. After the column norm division in
`operators`, this is the normalized original-pencil polynomial at
`(a,b,b,a)`.

If `e=sum_j e_j B_j` in the original coefficient basis, then

```math
\operatorname{rank}(e)=operatorname{trace}(e)
=\sum_{a,b}\sum_j e_j\,[z^{d_j}]P_{a,b,b,a}(z).      \tag{1}
```

The checker evaluates (1) at length ten using the seven published
idempotent coefficient vectors and obtains the primary dimensions
modulo five. This is a checksum on the published rank model; it does
not replace the exact rank computation.

## 2. Ten-site control

The trace values, reduced modulo five, are

```text
4, 1, 0, 0, 0, 1, 0.
```

They match the seven primary dimensions

```text
34, 1, 225, 300, 35, 406, 155
```

modulo five. The exact dimensions come from
[the ten-site algebra model](../evidence/bounded-spin-ten-algebra-model.json);
the trace calculation is an independent diagonal check.

## 3. Twelve-site trace forms

At length twelve there are `89^2=7921` path pairs and 35 distinct
normalized trace forms. Their multiplicities begin

```text
800, 624, 616, 544, 544, 476, 384, 368, 320, 304, 304, ...
```

and have class-size sum `7921`. The identity trace is one modulo five
for every pair, as it must be at the crossing endpoint. The
square-root diagnostic

```math
\sum_{C} \sqrt{|C|}=467.988540814\ldots
```

is below the twelve-site label target `89^sqrt(2)=571.286269...`.
**This is not a rank lower bound.** The same diagnostic at ten sites
is `137.673...`, below the `146.500...` target, while the actual
central idempotent rank sum is much larger. The number is recorded only
as a diagnostic for the next primary computation.

The ordinary diagonal weights at `(a,b,a,b)` in
[the published table](bounded-spin-diagonal-weights.md) are different:
they have 124 classes and square-root sum `808.241...`. The trace forms
here use the swapped pair `(a,b,b,a)` because they are the diagonal
entries of the cut matrices, not the endpoint evaluations.

## 4. What is reduced and what remains

The twelve-site primary dimension problem now separates into two
objects:

1. the twelve-site coefficient idempotents in the original algebra;
2. the trace forms already computed here.

The forms can be evaluated in about two minutes with the standard
library constructor. They give a modulo-five checksum once the
idempotents are known. They do not determine the integer ranks, and a
rational lift or an exact rank computation is still required before any
comparison with the `571.286...` target.

The next concrete target is therefore the twelve-site original algebra
model, not a larger fixed-block sweep. No factoring exponent changes;
the complete bound remains `1/2`.

## Reproduction

```sh
python3 tools/check_bounded_spin_diagonal_trace.py \
  --output evidence/bounded-spin-diagonal-trace.json
make verify-bounded-spin-diagonal-trace
```

The checker uses the published path and monomial-entry constructors and
does not build a cut operator or a dense tensor.
