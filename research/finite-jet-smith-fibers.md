# Finite crossing jets inside regular Smith-frame fibers

## Constructive result

Finite Hasse constraints can be imposed directly inside the complete-entry
reflection fibers used for the two original singular places. This retains
the specified regular-frame bounds while introducing a new finite value or
jet. The resulting tensor is explicit over `F4`, with its complete degree
and strict endpoint checked.

The construction does not require a physical-frame word or membership in
the paired-trace affine orbit. Its exact statement is feasibility inside
the chosen regular frames. It does not parameterize every tensor with the
same Smith profiles, and a feasible finite rank target is not by itself a
pressure improvement.

Two components are available:

- [F4 regular frames](../tools/f4_regular_smith_frames.py) support a reference
  whose coefficients are not binary.
- [Finite-jet fibers](../tools/finite_jet_smith_fibers.py) add full `F4` Hasse
  equations to the existing binary reference fibers without changing their
  construction or coefficient conventions.

## 1. Exact finite-precision regular frames

For a square polynomial matrix `M(t)`, choose a precision `P` above every
requested row target. The frame algorithm chooses a nonzero entry of least
order, moves it to the next diagonal position, and writes it as `t^a u(t)`
with `u(0)!=0`.

An entry `t^a q(t)` below that pivot is eliminated by the row operation

```math
R_i\longmapsto uR_i-qR_k.
```

The corresponding column operation eliminates the rest of the pivot row.
All computations are exact modulo `t^P`. In characteristic two, subtraction
is addition. Multiplication of the two coefficient bitplanes uses
`omega^2=omega+1`.

Every operation has an invertible constant determinant: its diagonal
multiplier is the unit `u`. Thus the returned polynomial representatives
`U,V` are invertible over `F4[[t]]`. The implementation replays

```math
U M V\equiv D\pmod{t^P}
```

using the full, untruncated source. It returns only when all diagonal
channels are exposed and their orders are below `P`. At that precision,
the remaining error cannot alter any Smith order. Insufficient precision
raises an error.

For a new tensor cut `N`, impose that row `i` of `U N V` is divisible by
`t^(tau_i)`, with `tau_i<=P`. These are finite linear coefficient
equations. They imply

```math
\Gamma_k(N)\ge\sum_{i=1}^k\tau_{(i)},
```

where the targets on the right are sorted and `Gamma_k` is the `k`th
determinantal order. Indeed, every minor of the transformed matrix has the
corresponding row factors, while the two frames preserve its Smith
invariants. The implication is a lower bound, not an equality claim.

## 2. Complete endpoint-preserving coordinates

Let `T` be a binary reference crossing of degree at most `d`. Every
`F4` tensor `S` of degree at most `d` with the same strict endpoint has
the unique form

```math
S_e(w)=T_e(w)+(w+1)f_e(w),\qquad\deg f_e<d                 \tag{1}
```

at each of its `q^4` entries. There are exactly `q^4*d` free field
coefficients before the regular-frame constraints.

The existing binary frame elimination supplies affine pivot equations
for these coefficients. Since their matrix coefficients are binary,
the same pivots parameterize their full `F4` solution space. Its two
binary coefficient planes use the same homogeneous reduction; the
particular affine constants occur in the first plane.

For the four-letter reflection source, natural heights `(0,6,2,8)` and
a six-order relaxation of the specified row floors leave 394 `F4`
coordinates. Fixing all original column-infinity row targets inside
that relaxed fiber leaves 92 coordinates. These dimensions refer to
the declared frames and floors.

## 3. Reduce the finite Hasse equations into that fiber

At a field point `a`, the order-`j` Hasse functional on the coefficient
of `w^k` in `f_e` is

```math
\binom{k}{j}a^{k-j}+\binom{k+1}{j}a^{k+1-j}.              \tag{2}
```

Equation (2) follows directly from (1). Terms with negative exponents
are absent through their zero binomial coefficient. The requested
right side is the target jet minus the source jet.

Each such equation is reduced through the already fixed binary affine
pivots. When its leading coefficient is `omega` or `omega^2`, that
scalar multiplies both the pivot row and its affine constant. The result
is an exact `F4` equation on the free coordinates.

Ordinary field elimination then decides simultaneous feasibility with
the requested row targets. Increasing one row target by one adds only
the next coefficient equations for that row. The full coefficient
tensor is reconstructed by back-substitution. Every retained tensor
is replayed for its endpoint, all requested Hasse entries, and the
complete regular-frame conditions.

The method is complete for this fixed linear system: the starting
coordinates cover every endpoint-preserving tensor, the first
elimination imposes exactly the chosen regular-frame equations, and
the second imposes exactly the additional finite jets.

## 4. A determinant budget for preserving the original infinity side

The original four-letter column-cut determinant is a nonzero constant.
The sum of its natural infinity profile, after subtracting the physical
height shifts, is therefore zero.

If all original column-infinity row floors are retained, the same
determinantal lower bound holds for a new crossing. A polynomial
determinant has infinity order equal to minus its degree, which is at
most zero. The strict endpoint makes it nonzero. Its infinity order
must therefore equal zero, and the new column determinant is also
constant. Such a candidate has no finite singular point on that cut.

If the original row-zero floors are also retained, its row determinant
has zero order at least 256. Its degree is at most `16*24=384`, leaving
at most 128 total multiplicity at other finite points. A new complete
rank-one row jet of length `L` uses at least `15L` of that budget, so
`L<=8`.

These bounds have the stated full-floor hypotheses. Relaxing some
infinity or zero targets changes the available determinant budget;
the search records that relaxation explicitly.

## 5. Search and normalization

The [bounded search](../tools/search_finite_jet_smith_fibers.py) first
tests full Bell jets and then complementary rectangles. Its default
keeps the original infinity targets. A declared `infinity-window`
allows a separately charged relaxation. Each feasible target is
tightened within the fixed frame, and its complete profiles are
recomputed.

The quick mixed-place objective retains the natural physical heights
at zero and infinity, a zero height at the proposed finite place,
and one common Gibbs type parameter. It is a valid restricted search
objective, not an optimization over all physical heights and paired
types. Full mixed-place refinement can vary the two paired types
independently and charge a separate physical height at the finite
place.

Candidates may have a nonconstant common scalar factor because this
fiber does not impose the old paired-trace identities. The search
computes the full polynomial gcd, divides it exactly, restores the
strict endpoint, and stores the actual degree. The initial frame and
quick pressure records concern the original degree-cap representative;
the stored normalized seed is used for subsequent complete-oracle
refinement. The common factor and original hash preserve the link
between the two representatives.

## 6. Independent controls

```sh
python3 tools/check_f4_regular_smith_frames.py \
  --output evidence/f4-regular-smith-frame-controls.json
python3 tools/check_finite_jet_smith_fibers.py \
  --output evidence/finite-jet-smith-fiber-controls.json
```

The [frame checker](../tools/check_f4_regular_smith_frames.py) reconstructs
random matrices with known complete Smith profiles, replays both frame
matrices using independent coefficient dictionaries, checks their
constant invertibility, and compares four literal CRT cuts with a separate
fraction-free oracle. Its [record](../evidence/f4-regular-smith-frame-controls.json)
also includes insufficient-precision rejection.

The [fiber checker](../tools/check_finite_jet_smith_fibers.py) compares small
field systems with exhaustive assignment lists, reconstructs eight full
four-letter tensors, and uses integer binomial coefficients for independent
finite-jet and endpoint checks. Two quick pressure values are matched with
the complete mixed-place oracle. Its
[record](../evidence/finite-jet-smith-fiber-controls.json) records each
scope and source hash separately.
