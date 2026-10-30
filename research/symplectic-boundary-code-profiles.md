# Exact rank profiles for symplectic boundary codes

The sparse symplectic construction produces large label spaces with exact
crossing endpoints. This note records complete finite coefficient-space
tests for its boundary, radical, and oriented-radical variants. It also
isolates a six-dimensional polynomial family on 26 labels for further
tensor-power analysis.

The finite tests below do not improve the factoring exponent. They exclude
all decompositions inside each stated coefficient space over the stated
field. They do not exclude larger codes, other output projections,
extension fields beyond the one tested, or decompositions available only
after taking tensor powers.

## 1. Construction and exact endpoints

Use the universal ternary Pfaffian tensor from
[the edge-family note](pfaffian-ternary-edge-family.md). There are `P`
complementary mode pairs, and all four external labels have exterior
grade `L`. Let `J` swap the two physical modes in every pair. It is not
the operation that complements occupied and unoccupied modes.

Write
\[
 B=\operatorname{im}\bigl(\Omega:\Lambda^{L-2}F^{2P}
                                    \longrightarrow\Lambda^LF^{2P}\bigr),
 \qquad R=\operatorname{Rad}(B).
\]
The [symplectic radical theorem](noncommutative-crossing-symplectic-radicals.md)
gives exact dimensions, a sparse basis, and a coordinate projection for
these spaces. The checker also constructs the spaces independently by
binary elimination, grouped by the single-occupation syndrome.

Let `I` be the input injection and `Π` the coordinate left inverse. The
two input legs are `y,v`; the output legs are `x,u`. The three vertices
are
\[
 [p^L]G=(\Pi\Pi^{\mathsf T})_{xu}(I^{\mathsf T}I)_{yv},\quad
 [q^L]G=(\Pi JI)_{xy}(\Pi JI)_{uv},\quad
 [r^L]G=\delta_{xv}\delta_{yu}.
 \tag{1}
\]
Every entry of all three identities is checked during tensor construction.
The two complete matrix cuts use rows and columns
\[
 (x,y)\mid(v,u),\qquad (x,u)\mid(v,y).
 \tag{2}
\]

For the full radical, the `p` vertex vanishes. The space is `J`-stable,
so the `q` vertex has ranks `(1,d²)`. The oriented half-radical retains
the low mode at the least singly occupied pair and omits the syndrome
with no singles. Its coordinate projection kills `JI`, so both the
`p` and `q` vertices vanish. The `r` vertex is the crossing in every case.

The literal polynomial evaluator sums the sparse input supports before
computing either cut rank. It retains all Pfaffian contractions, including
the terms created by occupied complementary pairs. It uses the exact
parity block decomposition to reduce the matrix sizes.

## 2. Complete binary coefficient-space exclusions

For every monomial `p^a q^b r^c`, the checker computes both matrix ranks
and its coordinates in an independent coefficient basis. It then tests
all nonzero binary combinations of that basis.

Let
\[
 S=\sum_i\sqrt{r_i c_i}.
\]
A winning finite seed must satisfy `S < d^{√2}`. The native screen uses
an integer gate strictly above this threshold and lists every vector
with individual cost below that gate. If these eligible vectors do not
span the crossing, no winning decomposition in the coefficient space
exists. Every eligible vector has its two ranks independently recomputed
by the Python checker.

| Code | `P,L` | Labels `d` | Coefficient dimension | Nonzero vectors tested | Eligible vectors | Eligible span dimension |
|---|---:|---:|---:|---:|---:|---:|
| Full boundary | 4,3 | 8 | 5 | 31 | 4 | 4 |
| Full boundary | 4,4 | 27 | 10 | 1,023 | 28 | 8 |
| Full boundary | 5,4 | 44 | 12 | 4,095 | 10 | 9 |
| Radical | 4,4 | 26 | 6 | 63 | 4 | 4 |
| Radical | 5,3 | 10 | 4 | 15 | 3 | 3 |
| Radical | 6,4 | 64 | 8 | 255 | 5 | 5 |
| Radical | 5,5 | 100 | 9 | 511 | 13 | 7 |
| Oriented half-radical | 4,4 | 12 | 5 | 31 | 3 | 3 |
| Oriented half-radical | 5,5 | 50 | 8 | 255 | 7 | 5 |

In all nine cases, the crossing lies outside the eligible span. The
reports retain the exact code injections and output coordinates, all
coefficient ranks, the independent scalar coefficient polynomials, and
the native screen hashes. The screen logs contain a completion record
and the exact eligible vectors.

## 3. Complete four-element-field exclusions

The six spaces of coefficient dimension at most eight were also screened
over `F₄ = F₂[ω]/(ω²+ω+1)`. Scaling a nonzero tensor does not change either
rank, so it suffices to enumerate `(4^k−1)/3` projective coefficient
vectors. The ranks here are ranks over `F₄`, not ranks of a binary
expansion.

| Code | Labels | Coefficient dimension | Projective vectors tested | Eligible vectors | Eligible span dimension |
|---|---:|---:|---:|---:|---:|
| Radical, `P5,L3` | 10 | 4 | 85 | 3 | 3 |
| Boundary, `P4,L3` | 8 | 5 | 341 | 4 | 4 |
| Half-radical, `P4,L4` | 12 | 5 | 341 | 3 | 3 |
| Radical, `P4,L4` | 26 | 6 | 1,365 | 4 | 4 |
| Radical, `P6,L4` | 64 | 8 | 21,845 | 5 | 5 |
| Half-radical, `P5,L5` | 50 | 8 | 21,845 | 11 | 5 |

Again the crossing is outside every eligible span. A separate scalar
field implementation verifies every eligible vector and forty additional
deterministically selected projective vectors per case, as well as the
crossing. The complete exclusion uses the exhaustive native enumeration;
the scalar audit is an independent check on its rank calculations.

The four-element-field audit originally used exports made before the
experiment source was frozen. The saved audit explicitly verifies exact
equality of both coefficient matrix lists, block sizes, endpoint
coordinates, and native input bytes against the frozen exports. It
records both the historical audited hash and the current frozen basis
hash. The replay command below computes a fresh audit directly from the
frozen exports.

## 4. A small scalar family inside the 26-label example

For `P4,L4`, put
\[
 u=q+r,\qquad v=p+q,\qquad s=qr.
\]
The entire tensor is
\[
 G=u^4E+u^3vA+u^2v^2B+uv^3C+s u^2D+suvF.
 \tag{3}
\]
Here `E` is the crossing. If `B₀,…,B₅` are the selected independent
coefficients at monomials
\[
 r^4,\ pr^3,\ p^2r^2,\ p^3r,\ qr^3,\ pqr^2,
\]
then
\[
 E=B_0,\quad A=B_1,\quad B=B_2,\quad C=B_3,\quad
 D=B_1+B_4,\quad F=B_1+B_5.
 \tag{4}
\]
The scalar polynomials in (3) are recovered from an invertible entry
minor and checked against every original coefficient matrix. Thus (3)
is a tensor identity, rather than a fit to selected entries.

Every tested radical and half-radical tensor has exact common factor
`q+r` and no common factor `p+r`. The
[Ward identity](noncommutative-crossing-radical-ward.md) proves the first
factor for all sizes before the output projection is applied. Further
nesting of the same global-pivot construction adds only a scalar factor,
as proved in that note. The six-generator family in (3) remains a live
source for common-parameter tensor powers: the finite coefficient-space
exclusions above do not settle their asymptotic merged rank costs.

## 5. Reproduction and evidence

The construction and profile checker is
[`explore_symplectic_boundary_codes.py`](../tools/explore_symplectic_boundary_codes.py).
Its frozen SHA-256 is
`33995930c5d2181dc13a9ac9da661be201c2f1cfbd6be07b1411c08f3f5d8720`.

From the repository root, compile the native field screen and generate
all nine exports:

```sh
c++ -O3 -std=c++20 tools/screen_pfaffian_span_f4.cpp -o /tmp/omega-span-screen
python3 tools/explore_symplectic_boundary_codes.py --cases 4:4,5:3,6:4,5:5 --native /tmp/omega-span-screen --export-dir /tmp/omega-codes --output /tmp/radical-profiles.json
python3 tools/explore_symplectic_boundary_codes.py --cases 4:3,4:4,5:4 --full-boundary --native /tmp/omega-span-screen --export-dir /tmp/omega-codes --output /tmp/boundary-profiles.json
python3 tools/explore_symplectic_boundary_codes.py --cases 4:4,5:5 --half-orientations --native /tmp/omega-span-screen --export-dir /tmp/omega-codes --output /tmp/half-radical-profiles.json
python3 tools/check_symplectic_code_f4.py --basis-dir /tmp/omega-codes --output /tmp/omega-F4-audits.json
```

The last command independently audits the six saved complete `F₄` logs.
To repeat an exhaustive native `F₄` enumeration itself, run, for example:

```sh
/tmp/omega-span-screen < /tmp/omega-codes/omega_radical_P4_L4.txt > /tmp/omega-radical-P4-L4-F4.jsonl
python3 tools/check_pfaffian_span_f4.py --basis /tmp/omega-codes/omega_radical_P4_L4.json --screen /tmp/omega-radical-P4-L4-F4.jsonl --output /tmp/omega-radical-P4-L4-F4-audit.json
```

Durable evidence:

- [Radical profiles](../evidence/symplectic-radical-code-profiles.json).
- [Boundary profiles](../evidence/symplectic-boundary-code-profiles.json).
- [Oriented half-radical profiles](../evidence/symplectic-half-radical-code-profiles.json).
- [Four-element-field audits](../evidence/symplectic-code-F4-audits.json).
- [All fifteen native completion logs](../evidence/symplectic-code-screens).
- [Complete 26-label coefficient matrices](../evidence/symplectic-radical-P4-L4-basis.json).

All larger coefficient exports are deterministically regenerated by the
commands above. No formalization is needed for these ordinary finite
linear-algebra certificates.
