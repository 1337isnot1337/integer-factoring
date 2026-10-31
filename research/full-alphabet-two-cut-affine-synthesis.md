# Complete-entry crossing synthesis on larger physical alphabets

## Result and scope

Every entry of a degree-bounded strict crossing can be varied independently
over F4 while imposing simultaneous row-lattice conditions at the two cost
cuts. This note extends the earlier binary implementation to an arbitrary
finite physical alphabet and to arbitrary integer heights for each actual
symbol. The coefficient equations and their compressed solution fibers are
exact. The selected local frames remain sufficient constraints; they do not
represent every tensor having a given pair of Smith profiles.

Two completed four-symbol searches each performed 20,000 target walks. They
found 45 and 1,001 full-profile classes, respectively. Neither supplied a
critical crossing certificate. The established factoring score remains 1/2.

## 1. All endpoint-preserving polynomial coefficients

Let the physical alphabet have size q and let B(w;x,y,u,v) be a binary
polynomial tensor of degree at most D with

\[
B(1;x,y,u,v)=[x=v][y=u].
\]

Every F4 polynomial tensor T of degree at most D with the same endpoint has
the unique form

\[
T=B+(1+w)Q,\qquad \deg Q<D.
\tag{1}
\]

The divisibility follows entry by entry from T(1)=B(1). Thus (1) has exactly
Dq^4 free F4 coefficient variables. No sparsity, Gaussian structure, or
occupation-preserving condition is imposed on Q.

Use the cuts (x,u)|(y,v) at infinity and (x,y)|(u,v) at zero. At either
place, allow four independent integer height functions on the physical
alphabet. Negative heights are permitted. A common scalar monomial clears
every entry before the local coefficient equations are formed; its shift is
removed again from the reported Smith orders.

Apply the exact fraction-free regular Smith reduction to the weighted source
B. It supplies polynomial matrices L and R, each invertible at the local
origin, for which LBR is diagonal. The endpoint guarantees full generic rank
q^2. For a proposed row target tau_i, impose

\[
\operatorname{coeff}_{z^k}(LTR)_{ij}=0\qquad(k<\tau_i,\ \text{all }j).
\tag{2}
\]

These are affine equations in the coefficients of Q. Since L and R are
regular units, (2) gives a valid row factorization of the original cut. In
particular the first k Smith orders have sum at least the sum of the k
smallest row targets, after the common scalar shift is removed.

The implementation constructs every equation below an explicit row cap.
Requests beyond that cap are rejected. It never treats unconstructed high
coefficients as zero. The two places use the same Q variables, including
the exact reversal and the factor (1+w) at infinity.

## 2. Exact compression of repeated feasibility queries

Choose common row floors and solve their full affine system once. Eliminate
its pivots from every later equation. The remaining coordinates parameterize
the entire floor fiber, and a later target adds only equations in these free
coordinates. Feasibility and dimension agree with solving the original full
system directly.

Because the frame source and equations have binary coefficients, writing
Q=Q_0+omega Q_1 splits the F4 system into its binary affine system for Q_0
and the corresponding homogeneous system for Q_1. Both planes are sampled,
then lifted through the common floor system with the correct affine right
hand side. Every retained witness is a complete literal tensor, not merely
a compatible list of desired orders.

[`full_alphabet_two_cut_lattices.py`](../tools/full_alphabet_two_cut_lattices.py)
contains this construction. The reference must currently be binary, although
every synthesized coefficient is free over F4.

## 3. Completed four-symbol searches

The reference is the interval-four reflection crossing, of degree 24. Each
search therefore starts with 6,144 F4 coefficient variables. The physical
height lists are zero on X,Y at infinity and U,V at zero, and half the stated
paired energy list on the other four legs. Both paired types use that list.

| Paired energies | Floor-fiber variables | Target walks | Distinct targets | Full profile classes | Best rho |
|---|---:|---:|---:|---:|---:|
| (0,12,4,16) | 394 | 20,000 | 42 | 45 | 0.7310809359736007 |
| (0,0,12,16) | 348 | 20,000 | 695 | 1,001 | 0.9130262562369413 |

Each walk lowers a selected set of row targets, then greedily raises targets
in a seeded random order subject to exact feasibility. It uses three F4
samples per new target and computes both complete Smith profiles. The
coefficient window is six orders below and eight orders above the reference
row bounds. This is a completed finite stochastic search, not an exhaustive
enumeration of its target box or coefficient fiber.

The natural-energy run retained the reference rate. In the degenerate-energy
run, the starting rate was 0.916266732101722, so the reported rate is a small
improvement within that particular family. It remains far above the critical
value 1/sqrt(2), and is also worse than the natural-energy reference.
All reported rate optimizations are numerical. The exact coefficient and
endpoint claims do not depend on those optimizations.

The two records preserve their source hashes, complete search arguments,
profile multiplicities, and forty literal witnesses each:

- [`full-alphabet-q4-natural-affine-search.json`](../evidence/full-alphabet-q4-natural-affine-search.json)
- [`full-alphabet-q4-degenerate-affine-search.json`](../evidence/full-alphabet-q4-degenerate-affine-search.json)

## 4. Independent controls and reproduction

The independent checker compares nine binary systems with the earlier solver
and 24 full-system/fiber feasibility queries at three and four symbols. It
checks 1,348 literal endpoint entries, 32 complete native or independent
profiles, and 58 block-Toeplitz ranks. Restriction from F4 to F2 duplicates
each Smith order, providing an independent polynomial calculation. Three
invalid endpoint or coefficient-window requests are rejected. It also
reconstructs all eighty retained search tensors and verifies their exact
row conditions and complete profiles.

```sh
make verify-full-alphabet-affine
OPENBLAS_NUM_THREADS=1 python3 tools/search_full_alphabet_two_cut_lattices.py \
  --energies 0,12,4,16 --trials 20000 \
  --output evidence/full-alphabet-q4-natural-affine-search.json
OPENBLAS_NUM_THREADS=1 python3 tools/search_full_alphabet_two_cut_lattices.py \
  --energies 0,0,12,16 --trials 20000 \
  --output evidence/full-alphabet-q4-degenerate-affine-search.json
```

The checker output is
[`full-alphabet-two-cut-lattice-controls.json`](../evidence/full-alphabet-two-cut-lattice-controls.json).
Independent symbol-height optimization and subsequent changes of reference
frame are separate extensions; the two records above keep their heights and
frames fixed throughout.
