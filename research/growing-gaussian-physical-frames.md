# Growing Gaussian changes of physical coordinates

## Status and scope

This note gives a constructive family and a complete profile interface. It
does not improve the achieved factoring bound. The number of changed physical
modes can grow with the interval parameter. The argument applies to arbitrary
independent mode heights, including heights outside the old regularity cone.

The exact coefficient implementation is
[`reflection_gaussian_physical_frames.py`](../tools/reflection_gaussian_physical_frames.py).
It uses the unrestricted Clifford theorem in
[`clifford-gaussian-smith-profiles.md`](clifford-gaussian-smith-profiles.md).
The original Fock tensor remains the object whose two Smith profiles are used.

## 1. The family

Let M be a power of two, M at least 4, and put m=M/2. Each of the four
physical letters x,y,u,v has m binary occupation modes. On one physical letter,
choose a set E of unordered pairs of distinct modes. A creation block is

\[
 U_E(g)=\prod_{\{i,j\}\in E}(1+g C_iC_j),
 \qquad g(w)=\gamma(w^a+w^b),\quad a\ne b,
 \quad\gamma\in\mathbb F_4^\times.                 \tag{1}
\]

An annihilation block replaces both C operators by A operators. All factors
within one block commute. A word may contain blocks on any physical letters,
with creation and annihilation blocks in any order. Examples supported by the
implementation are adjacent pairs, reversal pairs, and every pair.

Every factor in (1) has inverse itself: its nonconstant part squares to zero.
It has determinant one on the physical Fock space. Since g(1)=0, the word is
the identity at w=1. Applying the word to the reflection tensor therefore
preserves the strict crossing endpoint

\[
 T_{xyuv}(1)=\mathbf1_{x=v}\mathbf1_{y=u}.             \tag{2}
\]

The word can change a positive fraction of all modes. It is not restricted to
a fixed finite group of modes.

## 2. Exact action on both cuts

Write the Clifford vector basis as (C_1,...,C_M,A_1,...,A_M), and let B be the
conjugation matrix of a Fock cut F. A creation block on row modes left multiplies
B by the alternating upper shear. For one pair, its changed rows are

\[
 B_{C_i,*}\leftarrow B_{C_i,*}+gB_{A_j,*},\qquad
 B_{C_j,*}\leftarrow B_{C_j,*}+gB_{A_i,*}.             \tag{3}
\]

An annihilation block exchanges the C and A roles. The source rows are
unchanged, so (3) also handles overlapping pairs in one block.

A physical operation on a column letter acts on the cut as F U^T. Thus
creation on a column letter becomes annihilation after transposition. The
resulting right lower shear changes the creation columns:

\[
 B_{*,C_i}\leftarrow B_{*,C_i}+gB_{*,A_j},\qquad
 B_{*,C_j}\leftarrow B_{*,C_j}+gB_{*,A_i}.             \tag{4}
\]

Column annihilation changes the annihilation columns from creation columns.
Equations (3) and (4) specify the action convention completely. They apply
separately to the cuts (x,y)|(u,v) and (x,u)|(y,v).

All physical frame factors have determinant one on either cut. Consequently
the ungauged mean Smith orders remain M^2 at zero and 0 at infinity, in the
normalization of the reflection tensor. If h is the vector of all 2M physical
mode heights at one place, its additional mean is (sum h)/2. If the positive
half of the resulting Bogoliubov Smith spectrum is a_1,...,a_M, the complete
Fock spectrum is therefore

\[
 \left\{\mu+\sum_{i\in S}a_i:S\subseteq\{1,\ldots,M\}\right\},
 \qquad
 \mu=\operatorname{mean}(F)-\frac12\sum_i a_i.        \tag{5}
\]

The Clifford theorem proves (5) in characteristic two. In particular, this
calculation does not require an integral Gaussian kernel factorization after
the physical heights are applied.

## 3. All degree and pole costs

In one creation block, a nonzero term uses at most m/2 quadratic factors,
because every mode can occur at most once. The same argument applies to an
annihilation block. If block j has powers a_j<b_j, the word can lower the
minimum Laurent exponent by at most

\[
 C=\frac m2\sum_j\max(0,-a_j),                       \tag{6}
\]

and raise the maximum exponent by at most
\(\frac m2\sum_j\max(0,b_j)\).

The original reflection tensor has degree at most 3M^2/2. One way to see this
bound is to use its natural infinity heights: every kernel edge exponent is
at most the sum of the two mode heights, their total is M^2, and the common
scalar has exponent M^2/2. After multiplying the new tensor by w^C, its degree
is at most

\[
 \frac{3M^2}{2}+\frac m2\sum_j
       \bigl(\max(0,-a_j)+\max(0,b_j)\bigr).          \tag{7}
\]

The clearing scalar equals one at the endpoint. It adds C to every zero
Smith order and subtracts C from every infinity Smith order. The product of
the two profile sums is unchanged, so neither pole cost is omitted.

For a bounded number of blocks with powers O(M), (7) is O(M^2). The original
coefficient table and the framed table can be constructed in 2^O(M)
operations. Thus this family stays within the effective resource scale of
the existing reflection construction.

## 4. Verification and search scope

The checker constructs the original tensor coefficients independently, applies
the occupation operators directly, and verifies every endpoint coefficient.
It compares each complete Fock Smith spectrum with (5), using independently
chosen integer heights at both places. The controls include all four physical
letters, both operator kinds, nontrivial F4 coefficients, overlapping pairs,
and noncommuting words. Small F2 cases also use rational-function elimination
without a truncated series window, and compare the full conjugation matrix
against the original Fock operator.

The native matrix engine increases its series window if a promised pivot is
unresolved. A failure at the supported maximum remains an error; it is never
interpreted as a vanishing pivot or a favorable profile.

An exploratory screen at the natural heights is recorded separately. Its
finite list of tested words does not exclude other words or other heights.
The pressure calculation uses the complete Gibbs root and restores both
means in (5). It does not identify a small positive critical pressure with a
factoring improvement.

Run the exact controls with:

```sh
python tools/check_reflection_gaussian_physical_frames.py \
  --output evidence/reflection-gaussian-physical-frames.json
```
