# Exact finite cutoff 25 for the joint antipodal-sector span

## Result and scope

At binary label width two, the complete joint coefficient symmetry span
of the one-pair and two-pair constructions has dimension 33 over `F_2`.
Its optimal maximum row-rank times column-rank cutoff is **exactly 25**.

The [six-piece witness](../evidence/hamming-joint.json) gives the upper
bound. A separating functional and a retained clause refutation give the
lower bound. The refutation is replayed by a small forward unit-propagation
checker; checking requires no SAT solver, Python package, or Lean installation.
A second, independently obtained RUP checker has also accepted the same
CNF and proof.

This is a finite theorem about the declared 33-dimensional `F_2` span.
It does not apply to arbitrary tensors, other label selections, or extension
fields. The maximum rounded binary channel-width sum remains six. Thus it
is neither a new width-rate theorem nor an exponent above one half.
The construction, encoding, and proof-checking work used substantial OpenAI
assistant support.

## 1. The exact source space

Use the binary width-two label selection in
[hamming-sectors.md](hamming-sectors.md), with all bit weights one:

- One pair per label: free width two and all four pair labels.
- Two pairs per label: free width two and the first four two-subsets of
  the four basic pairs, in lexicographic order.

Take every polynomial coefficient of each family and all eight
crossing-preserving coordinate permutations. Their joint span `S` is
33-dimensional. The source checker constructs a deterministic echelon
basis `G_0,...,G_32` from these actual tensors; the first family spans
13 dimensions. No candidate truncation or heuristic low-rank list enters
the definition of `S`.

The upper-bound tensors have rank pairs

```text
(16,1), (1,16), (4,4), (4,4), (3,4), (5,5).
```

They sum to the crossing. Their complete source expressions use only
14 coefficient images in total, at most four in one tensor. Independent
dense elimination checks both ranks of every tensor.

## 2. A three-entry separating functional

For a tensor on four labels `0,1,2,3`, define

```math
ell(F)=F(0,1,1,0)+F(0,1,2,0)+F(1,0,0,2).
\tag{1}
```

It has `ell(E)=1` on the crossing. On the declared source space it is
invariant under every crossing-preserving coordinate symmetry. This is a
finite linear identity, checked exactly on all 33 generators for each
of the eight permutations. In canonical source coordinates it is the
parity functional with mask `0x686610c0`. The three-entry expression is
also checked against that mask on every generator.

The retained certificate proves the following finite statement:

> No tensor `F` in `S` with `ell(F)=1` has row rank at most four.

Applying the row/column exchange, and using the verified invariance of
`ell` on `S`, proves the same assertion for column rank. Therefore every
`F in S` with `ell(F)=1` has both ranks at least five.

In any `F_2` decomposition `sum_i F_i=E`, linearity gives
`sum_i ell(F_i)=1`. Some summand has `ell(F_i)=1`, hence rank product at
least `5*5=25`. The displayed witness attains 25, proving optimality.

## 3. Exact Boolean encoding of the row-rank statement

There are 33 Boolean source coefficients `c_j`, and the candidate tensor
is `F=sum_j c_j G_j`. Its row arrangement is a 16-by-16 binary matrix.
The formula imposes `ell(F)=1` and a factorization

```math
F=AB,\qquad A\in F_2^{16\times4},\quad B\in F_2^{4\times16}.
\tag{2}
```

Each product is an AND gate and each sum is an XOR gate. Tseitin clauses
encode these gates with their exact two-input truth tables. Equality
clauses tie every entry of (2) to the corresponding source linear form.
Thus every satisfying assignment gives a tensor of row rank at most four.

Conversely, every matrix of row rank at most four has such a factorization.
For a finite symmetry reduction, require `B` to be the reduced row-echelon
basis of its row space, padded by zero rows. One-hot pivot variables select
a pivot in columns 0 through 15, or sentinel 16 for an inactive row.
Clauses impose strictly increasing active pivots, inactive rows last,
zeros before a row's pivot, and zero entries in every other row at a pivot
column. Inactive rows are zero. Any row space of dimension at most four
has exactly this type of basis, and its original rows give the required
coordinates in `A`. Hence the additional constraints remove no relevant
matrix.

The resulting deterministic formula has 2,357 Boolean variables and 9,802
clauses. Source generation and the gate/row-space argument above prove the
encoding's meaning; an UNSAT answer alone would not establish that meaning.
As implementation controls, exact numeric row reduction supplies satisfying
assignments for five tensors with `ell=0`, with ranks at most four,
and for a rank-five tensor with `ell=1` in the rank-five version of the same
encoding. The default checker verifies every clause on those assignments.

## 4. A solver-independent clause certificate

Glucose 4, through `python-sat 1.9.dev17`, produced an UNSAT proof.
Deletion lines were omitted: retaining a previously justified clause can
only strengthen unit propagation. The stored
[compressed RUP proof](../evidence/hamming-joint-rank4.rup.gz) has 8,532
addition steps, ending in the empty clause.

For each proposed clause `C`, the checker assumes the negation of every
literal in `C` and repeatedly scans all earlier clauses for forced unit
literals. It accepts `C` only if this propagation reaches a contradiction.
Every accepted clause is consequently a logical consequence of the original
CNF and the previously accepted consequences. Induction proves that an
accepted final empty clause certifies unsatisfiability of the original
formula. This argument justifies the checking rule; it does not trust the
search solver's status or heuristics.

The checking program is [check_rup.c](../tools/check_rup.c), an independent
forward scanner written for this package. It accepts only addition steps,
checks literal bounds and complete clauses, and requires a checked empty
clause. Negative controls reject an unproved empty clause, truncated proof
input, an out-of-range literal, and an empty-clause claim for a satisfiable
small formula. A small two-step refutation supplies a positive control.
The default replay then checks every step of the retained proof.

The method and public comparison implementation are described by
[Heule, Hunt, and Wetzler's RUP/DRUP page](https://www.cs.utexas.edu/~marijn/drup/),
which cites *Trimming while checking clausal proofs*, FMCAD 2013.
The unmodified public `rup-forward.c` also reports `s VERIFIED` on the exact
same regenerated CNF and addition-only proof. Its source hash, commands,
formula/proof hashes, and result are recorded in the
[independent audit receipt](../evidence/hamming-joint-rup-independent.json).
That historical audit is distinct from the repo's default forward replay.

## 5. Reproduction and limits

```sh
python3 tools/check_hamming_joint.py
python3 tools/check_hamming_joint_lower.py
make verify
```

Default replay regenerates every source tensor and clause, verifies the
separating functional and its symmetries, checks the numeric encoding
controls, compiles the C checker with fortified warnings, and checks the
retained proof. Python's standard library and a C11 compiler suffice.
The evidence records hashes of the regenerated CNF, uncompressed proof,
checker sources, and compressed artifact.

Regenerating the proof is optional and requires an external search solver:

```sh
python3 -m venv /tmp/factoring-proof-venv
/tmp/factoring-proof-venv/bin/pip install python-sat==1.9.dev17
/tmp/factoring-proof-venv/bin/python tools/check_hamming_joint_lower.py --regenerate-proof
```

The solver is needed only for that regeneration command. A different
solver version may produce a different valid trace, which the independent
checking rule must still accept. It is not a dependency of `make verify`.

All ranks and coefficients in this lower-bound certificate are over `F_2`.
Fixed extension-field optimization remains a separate direction. The
all-size joint-family construction remains effective, but a uniform rank
saving is open. This exact finite optimum provides a reliable test case
for that further analysis.
