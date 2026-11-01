# Graph rewrites of path B: relation cycles and multiplier treewidth

Base: local checkout of `1337isnot1337/integer-factoring` at `7cb4576`. This
follows the path-B note in PR #68. All files are under
`tools/multiplier_treewidth/`. AI-assisted; see disclosure below.

Labels: **PROVEN** means proved here, with the full argument given.
**VERIFIED** means checked by a script in this directory, with its output
saved. **KNOWN** means cited. **CONJECTURE** means unproved.

## Summary

| Angle | Verdict |
|---|---|
| 1. Relations as cycles in a factor-base graph | The cycle *guarantee* is trivial. Graph theory gives existence, including short ones, but not construction. The construction step reduces back to the existing inputs: connectivity needs B1b, the least-nonresidue / generation bound, and finding *useful* cycles needs U_beta, the smooth-value localization. No door opens. |
| 2. Factoring as SAT on a multiplier network | **Ruled out.** Sub-`ell/sqrt(log ell)` treewidth is impossible for any network in which N enters only by pinning product bits and x and y appear as bit variables. The treewidth is in fact `Theta(ell)` (PROVEN, Theorem T1 plus the row-by-row upper bound). The only escape is to precompute a global, "spread" code of N into the tables. For that model, communication arguments provably cannot give lower bounds, and a polynomial-table lower bound would separate factoring from P. |

Neither angle changes the repository score of `1/2`.

---

## Angle 2: treewidth of multiplication networks

### Setting

A *network* F is a set of local constraints, such as CNF clauses or tables,
over variables `x_0..x_{n-1}`, `y_0..y_{n-1}`, `z_0..z_{2n-1}` and arbitrary
auxiliary variables. F *encodes multiplication on a family Y* if, for every
`y` in Y and all x and z, `∃aux F(x,y,z,aux)` holds iff `x*y = z`. The primal
graph `G_F` joins two variables that share a constraint. There is no
restriction on the size of F or on its auxiliary variables.

### Theorem T1 (PROVEN)

Let `n' = 6*floor(n/6)`. If F encodes multiplication on
`Y = {2^j : 0 <= j < n'}`, then

```
tw(G_F) >= n'/120 - 1.
```

The same bound holds, minus 2, if Y is only
`{2^j + 1 : n'/2 <= j < n}` and x is restricted to odd values
`< 2^(n'/2)`. This covers encodings that build in "x, y odd" and
"1 < x <= y". (For odd-only encodings, read `n'/2` in place of `n'` in the
proof below.)

**Proof.**

*Step 1: rectangles from a separator.*
- Take the vertex set W of the 2n' product bits `z_0..z_{2n'-1}`. Fix the
  higher x and y bits to 0, which is a restriction and does not change the
  graph.
- Let k = tw(G_F). Standard balanced-separator lemma (Robertson–Seymour; see
  Cygan et al., *Parameterized Algorithms*, Ch. 7): there is a separation
  (A, B) with `|A∩B| <= k+1` and `|W ∩ (A\B)|, |W ∩ (B\A)| <= 2|W|/3`. Put
  `S = A∩B`.
- Every constraint lies inside A or inside B. Assign the constraints on S to
  A. Then, for each assignment σ to S,

  ```
  f = OR_σ [∃aux_A F_A(σ,·)] AND [∃aux_B F_B(σ,·)].
  ```

- Each disjunct is a combinatorial rectangle of 1-inputs with respect to the
  bipartition (visible variables in A, visible variables in B\A). Visible
  variables in S go to the A side, together with the condition that they
  agree with σ. So f has a 1-rectangle cover of size `2^(k+1)`.
- This step is the known communication-complexity route to width lower
  bounds for encodings with auxiliary variables (KNOWN: Briquel–Koiran–Meer
  2011; Mengel–Wallon, JAIR 2020, arXiv 1905.05290, Thm 9).

*Step 2: both sides hold many product bits.*
- If `k+1 > n'/30`, the theorem holds.
- Otherwise each side owns at least `2n'/3 - (k+1) >= 19n'/30` of the z bits.
  Colour each variable A or B by its side.

*Step 3: the shift lemma (PROVEN).*
- Define `Cross(j) = #{i < n' : colour(x_i) != colour(z_{i+j})}`. Then
  `Σ_j Cross(j) = Σ_i O_i`, where `O_i` counts the z bits of the opposite
  colour in the window `[i, i+n')`.
- Split the indices i into six blocks of n'/6. All windows in block t contain
  the core `C_t`, which has size `5n'/6 + 1`.
- Suppose every core had minority count `< n'/20`. Adjacent cores share
  `>= 4n'/6` elements, so all cores have the same majority colour. The union
  of the cores misses only `n'/3 - 1` z bits. So the other colour has
  `< n'/3 + 6n'/20 = 19n'/30` z bits, which contradicts Step 2.
- So some block contributes at least `(n'/6)(n'/20)`. Hence
  `Σ_j Cross(j) >= n'^2/120`, and some j has `m := Cross(j) >= n'/120`.

*Step 4: fooling set.*
- Restrict to `y = 2^j`. Then `z = x << j` with no carries (VERIFIED, L4).
- Fix the x bits outside the crossing pairs to 0, and the z bits to the
  forced values. What remains of f is `AND_i [x_i = z_{i+j}]` over m pairs,
  and each pair is split between the two sides.
- Its 2^m solutions form a fooling set: mixing two of them breaks some split
  pair. So `2^(k+1) >= 2^m`, which gives `k+1 >= n'/120`.

*Odd variant.*
- Use x bits below `n'/2`, the shifts `j in [n'/2, n')`, and W equal to the
  middle n' product bits. Then `x*(2^j+1) = x + (x << j)` has no carries
  (VERIFIED, L4).
- Each crossing gives a triple `x_i = z_i = z_{i+j}` that is not all one
  colour, and fooling sets work for triples too. Dropping `x_0 = 1` and the
  excluded solution x = 1 costs at most 2. ∎

**Provenance.** The shift / information-transfer idea is KNOWN. It is
Brent–Kung, "The area-time complexity of binary multiplication", JACM 28
(1981), `AT^2 = Omega(n^2)` via `b = 2^j`. Related work: Abelson–Andreae
(CACM 1980); Ja'Ja'–Prasanna Kumar (information transfer, JACM 1984);
Paterson–Fischer–Meyer (on-line multiplication overlap argument, 1974). The
transfer of rectangle covers to treewidth with auxiliary variables is also
KNOWN (above).

I did not find T1 for multiplication stated in print. It is a direct
combination of known pieces, so treat it as **probably folklore**, not as a
new result. The constant 1/120 is weak. The exact minimum of
`Σ_i min(#A, #B)` is `n^2/9` at n = 6 and n = 12 (VERIFIED, L1), against the
proved `n^2/120`.

### Corollary T1.1: hard-wired N (PROVEN)

Suppose a factoring network G0 for a given N is obtained from such an F by
pinning the product bits to N. There are two cases.

- **Pinned bits stay as vertices** (unary constraints). Then `G0 = G_F` and
  T1 applies directly.
- **Pinned bits are propagated away** (`G0 = G_F - Z`). Suppose each z_k is
  simplicial in G_F, for example because it occurs in one constraint of arity
  at most w. Then `tw(G_F) <= max(tw(G0), w-1)`: hang a leaf bag holding z_k's
  scope off a bag that contains its neighbour clique.

  So `tw(G0) >= n'/120 - 1` whenever `w - 1` is below that. This covers the
  standard Tseitin "multiplier + output unit clauses" SAT encoding (w = 3),
  and any bit-level circuit that ends by comparing recovered product bits
  with N.

With `ell = 2n`, the bound is `tw >= ell/240 - O(1)`. That is
`omega(ell/sqrt(log ell))`, so the repository's open question
(`sanftenberg-transform-treewidth-floor.md`, last section;
`sanftenberg-root-four-bit-circuits.md`, "global separator") is answered
**negatively** for every whole bit circuit of that shape. This holds whatever
the transforms, twiddles, pointwise products and carry recovery look like
inside.

Scope caveat: I did not inspect the upstream Lean sources. The claim applies
if x and y occur as binary variables and N enters only as pinned product
bits, which is how the repository notes describe the circuit ("canonical-word
checks", "coefficient/carry recovery"). If x or y exist only in transform
coordinates, T1 needs a re-check (CONJECTURE: the same bound holds whenever
the binary x is recoverable by local constraints).

### Upper bound and growth (PROVEN / VERIFIED)

- **PROVEN:** an array (row-by-row ripple) multiplier has pathwidth
  `<= 2n + O(1)`. Process the rows in order and, inside a row, the ripple
  positions in order. A bag holds all of x, y_j, the live accumulator window
  of n+1 bits, and O(1) adder wires. Each pinned z_k hangs off as a leaf.
  Together with T1, the treewidth is `Theta(n) = Theta(ell)`.
  (`sanftenberg-root-four-bit-circuits.md` already notes an `O(T)` path
  decomposition per product gadget.)
- **VERIFIED** (`tw_multipliers.py`, output in `output_tw_multipliers.txt`):
  certified brackets `MMD+ lower bound <= tw <= min-fill upper bound` on the
  Tseitin primal graphs, n = 4..24 (min-degree only at n = 32).

| Circuit | Fit of the upper bound (n ≤ 24) | Fit of the lower bound | log-log slope (UB), n = 20 → 24 |
|---|---|---|---|
| array | `3.05 n − 0.5` (R² = 0.996) | `1.08 n + 2.5` | 0.91 |
| Wallace | `3.66 n − 1.4` (R² = 0.987) | `1.33 n + 1.3` | 1.10 |

  Both bounds grow linearly. The gap between them is heuristic slack. Wallace
  trees come out somewhat *worse* than arrays: their column compression mixes
  far-apart partial products. Removing Z leaves the upper bound unchanged.
  Small n cannot distinguish `n` from `n/sqrt(log n)`, so the asymptotic
  verdict rests on the proof, not on the fit.

### What is still open (the only door)

**PROVEN (trivial):** suppose the tables may depend on N arbitrarily, with
polynomial-time construction. If factoring were in P, unit constraints
`x_i = bit_i(p)` would give width 0. So no superconstant treewidth lower
bound for that model can be proved without proving factoring ∉ P.

Rectangle methods also cannot reach it, because `[xy = N]` has at most d(N)
solutions and therefore a cover of size d(N). So any encoding that beats T1
must read N in a **spread** way: every balanced separator must see
information about N that is precomputed globally.

**CONJECTURE:** the natural candidate, a residue-number-system encoding
(local `r_i * s_i ≡ N mod m_i`, plus the global range checks `1 < x, y < 2^n`),
moves the whole difficulty into the range check.

- Supporting fact (PROVEN, elementary): with balanced moduli
  `M_A, M_B >= 2^n`, the predicate `[CRT(r) < 2^n]` alone has 1-rectangles of
  size at most `(2^n/M_A + 1)(2^n/M_B + 1) <= 4`, hence cover size
  `>= 2^(n-2)`.
- This does not bound the conjunction, because of the few-solutions
  objection above. So the RNS route is **open, and I would rate it
  implausible**. RNS magnitude comparison is the classical hard RNS
  operation.

## Angle 1: relations as cycles

### Class-group Cayley graph

**Construction.** Vertices are the classes of `Cl(D)`, with `D = -4N` or a
fundamental D, and `h = |Cl(D)|`. Edges are `f -- f·[p]` for split primes
`p <= B`, so the edges are explicitly constructible by composition. A closed
walk with exponent vector `e ∈ Z^r` is a relation `prod [p_i]^{e_i} = 1`.

**Facts:**

1. **KNOWN (Schreier's lemma).** If the graph is connected, the fundamental
   cycles of any BFS tree generate the full relation lattice
   `Λ = ker(Z^r → Cl(D))`, which has rank r and determinant h.
   VERIFIED for 23 discriminants (`check_cayley_relations.py`): the HNF
   determinant equals h in every case.
2. **PROVEN (trivial).** The cycle space has rank `(r−1)h + 1`, but almost
   all of it is commutation cycles `g h g^-1 h^-1` with exponent vector 0. So
   the girth is ≤ 4 as soon as r ≥ 2, and girth/density (Moore-bound)
   reasoning only produces useless cycles. The right object is the lattice Λ,
   not the cycle space. Example: D = −8084 has cycle-space rank 409 but
   lattice rank 7.
3. **KNOWN (Minkowski / pigeonhole, the abelian analogue of the Moore
   bound).** A nonzero relation with `|e|_1 <= 2L` exists once
   `#{|e|_1 <= L} > h`. VERIFIED on the same set: the shortest relation found
   is always within the pigeonhole bound. Abelian Cayley graphs grow
   polynomially, so this gives only `|e|_∞ ≈ h^{1/r}`.
4. **Construction is the hard part (KNOWN).** For D = −4N, knowing Λ gives
   the group structure, and from it ambiguous forms and a factor
   (Schnorr–Lenstra 1984; Lenstra–Pomerance 1992; Hafner–McCurley 1989 for
   the subexponential class-group algorithm under GRH). With generic group
   access, relation and discrete-log problems need about `Omega(sqrt(p))`
   operations, where p is the largest prime divisor of h (KNOWN for DL:
   Shoup 1997; for order computation: Sutherland's 2007 thesis). I have not
   checked that the exact "nonzero cycle" form is in print, so treat that
   form as CONJECTURE. Graph structure alone therefore gives nothing
   subexponential.
5. **The non-generic input is smoothness.** A reduced form `(a,b,c)` has
   `a <= sqrt(|D|/3)`. If a is B-smooth, its class is a short path. This is
   the Hafner–McCurley relation step. Producing such forms deterministically
   is exactly U_beta from the path-B note.
6. **Connectivity is B1b.** The graph is connected iff the split primes
   ≤ B generate `Cl(D)`. Under GRH, `B = O(log^2 |D|)` suffices (KNOWN:
   Bach 1990). Unconditionally the only bounds are Burgess-type,
   least-nonresidue-type bounds. That is the B1b input again. VERIFIED: the
   smallest connecting B is small in all 23 cases (column `Bconn`), but the
   cases are tiny.

### Residues mod N and large-prime graphs

**Residues mod N.** Use `Cay((Z/N)^*, {p <= B})`. A cycle gives
`prod p^e ≡ 1 (mod N)`, i.e. order information. That factors N only through
Miller-type splitting, which is randomized or ERH-based (KNOWN: Miller 1976).
Cycles exist trivially, finding them is order finding, and connectivity is
again a generation / least-nonresidue bound (Bach, under ERH).

**Large-prime-variation graphs.** Vertices are large primes and edges are
partial relations. A cycle is a full relation, and only a fundamental cycle
basis matters (KNOWN: Lenstra–Manasse, *Math. Comp.* 63 (1994)). The
guarantee is graph-theoretic and trivial: once `|E| > |V| − c`, a cycle
exists. But edges come from smooth-up-to-cofactor values, so deterministic
edge production is U_beta with Y replaced by a larger cofactor bound. The
large-prime variation changes only the `o(1)` in the L-exponent (KNOWN:
Lenstra–Manasse 1994, §1).

**Verdict:** every graph-theoretic guarantee here is either trivial
(existence) or reduces to B1b (connectivity) or U_beta (edge production).
There is no new door.

## Literature searched

- Treewidth/CNF encodings of multiplication.
- Mengel–Wallon (JAIR 2020).
- Briquel–Koiran–Meer (2011).
- Bova–Capelli–Mengel–Slivovsky (IJCAI 2016, DNNF ↔ rectangle covers).
- Pipatsrisawat–Darwiche (2010).
- Bryant (1991, OBDD for the middle bit, `2^{n/8}`).
- Woelfel; Bollig–Woelfel; Sauerhoff–Woelfel (branching programs for
  multiplication).
- Beame et al. ("Towards verifying nonlinear integer arithmetic").
- Brent–Kung (1981); Abelson–Andreae (1980); Ja'Ja'–Prasanna Kumar (1984);
  Paterson–Fischer–Meyer (1974).
- Lenstra–Manasse (1994); Hafner–McCurley (1989); Shoup (1997); Bach (1990).
- Repo: `research/` grep for "treewidth" and INDEX line 123.

I did not find T1 itself in print. Given how close it is to Brent–Kung plus
Briquel–Koiran–Meer, I still assume someone has written it.

## Rerun

```
cd tools/multiplier_treewidth
python3 check_shift_lemma.py      > output_check_shift_lemma.txt       # ~20 s, PASS
python3 check_cayley_relations.py > output_check_cayley_relations.txt  # <1 s, PASS
python3 tw_multipliers.py         > output_tw_multipliers.txt          # ~2.5 min
```

These use pure Python 3 with no dependencies. `check_cayley_relations.py`
imports `../smooth_relations_scope/bqf.py`.

This work used substantial Grok Bot assistance.
