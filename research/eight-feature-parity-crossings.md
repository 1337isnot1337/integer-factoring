# Eight local parity features in a constant-bond crossing

## Result and scope

This note gives an explicit order-dependent crossing family over characteristic
two. Its local tensors use all eight even-parity bit patterns, including a feature
that the earlier XOR-circulant family did not allow. A transition delay preserves
the exact endpoint and leaves only linearly many coefficient pieces.

The accompanying native engine computes both flattening ranks without assuming
XOR-circulant structure. Structural screens selected eleven models for exact
profiles through eight positions. Those profiles do not improve the rank-sum
criterion. The best rank-sum branch has a stronger explanation: an exact
multiplicative law makes its coefficients orthogonal idempotents, and both ranks
are `2^t binom(t,k)`. Its rank sum is exactly `4^t` for every length.

The new construction and exchange identity are positive algebraic results. No
better integer-factoring exponent follows from this checkpoint. All screen ranks
are over `F_2`; their extension-field ranks are unchanged. The all-length proofs
below apply over every characteristic-two field.

## 1. A local matrix product with an exact crossing endpoint

For one bit quadruple `(x,y,u,v)`, the local matrix is zero unless
`x+y+u+v=0`. On even parity, define

```math
g=x+y,\qquad h=x+u,\qquad q=4x+2g+h.
```

Choose eight matrices `B_q` in `M_2(F_2)` and an integer delay `a>=0`.
Write `D(B)` for the diagonal part and `O(B)` for the off-diagonal part. Put

```math
K_q(z)=[g=h]I+(1+z)\bigl(D(B_q)+z^aO(B_q)\bigr).
\tag{1}
```

For words of length `t`, multiply the local matrices in position order and apply
a fixed linear functional `ell` with `ell(I)=1`:

```math
F_t(z;x,y,u,v)=\ell\!\left(\prod_{i=1}^t K_{q_i}(z)\right).
\tag{2}
```

The supplied screens use `ell_0(M)=M_00` and
`ell_1(M)=M_01+M_10+M_11`. They are two choices satisfying the endpoint condition;
no exhaustiveness claim about delayed boundary forms is needed.

At `z=1`, every perturbation in (1) vanishes. The remaining product is `I`
exactly when `g_i=h_i` for every position. Together with even parity, that is
`x=v` and `y=u`. Thus

```math
F_t(1)=E_t,\qquad
E_t(x,y,u,v)=[x=v][y=u].
\tag{3}
```

The degree is at most `(a+1)t`, so writing `F_t(z)=sum_k z^k A_k` gives an exact
crossing decomposition with at most `(a+1)t+1` pieces. For fixed delay and fixed
local matrices, entry generation is polynomial in `t`. Explicitly writing every
piece takes `2^{O(t)}` operations and bits. This is a uniform constant-bond family
with polynomially many pieces; it can be used by the existing score-transfer
framework if a suitable asymptotic rank estimate is eventually proved.

When `a=0` and `B_(0,g,h)=B_(1,g,h)`, equation (2) reduces to the four-feature
family. Allowing the two matrices to differ removes its independence from `x`.
The new engine does not infer ranks from the former circulant reduction.

### Persistent witnesses

An identity letter is a local label `q` with `g=h` and `B_q=0`. It contributes
`I` at every `z`, so appending it preserves any nonzero polynomial difference.
The engine records explicit words of length at most three witnessing:

- dependence on a simultaneous four-leg XOR;
- failure of an adjacent position interchange;
- a nonzero entry outside the crossing support.

An identity letter pads each witness to every larger length. For the XOR
witness, only the original positions are flipped; the new padding positions
remain fixed. A paired identity letter is unnecessary for this conclusion.
Absence of a short witness is not asserted to prove an invariance.

## 2. Exact ranks from parity blocks

Let `d=2^t`, and use the two cuts `(x,y)|(u,v)` and `(x,u)|(y,v)`.
Even parity splits the first cut into blocks indexed by
`g=x+y=u+v`. Each block has row index `x` and column index `u`.
The second cut splits into blocks indexed by `h=x+u=y+v`, with row index `x`
and column index `y`.

For each triple `(g,h,x)`, evaluate (2). Its coefficient of degree `k` is placed
at

```math
(R_{k,g})_{x,x+h},\qquad (C_{k,h})_{x,x+g}.
```

The exact ranks are the sums of the ranks of these disjoint blocks. General
packed Gaussian elimination is used on every nonzero block. No circulant,
quadratic, or symmetry rank shortcut is assumed.

The native implementation supports `1<=t<=8`, polynomial degree below 64, and
short-witness degree `3(a+1)<64`. A block has at most 256 rows, stored in four
64-bit words. The checker compares independent full four-leg tensors in 32
cases, checks 131,072 endpoint entries, compares six larger profiles with the
independent circulant engine, and checks a saved eight-position baseline.
It also independently evaluates every reported short witness in those profiles.

## 3. Screen selection and longer profiles

The compact screen evidence audits five reproducible files. The pools are
structural subsets of the full eight-matrix space. They overlap; the case count
must not be read as a count of distinct parameter choices.

The broad pool uses eight selected diagonal truth tables and twelve transition
truth tables. The efficient pool uses all 52 scalar diagonal pencils whose
one-position rank sum is four, together with singleton-supported upper and
lower transitions. Its local criterion is computed directly by four-by-four
Gaussian elimination. A smaller sparse pool supplies the delay-two comparison.
Every listed input is regenerated in a fixed order and checked against its raw
record. The audit checks record counts, all rank products, and integer lower
bounds for every rank sum. It also hashes a normalized exact-rank record, so a
metadata change does not alter the semantic reproduction comparison.

At four positions, the minimum integer rank-sum lower bound is at least 256 in
every pool except the delay-two pool, where it is at least 343. These are much
larger than `2^(sqrt(2)*4)<64`. These short screens serve model selection.

After discarding locally crossing-supported models, a shortlist was checked for
actual order, XOR, and off-crossing witnesses. Eleven models represent the
rank-sum/max-product tradeoff, both displayed boundaries, and delays zero through
two. Full integer ranks for lengths 4, 6, and 8 are retained in
`evidence/eight-feature-extensions.json`.

For a profile, define

```math
P_t=\max_k r_kc_k,\quad S_t=\sum_k\sqrt{r_kc_k},\quad
\rho_P(t)=\frac{\log_2P_t}{4t},\quad
\rho_S(t)=\frac{\log_2S_t}{2t}.
```

The last column measures `log2(P_8/P_6)/8`; it is an observed growth increment,
not an asymptotic exponent.

| Model | Boundary | Delay | rho_P(4) | rho_P(6) | rho_P(8) | rho_S(8) | P increment 6 to 8 |
|---|---:|---:|---:|---:|---:|---:|---:|
| model-00 | 0 | 0 | 0.778602 | 0.840511 | 0.879225 | 1.037842 | 0.995367 |
| model-01 | 0 | 0 | 0.780991 | 0.830254 | 0.869591 | 1.038284 | 0.987600 |
| model-02 | 0 | 1 | 0.787973 | 0.822000 | 0.864269 | 1.026921 | 0.991077 |
| model-03 | 0 | 0 | 0.823120 | 0.860161 | 0.883080 | 1.000000 | 0.951839 |
| model-04 | 0 | 0 | 0.955030 | 0.944411 | 0.953632 | 1.021271 | 0.981293 |
| model-05 | 0 | 0 | 0.778602 | 0.840511 | 0.879225 | 1.037842 | 0.995367 |
| model-06 | 0 | 1 | 0.787973 | 0.822000 | 0.864269 | 1.026921 | 0.991077 |
| model-07 | 0 | 2 | 0.823120 | 0.860161 | 0.884315 | 1.048692 | 0.956776 |
| model-08 | 0 | 2 | 0.832621 | 0.916549 | 0.929236 | 1.069077 | 0.967297 |
| model-09 | 1 | 0 | 0.809637 | 0.862578 | 0.875463 | 1.031675 | 0.914116 |
| model-10 | 1 | 0 | 0.918016 | 0.918369 | 0.943213 | 1.025589 | 1.017743 |

The best displayed eight-position `rho_P` is approximately 0.864269. Every
selected eight-position rank-sum exponent is at least one. The profiles provide
no finite seed satisfying `S_t<2^(sqrt(2)t)`. They also do not establish a lower
bound for the unrestricted eight-feature family at other lengths.

## 4. An order-dependent monoid with exact coefficient ranks

One selected model admits a complete all-length analysis. Take boundary `ell_0`,
delay zero, and the matrices encoded by

```text
[9,2,4,9,0,0,0,0].
```

Here bits `1,2,4,8` encode matrix positions `00,01,10,11`. Its local matrices are
`zI` for `(x,g,h)=(0,0,0),(0,1,1)`, `(1+z)E_01` for `(0,0,1)`,
`(1+z)E_10` for `(0,1,0)`, and `[g=h]I` when `x=1`.

The one-letter values at `q=0` and `q=4` are `z` and `1`. The two-letter words
`[1,2]` and `[2,1]` have values `(1+z)^2` and zero. Thus the family depends on
XOR translations and on position order, and has off-crossing entries. Padding
with `q=4` proves these properties at every length at least two.

Read its row matrix in order `(x,y)|(v,u)`, so the crossing is the identity.
At physical position `i`, let `D_i(z)` be diagonal with value `z` when `x_i=0`
and value one when `x_i=1`. Let

```math
e_i=|00\rangle\langle11|,\qquad
f_i=|01\rangle\langle10|,
```

where matrix rows precede columns. All pairwise products of `e_i,f_i` vanish.
Also `D_i(z)e_i=ze_i` and `e_iD_i(w)=e_i`, with the same identities for `f_i`.
The local matrix is therefore

```math
L_i(z)=D_i(z)I+(1+z)(e_iU+f_iV),\qquad
U=E_{01},\quad V=E_{10}.
\tag{4}
```

### The multiplicative law

Stack the row operators `F_t(z)` and `F_t(w)` using independent auxiliary bits
`A,B`. Factors at different physical positions commute across the two copies.
At a fixed position, the pairwise-zero identities give

```math
L_i^A(z)L_i^B(w)
=D_i(zw)I+e_iU'+f_iV',
```

where

```math
U'=aU_A+bU_B,\quad V'=aV_A+bV_B,\quad
a=1+z,\quad b=z(1+w),\quad c=a+b=1+zw.
```

Starting at the auxiliary row `v_0=<00|`, put `v_1=v_0U'`. Direct multiplication
in characteristic two gives

```math
v_0V'=0,\qquad v_1U'=0,\qquad v_1V'=c^2v_0.
```

The `00`-to-`00` matrix element of any word in `U',V'` is consequently the same
as the `0`-to-`0` matrix element of that word in `cU,cV`. This includes `c=0`;
no division by `c` is used. Applying the equality to every physical monomial
in the ordered product proves

```math
F_t(z)F_t(w)=F_t(zw).
\tag{5}
```

Writing `F_t(z)=sum_k z^k A_k` and comparing coefficients in (5) gives

```math
A_iA_j=[i=j]A_i,\qquad \sum_kA_k=I.
\tag{6}
```

Thus the coefficients are pairwise orthogonal idempotents.

### Exact ranks in both cuts

Every off-diagonal local term has first physical bit zero in its row and
first physical bit one in its column.
The normalized row matrix is upper triangular when the first word is ordered
numerically. Its diagonal entry on `(x,y)` is `z^(t-|x|)`.
The diagonal of `A_k` therefore contains exactly `2^t binom(t,k)` ones.
An idempotent matrix has only eigenvalues zero and one and is diagonalizable,
so its rank equals this integer multiplicity. Hence

```math
\operatorname{rank}_{row}A_k=2^t\binom tk.
\tag{7}
```

Swapping `y,u` exchanges `g,h` and transposes each auxiliary matrix. The boundary
`M_00` turns a product of transposed matrices into the reversed product.
Reversal of all physical positions only permutes rows and columns. Thus the two
cut ranks are equal for every coefficient:

```math
\operatorname{rank}_{col}A_k=2^t\binom tk.
\tag{8}
```

It follows exactly that

```math
S_t=4^t,\qquad
P_t=4^t\binom{t}{\lfloor t/2\rfloor}^{\!2}.
\tag{9}
```

In particular `rho_S(t)=1` for every positive `t`, while `rho_P(t)` tends to one.
This conclusion uses the exchange law, not an extrapolation from short profiles.
The checker verifies 54 complete coefficient products through length four,
511 formal two-auxiliary vacuum moments, and the independently generated
lengths 4, 6, and 8 profiles.

## 5. Reproduction

Run from the repository root. Native checks require a C++17 compiler.

```bash
python3 tools/check_eight_feature_crossings.py --output /tmp/eight-feature-verification.json
python3 tools/reproduce_eight_feature_extensions.py --directory /tmp/eight-feature-extensions --output /tmp/eight-feature-extensions.json
python3 tools/check_eight_feature_monoid.py --output /tmp/eight-feature-monoid-verification.json
```

The monoid checker reads the committed extension evidence. To rebuild the five
raw screens and their audit, use a separate temporary output directory:

```bash
python3 tools/reproduce_eight_feature_screen.py --directory /tmp/eight-feature-screen --output /tmp/eight-feature-screen.json
```

Raw JSONL files are deliberately not committed. The normalized exact-rank hashes
in the compact audit are the reproduction comparison; raw-file hashes can differ
when newer metadata fields are present. Both reproduction drivers accept
`--resume`; the screen driver checks complete record counts and parameter order
before reusing a file. Native output is explicitly flushed with stream errors
enabled, and each completed larger profile is checked before its aggregate
checkpoint is replaced.

The native profile interface is

```text
eight_feature_crossings profile t output.json boundary B0 B1 B2 B3 B4 B5 B6 B7 [delay]
```

The chosen matrices and boundaries are embedded in the extension reproducer and
also appear in every retained profile. The current evidence supports further
targeted perturbations that break the idempotent or triangular structure; it
does not support a score improvement for the screened models.
