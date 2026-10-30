# Regular roots with a rank-two terminal power

## Scope

This note supplies an exact family of nilpotent crossing paths. For every
`q=2^a`, `a>=1`, it attains the ideal column Smith spectrum and has row
identity at parameter zero. Its underlying unipotent operator also has a
specified high power of row rank two. The ranks of its early powers remain
an independent obligation; the displayed larger examples do not improve
the factoring exponent.

A separate theorem excludes the Smith-profile certificate for **every
full-cycle Feistel function** on `Z/qZ` when `q>=8`. This theorem does not
exclude Feistel permutations with several cycles or arbitrary coefficient
recombinations outside that sufficient certificate.

All construction identities hold over every field of characteristic two.
The matrices and the supplied witnesses have binary entries. There is no
new Lean formalization.

## 1. A modified transpose with only two-dimensional Jordan blocks

Let `q=2t`, `n=q^2`, and let `N` be block diagonal with `t` copies of

\[
N_2=\begin{pmatrix}1&1\\1&1\end{pmatrix}.
\]

Thus `N=N^T`, `N^2=0`, and `rank N=q/2`. Define a linear map on `Mat_q` by

\[
V(X)=X^{\mathsf T}+NX^{\mathsf T}N,
\qquad T=I-V=I+V.
\tag{1}
\]

Use the reshuffling convention

\[
\Phi(M)(x,y,u,v)=M(E_{vy})_{xu}
\tag{2}
\]

from the [nilpotent crossing construction](nilpotent-flag-crossings.md).
Then

\[
\boxed{\quad V^2=I,\qquad \operatorname{rank}T=n/2,
       \qquad \operatorname{rank}\Phi(V)=2.\quad}
\tag{3}
\]

For the first identity, write `V=tau (I+L_N R_N)`. Transposition commutes
with `L_N R_N` because `N` is symmetric. The square is the identity because
`N^2=0` and the characteristic is two.

For the row rank, the two terms in (1) give

\[
\Phi(V)=\operatorname{vec}(I)\operatorname{vec}(I)^{\mathsf T}
       +\operatorname{vec}(N)\operatorname{vec}(N)^{\mathsf T}.
\]

The two column vectors, and the corresponding two row vectors, are
independent. This proves rank exactly two.

To prove the middle equality of (3), identify
`Mat_(2t)=Mat_t tensor Mat_2`. Under this identification

\[
V=\tau_t\otimes V_2.
\]

On `Mat_2`, the image of `V_2-I` is spanned by the identity and the
all-one matrix. Its rank is two and its square is zero. Hence `V_2` is a
direct sum of two size-two unipotent Jordan blocks.

If `W^2=I` on a vector space of dimension `d`, then `W tensor J_2(1)`
has a fixed space of dimension exactly `d`. Indeed its fixed-point
equations, in block form, determine the second component uniquely as
`(W-I)` times the first; the remaining equation follows from `(W-I)^2=0`.
Apply this to `W=tau_t` and to the two blocks of `V_2`. The fixed space
of `V` has dimension `2t^2=n/2`, proving (3). Consequently `T` has exactly
`n/2` nilpotent Jordan blocks, all of size two.

Both `V` and `T` have the required normalizations:

\[
V(I)=I,\quad \operatorname{Tr}(VX)=\operatorname{Tr}X,
\quad T(I)=0,\quad\operatorname{Tr}(TX)=0.
\tag{4}
\]

The trace identity uses `Tr(N X^T N)=Tr(N^2 X^T)=0`.

## 2. Explicit paired-chain construction of a regular root

Now assume `q` is a power of two and put `h=n/2`. We construct a regular
nilpotent `S` satisfying

\[
\boxed{\quad S^h=T,\quad S(I)=0,\quad
        \operatorname{Tr}(SX)=0,\quad
        S^{n-1}(X)=\operatorname{Tr}(X)I.\quad}
\tag{5}
\]

The proof is constructive linear algebra. Since `T^2=0` and `rank T=h`,
its image equals its kernel. Equation (4) implies that this common space
consists of trace-zero matrices. Thus the functional

\[
\ell:\operatorname{im}T\longrightarrow F,
\qquad \ell(TX)=\operatorname{Tr}X
\tag{6}
\]

is well-defined. It is nonzero because trace is a nonzero functional on
the full matrix algebra.

Let

\[
B_0=I_t\otimes E_{01}.
\]

A direct two-by-two calculation gives `T(B_0)=I`, and `Tr(B_0)=0`.
Therefore `ell(I)=0`. Choose a basis

\[
A_0=I,A_1,\ldots,A_{h-2}
\]

of `ker ell`, and choose `A_(h-1)` with `ell(A_(h-1))=1`.
For every `j`, choose `B_j` with `T(B_j)=A_j`, retaining the displayed
choice of `B_0`. By (6), the traces of these preimages are forced:
`Tr(B_j)=0` for `j<h-1`, and `Tr(B_(h-1))=1`.

The ordered list

\[
A_0,\ldots,A_{h-1},B_0,\ldots,B_{h-1}
\tag{7}
\]

is a basis: applying `T` to a dependence first eliminates all preimage
coefficients, after which the image basis eliminates the remainder.
Define `S` to send each vector of (7) to its predecessor and to kill
`A_0`. This is one Jordan chain of length `n`. Its `h`-th power sends
`B_j` to `A_j` and kills every `A_j`, so it equals `T`. All basis vectors
except the last have trace zero, and the first is `I`. These facts prove
all remaining assertions in (5).

The checker makes every choice deterministic. It starts with `A_0=I`,
uses the standard trace-zero matrix basis and exact greedy independence
tests to extend `ker ell`, and takes `B_(h-1)=E_00`. This last choice is
valid because `ell(T(E_00))=1`. Gaussian elimination then constructs `S`
in the physical matrix-unit basis.

Set `U=I+S`. Since `h` is a power of two, the characteristic-two binomial
identity gives

\[
U^h=I+S^h=V.
\tag{8}
\]

Thus `U` is regular unipotent of order `n`, while its terminal nontrivial
dyadic power has row rank two under `Phi`. The construction of `S`
itself takes polynomially many field operations in `q`.

## 3. The resulting crossing paths

Apply the general nilpotent construction at index `n`:

\[
K(z)=\Phi\left(\sum_{j=0}^{n-1}
                    z^{n-1-j}S^j(I-S)\right).
\tag{9}
\]

Then

\[
K(1)=E,\qquad K(0)=I_{q^2}\text{ in the row cut},
\qquad \nu_0=(0,n,\ldots,n)\text{ in the column cut}.
\tag{10}
\]

The crossing endpoint and column spectrum follow from the general
construction. At zero, (5) gives the operator `X -> Tr(X)I`, whose row
reshuffling is exactly the identity matrix.

Equation (8) adds a condition not supplied by arbitrary regular flags:

\[
\operatorname{rank}\Phi(U^{n/2})=2.
\tag{11}
\]

It does not bound `rank Phi(U)`, `rank Phi(U^2)`, the intermediate powers,
or the complete row Smith spectrum. Those remain necessary parts of a
successful rank-cost construction.

The canonical paired-chain witnesses have the following exact row ranks
at successive dyadic powers `U,U^2,U^4,...,U^(n/2)`:

| `q` | Ranks |
| --- | --- |
| 2 | `2, 2` |
| 4 | `11, 16, 15, 2` |
| 8 | `58, 63, 63, 64, 63, 2` |
| 16 | `250, 254, 255, 255, 255, 255, 254, 2` |

The larger canonical witnesses therefore have large early ranks. No
improved factoring exponent follows from (10)--(11).

For a simple two-label control, let `A=E_01` and

\[
U_2(X)=X^{\mathsf T}(I+A)+AX^{\mathsf T}.
\tag{12}
\]

Its matrix-unit columns are `(3,13,2,10)` in row-major binary masks.
Its square is `V_2`, its row rank is two, and `I+U_2` is regular nilpotent.
Exact polynomial elimination gives row determinant one for (9), row
orders `(0,0,1,11)`, and column orders `(0,4,4,4)`. This agrees with the
existing two-label control; it is not a better numerical seed.

## 4. Full-cycle Feistel maps have an exact Smith barrier

Let `q=2^a`, `n=q^2`, and let `f:Z/qZ -> Z/qZ` be any function. Consider
the permutation of matrix units induced by

\[
U(i,j)=(j,i+f(j)).
\tag{13}
\]

Assume this permutation is a single cycle of length `n`. In characteristic
two, `S=I+U` is then regular nilpotent, so its column spectrum is the
ideal spectrum in (10). Put

\[
r=|f(Z/qZ)|,\qquad R_2=\operatorname{rank}\Phi(U^2).
\]

The first row reshuffling has rank exactly `r`. Its nonzero rows are
cyclic-shift permutation matrices indexed by the distinct values of `f`;
these matrices have disjoint supports and are independent over every
field.

For the squared map, every row `(x,y)` of `Phi(U^2)` has one nonzero
entry, at the column

\[
\Psi(x,y)=(y+f(x),x-f(y)).
\tag{14}
\]

Consequently `R_2=|im Psi|`. Every fiber of `Psi` has at most `r` points:
for a prescribed value `(a,b)`, choosing `f(x)` fixes `y=a-f(x)`, and
then fixes `x=b+f(y)`. Therefore

\[
R_2\ge n/r.
\tag{15}
\]

Write the reversed row polynomial as

\[
H(w)=\Phi\bigl(U(I-wS)^{-1}\bigr).
\]

Its constant coefficient has rank `r`, and its linear coefficient is
`Phi(U+U^2)`. Eliminate the constant pivots. The first Schur coefficient
is the compression of `Phi(U^2)` between a kernel and a quotient of
codimension `r`. It has rank at least `R_2-2r`. Hence the row Smith
orders contain exactly `r` zero orders and at least

\[
k=\max(0,n/r-2r)
\tag{16}
\]

orders equal to one; using a real lower bound in (16) only weakens the
integer count.

For the Smith-profile sufficient certificate, let

\[
B(s)=\left(\sum_{i=1}^{n}s^{n-1-\mu_i}\right)
     \left(1+(n-1)s^{-n}\right),\qquad s\ge1.
\tag{17}
\]

Every zero order contributes

\[
a_0(s)=s^{n-1}+(n-1)s^{-1}\ge n.
\tag{18}
\]

Its derivative is nonnegative for `s>=1`, so the minimum occurs at one.
Every order-one pivot contributes

\[
a_1(s)=s^{n-2}+(n-1)s^{-2}\ge n-1\qquad(n\ge4).
\tag{19}
\]

For completeness, its critical point satisfies
`s^n=2(n-1)/(n-2)`. The minimum is

\[
\frac{n(n-1)}{n-2}
\left(\frac{2(n-1)}{n-2}\right)^{-2/n}.
\]

This is at least `n-1` because Bernoulli's inequality gives

\[
\left(\frac{n}{n-2}\right)^{n/2}
=\left(1+\frac{2}{n-2}\right)^{n/2}
\ge1+\frac{n}{n-2}
=\frac{2(n-1)}{n-2}.
\]

Using only these zero and order-one contributions proves

\[
B(s)\ge nr+(n-1)\max(0,n/r-2r)
       \ge\frac{n^{3/2}}{\sqrt2}.
\tag{20}
\]

For the last minimization, if `r>=sqrt(n/2)`, the first term alone
suffices. Below that point, the expression is
`n(n-1)/r-(n-2)r`, which decreases with `r`; it has the displayed value
at the boundary.

Since `n=q^2`, the right side of (20) is `q^3/sqrt(2)`. For every `q>=8`,

\[
\frac{q^3}{\sqrt2}>q^{2\sqrt2}.
\tag{21}
\]

It suffices to check `q=8`, since `3-2sqrt(2)>0`. At eight, (21) is the
exact comparison `17/2>6sqrt(2)`, equivalent to `289>288`. Equivalently,
the normalized Smith exponent obeys

\[
\rho\ge\frac34-\frac{\log2}{8\log q}
       \ge\frac{17}{24}>\frac1{\sqrt2}.
\tag{22}
\]

Thus no full-cycle Feistel function at these sizes passes this sufficient
Smith certificate. This statement covers all functions `f`, without a
finite enumeration and over every characteristic-two extension field.
It is specifically a statement about (17), not a lower bound for every
possible decomposition of the same tensor powers.

## 5. Reproduction and remaining direction

Run the standard-library checker:

```sh
python3 tools/check_modified_transpose_roots.py \
  --output evidence/modified-transpose-roots.json
```

It constructs the paired chains at `q=2,4,8,16`, checks nilpotence,
regularity, the exact terminal power, both normalizations in (5), and all
reported dyadic ranks. It also checks the two-label determinant and Smith
orders, and independently verifies the two Feistel rank formulas on a
full-cycle four-label control. The all-size arguments are the proofs
above; the finite tests are controls, not a substitute for those proofs.

The constructive open target is a different regular root of (1) with
simultaneously small early row ranks, followed by a complete rank profile
or an explicit tensor-power decomposition. The paired-chain method
provides all roots obtained by changing the paired basis. The canonical
choices recorded here do not meet that target.
