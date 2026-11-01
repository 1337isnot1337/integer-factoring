# Four-feature matrix-product crossing identities

## Construction and scope

This is a constant-auxiliary-dimension family over `F_2` that adds one
local feature to the `R,H,C` constructions in
[finite-state-pairing-crossings.md](finite-state-pairing-crossings.md).
Its coefficients can depend on the order of the bit positions, while
their sum is always the exact crossing tensor. The construction has
degree at most the word length, so its coefficient count is linear.

The complete finite screen tests **131,072** parameter choices at six
positions. It finds no strict finite rank-sum certificate and makes no
asymptotic rank or factoring-score claim. In particular, the finite
results do not prove an obstruction for all constant-state constructions.

For a local bit quadruple `q=(x,y,u,v)`, set the local matrix to zero
when `x+y+u+v` is odd. Otherwise put

\[
 g=x+y,\qquad h=x+u,
\]

with addition in `F_2`. Choose four arbitrary matrices
`B_{00},B_{01},B_{10},B_{11}` in `Mat_2(F_2)` and define

\[
 K_{g,h}(z)=[g=h]I+(1+z)B_{g,h}.                         \tag{1}
\]

For a boundary functional `ell` with `ell(I)=1`, define

\[
 F_t(z;q_1,\ldots,q_t)
   =\ell\!\left(K_{g_1,h_1}(z)\cdots K_{g_t,h_t}(z)\right).
                                                               \tag{2}
\]

An odd-parity local quadruple makes the whole entry zero. At `z=1`,
each remaining factor is `[g_i=h_i]I`. The conditions of even parity
and `g_i=h_i` are equivalent to `x_i=v_i` and `y_i=u_i`. Therefore

\[
 \boxed{F_t(1)=E_t.}                                    \tag{3}
\]

The degree is at most `t`. Thus writing `F_t=sum_{d=0}^t z^d A_d`
gives an exact crossing decomposition with at most `t+1` pieces.
The auxiliary matrices have dimension two, independent of `t`.

## The additional local feature

On even-parity inputs, the values of `R,H,C` are respectively

\[
 1-g,\quad1,\quad1-h.
\]

Their span consists of affine functions of `g,h`. A matrix-valued
function belongs to this span exactly when its mixed difference

\[
 B_{00}+B_{01}+B_{10}+B_{11}
\]

vanishes. The four independent matrices in (1) allow this mixed
difference to be nonzero. Equivalently the construction allows a
nonzero `gh` coefficient. It therefore enlarges the local tensor
space, rather than only changing the transition costs of the earlier
three-symbol family.

## Order dependence at every length

Take

\[
 B_{00}=B_{11}=0,\qquad B_{01}=E_{01},\qquad B_{10}=E_{10},
 \qquad \ell(M)=M_{00}.                                 \tag{4}
\]

The two off-diagonal letters satisfy

\[
 \ell(K_{01}K_{10})=(1+z)^2=1+z^2,
 \qquad
 \ell(K_{10}K_{01})=0.
\]

They are related by interchanging the positions, so at least one
coefficient is not position invariant. Here `K_{00}=I`, so padding
both words by `00` letters preserves the difference at every length
`t>=2`. The mixed feature in (4) is `E_{01}+E_{10}`, which is nonzero.
Thus this is a genuine all-length example outside full position
invariance, with the exact endpoint (3) and only linearly many pieces.

The example is not itself a good rank certificate. Its role is to
verify that the extra local feature and order dependence are both
present in a scalable construction.

## Exact rank computation and independent checks

The tensor still has local even-parity support and is unchanged by
simultaneously flipping the four local bits. Consequently the exact
XOR-circulant reduction from the preceding note applies. In the row
cut, each block at `g=x+y=u+v` has entries `K_d(g,x+u)`, where the
kernel polynomial is evaluated by the two-state matrix product (2).
In the column cut the corresponding kernel grid is transposed. Both
cut ranks are sums of exact binary circulant ranks.

The native implementation is
[`parity_matrix_crossings.cpp`](../tools/parity_matrix_crossings.cpp).
Its fast rank routines use the already proved constant/linear
square-zero-algebra shortcuts and exact binary Gaussian elimination.

[`check_parity_matrix_crossings.py`](../tools/check_parity_matrix_crossings.py)
independently evaluates every entry of 24 three-position tensors, then
forms both full cut matrices and compares their exact ranks with the
native output. It verifies the endpoint on 98,304 entries. It also
checks an explicit order witness and its identity-padding extension
through six positions. These computations use the original four-leg
tensor entries, not the circulant block formula. Evidence is recorded
in `evidence/parity-matrix-verification.json`.

## Exhaustive six-position screen

There are `16^4=65,536` four-matrix choices. The screen tests every
choice with each of

\[
 \ell_0(M)=M_{00},\qquad
 \ell_1(M)=M_{01}+M_{10}+M_{11}.
\]

As explained in the preceding note, these represent all conjugacy
classes of binary boundary functionals with `ell(I)=1`. Simultaneous
conjugation of arbitrary local matrices preserves the exhaustive
parameter set. There are thus 131,072 recorded choices.

Let `P=max_d r_d c_d`, `S=sum_d sqrt(r_d c_d)`, and
`rho_max=log2(P)/(4t)`. The global optima at `t=6` are

\[
 \min P=132496,\qquad \min S\simeq1364.59716939.
\]

They coincide with values already found in the earlier auxiliary
family. The minimum exact integer lower bound
`sum_d floor(sqrt(r_d c_d))` is 1,363. Since

\[
 1363>512=2^9>2^{6\sqrt2},
\]

every screened model fails the strict finite seed criterion. This
exclusion uses integers; decimal square roots only order candidates.

The compact evidence also retains separate Pareto frontiers for
models with a nonzero mixed feature and for models with a certified
order witness. A witness consists of two permuted words of length
two or three whose highest polynomial coefficients differ. The
highest coefficient is `ell(B_{q_1}...B_{q_s})`, so this is an exact
binary certificate. Models with an identity diagonal letter receive
an additional flag certifying that the witness extends to every
larger length. These classes overlap and are not claimed to classify
every possible order-dependent model.

| Boundary | All choices | Nonzero mixed feature | Mixed feature and order witness | Best `P` in the last class |
|---|---:|---:|---:|---:|
| `ell_0` | 65,536 | 61,440 | 54,096 | 1,638,400 |
| `ell_1` | 65,536 | 61,440 | 60,144 | 1,638,400 |

The last value gives `rho_max approximately 0.86016067`. It supplies
finite search information about this degree-one local pencil, rather
than an asymptotic lower bound or a conclusion about larger auxiliary
dimensions, higher local degrees, or periodically varying matrices.

## Reproduction from the repository root

The native code requires a C++17 compiler. Keep raw screens outside
the repository:

```sh
mkdir -p /tmp/parity-matrix-screen
g++ -O3 -std=c++17 tools/parity_matrix_crossings.cpp -o /tmp/parity-matrix-screen/search
/tmp/parity-matrix-screen/search screen 6 /tmp/parity-matrix-screen/parity-matrix-screen-b0-t6.jsonl 0
/tmp/parity-matrix-screen/search screen 6 /tmp/parity-matrix-screen/parity-matrix-screen-b1-t6.jsonl 1
python3 tools/audit_parity_matrix_screen.py --directory /tmp/parity-matrix-screen --output /tmp/parity-matrix-screen/parity-matrix-screen.json
python3 tools/check_parity_matrix_crossings.py --output /tmp/parity-matrix-screen/parity-matrix-verification.json
```

The compact screen retains full integer ranks for its Pareto
representatives, explicit order witnesses, and hashes identifying the
two raw runs. A single model can be reproduced with the `profile`
mode, followed by its boundary index and four matrix masks. Each
matrix uses four row-major bits, with `E_00,E_01,E_10,E_11` encoded
as `1,2,4,8`, respectively.
