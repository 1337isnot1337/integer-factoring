# A copied-wire barrier for every full-cycle shift register

## Statement and scope

Every reversible shift register that is one full cycle has a structural
bare-profile obstruction, regardless of its feedback function. Growing
its length and optimizing a balanced coordinate partition gives

```math
\liminf_{m\to\infty}\rho_{\rm bare}\ge\frac34.                \tag{1}
```

Here `m` is an even register length and the physical alphabet is
`q=d^(m/2)`, where the register's digit alphabet `d` is any prime power.
The coefficient field has the same prime characteristic. The binary
case includes every de Bruijn feedback rule, including the
[linear](shift-register-crossings.md) and
[quadratic](quadratic-register-crossings.md) catalogues. All bounds are
uniform in the feedback and coordinate partition. For every
integer `1<=J<m`, put

```math
n=d^m,\qquad
E=\frac m2-\left\lceil\frac{m}{2J}\right\rceil-J,
```

the complete bare product obeys

```math
\boxed{\quad B(s)\ge\frac{d^{m+E}}{J^3}
                    \quad\hbox{for every }s\ge1.\quad}        \tag{2}
```

In particular, `rho_bare >= 17/24 > 1/sqrt(2)` for every even
`m>=576`. The explicit finite cutoff is large: in binary it corresponds
to `q>=2^288`. This leaves a gap between the completed small catalogues
and this uniform large-register result. It also leaves physical word
types and gauges, arbitrary physical bases, nonorthogonal operators,
other polynomial generator transforms, and other interpolation schemes
outside the theorem. This is a lower bound for a sufficient crossing
criterion, not for factoring algorithms. The factoring score stays `1/2`.

## 1. Why the first two power ranks do not settle the question

Let `T` swap the two physical factors. The row rank of `Phi(U^j)` is
the operator Schmidt rank of `V_j=U^j T`: its least number of terms
`A_i tensor B_i`. These ranks are submultiplicative. If `U` is a
permutation matrix, `V_j` is orthogonal and its inverse has the same
operator Schmidt rank, since taking the full transpose preserves it.
The identity

```math
T=V_j^{-1}V_{j+k}V_k^{-1}
```

therefore gives the field-independent necessary inequality

```math
q^2\le r_jr_kr_{j+k},\qquad r_j=\operatorname{rank}\Phi(U^j).
                                                               \tag{3}
```

One must not replace the right-hand side by `r_j r_(2j)`. A full
4096-cycle from the preceding construction has `r_1=17`, `r_2=18`:
its tap mask is `2100` and its first physical coordinate set is
`{0,1,3,6,7,10}`. Thus `17*18=306<4096`, whereas
`17^2*18=5202>=4096`. This is an exact counterexample within regular
full cycles, not only within smaller-order wire permutations.

Equation (3) allows first power ranks on the scale `n^(1/3)`.
The uniform obstruction below instead uses a growing range of powers,
then relates one of those ranks to the Smith spectrum with every
coefficient-expansion cost retained.

## 2. A colour-correlation lemma

Colour each bit position according to its physical half. Write the
colours as `a_i in {-1,1}`, cyclically indexed modulo `m`, and let

```math
c_j=\#\{i:a_i=a_{i+j}\},\qquad
R_j=\sum_i a_i a_{i+j}=2c_j-m.
```

For any `J>=1`, expanding squares gives

```math
0\le\sum_i\left(\sum_{t=0}^J a_{i+t}\right)^2
 =m(J+1)+2\sum_{j=1}^J(J+1-j)R_j.
```

The positive weights in the last sum total `J(J+1)/2`. Consequently
some `1<=j<=J` satisfies

```math
c_j\ge\frac m2-\frac{m}{2J}.                              \tag{4}
```

This elementary real inequality does not depend on the finite field
used for crossing ranks. It holds for every colouring, including the
balanced colourings required here.

## 3. Copied wires force disjoint nonzero blocks

A shift register on `d` digits has the form

```math
U(s)=d s\bmod d^m+f(s),\qquad 0\le f(s)<d.
```

Assume this map is a permutation; no linearity condition is imposed
on `f`. After `j<m` steps, every output digit at position `i>=j`
is exactly the input digit at position `i-j`. Let `w_j` count these
copied wires whose two positions have the same physical colour.
The cyclic agreement count `c_j` includes at most `j` additional
wrapped positions, so

```math
w_j\ge c_j-j.                                             \tag{5}
```

In the row reshuffling, the variables are `(output A,input B)` versus
`(output B,input A)`. Each same-colour copied wire therefore equates
one row digit with one column digit. These are `w_j` disjoint equality
constraints. Reorder rows and columns by their values: the matrix is
block diagonal with `d^w_j` blocks. Every block is nonzero, because
every assignment to the copied input digits extends to an input state
and the permutation supplies an edge. Thus, over every field,

```math
\operatorname{rank}\Phi(U^j)\ge d^{w_j}\ge d^{c_j-j}.        \tag{6}
```

Combining (4)--(6) yields some `1<=j<=J` with

```math
\operatorname{rank}\Phi(U^j)\ge d^E.                       \tag{7}
```

This block argument covers arbitrary nonlinear feedback. It uses
literal copied coordinates, rather than a linear approximation to
the register transition.

## 4. A power rank forces early Smith channels

For any prime-power digit alphabet, a full cycle has prime-power
length. In the corresponding characteristic, `S=I-U` is regular
nilpotent, and the complete crossing is
`K(z)=Phi(sum_(k=0)^(n-1) z^(n-1-k) U S^k)`. Its endpoint is `Phi(I)`
because `U sum_k S^k=I`; its column zero orders are `(0,n,...,n)`.
Use the normalized reversed row matrix

```math
H(w)=\Phi\bigl(U(I-wS)^{-1}\bigr)=\sum_{k=0}^{n-1}w^k A_k.
```

Let its nonnegative Smith orders be `mu_i`, and let `N_j` count those
strictly smaller than `j`. A Smith factorization over the power-series
ring expresses `A_k` as a sum of at most
`sum_(mu_i<=k)(k-mu_i+1)` rank-one matrices. This explicitly counts the
ways to divide the coefficient order between both unit factors.

Since `U^j=U(I-S)^(j-1)`, the binomial formula in any characteristic gives

```math
\operatorname{rank}\Phi(U^j)
\le\sum_{k=0}^{j-1}\operatorname{rank}A_k
\le\frac{j(j+1)}2N_j\le j^2N_j.                           \tag{8}
```

No coefficientwise rank monotonicity is assumed. By (7)--(8),
`N_j >= d^E/j^2` for one such `j`.

Each of these channels contributes at least

```math
s^{n-j}+(n-1)s^{-j}
```

to the bare product. Its minimum over `s>=1` is at least `n/j`.
Indeed, putting `x=s^n`, the unconstrained positive minimum is

```math
\frac{n(n-1)^{(n-j)/n}}
 {(n-j)^{(n-j)/n}j^{j/n}}\ge\frac nj,
```

using `(n-1)/(n-j)>=1` and `j^(j/n)<=j`. Restriction to `x>=1`
can only raise that minimum. Therefore

```math
B(s)\ge\frac njN_j
      \ge\frac{nd^E}{J^3},
```

which proves (2). The opposite bare orientation is already at least
`q^4`: the column leading coefficient is invertible and the row
determinant degree is at most `n(n-1)`. This is the same
degree/AM--GM argument as in the
[controlled-monomial predecessor](controlled-permutation-smith-floor.md).

## 5. Rate and explicit large-register cutoff

Choose `J=floor(sqrt(m))`. Then `E=m/2-O(sqrt(m))` and
`log J=O(log m)`. Equation (2) gives

```math
\log_d B(s)\ge\frac{3m}2-O(\sqrt m),\qquad
\rho_{\rm bare}\ge\frac34-O(m^{-1/2}),
```

uniformly, proving (1).

Here are explicit constants. Put `r=sqrt(m)` and suppose `m>=1296`.
Then `J>=r-1`, `J<=r`, and
`ceil(m/(2J))+J<=3r/2+2`. Thus `E>=m/2-3r/2-2`.
Since `d>=2`, equation (2) implies

```math
\log_d B(s)\ge\frac{3m}2-\frac{3r}2-3\log_2r-2
             \ge\frac{3m}2-2r.
```

For the last step use `2+3 log_2 r<=r/2` for `r>=36`: at thirty-six,
`log_2 36<16/3` proves it, and the derivative inequality follows from
`ln 2>=1/2`. Finally, `r>=36` gives `2r<=m/12`. Hence
`log_d B(s)>=17m/12`, and division by `4 log_d q=2m` proves
`rho_bare>=17/24` for this infinite tail.

For every even `576<=m<1296`, use `J=floor(sqrt(m/2))`. The 360
exact integer comparisons

```math
2^{12E-5m}\ge J^{36}
```

are all true and are replayed by the checker. They prove (2) is at
least `q^(17/6)` throughout this complete finite prefix. The same
comparison holds for every `d>=2`, because its exponent is nonnegative
and `d^(12E-5m)>=2^(12E-5m)`. Combined with the analytic tail, this
establishes the explicit cutoff `m>=576`. The comparison with
`1/sqrt(2)` is the exact integer
inequality `2*17^2>24^2`.

## 6. Verification and remaining directions

The [checker](../tools/check_shift_register_copy_wire_barrier.py) verifies
the colour lemma exhaustively on small balanced colourings, literal
copied-wire blocks for linear and nonlinear registers, coefficient rank
accounting on complete polynomial tensors, the orthogonal rank triple,
and the 4096-cycle counterexample. Ternary full-cycle controls check
the field and digit extension. It also replays all 360 finite cutoff
comparisons and selected larger sizes without materializing their tensors.
The all-size proof is Sections 2--5; finite matrix checks are its controls.

```sh
make verify-register-barrier
```

Saved record: [`shift-register-copy-wire-barrier.json`](../evidence/shift-register-copy-wire-barrier.json).
This work uses the preceding full-cycle constructor and source transfer,
retaining their credit. The abstract argument and controls used substantial
OpenAI Codex assistance. No Lean formalization or improved factoring
exponent is claimed.

The obstruction identifies what a new route must escape: the persistent
copied-coordinate wires, the physical coordinate partition, or the bare
Smith pressure estimate. Nonorthogonal regular operators and more general
physical encodings remain constructive possibilities.
