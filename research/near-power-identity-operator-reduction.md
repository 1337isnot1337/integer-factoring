# Identity-operator cores just below a power of two

## Result and scope

Let \(n=2^k-1\). An arbitrary occupation-preserving operator on
\(r\) selected polynomial identity directions has a regular zero
reduction to \(2r+1\) active Gaussian modes in either of these cases:

| Actual selected identity directions | Retained auxiliary \(X,V\) indices |
|---|---|
| \(n-r,\ldots,n-1\) | \(0,\ldots,r\) |
| \(0,n-r+1,\ldots,n-1\) | \(0,\ldots,r-1,n-1\) |

Assume \(1\le r<n\); the second case requires \(r\ge2\).
The statement holds over any field of characteristic two. The exact
implementation below uses \(\mathbb F_4\).

One formal auxiliary identity direction permits a further local Weyl
reduction to a matrix of order \(2^{r+1}\). An exact multiplicity
division removes its known Gaussian factor afterward. Thus a
three-direction zero profile of order 128 is obtained from a matrix
of order 16. The formal direction adds no mode to the physical
crossing, no endpoint condition, and no global degree cost.

This extends the operator search beyond prime-power interval lengths,
where the contraction matrix is scalar. It is a construction and
profile theorem; it does not assert an improved factoring exponent.

## 1. A rectangular simultaneous concentration

Write

\[
N=n+1=2^k,\qquad q=z^2,\qquad \tau=q^P,\qquad c=1+\tau,
\tag{1}
\]

where \(P\) is a positive odd integer. Use the existing polynomial
forms

\[
f_j(i;c)=\sum_h\binom jh\binom ih c^h.
\tag{2}
\]

Append the formal column with index \(j=n\) to the actual selected
columns. This column is \(f_n(i;c)=\tau^i\). It is an interpolation
device at this stage, not a physical identity summand. Let
\(\widehat S=S\cup\{n\}\), let \(T\) be the retained set in
the table, and define

\[
E_{ij}=f_j(i;c),\quad 0\le i<n,\ j\in\widehat S,
\qquad F=E_T.
\tag{3}
\]

There are \(r+1\) columns and retained rows. The parent interval
length \(N\) gives the reversal identity

\[
f_{n-h}(i;c)=\tau^i f_h(i;c/\tau).
\tag{4}
\]

It is the characteristic-two specialization of the prime-power
reversal identity in the
[mixed-end reduction](split-identity-mode-reduction.md). Its use here
is justified by the parent length \(N\), even though the physical
interval contains only the first \(n\) rows.

### Terminal directions

In the first case, the reversed column degrees \(h\) run from
zero to \(r\). Factoring (4) through the constant binomial basis
\(\binom{i}{h}\) shows that interpolation on rows \(0,\ldots,r\)
has the exact form

\[
E_i=B_iF,\qquad (B_i)_t=\tau^{i-t}(B_0)_{it},
\quad 0\le t\le r,
\tag{5}
\]

where \(B_0\) is a constant matrix over the base field. The retained
binomial evaluation matrix is triangular with diagonal one. The
retained terminal matrix has determinant equal to a power of \(c\),
up to a nonzero constant and a column permutation. In particular,
\(F\) is a zero unit.

For every discarded row \(i>r\), (5) gives

\[
\operatorname{ord}_\tau(B_i)_t\ge i-t.
\tag{6}
\]

### One initial and the remaining terminal directions

Put \(b=r-1\). The augmented terminal columns are
\(n-b,\ldots,n\), with low retained rows \(L=\{0,\ldots,b\}\).
The extra initial column is the constant column \(f_0(i;c)=1\),
and the high retained row is \(H=n-1\).

Let \(r_i\) be terminal-only interpolation on \(L\). Formula (5)
gives \((r_i)_t=\tau^{i-t}(B_0)_{it}\). The full interpolation
coefficients satisfy

\[
B_{i,H}=\frac{1-\sum_{t\in L}(r_i)_t}
                   {1-\sum_{t\in L}(r_H)_t},
\qquad
B_{i,L}=r_i-B_{i,H}r_H.
\tag{7}
\]

The denominator in (7) is a zero unit because every entry of
\(r_H\) has positive order. For a discarded row \(b<i<H\),
\(B_{i,H}\) is regular, and each low component has order at least
\(i-t\). The second term has order at least \(H-t>i-t\).
The retained matrix is also a unit: modulo \(\tau\), its terminal
low block is invertible, the terminal high row is zero, and the
initial high entry is one.

Both cases therefore satisfy the same complete weighted certificate:

\[
\boxed{\operatorname{ord}_\tau B_{i,t}\ge\max(0,i-t).}
\tag{8}
\]

The auxiliary diagonal orders are \(2P(i+1)\) at aggregate gauge
\(J=1\). Thus (8) makes both interpolation and its diagonal Gaussian
companion regular. The weighted concentration theorem from the
[coupled reduction](coupled-identity-filter-reduction.md) applies.
It retains \(r+1\) auxiliary modes and the \(r\) actual identity
modes, for a total of \(2r+1\).

## 2. Complete polynomial data for the rectangular core

The implementation uses exact adjugates to avoid series truncation.
Put

\[
\delta=\det F,\qquad N_i=E_i\operatorname{adj}(F),
\qquad B_i=N_i/\delta.
\tag{9}
\]

After scaling retained \(V\) variables by the zero unit \(\delta\),
the augmented form matrices are

\[
F_X=F,\qquad
(F_V)_{tj}=\delta F_{tj}
 +\sum_{i\notin T}\tau^{t-i}(N_i)_t E_{ij}.
\tag{10}
\]

All exponents in (10) are nonnegative by (8). These are local basis
changes; no denominator is removed from the global polynomial tensor.
For completeness, \(F_V\) is a zero unit uniformly. With
\(D_i=\tau^{i+1}\), write

\[
F_V=\delta W F,\qquad
W=I+D_T B^{\mathsf T}D_{T^c}^{-1}B.
\tag{11}
\]

Each diagonal correction in \(W\) has positive order. For retained
indices \(t>u\), each correction has order at least \(t-u>0\).
Thus \(W\) reduces to an upper triangular matrix with diagonal one.
The weighted conjugate and inverse are regular as well.

Use only the columns in \(S\) to construct the actual rectangular
core. At \(J=1\), its linear forms are

\[
L_j=\tau^jY_j+c\sum_{t\in T}(F_V)_{tj}V_t,
\qquad
R_j=U_j+\tau^{-j}\sum_{t\in T}(F_X)_{tj}X_t.
\tag{12}
\]

Its diagonal Gaussian is

\[
G_{\rm act}=\prod_{j\in S}(1+\tau^jY_jU_j)
             \prod_{t\in T}(1+\delta\tau^{t+1}V_tX_t).
\tag{13}
\]

The standard identity-operator construction applies to (12)--(13).
The row/column charge is \(|Y|-|X|\), ranging from \(-r-1\) to
\(r\). General aggregate gauge \(J\) adds
\(P(1-J)\) times this charge to each Smith order.

## 3. A local dummy direction and exact profile division

Temporarily append one independent Gaussian pair
\(1+\tau^nY_*U_*\). Use the formal column \(n\) of (10) to
define \(L_*\) and \(R_*\), and extend the occupation kernel from
\(F_{ab}\) to \(F\otimes I_2\).

The extra identity factor in the occupation sum is \(1+L_*R_*\).
It cancels the corresponding removed identity factor exactly:

\[
(1-L_*R_*)(1+L_*R_*)=1.
\tag{14}
\]

Indeed, the linear forms have square zero in characteristic two.
The augmented coefficient tensor is therefore the original
rectangular tensor times the independent Gaussian pair, with no
remaining formal-summand perturbation.

The square [weighted Weyl theorem](identity-operator-weyl-reduction.md)
now uses \(r+1\) identity modes and \(r+1\) auxiliary modes.
Its differential-operator matrix has order \(2^{r+1}\). All rational
denominators are cleared by a proved zero unit, using complete
polynomials. Section 4 below proves its weighted regularity uniformly.

Let \(\mathcal Z_c(x)\) record Smith-order multiplicities of the
rectangular charge block, with integer coefficients, and let
\(\widehat{\mathcal Z}_c(x)\) denote the augmented block. At \(J=1\),

\[
\widehat{\mathcal Z}(u,x)
=\mathcal Z(u,x)(1+u x^{2Pn}).
\tag{15}
\]

Consequently the exact recursion

\[
\mathcal Z_c(x)=\widehat{\mathcal Z}_c(x)
                      -x^{2Pn}\mathcal Z_{c-1}(x)
\tag{16}
\]

recovers every rectangular block in ascending charge order. Start
below the minimum charge with zero. Every resulting multiplicity
must be nonnegative, and the final extra block must vanish. The
implementation checks both conditions and the resulting total of
\(2^{2r+1}\) orders. At general \(J\), the dummy order would be
\(2Pn+P(1-J)\); the code instead divides at \(J=1\) and restores
charge shifts afterward.

Equation (15) concerns only local profile arithmetic. The dummy pair
is removed before evaluating the physical construction. In particular,
it is never counted as a physical letter, an infinity mode, or a
contribution to the global polynomial degree.

## 4. The near-power contraction is scalar plus rank one

This section also supplies a compact independent formula for the
parameter used by the adjugate implementation. Let the augmented
indices be \(j_1,\ldots,j_{r+1}\), and put

\[
v=(\tau^{j_1},\ldots,\tau^{j_{r+1}})^{\mathsf T},\quad
s=\mathbf1^{\mathsf T}v,\quad d=1+\tau^N,\quad e=d+cs.
\tag{17}
\]

The full parent Pascal Gram identity from the
[prime-power theorem](prime-power-joint-identity-filters.md) is
diagonal. Removing its final row, whose entries are \(\tau^j\),
gives

\[
H=E^{\mathsf T}\operatorname{diag}_{0\le i<n}(\tau^{-i-1})E
 =\tau^{-N}\left[c^{N-1}\operatorname{diag}(v)+vv^{\mathsf T}\right].
\tag{18}
\]

Since \(c^N=1+\tau^N=d\), inversion of this diagonal plus rank-one
matrix yields the Weyl parameter. More explicitly, for
\(A=\operatorname{diag}(v)\), the normal-form contraction gives
\(C^{-1}=I+A(cH)^{-1}\), and the rank-one inversion gives

\[
\boxed{
C^{-1}=d^{-1}I+\frac{\tau^N c}{d e}v\mathbf1^{\mathsf T}.
}
\tag{19}
\]

For terminal-only selection, \(e=1+\tau^{n-r}\). With one initial
and \(b=r-1\) terminal directions, \(e=\tau+\tau^{n-b}\).
Write \(e=\tau^h e_0\). In these cases \(h=0\) and \(h=1\),
respectively, and \(e_0\) is a zero unit. A complete polynomial
numerator and zero-unit denominator for (19) are

\[
\Delta=d e_0,\qquad
\mathcal N=e_0I+\tau^{N-h}c\,v\mathbf1^{\mathsf T}.
\tag{20}
\]

Both \(C^{-1}\) and its conjugate by
\(A=\operatorname{diag}(\tau^{j_1},\ldots,\tau^{j_{r+1}})\)
are the identity modulo \(\tau\). For the conjugate
\(A^{-1}C^{-1}A\), every rank-one entry has order at least
\(N-h+j_k>0\). Hence both matrices and their inverses are regular.
This is precisely the weighted assumption needed by the Weyl theorem.

## 5. Endpoint, degree, infinity, and spectators

The physical crossing uses only its original \(r\) polynomial
identity summands. An invertible endpoint preserving the physical
energy filtration is repaired by the same constant physical \(Y\)
map as in the [operator construction](coupled-identity-operator-filters.md).
The map and its inverse are regular at both cuts. The original
coefficientwise endpoint and ordinary degree bound therefore apply
without alteration.

At infinity the active matrix is still

\[
\left(F_{ab}(z^{-2})\,
z^{P\sum_{j\in S}(2j+1-2n)b_j}\right)_{a,b}.
\tag{21}
\]

The zero spectators have energies \(P(2i+1-J)\) for
\(i\notin S\) and \(P(2i+1+J)\) for \(i\notin T\).
The infinity spectators retain the ordinary Gaussian spectrum with
only the actual selected directions removed. Both paired physical
partitions remain explicit. No pressure calculation omits the extra
retained auxiliary mode or counts the formal dummy after division.

## 6. Implementation and verification

The rectangular constructor and original physical tensor builder are
[`near_power_terminal_operators.py`](../tools/near_power_terminal_operators.py).
The dummy Weyl engine and complete F4 affine compiler are
[`near_power_terminal_weyl.py`](../tools/near_power_terminal_weyl.py).
The optional parameter `initial_modes=1` selects the second geometry;
the default is the terminal-only geometry.

The independent
[`checker`](../tools/check_near_power_terminal_operators.py) verifies
Pascal reversal, exact interpolation, weighted orders, Gaussian
recovery, complete rectangular versus Weyl profiles, compiled versus
literal polynomial profiles, and the original physical endpoints and
both cuts. A degree-600 coefficient in the second F4 bit plane checks
that complete sources survive the compiler. The saved
[`extended evidence`](../evidence/near-power-terminal-operators.json)
also compares the pressure evaluator with independently expanded
spectator and physical partitions.

The extended replay includes nine complete rectangular/Weyl comparisons,
68 compiled two-cut comparisons, ten original physical cut profiles,
3,964 endpoint coefficients, and 81 independent pressure evaluations.
It checks 1,039 reversal identities and 878 interpolation identities.
The closed rank-one formula is cross-multiplied against the separate
adjugate construction in 31 geometries, covering 438 matrix entries.

```bash
OPENBLAS_NUM_THREADS=1 python3 tools/check_near_power_terminal_operators.py --extended
```

The finite
[`projector driver`](../tools/search_near_power_terminal_operators.py)
enumerates the same 42,560 F4 occupation-sector kernels as the earlier
three-direction search. Its profiles are exact; scale and aggregate
charge optimization remain numerical search evidence. A completed
finite run does not certify optimality over all Laurent operators or
all continuous pressure parameters.

Two complete \(n=31\), \(P=1\) runs are retained in
[`near-power-operator-projectors-n31.json`](../evidence/near-power-operator-projectors-n31.json).
The terminal-only selection \((28,29,30)\) gives 973 exact profile
classes. The split selection \((0,29,30)\) gives 1,153. Each run
enumerates all 42,560 specified kernels and uses a 72-point scale/charge
grid, followed by 30 local refinements. Both attain the Gaussian
critical pressure, approximately \(0.033915091865694\), without a
detected improvement. Relative differences of magnitude below
\(6\times10^{-15}\) are floating-point noise, not negative certificates.

Both saved runs and the extended checker were completed again after
the execution-environment interruption. No incomplete pre-interruption
run is used as a complete result here.
