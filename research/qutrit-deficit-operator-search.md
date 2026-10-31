# Exact qutrit operators just below a power of three

## Result and scope

The [fixed-deficit theorem](fixed-deficit-identity-cores.md) now has a
working prime-field operator evaluator. In characteristic three, two
selected directions at \(n=26=3^3-1\) use a complete 27-by-27 Weyl
matrix, followed by exact removal of one three-digit Gaussian factor.
The resulting zero profile has all 243 active orders.

Two bounded complete projector families were tested, one on the terminal
indices \((24,25)\) and one on the split indices \((0,25)\).
Each family contains 1,129 distinct kernels with the exact identity
endpoint. Their profiles are exact. The numerical pressure search
found no improvement: both runs retain the Gaussian critical pressure
of approximately \(0.07235825135226\).

This is a new constructive evaluator and a finite negative search
result. It does not give an improved factoring exponent or a general
obstruction for qutrit operators.

## 1. The small exact evaluator

For an actual \(r\)-direction kernel \(F\), the evaluator uses the
polynomial Weyl operator for \(F\otimes I_{p^d}\). It includes the
general-prime signs and factorial units from the
[weighted Weyl theorem](identity-operator-weyl-reduction.md). Every
coefficient is represented as a complete finite-field Laurent
polynomial.

In the one-deficit case, a short rational parameter improves the
coefficient sizes. Put

\[
N=n+1=p^k,\quad c=1-\tau,\quad D=1-\tau^N,
\quad v_j=\tau^j,\quad e=D-c\sum_jv_j=\tau^h e_0.
\tag{1}
\]

The indices include the formal dummy direction. The parameter is

\[
C^{-1}=
\frac{e_0I+\tau^{N-h}c\,v\mathbf1^{\mathsf T}}{D e_0}.
\tag{2}
\]

Here \(e_0\) is a zero unit. The implementation cross-multiplies
every entry of (2) against the independent complete adjugate
construction before using it. It clears only a local zero unit; it
never divides the global physical tensor by the denominator in (2).

For \(n=26\), both geometries have \(e_0=1-\tau^{24}\).
The terminal geometry has \(h=0\), while the split geometry has
\(h=1\). In each case there are three augmented modes, hence 27
Weyl basis states. Occupation splits this matrix into blocks of sizes
\(1,3,6,7,6,3,1\). The auxiliary shifts and exact Gaussian-factor
division recover the 243 physical active zero orders. Infinity uses
only the actual two selected directions, with block sizes
\(1,2,3,2,1\).

All unselected Gaussian spectators and both paired physical
partitions remain in the pressure evaluator. The recorded relative
cost compares with the Gaussian at the same scale and aggregate
charge. Common Laurent clearing shifts the two local profiles by
opposite scalar orders and is charged in the polynomial degree.

## 2. A bounded complete off-diagonal family

For each occupation sector \(s=1,2,3\), take every rank-one
idempotent matrix \(P_s\) over \(\mathbb F_3\). Extend it by zero
on all other occupation sectors and set

\[
F(q)=I+\lambda(1-q^h)P_s,\qquad
\lambda\in\{1,2\},\quad
h\in\{-nP,-1,1,nP\}.
\tag{3}
\]

The index \(P\) in the exponent set is the ramification, whereas
\(P_s\) denotes the projector. All saved runs use ramification one.
Every kernel has \(F(1)=I\), so no endpoint repair is needed for
this search family. Each kernel preserves total occupation. Its
off-diagonal entries are retained in both cuts.

A rank-one idempotent in dimension \(m\) has the form \(uv^{\mathsf T}\)
with \(v^{\mathsf T}u=1\). Normalize the first nonzero coordinate of
\(u\) to one. There are \((3^m-1)/2\) choices for that projective
line and \(3^{m-1}\) covectors with pairing one. Thus the exact
counts are:

| Occupation | Sector dimension | Rank-one projectors |
|---|---:|---:|
| 1 | 2 | 12 |
| 2 | 3 | 117 |
| 3 | 2 | 12 |

There are \((12+117+12)\cdot2\cdot4=1,128\) nonidentity kernels.
Adding the Gaussian gives 1,129. The code verifies idempotence and
deduplicates the complete coefficient tables to check this count.

### Explicit global degree accounting

Every original one-particle kernel entry has degree at most \(2n\).
The truncated symmetric Gaussian uses at most \(2n(p-1)\) such
factors, so

\[
\deg G\le4n^2(p-1).
\tag{4}
\]

The existing operator construction gives the sufficient bound

\[
\deg T_F\le
P\bigl(4n^2(p-1)+2nr(p-1)\bigr)+D_F,
\tag{5}
\]

where \(D_F\) is the maximum coefficient degree after one common
Laurent clearing. For (3), that clearing has degree at most
\(nP\), and \(D_F\le nP\). At \(p=3,n=26,r=2,P=1\),
equation (5) is at most 5,642. The bound is conservative, but it is
explicit and includes every added coefficient power. Each stored
representative records its actual coefficient support, clearing,
and the corresponding bound (5).

## 3. Exact controls

The evaluator and original physical tensor constructor are
[`prime_deficit_identity_operators.py`](../tools/prime_deficit_identity_operators.py).
The independent
[`checker`](../tools/check_prime_deficit_identity_operators.py) has a
base replay and an extended replay. The saved
[`extended evidence`](../evidence/prime-deficit-identity-operators.json)
contains:

- Four complete rectangular-versus-Weyl controls, including two
  nontrivial 243-order qutrit operator profiles.
- Six original 81-by-81 physical cut profiles at \(n=2\), with
  ramification one or three and identity or filtered endpoint repair.
- 2,937 original endpoint coefficients and the ordinary polynomial
  degree assertions on those physical tensors.
- Two independent prime-power specializations and Gaussian recovery
  at both \(n=26\) geometries.
- Eighteen independent evaluations using the complete spectator and
  paired physical partitions, with maximum discrepancy below
  \(3.5\times10^{-13}\).
- A degree-601 source control that retains its complete infinity
  profile, and exact checks of the finite projector counts.

```bash
OPENBLAS_NUM_THREADS=1 python3 tools/check_prime_deficit_identity_operators.py --extended
```

## 4. Completed finite searches

The [driver](../tools/search_qutrit_deficit_operators.py) enumerates the
entire specified family. The
[`saved n26 evidence`](../evidence/qutrit-deficit-operator-projectors-n26.json)
retains every exact profile class, its multiplicity, representative
kernel, degree accounting, and numerical pressure result.

| Actual selected indices | Complete kernels | Exact profile classes | Best refined critical pressure |
|---|---:|---:|---:|
| \((24,25)\) | 1,129 | 101 | \(0.0723582513522645\) |
| \((0,25)\) | 1,129 | 207 | \(0.0723582513522599\) |

Each run uses a 72-point scale/charge grid and refines the best 30
profile classes numerically. The last-digit variation in the two
displayed minima is numerical. Relative values near
\(-10^{-14}\) are rounding noise and are not negative pressure
certificates. No global optimality over continuous parameters, other
Laurent supports, or other occupation operators follows from these
finite runs.
