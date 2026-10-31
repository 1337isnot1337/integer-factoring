# Positive-energy quartics and removable endpoint poles

A quartic perturbation need not vanish locally at the crossing endpoint.
It is enough that every perturbed endpoint term strictly increases an
index sum that is fixed by the physical word type. The same argument also
permits certain rational coefficients whose poles cancel only after the
type restriction.

The two concrete consequences are:

- \(f(w)=w^{-2}\) is allowed for any quartic with strictly positive index
  increase.
- \(f(w)=1/[w(1+w)]\) is allowed when that increase exceeds the width of
  the mode interval.

These are endpoint and polynomiality statements. They do not establish a
better pressure bound by themselves.

## 1. Kernel, energy, and fixed types

Work over a field of characteristic two.

Use the one-species polynomial Gaussian kernel

\[
\mathcal D(w)=
\begin{pmatrix}
\chi QP^{\mathsf T}&w^{2M}P\\
I&\chi H
\end{pmatrix},
\qquad \chi=1+w,\quad z=w^2,
\tag{1}
\]

from [the one-species construction](schur-one-species-interval.md).
The finite mode set \(S\) lies in \([0,M)\), where \(M\) is a power of two.
The rows are the plus variables on physical boundaries \(V,Y\), and the
columns are the minus variables on \(U,X\).

Write \(F(w)=\Gamma(\mathcal D(w))\) for its Gaussian tensor. A coefficient
of \(F\) is a minor of \(\mathcal D\), hence a sum of matching terms.
For a physical monomial, define the integer energy

\[
E=
\sum_{\text{occupied modes in }Y,V} i
-
\sum_{\text{occupied modes in }X,U} i.
\tag{2}
\]

An edge matching a row mode \(i\) to a column mode \(j\) contributes
\(i-j\).

The two cross-edge families in (1) have nonnegative energy: the identity
edge has energy zero, and the Pascal edge is supported on \(j\subseteq i\),
so \(i\ge j\). Each of the two internal-edge families has a factor
\(\chi\). If

\[
L=\max S-\min S,
\]

an internal edge has energy at least \(-L\).

In a tensor power, prescribe the same occupation count for every mode on
all four physical boundaries. Every retained entry has total energy zero.
The argument below also works for any finer common word type.

## 2. Positive-energy endpoint lemma

Let

\[
Q_4=x_a y_b u_c v_d,\qquad
\Delta=b+d-a-c>0,
\tag{3}
\]

where the variables lie on their respective physical boundaries. Let
\(f(w)\) be Laurent polynomial in \(w\), regular at \(w=1\), and set

\[
T(w)=(1+f(w)Q_4)F(w).
\tag{4}
\]

**Lemma.** Every tensor power of (4), restricted to a common fixed-mode
type, has the same crossing endpoint as the unperturbed kernel.

**Proof.** At \(w=1\), the two internal edge families vanish because
they contain \(\chi\). Every surviving unperturbed matching therefore
has nonnegative energy. A term using \(k\ge1\) quartic factors has energy
at least \(k\Delta>0\). It cannot survive a restriction whose total energy
is zero.

The surviving terms use no quartic factor. In those terms every Pascal
matching must also have zero index increase. Its subset inclusions are
therefore equalities, so the endpoint is exactly the crossing. \(\square\)

In particular, \(f(w)=w^{-2}\) is valid even though \(f(1)=1\).

For \(d=M-1\), \(b=M-1-a\),

\[
\Delta=2b-c.
\tag{5}
\]

At \(M=4\), \(a=1,b=2,d=3\), all \(c=0,1,2,3\) give positive
energy. The localized choice \(a=b=c=0,d=3\) also has positive energy.

## 3. Rational cancellation after type restriction

The endpoint restriction can remove poles at \(w=1\) as well. Let

\[
f(w)=w^{-r}\chi^{-s}h(w),
\qquad r,s\in\mathbb Z_{\ge0},\quad h(w)\in F[w].
\tag{6}
\]

**Rational endpoint lemma.** If

\[
\Delta>sL,
\tag{7}
\]

the fixed-mode restriction of every tensor power of (4) is a Laurent
polynomial in \(w\), and its endpoint is exactly the crossing.

**Proof.** Expand a retained coefficient into quartic choices and Gaussian
matching terms. Consider a term with \(k\ge1\) quartic factors and \(j\)
internal Gaussian edges. Cross edges have nonnegative energy, each
internal edge has energy at least \(-L\), and the total energy is zero.
Thus

\[
0\ge k\Delta-jL,
\qquad
j\ge \left\lceil\frac{k\Delta}{L}\right\rceil.
\tag{8}
\]

Condition (7) implies \(j\ge sk+1\). The matching term contains
\(\chi^j\), while the quartic coefficients contribute only
\(\chi^{-sk}\). The denominator therefore cancels term by term, leaving
at least one factor \(\chi\). Every term using a quartic consequently
vanishes at \(w=1\).

The remaining denominator is a power of \(w\). Its order is at most
\(rm\) in the \(m\)-fold tensor power. The upper polynomial degree is also
linear in \(m\), since the mode set and local polynomials are fixed.
Hence the projected tensor is a Laurent polynomial with degree range
linear in \(m\). Its unperturbed endpoint is the crossing by the previous
argument. \(\square\)

The case \(L=0\) has no positive-energy quartic on \(S\), so the
nontrivial statement concerns \(L>0\).

For \(s=0\), this recovers the positive-energy lemma. For \(s=1\), the
coefficient

\[
f(w)=\frac1{w(1+w)}
\]

is valid whenever \(\Delta>L\). A concrete full-interval example is

\[
a=0,\qquad b=1,\qquad c=0,\qquad d=M-1,
\]

which has \(\Delta=M>M-1=L\).

The strict inequality matters. At \(\Delta=sL\), this proof only
guarantees cancellation of the denominators. Terms with exactly \(sk\)
internal edges can survive at the endpoint.

### Several quartic factors

The same proof applies to a finite product or sum of perturbations with
quartics \(Q_j\), energy increases \(\Delta_j\), and endpoint-pole orders
\(s_j\), provided every selected factor satisfies
\(\Delta_j>s_jL\). For any nonempty selection, the total energy excess
over the required denominator order is strictly positive. Nilpotent
products that vanish cause no difficulty.

## 4. Relation to the determinant calculation

For the localized quartics whose two-body contraction is
\(\mu=1+w^2=\chi^2\), the scalar factors in the column-cut determinant
satisfy

\[
\begin{aligned}
1+w^{-2}\mu&=w^{-2},\\
1+\frac{\mu}{w(1+w)}&=w^{-1},\\
1+(1+w^2)w^{-4}\mu&=w^{-4}.
\end{aligned} \tag{9}
\]

These identities explain why allowing nonvanishing endpoint coefficients
and removable endpoint poles enlarges the candidate family. They do not
control the full weighted Smith profile. A determinant gain can coexist
with a worse pressure bound.

For a rational local tensor satisfying the lemma, the zero and infinity
rank bounds can be applied to the projected Laurent polynomial after
clearing its bounded \(w\)-pole. Polynomiality here is a property of the
projected tensor powers, not of the unprojected local tensor.

## 5. Exact four-mode controls

The checker is
[check_positive_energy_quartic_endpoints.py](../tools/check_positive_energy_quartic_endpoints.py);
its output is
[positive-energy-quartic-endpoints.json](../evidence/positive-energy-quartic-endpoints.json).

It reconstructs every Gaussian coefficient as an exact polynomial minor
over \(\mathbb F_2\). For six positive-energy quartics, it checks all
\(16^4=65\,536\) local endpoint entries. It then checks all 81 fixed-mode
types in two copies, covering \(104\,976\) typed entries per quartic.

For \(Q_4=x_0y_1u_0v_3\) and \(f=1/[w(1+w)]\), it additionally checks
the denominator cancellation before summing terms: 1744 nonzero terms
with one quartic and six nonzero terms with two quartics. Each has the
extra factor \(1+w\) required to make its endpoint vanish. The complete
two-copy projected endpoint is also checked.

A zero-energy negative control, \(Q_4=x_0y_0u_0v_0\), cancels a required
diagonal endpoint entry when its coefficient equals one at \(w=1\).
Thus nonvanishing perturbations cannot be admitted without an endpoint
argument.

Reproduce:

~~~sh
python tools/check_positive_energy_quartic_endpoints.py \
  --output evidence/positive-energy-quartic-endpoints.json
~~~
