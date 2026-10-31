# Complete finite Laurent resonance catalogues

## Result and limits

The [four-mode resonance reduction](binary-four-mode-resonance.md) gives
exact two-cut profiles for an arbitrary two-by-two Laurent block over
\(\mathbb F_4\). The finite catalogues below include every invertible
endpoint, both off-diagonal entries, and every coefficient in two
specified three-power windows. They therefore test branches beyond the
polynomial-block obstruction in the preceding note.

At \(n=32\), \(P=1\), all 2,360 resonant index quadruples were included.
The complete coefficient enumeration and each profile signature are
exact. Pressure evaluation uses floating-point grids and local
refinements; it is search evidence. These runs found no improvement.
They do not prove a general Laurent obstruction or an improved
factoring exponent.

| Coefficient family for the resonant block \(K\) | Valid kernels | Exact signatures | Kernel-selection pairs |
|---|---:|---:|---:|
| \(K=E+(1+q^{-1})N\), \(E\in\mathrm{GL}_2(\mathbb F_4)\), arbitrary \(N\) | 46,080 | 105 | 108,748,800 |
| Every entry supported on \(q^{-2},q^{-1},1\) | 11,796,480 | 556 | 27,839,692,800 |
| Every entry supported on \(q^{-1},1,q\) | 11,796,480 | 575 | 27,839,692,800 |

The two large windows overlap; their counts must not be added as a
count of distinct kernels. The best refined critical pressure in all
three runs was the Gaussian value, approximately
\(0.0330679625084472\). Every catalogue was evaluated on 84 scale/charge
points for every resonant selection. The best 30 signature records
received numerical scale/charge refinements at their best grid selection.

## 1. Exact grouping of all coefficient tables

Write

\[
K=\begin{pmatrix}A&C\\D&B\end{pmatrix},\qquad
\delta_A=A-1,\quad\delta_B=B-1,\quad m=nP.
\]

An entry with support in three prescribed powers has \(4^3=64\)
possibilities. Thus there are \(64^4=16,777,216\) raw matrices.
Evaluation at one has \(4^2\) preimages for each entry value, and
\(|\mathrm{GL}_2(\mathbb F_4)|=180\). Exactly

\[
180\,4^8=11,796,480
\]

of the raw tables have an admissible endpoint.

For fixed diagonal entries, the off-diagonal pair affects both profiles
only through

\[
\left(CD,\ \min(\operatorname{ord}C,\operatorname{ord}D),\
\max(\deg C,\deg D)\right).
\tag{1}
\]

Zero polynomials use the usual infinite valuation and negative infinite
degree conventions. All 4,096 off-diagonal pairs split into exactly
607 classes under (1). The implementation retains each class's exact
multiplicity and one literal representative. This reduces the complete
enumeration to 4,096 diagonal pairs times 607 off-diagonal classes,
without omitting coefficient cancellation.

The sufficient seven-integer profile signature is

\[
(v_0,v_A,v_B,\lambda_1,\lambda_2,\eta_1,\eta_2),
\tag{2}
\]

where the first three entries are the vacuum and two singleton
\(q\)-orders, \(\lambda\) are the two \(q\)-orders of
\(R=(q^{2m}-1)I+K\), and \(\eta\) are the two \(z\)-orders
of \(K(z^{-2})\). The selected indices enter only through the explicit
physical row/column weights; they do not change (2).

For either listed window and \(m\ge4\), the low coefficient bands,
the \(q^{2m}\) band, and the \(q^{4m}\) band in

\[
\det R=q^{4m}+q^{2m}(\delta_A+\delta_B)
             +\delta_A\delta_B+CD
\tag{3}
\]

are disjoint. Thus its order is the order of the first nonzero band.
The first Smith order is the minimum entry order, and the second is
the determinant order minus the first. The corresponding infinity sum
is \(-2\deg\det K\). The scalar formulas in the resonance note give
the remaining three entries. All of these computations use complete
finite-field polynomials or exact finite coefficient tables.

## 2. Reproducible evidence

The grouping implementation is
[`resonant_laurent_catalogue.py`](../tools/resonant_laurent_catalogue.py).
The independent
[`checker`](../tools/check_resonant_laurent_catalogue.py) compares grouped
counts with 131,072 literal matrix endpoint decisions, including 91,136
admissible literal profiles. It then reconstructs and checks every
representative of both complete catalogues. The saved
[`control record`](../evidence/resonant-laurent-catalogue-controls.json)
contains these exact counts.

The pressure driver is
[`search_resonant_identity_operators.py`](../tools/search_resonant_identity_operators.py).
All three complete run records, including every signature, multiplicity,
representative, resonant index selection, grid minimum, and refinement,
are retained in
[`resonant-laurent-catalogues-n32.json`](../evidence/resonant-laurent-catalogues-n32.json).

```bash
python3 tools/check_resonant_laurent_catalogue.py
python3 tools/search_resonant_identity_operators.py --n 32 --window-lower -2 --output /tmp/resonance-minus2.json
python3 tools/search_resonant_identity_operators.py --n 32 --window-lower -1 --output /tmp/resonance-minus1.json
```

The original physical endpoint and both complete physical cut profiles
of the resonance construction were checked separately in the
[four-mode reduction](binary-four-mode-resonance.md). Every pressure
comparison retains all Gaussian spectators and both paired physical
partitions.
