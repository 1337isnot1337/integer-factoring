# A nonbinary output section for the 24-label code

The full binary output-section scan can be extended by changing every
output row over \(\mathbb F_4\). This gives a coefficient space outside
the fixed 24-label span covered by the earlier all-power theorem.
The uniform nonbinary section below still fails the finite-seed
criterion over \(\mathbb F_4\).

Let \(\omega^2=\omega+1\). For the logical input
\(m(D_a+D_b)\), with \(a<b\), define

\[
\Pi_\omega=(1+\omega)\Pi_{\rm lower}+\omega\Pi_{\rm upper}.
\]

The two coefficients sum to one, so \(\Pi_\omega I=I_{24}\).
This changes all 24 output rows and retains the crossing endpoint.
Both unwanted endpoint conditions are inherited from the common
input code and the doubled-output parent tensor.

There are two output slots in the four-index tensor. Their scalar
products are

\[
(1+\omega)^2=\omega,\qquad
\omega(1+\omega)=1,\qquad
\omega^2=1+\omega.
\]

Consequently the new tensor is assembled exactly from the four
parent output choices as

\[
G_\omega=\omega G_{00}+G_{01}+G_{10}+(1+\omega)G_{11}.
\]

The checker reconstructs these matrices from the same literal
doubled-output parent used for the binary section search. It does
not interpolate from binary parameter values.

## Complete finite result

The new coefficient space has dimension six over \(\mathbb F_4\).
All \(4^6-1=4095\) nonzero vectors were checked, representing 1,365
projective classes. Fifteen affordable vectors, representing five
projective classes, span only four dimensions over \(\mathbb F_4\),
and this span misses the crossing.

As in the other finite gates, a hypothetical decomposition of cost
less than \(24^{\sqrt2}\) could use only vectors below the conservative
integer rank-product threshold. Since those vectors do not span the
target, the criterion fails throughout this complete coefficient
space over \(\mathbb F_4\).

The result applies to this uniform section and this field. It is not
an all-field or all-power statement about altered output sections.

## Restriction of scalars and independent audit

For a matrix \(A+\omega B\), the binary multiplication representation
is

\[
\begin{pmatrix}A&B\\B&A+B\end{pmatrix}.
\]

Its binary rank is twice the rank over \(\mathbb F_4\). Including each
coefficient matrix and its \(\omega\)-multiple therefore gives twelve
independent binary directions. Scanning that binary span is exactly
the full six-dimensional \(\mathbb F_4\) scan.

The native file uses doubled ranks and a doubled cost threshold. It
also includes explicit zero padding to meet the native square-size
input convention. The padding changes no rank. These bookkeeping
dimensions do not change the actual logical label count, which is 24.

Every one of the fifteen affordable rank pairs was independently
recomputed by direct \(\mathbb F_4\) elimination. The checker also
verifies closure of the affordable set under multiplication by
\(\omega\), the twelve-dimensional binary independence certificate,
and the exact crossing remainder.

Files:

- [check_uniform_f4_two_even_sections.py](../tools/check_uniform_f4_two_even_sections.py);
- [uniform-f4-two-even-sections.json](../evidence/uniform-f4-two-even-sections.json);
- [uniform-f4-two-even-sections-lift.jsonl](../evidence/uniform-f4-two-even-sections-lift.jsonl).

Reproduce after compiling the existing weighted native screen:

```sh
python tools/check_uniform_f4_two_even_sections.py \
  --native /tmp/weighted-screen \
  --screen-output /tmp/uniform-f4-section-lift.jsonl \
  --output /tmp/uniform-f4-section.json
```
