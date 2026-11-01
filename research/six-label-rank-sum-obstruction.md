# Six labels cannot pass the rank-sum amplification criterion

## Result and scope

For every field and every exact decomposition of the six-label crossing

\[
E_6(x,y,u,v)=[x=v][y=u]=\sum_i A_i(x,y,u,v),
\]

let \(r_i,c_i\) be the ranks in the cuts
\((x,y)\mid(u,v)\) and \((x,u)\mid(y,v)\), respectively. Then

\[
\boxed{\quad \sum_i\sqrt{r_ic_i}>6^{\sqrt2}.\quad}       \tag{1}
\]

The [rank-sum amplification theorem](finite-seed-amplification.md)
would give an exponent beyond one half from a finite seed with
\(\sum_i\sqrt{r_ic_i}<d^{\sqrt2}\). Equation (1) shows that such a
seed needs **at least seven labels**. The earlier universal bound in
that note already excludes \(d\le5\).

This only limits that sufficient amplification criterion. It does not
exclude a six-label seed amplified through a different construction,
or a favorable finite maximum rank product.

## 1. Three or more pieces have too much rank-sum cost

Discard zero pieces, and write \(k\) for the number left. Their ranks
are positive integers at most \(36\), and rank subadditivity gives

\[
\sum_i r_i\ge36,\qquad \sum_i c_i\ge36.                \tag{2}
\]

For \(2\le k\le36\), the fixed-piece-count version of the elementary
concavity bound is

\[
\sum_i\sqrt{r_ic_i}
\ge2\sqrt{37-k}+k-2.                                 \tag{3}
\]

For completeness, first decrease the ranks as real numbers, keeping
each at least one, until each sum in (2) is exactly 36. This cannot
increase the objective. Minimizing a concave function of the row
ranks on their simplex moves their surplus to one coordinate; then
do the same with the column ranks. If the surplus coordinates agree,
the objective is 36. If they differ, it is the right side of (3).
The latter is the smaller value. This argument relaxes the integer
constraints, so its lower bound applies to the original ranks.

The expression in (3) is nondecreasing in \(k\). Thus for \(k\ge3\)
it is at least \(1+2\sqrt{34}\), unless \(k>36\), when the simpler
bound \(\sum_i\sqrt{r_ic_i}\ge k\) is enough. A single piece has cost
36.

The following integer checks give a rational separation:

\[
6^{\sqrt2}<\frac{1261}{100}<1+2\sqrt{34}.               \tag{4}
\]

Indeed \(\sqrt2<99/70\), and

\[
1261^{70}>6^{99}100^{70},\qquad
34\cdot200^2>1161^2.
\]

It remains to exclude a two-piece decomposition with small cost.

## 2. Small two-piece cost forces two rank-one cuts

Suppose there are two nonzero pieces and their total cost is below
13. The positive integer ranks lie in \(1,\ldots,36\) and obey

\[
r_1+r_2\ge36,\qquad c_1+c_2\ge36.
\]

The only rank patterns compatible with these inequalities and cost
below 13 are

\[
(r_1,c_1,r_2,c_2)=(1,a,b,1)
\quad\text{or}\quad (a,1,1,b),
\qquad a,b\in\{35,36\}.                              \tag{5}
\]

This tiny finite arithmetic claim is checked using integers alone in
[check_six_label_rank_sum.py](../tools/check_six_label_rank_sum.py).
For each admissible quadruple, it tests the necessary inequality

\[
\left\lfloor100\sqrt{r_1c_1}\right\rfloor+
\left\lfloor100\sqrt{r_2c_2}\right\rfloor<1300.
\]

Exactly the eight patterns in (5) pass even this weaker test. This
avoids comparisons of floating-point square roots. Since the threshold
in (4) is below 13, (5) applies to any putative violation of (1).

If a four-tensor has column rank one, it has the form
\(P(x,u)Q(y,v)\). Its row matrix is \(P\otimes Q\), so its row rank
is \(\operatorname{rank}P\operatorname{rank}Q\). Both factors have
rank at most six. Their product is at least 35 only when both ranks
are six; the product is then 36. The same argument with the cuts
exchanged applies to the other piece.

Consequently a hypothetical low-cost two-piece decomposition would
have one piece of ranks \((36,1)\), and the other of ranks \((1,36)\).

## 3. A swap cannot be a product plus one rank-one matrix

The preceding pattern would require invertible six-by-six matrices
\(P,Q\) such that

\[
\operatorname{rank}(\operatorname{Swap}-P\otimes Q)=1. \tag{6}
\]

More generally, for any dimension \(d\ge2\) and invertible \(P,Q\),

\[
\operatorname{rank}(\operatorname{Swap}-P\otimes Q)
\ge d-1.                                             \tag{7}
\]

Fix a nonzero vector \(u\), and restrict the input to \(u\otimes V\).
The swap has image \(V\otimes u\), whereas the product has image
\(Pu\otimes V\). These are both \(d\)-dimensional and intersect in
the one-dimensional space spanned by \(Pu\otimes u\).
Pass to the quotient by \(Pu\otimes V\). The product disappears,
and the swap still has rank \(d-1\). This proves (7).

At \(d=6\), it contradicts (6). Together with Sections 1--2, this
proves (1) over every field.

## Reproduction

The checker verifies the exact threshold inequalities and every
admissible integer rank quadruple. It also tests the swap/product
inequality over several small finite fields as a negative control;
the all-field proof is (7).

```sh
python tools/check_six_label_rank_sum.py \
  --output evidence/six-label-rank-sum-obstruction.json
```
