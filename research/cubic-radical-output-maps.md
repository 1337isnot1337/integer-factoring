# Cubic radical codes with syndrome-preserving output maps

For every odd number \(d\ge3\) of complementary mode pairs, the
one-single-mode cubic radical code has \(d\) logical labels. This note
allows every linear output functional supported on its own
one-single-low-mode syndrome.

The resulting complete coefficient spaces satisfy
\[
 \boxed{\quad
 \sum_\nu\sqrt{r(T_\nu)c(T_\nu)}
 \ \ge\ \frac{d(d-1)}{\sqrt2}
 \quad}
 \tag{1}
\]
for every exact crossing decomposition drawn from them, over every
field of characteristic two. For odd \(d\ge5\), this is at least
\(d^{3/2}>d^{\sqrt2}\). The three-label case is excluded by the
universal finite gate. Thus no member of this output-map family gives
a strict finite seed.

The scope is **syndrome-preserving outputs**. Arbitrary left inverses
that add other single-mode syndromes, high-mode orientations, or
three-single-mode components are not covered.

## 1. Inputs and allowed outputs

Let the complementary modes at pair \(i\) be \(e_i,f_i\), and write
\(D_i=e_if_i\). The input vector for logical label \(i\) is
\[
 I_i=e_i\sum_{a\ne i}D_a.
 \tag{2}
\]
It lies in exterior grade three. Since \(d\) is odd,
\[
 \Omega I_i=\Lambda I_i=0,
\]
so this is the \(L=3\) half-oriented radical family from
[the symplectic radical construction](noncommutative-crossing-symplectic-radicals.md).

Choose arbitrary weights \(w_{ia}\in F\) satisfying
\[
 w_{ii}=0,\qquad\sum_{a\ne i}w_{ia}=1.
 \tag{3}
\]
The \(i\)-th output functional is
\[
 \Pi_i=\sum_{a\ne i}w_{ia}\,[e_iD_a].
 \tag{4}
\]
Here \([e_iD_a]\) extracts that monomial coefficient. These are exactly
the linear outputs supported on the specified syndrome of label \(i\)
and satisfying \(\Pi_i I_i=1\).

Distinct labels have distinct single-mode syndromes. Consequently
\[
 \Pi I=I_d,\qquad I^{\mathsf T}I=0,\qquad\Pi JI=0.
\]
The second identity uses the even number \(d-1\) of monomials in each
input. The third uses the low-mode output orientation, because \(J\)
changes \(e_i\) to \(f_i\). The \(p^3,q^3\) vertices vanish, and the
\(r^3\) vertex is the crossing.

## 2. Five scalar polynomials and the pair values

Put \(u=q+r\), and define
\[
 \begin{aligned}
 H&=(p+r)u^2,\\
 K&=pu(p+u),\\
 A&=u(r^2+pu),\\
 B&=ru^2,\\
 C&=B+p^2u.
 \end{aligned}
 \tag{5}
\]
In particular \(K=H+C\), so these polynomials span at most four scalar
directions.

For every unordered pair \(\{i,j\}\), define
\[
 \begin{aligned}
 A_{ij}={}&B+(A+B)(w_{ij}+w_{ji})\\
 &+(C+B)\left(w_{ij}w_{ji}
                  +\sum_h w_{ih}w_{jh}\right).
 \end{aligned}
 \tag{6}
\]
This is symmetric in \(i,j\).

**Theorem 1.** The restricted tensor has the following exact entries:
\[
 T(x,y,u,v)=
 \begin{cases}
 H,&x=y=u=v,\\
 A_{ij},&(x,y,u,v)=(i,j,j,i),\quad i\ne j,\\
 K,&(x,y,u,v)=(i,j,i,j),\quad i\ne j,\\
 0,&\text{otherwise}.
 \end{cases}
 \tag{7}
\]
In the second and third lines, both orders of an unordered pair are
allowed.

### Proof of the entry formula

Local occupation parity forces
\[
 e_x+e_y+e_u+e_v=0
\]
in the syndrome space. Thus the only possible patterns are all equal,
or two pairs of equal indices.

First fix individual output monomials \(e_xD_a,e_uD_b\), before
summing them with the weights in (4). Expand the two input vectors (2)
literally. The local Pfaffian identities give:

| Single-index pattern | Contracted scalar before output weighting |
|---|---|
| \(x=y=u=v=i\) | \(H\), for every \(a,b\ne i\) |
| \(x=y=i,\ u=v=j,\ i\ne j\) | \(0\) |
| \(x=u=i,\ y=v=j,\ i\ne j\) | A symmetric function of \(a,b\), equal to \(K\) when \(a=b\) |
| \(x=v=i,\ y=u=j,\ i\ne j\) | The case distinction below |

For the last line, \(a\ne i,\ b\ne j\), and the scalar is
\[
 \begin{cases}
 C,&(a=j\text{ and }b=i)\text{ or }a=b,\\
 A,&\text{exactly one of }a=j,\ b=i\text{ holds},\\
 B,&\text{otherwise}.
 \end{cases}
 \tag{8}
\]

These are all-size polynomial identities. A direct finite verification
of the table suffices for the following reason. Mark the single-index
sites and the two fixed output-double sites; there are at most four.
Outside the marked sites, an input-double site must occur in both
inputs, or its local Pfaffian is zero. A common unmarked input-double
site contributes the same scalar regardless of its position. Summing
such choices depends only on the number of unmarked sites modulo two.
That parity is fixed by odd \(d\) and the marked-site coincidence
pattern. Hence every case reduces to a local perfect-matching
calculation on at most five sites. The exact checker enumerates all
these patterns, and also replays them at \(d=3,5,7,9\).

The all-equal contribution remains \(H\), by the two row-sum identities
in (3). In the third line of the table, the off-diagonal terms in
\(a,b\) cancel in pairs. The diagonal contribution is
\[
 K\sum_a w_{ia}^2
 =K\left(\sum_a w_{ia}\right)^2=K.
\]
This uses characteristic two but does not require binary weights.

Finally, (8) equals
\[
 B+(A+B)([a=j]+[b=i])
       +(C+B)([a=j][b=i]+[a=b]).
\]
Multiply by \(w_{ia}w_{jb}\) and sum. The row sums give exactly (6).
This proves (7). \(\square\)

## 3. The complete two-cut block matrices

Use row cut \((x,y)\mid(v,u)\) and column cut
\((x,u)\mid(v,y)\). Separate the diagonal pairs \((i,i)\) from the
two ordered versions of each unordered pair \(\{i,j\}\).

In the row cut, the complete matrix is the direct sum of
\[
 H I_d
 \quad\text{and}\quad
 \begin{pmatrix}A_{ij}&K\\K&A_{ij}\end{pmatrix}
 \quad(i<j).
 \tag{9}
\]
In the column cut, the complete matrix is the direct sum of
\[
 (H+K)I_d+KJ_d
 \quad\text{and}\quad
 A_{ij}I_2
 \quad(i<j),
 \tag{10}
\]
where \(J_d\) is the all-ones matrix.

Equations (9)--(10) hold after any scalar linear combination of
coefficient tensors. One simply applies that linear combination to
\(H,K,A_{ij}\). They therefore describe the entire coefficient space,
not only one polynomial specialization.

## 4. An all-field rank charge

Consider one coefficient-space tensor, and let \(N\) count the
unordered pairs for which its value \(A_{ij}\) is nonzero.

For each such pair, its row block in (9) has rank at least one, and its
column block in (10) has rank two. The blocks are disjoint, so
\[
 r(T)\ge N,\qquad c(T)\ge2N,\qquad
 \sqrt{r(T)c(T)}\ge\sqrt2\,N.
 \tag{11}
\]
No assumption about the ranks of the diagonal blocks is used.

At the crossing endpoint \(p=q=0,r=1\),
\[
 H=1,\qquad K=0,\qquad A=B=C=1,
\]
and hence every \(A_{ij}=1\). In any exact decomposition of the crossing,
each unordered pair must therefore have a nonzero \(A_{ij}\) in at
least one piece. Summing (11) gives
\[
 \sum_\nu\sqrt{r(T_\nu)c(T_\nu)}
 \ge\sqrt2\binom d2
 =\frac{d(d-1)}{\sqrt2}.
\]
This proves (1).

For \(d\ge5\),
\[
 (d-1)^2\ge2d,
\]
so (1) is at least \(d^{3/2}\). Since \(3/2>\sqrt2\), every odd
\(d\ge5\) misses the strict finite-seed target. For \(d=3\), use
\(S\ge2\sqrt{d^2-1}>d^{\sqrt2}\) from
[the universal finite gate](noncommutative-crossing-finite-gates.md).

The charging argument depends only on the block shapes. Thus it also
allows a decomposition whose pieces come from different choices of
the syndrome-preserving weight matrix, with the same logical labels
and input family.

## 5. Verification and precise scope

Run
\[
 \texttt{python tools/check\_cubic\_radical\_output\_maps.py}.
\]
The checker makes 600 exact polynomial comparisons covering the
marked-site patterns at \(d=3,5,7,9\). It independently constructs dense
matrices for five output maps and verifies every binary basis
combination against the compressed block-rank computation.

The polynomial encoding sets \(r=1\) and retains independent formal
\(p,q\). Because every entry is homogeneous of degree three, this
retains every coefficient of the full \(p,q,r\) polynomial. The
checks are not numerical evaluations at one field point.

Sources and evidence:

- [Polynomial and dense-rank checker](../tools/check_cubic_radical_output_maps.py).
- [Exact audit evidence](../evidence/cubic-radical-output-maps.json).

The ordinary proof extends from binary checker data to arbitrary
characteristic-two weights by the polynomial identities (5)--(8), the
row sums in (3), and the rank argument (11). It covers all scalar
coefficient recombinations and mixtures of the stated output maps.

It does not cover arbitrary linear left inverses. A left inverse may
contain components supported on another logical syndrome, on high-mode
orientations, or on three-single-mode monomials, without changing
\(\Pi I=I_d\). Those components can destroy the block shapes
(9)--(10), and require separate analysis.

## 6. Exact exploratory profiles

The computations that led to the block formula are retained as
reproducible evidence. The cyclic selector chooses `w[i,i+1]=1`, with
indices taken modulo `d`. Its diagonal cut blocks are circulant, so
their ranks can be computed exactly by polynomial greatest common
divisors over the binary field. The remaining blocks have size two.

The recorded cyclic cases are `d=3,5,7,9,15,31,63,127,255,511,1023`.
Every coefficient combination is included for each listed case.
For example, the integer cost lower bounds at `d=31,127,1023` are
`699,11486,740668`, respectively. None passes the seed threshold.

A second screen permits every binary circulant output row with odd
weight. It covers all 32 selectors at `d=7`, all 128 at `d=9`, and all
512 at `d=11`. At `d=15,31`, it records 256 sampled selectors each,
with deterministic seed `20261030`. The latter two selector searches
are not exhaustive. Each recorded selector still receives a complete
binary coefficient-space calculation. The all-field theorem above
now covers these selectors and the other permitted output weights.

```sh
python tools/explore_cyclic_radical_cubic.py \
  --output /tmp/cyclic-radical-cubic.json
python tools/explore_circulant_cubic_outputs.py --pairs 7,9,11 \
  --output /tmp/circulant-cubic-outputs-small.json
python tools/explore_circulant_cubic_outputs.py --pairs 15,31 --limit 256 \
  --output /tmp/circulant-cubic-outputs-large.json
```

The [cyclic rank evaluator](../tools/explore_cyclic_radical_cubic.py)
and [circulant output screen](../tools/explore_circulant_cubic_outputs.py)
produce the [cyclic profiles](../evidence/cyclic-radical-cubic.json),
[complete small selector profiles](../evidence/circulant-cubic-outputs-small.json),
and [sampled larger selector profiles](../evidence/circulant-cubic-outputs-large.json).

