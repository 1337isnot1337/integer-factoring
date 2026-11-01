# The original constant-weight three-label restriction

## Result and scope

Consider the twenty-piece three-label certificate from
[PR #24](https://github.com/1337isnot1337/integer-factoring/pull/24),
at revision
[e96ee50560284301428d72421330d8773155d25b](https://github.com/1337isnot1337/integer-factoring/tree/e96ee50560284301428d72421330d8773155d25b).
The exact matrices and their common flags are recorded in
[finite-seed-amplification.md](finite-seed-amplification.md).
Number its pieces \(A_0,\ldots,A_{19}\) in the pinned witness order.

Take \(t\) tensor copies, group the copies by their twenty piece counts,
and restrict each external word to have exactly \(m\) occurrences of
\(w=e_2\). Its other coordinates lie in
\(U=\operatorname{span}(e_0,e_1)\). The label dimension is

\[
D_{t,m}=\binom tm2^{t-m}.
\tag{1}
\]

For every fixed density \(m/t\longrightarrow\theta<1\), some literal
piece-count coefficient has ranks \(r,c\) satisfying

\[
\log_2 r,\ \log_2 c
\ \ge\ \frac{36}{25}\log_2 D_{t,m}-O(\log(t+2)).
\tag{2}
\]

Thus this restriction cannot have a joint-width coefficient
\(\rho<18/25=0.72>1/\sqrt2\). For \(\theta\ge1/3\), the proof gives
the stronger coefficient \(3/4\). The exact finite bound below is
valid without a limiting-density assumption.

This is a statement about the **original basis**, with the line
\(w=\mathbb F_2e_2\), and its literal twenty-variable type coefficients.
Changing the one-label basis or recombining different types can change
the problem. In particular the isotropic line
\(\mathbb F_2(e_0+e_1)\) is outside this theorem.
There is no lower bound on arbitrary crossing constructions or factoring.

The separate
[permutation-invariant defect theorem](fixed-defect-invariant-barrier.md)
applies to this family when \((t-m)/t\) tends to zero. It gives
\(\rho\ge1-o(1)\) there. Together with the uniform version of (2) on
every closed interval \(\theta\le1-\varepsilon\), it also excludes a
uniform \(\rho<18/25\) along any sequence of these original-basis
restrictions with \(D_{t,m}\) tending to infinity.

## 1. A degeneration that preserves every position pattern

For the row cut, use the matrix with rows \((x,y)\) and columns
\((u,v)\), and multiply it on the right by the crossing permutation.
Its input pair is then \((v,u)\).
For the column cut, use rows \((x,u)\), columns \((y,v)\), and do the
same normalization. Its input pair becomes \((v,y)\).
The normalization preserves both the unrestricted rank and the
constant-weight restriction.

The following ordered bases, written as nine-bit integers, are the
common flags from the finite certificate:

| Cut | Pair-space basis | Height vector |
| --- | --- | --- |
| Row | \(27,36,192,10,256,3,4,64,1\) | \(0,1,2,3,4,6,7,8,9\) |
| Column | \(27,192,36,10,3,256,64,4,1\) | \(0,2,3,4,8,6,10,11,12\) |

Every basis vector is homogeneous for both individual indicators
\([\,\text{first coordinate}=w\,]\) and
\([\,\text{second coordinate}=w\,]\).
Consequently these basis changes, their tensor powers, and arbitrary
diagonal rescalings in these bases preserve every position-by-position
pattern of occurrences of \(w\). This is stronger than preserving
only the total count.

In either basis let \(h_j\) be the indicated height. Conjugate an
operator by the diagonal matrix with entries \(z^{-h_j}\).
Its entry in row \(r\), column \(c\) is multiplied by \(z^{h_c-h_r}\).
Divide each of the following six matrices by its minimum valuation.
Their constant terms have exactly these nonzero entries; indices start
at zero.

| Piece | Row valuation | Row initial entries | Column valuation | Column initial entries |
| --- | ---: | --- | ---: | --- |
| \(4\) | \(1\) | \((2,3),(7,8)\) | \(2\) | \((1,3),(6,8)\) |
| \(5\) | \(2\) | \((0,2),(5,7)\) | \(2\) | \((0,1),(4,6)\) |
| \(12\) | \(0\) | \((4,4)\) | \(0\) | \((5,5)\) |
| \(15\) | \(3\) | \((3,5)\) | \(0\) | \((0,0),(3,3),(4,4),(8,8)\) |
| \(16\) | \(0\) | \((0,0),(3,3),(5,5),(8,8)\) | \(4\) | \((3,4)\) |
| \(17\) | \(0\) | \((1,1),(6,6)\) | \(0\) | \((2,2),(7,7)\) |

All other nonzero entries have strictly larger valuation. These are
exact comparisons of small integers, checked directly from the pinned
twenty matrices.

Write \(\overline A_i\) for these initial matrices, separately in each
cut. Their pair-position patterns are the same in the two cuts:

| Piece | Output pattern | Input pattern | Rank |
| --- | --- | --- | ---: |
| \(4\) | \(wU\) | \(UU\) | \(2\) |
| \(5\) | \(UU\) | \(wU\) | \(2\) |
| \(12\) | \(ww\) | \(ww\) | \(1\) |
| \(17\) | \(Uw\) | \(Uw\) | \(2\) |

Pieces \(15,16\) act only on \(UU\). In the row cut they are,
respectively, a rank-one square-zero operator \(N\) and the identity
of order four. In the column cut these roles are reversed.

For a fixed piece-count vector, its total extracted valuation is
independent of the order of the tensor positions. The initial matrix
of that type coefficient is therefore exactly the corresponding type
coefficient of the six initial matrices. All conjugations preserve the
restricted pair space. Rank over \(\mathbb F_2(z)\) equals the original
rank, and specialization at \(z=0\) cannot increase rank. It follows
that any lower bound obtained from these initial matrices is a lower
bound for the original restricted type.

This step is necessary. In the original matrices, piece \(15\) has
column-cut entries between \(ww\) and \(UU\). Those entries prevent an
unqualified direct-sum argument over common \(w\)-position sets.
The displayed degeneration removes that issue.

## 2. An exact core and a binary remainder

Choose integers \(q,k,b\ge0\) with

\[
k=m-q,\qquad b=t-m-2q.
\tag{3}
\]

Consider the literal type

\[
12^k\,4^q\,5^q\,17^q\,15^i\,16^{b-i},
\qquad 0\le i\le b.
\tag{4}
\]

After the degeneration, piece \(12\) is the projector onto \(ww\).
It is the only piece that can have \(ww\) on either side.
Its positions form a common set \(K\) of size \(k\) on all four legs.
Different sets \(K\) give disjoint row and column blocks in both cuts.

For each \(K\), choose any fixed set \(J_K\) of \(3q\) positions outside
\(K\). Restrict all remaining occurrences of \(w\) to \(J_K\).
Every copy of \(4,5,17\) has an occurrence of \(w\), while \(15,16\)
have none. Hence the \(3q\) positions of \(J_K\) contain exactly the
copies of \(4,5,17\), and all other positions contain \(15,16\).

The three core pieces have distinct output patterns and distinct input
patterns. The \(q,q,q\) allocations therefore form a direct sum of
\((3q)!/(q!)^3\) blocks. Each block is a tensor product of \(3q\)
rank-two maps, and has rank \(8^q\), in either cut.
The remaining \(b\) positions give the binary coefficient of
\((N+zI)^{\otimes b}\), with its degree reversed in the other cut.

For clarity, if \(R_j(b)\) is the binary rank from
[binary-coefficient-entropy.md](binary-coefficient-entropy.md), then
the coefficient containing \(i\) copies of \(N\) has rank
\(R_{b-i}(b)\). The other cut has rank \(R_i(b)\). The equality follows
from the same square-zero elementary-symmetric operator calculation:
the extra factor appearing for \(N+z(I+N)\) is an invertible multiplier.

Taking \(i=\lfloor b/2\rfloor\), the exact binary bound gives

\[
R_i(b),R_{b-i}(b)
\ \ge\ \frac{\binom{2b}{\lfloor b/2\rfloor}}{b+1}.
\tag{5}
\]

Combining the common-set blocks, the core blocks, and (5), we obtain
the promised finite bound:

\[
\boxed{\quad
r,c\ \ge\
L(t,m,q):=
\binom{t}{m-q}
\frac{(3q)!}{(q!)^3}\,8^q\,
\frac{\binom{2(t-m-2q)}{\lfloor(t-m-2q)/2\rfloor}}
     {t-m-2q+1}.
\quad}
\tag{6}
\]

This holds for every integer
\(0\le q\le\min(m,\lfloor(t-m)/2\rfloor)\).
The right side may be nonintegral; ranks satisfy the stated real
inequality, and can of course be rounded up.

In particular, at \(t=3q,m=q\), the type \(4^q5^q17^q\) already has
both ranks exactly

\[
\frac{(3q)!}{(q!)^3}8^q.
\tag{7}
\]

For the equality in this special case, one can also read the original
matrices directly. With \(v=e_0+e_1\) and \(J\) the swap on \(U\),
they are
\(A_4=w_xv_vJ_{yu}\),
\(A_5=v_xw_vJ_{yu}\), and
\(A_{17}=w_u(w_yI_{xv}+J_{xy}v_v)\).
The \(x,v,u\) occurrences of \(w\) force the three disjoint position
sets. The required \(y\)-count then permits only the first term of
\(A_{17}\). This gives exactly the blocks counted above.

Its first example has ranks \(48,48\). The rate in (7) gives
\(\rho\ge\log 6/(2\log3)=0.81546\ldots\).

## 3. A uniform exact entropy certificate below density one third

Put \(H=H_2(1/4)=2-\tfrac34\log_2 3\), where
\(H_2(u)=-u\log_2u-(1-u)\log_2(1-u)\).
Let \(\theta=m/t\) and \(\eta=q/t\).
Elementary binomial estimates applied to (6) give

\[
\frac1t\log_2 L
\ \ge\
H_2(\theta-\eta)+3\eta\log_2 6
+2H(1-\theta-2\eta)-O\!\left(\frac{\log(t+2)}t\right).
\tag{8}
\]

Also

\[
\frac1t\log_2D_{t,m}
=d(\theta)+O\!\left(\frac{\log(t+2)}t\right),
\qquad d(\theta)=H_2(\theta)+1-\theta.
\tag{9}
\]

These errors are uniform; changing any count by a bounded integer
changes the logarithmic estimates by \(O(\log(t+2))\).

Set

\[
c_0=\frac{32}{761},\qquad
A=4+\log_2 761-\frac{15}{2}\log_2 3,\qquad
B=-9+\frac{15}{2}\log_2 3.
\tag{10}
\]

For \(c_0\le\theta\le1/3\), choose
\(k=m-q\) to be a nearest feasible integer to \(c_0t\).
Then the main term of (8) is

\[
E(\theta)=A+B\theta.
\tag{11}
\]

Indeed the optimum common-set density is characterized by
\((1-c_0)/c_0=729/32\).

Concavity of \(d\), using its tangent at \(\theta=1/9\), gives

\[
d(\theta)\le\log_2(9/4)+2\theta.
\tag{12}
\]

The following two integer inequalities are exact certificates:

\[
3^{375}>2^{594},\qquad
2^{242776}761^{38050}>3^{382959}.
\tag{13}
\]

The first is equivalent to \(B-72/25>0\).
The second is equivalent to

\[
E(c_0)-\frac{36}{25}
       \bigl(\log_2(9/4)+2c_0\bigr)>0.
\]

Thus (11)--(13) prove
\(E(\theta)\ge(36/25)d(\theta)\) throughout this interval.
No numerical minimization or rounded logarithm is used in this
certificate.

For \(0\le\theta<c_0\), take \(q=0\). The main rank-to-label ratio is

\[
\frac{H_2(\theta)+2H(1-\theta)}{2d(\theta)}
=\frac12+
\frac{2H-1}{2\bigl(1+H_2(\theta)/(1-\theta)\bigr)}.
\tag{14}
\]

It decreases with \(\theta\), since
\(H_2(\theta)/(1-\theta)\) increases. At \(c_0\), it agrees with
the value certified above, so it is at least \(18/25\) throughout.
This proves (2) up to density \(1/3\).

## 4. Densities at least one third

For \(\theta\ge1/3\), take \(q=\lfloor(t-m)/2\rfloor\). Then
\(b\) is zero or one. Put \(a=(1-\theta)/2\), so
\(0<a\le1/3\).
In natural logarithms, the leading rank and label exponents are

\[
\begin{aligned}
R(a)&=-3a\log a-(1-3a)\log(1-3a)+3a\log2,\\
D(a)&=-2a\log a-(1-2a)\log(1-2a).
\end{aligned}
\tag{15}
\]

Their difference \(f(a)=R(a)-\tfrac32D(a)\) is

\[
f(a)=-(1-3a)\log(1-3a)
+\frac32(1-2a)\log(1-2a)+3a\log2.
\]

Here \(f(0)=0\), \(f(1/3)=\log(2/\sqrt3)>0\), and

\[
f'(a)=3\log\frac{2(1-3a)}{1-2a}.
\]

The derivative is positive before \(a=1/4\) and negative afterward.
Hence \(f(a)\ge0\) on the entire interval.
This gives
\(\log_2 r,\log_2 c\ge\tfrac32\log_2D_{t,m}-O(\log(t+2))\),
and therefore \(\rho\ge3/4\) whenever the logarithmic error is
negligible. In particular this holds at every fixed density below one.

For densities approaching one, use the separate invariant theorem
mentioned at the start. Its hypotheses hold because every literal
type is invariant under simultaneously permuting the tensor positions.
This final step is independent of the six-piece degeneration.

## 5. Reproducible verification

Run:

    python3 tools/check_balanced_three_label_barrier.py \
      --output evidence/balanced-three-label-barrier.json

The checker uses the pinned matrices from the finite-seed checker.
It verifies all \(972\) matrix cells needed for the two six-piece
degenerations, the individual \(w\)-pattern homogeneity, every initial
support pattern, and the rank-one square-zero binary operators.
It then constructs six restricted initial type matrices and checks
their direct-sum lower bounds, including nontrivial common-set,
core, and binary-remainder cases. It checks \(97\) exact middle binary
rank bounds and both integer inequalities in (13).

The finite computations support the displayed local identities.
The all-size statement follows from the degeneration, block
decomposition, and entropy proofs above; it is not an extrapolation
from those finite matrices.
