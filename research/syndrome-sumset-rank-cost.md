# A sumset lower bound from the two syndrome block decompositions

An additive conservation law can by itself exclude a crossing family,
even when its individual coefficient formulas and output maps vary.
The relevant quantity is the number of distinct syndrome sums.

## Theorem

Let the logical labels be `1,...,d`, with syndromes `sigma_i` in a binary
vector space. Put `Sigma={sigma_i}` and `m=|Sigma+Sigma|`. Suppose every
allowed four-leg tensor satisfies

\[
T_{xy,vu}=0
\quad\text{unless}\quad
\sigma_x+\sigma_y+\sigma_v+\sigma_u=0.
\tag{1}
\]

Over any coefficient field, every finite decomposition of the crossing
`E_(xy,vu)=delta_xv delta_yu` into such tensors obeys

\[
\boxed{\sum_T\sqrt{r(T)c(T)}\ge m.}
\tag{2}
\]

For the full tensor-product support class at power `n`, the corresponding
bound is `m^n`. It applies to arbitrary coefficient recombinations and
mixtures of different constructions that obey the same syndrome support
law. No common polynomial parametrization is required.

**Proof.** In the row cut, group the row `(x,y)` by `sigma_x+sigma_y`
and the column `(v,u)` by `sigma_v+sigma_u`. Equation (1) makes the
matrix block diagonal. In the column cut, group `(x,u)` and `(v,y)` by
their corresponding syndrome sums; that matrix is also block diagonal.

For each `s` in `Sigma+Sigma`, choose labels `i_s,j_s` with syndrome
sum `s`, and define the common entry functional

\[
\chi_s(T)=T_{i_s,j_s,i_s,j_s}.
\]

This entry belongs to block `s` in both cuts. If it is nonzero, both
blocks have rank at least one. Thus, writing `q(T)` for the number of
nonzero `chi_s(T)`, we have `r(T)>=q(T)` and `c(T)>=q(T)`. It follows
that `sqrt(r(T)c(T))>=q(T)`.

For every `s`, the crossing satisfies `chi_s(E)=1`. In a decomposition
of `E`, at least one piece must therefore have nonzero `chi_s`.
Summing the preceding inequality over pieces proves (2).

At tensor power `n`, use the ordered tuple of the `n` syndrome sums as
the block label. There are exactly `m^n` such labels, with a common
crossing entry in each block. The same argument proves the tensor-power
claim. It also proves the mixture claim, since only (1) was used.
\(\square\)

The zero sum contributes one block as well. Omitting it gives a slightly
weaker valid lower bound, but it is included in (2).

## Sidon syndromes

Suppose the syndromes are distinct and all sums `sigma_i+sigma_j` with
`i<j` are distinct and nonzero. This is the pair-sum property needed
here; linear independence of the syndromes is sufficient but unnecessary.
Then

\[
m=1+\binom d2,
\qquad
S_n\ge\left(1+\binom d2\right)^n.
\tag{3}
\]

For `d>=6`, even `binom(d,2)>d^(3/2)>d^sqrt(2)`: the first inequality
is equivalent to `(d-1)^2>4d`, which holds at six and increases
thereafter. At `d=5`, the exact comparison `11^70>5^99`, together with
`99^2>2*70^2`, gives the same strict-seed exclusion. These comparisons
exclude every tensor power and every allowed mixture in the stated
support class.

## Application to cubic output maps

The cubic inputs `Omega e_i` have independent single-site syndromes
`sigma_i=e_i`. Any output map preserving those logical syndromes obeys
(1), including outputs that use both low and high single modes at site
`i`. In particular, the scalar and diagonal complementary images in
[the five-direction cubic construction](cubic-complementary-output-maps.md)
remain in this support class.

For its seven-label scalar-image family, equation (3) gives

\[
S_n\ge22^n>7^{n\sqrt2}
\]

over every field and for every tensor power. This does not contradict
the improved finite cost 36: it establishes that the remaining gap to
the strict-seed target cannot be removed within the same support class.

A non-diagonal complementary image `C=Pi Omega f`, or an output
correction with triple-single syndrome, need not preserve the logical
syndrome. Those maps are outside this theorem. Likewise, codes with many
repeated syndromes or repeated pair sums can have a much smaller value
of `|Sigma+Sigma|`; (2) must be evaluated for the actual syndrome set.
