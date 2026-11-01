# Quadratic charges, residue syndromes, and label entropy

## Result and scope

In the characteristic-two regular Fock space, the joint eigenspaces of
all quadratic zero modes of parameter weight at most a power of two
have an exact description. Their characters record the parity of
occupation in opposite odd-residue pairs.

If the charge cutoff \(M\) grows faster than \(\sqrt h\), every such
joint eigenspace at state energy \(h\) has dimension
\(\exp(o(\sqrt h))\). Thus a construction that retains
\(\exp(c\sqrt h)\) state labels cannot impose a common eigencharacter
for **all** these charges with \(M/\sqrt h\to\infty\).
This restricts that particular strategy. It is not a barrier for
selected charges, parameter-dependent cancellation, other crossing
constructions, or integer factoring algorithms.

The result also supplies explicit homogeneous vacuum duals built
from complementary pairs. The more general fixed-length-shell
construction is proved in
[complement-shell-crossings.md](complement-shell-crossings.md).
Neither result establishes a factoring exponent above \(1/2\).

The underlying operators are those of Anthropic's Appendix B at
[16a5ce41304fcdaea4e0be7c132e1c361dc26da5](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
Their characteristic-two quotient and normalization are recorded in
[modular-crossing.md](modular-crossing.md). These calculations and
proofs used OpenAI assistant support. No new Lean formalization is
asserted.

## 1. Every quadratic zero mode is diagonal

Work in

\[
W=\mathbf F_2[t_1,t_3,\ldots]/(t_1^2,t_3^2,\ldots),
\qquad \deg t_k=k.
\]

For a finite set \(S\) of positive odd integers write
\(t^S=\prod_{k\in S}t_k\). A quadratic parameter is
\(t_{2a+1}t_{2b+1}\), where \(a\ne b\). Its zero mode is diagonal:

\[
\bigl(t_{2a+1}t_{2b+1}\bigr)[0]t^S
=\left(\sum_{2c+1\in S}d_{a,b}(c)\right)t^S,
\tag{1}
\]

where the sum is in \(\mathbf F_2\) and

\[
d_{a,b}(c)
=\binom ca\binom{c+b}b+\binom cb\binom{c+a}a
\pmod2.
\tag{2}
\]

A quadratic normal-order term of shift zero must contain one
creation and one annihilation of the same odd mode. The two orders
of the parameter parts give the two summands in (2). Terms with two
creations or two annihilations have nonzero energy shift.
Lucas's binary criterion rewrites (2) as

\[
d_{a,b}(c)
=[a\subseteq c,\ b\mathbin{\&}c=0]
 +[b\subseteq c,\ a\mathbin{\&}c=0].
\tag{3}
\]

Here \(a\subseteq c\) means containment of binary set bits.
If \(a\mathbin{\&}b\ne0\), both summands vanish. If
\(a\mathbin{\&}b=0\) and \(a\ne b\), putting \(A=a\lor b\) gives

\[
d_{a,b}(c)=[c\mathbin{\&}A\in\{a,b\}].
\tag{4}
\]

For completeness, the expression with \(a=b\) is zero, consistent
with the square-zero quotient.

## 2. The full low-charge character is a residue syndrome

Fix a power of two \(M=2^{r+1}\ge4\). For each odd
\(1\le v<M/2\), define

\[
\sigma_v(S)
=\#\{k\in S:k\equiv v\ \text{or }-v\pmod M\}\pmod2.
\tag{5}
\]

There are \(M/4\) coordinates. Let
\(\sigma_M(S)=(\sigma_v(S))_v\).

**Theorem 1.** Two regular monomials have the same eigenvalues for
every quadratic zero mode of parameter weight at most \(M\) if
and only if they have the same syndrome \(\sigma_M\).

**Proof.** The parameter weight is \(2(a+b+1)\). If it is at most
\(M\), then \(a,b<2^r\). Formula (3) depends only on the bottom
\(r\) bits of \(c\). Formula (4) is invariant under complementing
all those bits: this exchanges its two accepted patterns \(a,b\).
In mode coordinates \(k=2c+1\), that complement is
\(k\mapsto M-k\pmod M\). Every charge eigenvalue is therefore a
linear combination of the coordinates (5).

Conversely, fix an odd \(v<M/2\) and take

\[
a=(v-1)/2,\qquad b=(M-v-1)/2.
\]

Their binary masks are complementary within \(r\) bits, so
\(a\lor b=2^r-1\). The parameter has weight exactly \(M\), and
(4) is the indicator of the residue pair \(\{v,M-v\}\).
Its eigenvalue is exactly \(\sigma_v(S)\). This proves both
directions. \(\square\)

Some syndromes may be absent at a fixed energy. For example, the
sum of their coordinates equals \(|S|\bmod2=h\bmod2\).
The theorem does not assert that every character occurs.

## 3. A common character forces a short symmetric difference

Let \(V_h[\chi]\) be the span of weight-\(h\) monomials with
syndrome \(\chi\). Suppose it is nonzero and choose one such
monomial \(S_0\). For every other \(S\) in this character space,
the symmetric difference \(T=S\mathbin{\triangle}S_0\) has
syndrome zero and total weight at most \(2h\).

**Pairing lemma.** If a set \(T\) of distinct positive odd integers
has syndrome zero, its elements can be partitioned into pairs,
each with sum at least \(M\).

**Proof.** In each opposite residue pair \(\{v,M-v\}\), the
counts in the two individual residues have the same parity. If
both are odd, match one element from each residue. Their sum is
a positive multiple of \(M\). The remaining counts are even;
match within each individual residue. Distinct positive
integers in one residue differ by at least \(M\), so each such
pair also has sum at least \(M\). \(\square\)

The lemma gives

\[
|T|\le \frac{4h}{M}.
\tag{6}
\]

The map \(S\mapsto T\) is injective. Every part of \(T\) is at
most \(h\), so with \(n=\lceil h/2\rceil\) and
\(K=\lfloor4h/M\rfloor\) we obtain the finite bound

\[
\dim V_h[\chi]
\le \sum_{\ell=0}^{\min(n,K)}\binom n\ell.
\tag{7}
\]

The zero character satisfies the stronger length bound
\(|S|\le2h/M\), by applying the lemma directly to \(S\).

For the asymptotic conclusion, use also the total-weight bound.
The number of sets of \(\ell\) distinct positive integers with
sum at most \(2h\) is at most
\(\binom{2h}{\ell}/\ell!\). Each such set has \(\ell!\)
orderings, and there are \(\binom{2h}{\ell}\) ordered positive
\(\ell\)-tuples with sum at most \(2h\). Dropping the odd-part
restriction can only enlarge the count. Therefore

\[
\dim V_h[\chi]\le
\sum_{\ell=0}^{K}\frac{\binom{2h}{\ell}}{\ell!}.
\tag{8}
\]

Terms with \(\ell>2h\) are zero. For \(\ell\ge1\),

\[
\frac{\binom{2h}{\ell}}{\ell!}
\le \left(\frac{2e^2h}{\ell^2}\right)^\ell.
\tag{9}
\]

**Theorem 2.** If \(M/\sqrt h\to\infty\), then

\[
\log\max_\chi\dim V_h[\chi]=o(\sqrt h),
\tag{10}
\]

where the maximum is over nonempty character spaces.
If \(M\ge\delta h\) for a fixed \(\delta>0\), every such
dimension is polynomial in \(h\).

**Proof.** Now \(K=o(\sqrt h)\). If \(K=0\), (8) is one.
Otherwise the bound in (9) increases for \(1\le\ell\le K\)
once \(h\) is sufficiently large, and

\[
\log\dim V_h[\chi]
\le \log(K+1)+K\log(2e^2h/K^2)=o(\sqrt h).
\]

The last step uses \(x\log(1/x)\to0\) as \(x\downarrow0\),
with \(x=K/\sqrt h\). If \(M\ge\delta h\), then \(K\) is
bounded, and (7) gives the polynomial assertion. \(\square\)

For fixed \(M\), at least one character has dimension at least
\(\dim W_h/2^{M/4}\), since the character spaces partition the
monomial basis. This does not identify that character as zero.

## 4. Exact paired vacuum labels

There is a positive construction inside the zero character.
Choose \(R\) of the \(M/4\) pairs

\[
\{k,M-k\},\qquad k=1,3,\ldots,M/2-1,
\]

and let \(S\) be their union. The resulting
\(\binom{M/4}{R}\) monomials have weight \(h=RM\), length
\(2R\), and syndrome zero. Their raw vacuum Gram is the identity:

\[
[1]t^S[-h]t^T=[S=T].
\tag{11}
\]

Use the one-part contraction \([n\mathbin{\&}m=1]\).
Replacing each row part \(n\) by \(M-n\) turns this into
Boolean inclusion on half-indices. Each set \(S\) is invariant
under that complement. In a nonzero determinant matching,
every column half-index is at most its matched row half-index.
Equality of total weights forces equality in every pair.
This is possible exactly when \(S=T\), and then only the
identical matching survives. The general shell proof gives
the details without assuming complement invariance.

Quadratic annihilation does **not** imply annihilation by all
higher normal products. For example,

\[
\sigma_8(\{1,7\})=\sigma_8(\{3,5\})=0,
\qquad
[t_3t_5]\,(t_1t_3t_5t_7)[0]t_1t_7=1.
\tag{12}
\]

Higher Taylor identities still require their contact terms and
higher-degree charges.

## 5. Exact checks

Reproduce the evidence with:

    python3 tools/check_binary_quadratic_charge_shells.py \
      --output evidence/binary-quadratic-charge-shells.json

The checker verifies:

- \(66{,}560\) Lucas coefficient identities and the disjoint-mask
  form (4).
- \(2{,}548\) zero-mode entries against the original normal-order
  operator, independently of the residue formula.
- \(10{,}044\) residue-function identities and \(18{,}665\)
  symmetric-difference pairings.
- Every character dimension through energy \(64\), for
  \(M=4,8,16,32,64\), against the finite bound (7).
- Five full paired vacuum Grams, including \(70\) labels at
  \((M,R)=(32,4)\), and the quartic counterexample (12).

At energy \(32\), charges through weight \(8\) split the \(23\)
regular states into character dimensions \(22\) and \(1\);
the exceptional state is \((1,5,9,17)\).
At energy \(64\), the largest common character dimensions for
\(M=8,16,32,64\) are respectively \(231,113,44,16\).
These values illustrate the exact classification.
The asymptotic conclusion follows from (6)--(10), not from
extrapolating the finite data.
