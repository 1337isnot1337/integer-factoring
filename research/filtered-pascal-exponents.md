# Nonlinear exponent maps that preserve a complete Pascal cut

## Constructive result

The diagonal exponents in the [weighted Pascal constructor](weighted-pascal-flag-crossings.md)
can follow a much larger nonlinear family while retaining the entire
zero-cut spectrum. The permitted maps have an exact classification.
Their infinity spectra remain available from the actual polynomial
kernel and its Clifford matrix.

For `L=2^r`, write the binary coordinates of `i` as `i_0,...,i_(r-1)`.
Choose constants `a_j` in `F_2` and arbitrary Boolean functions `f_j`
of the lower `j` bits. Set

\[
 \phi_j(i)=a_ji_j+f_j(i_0,\ldots,i_{j-1}),\qquad
 \phi(i)=\sum_{j=0}^{r-1}2^j\phi_j(i).                 \tag{1}
\]

There are exactly `2^(L+r-1)=L 2^(L-1)` such maps. They include
`2^(L-1)` permutations, obtained by setting all `a_j=1`.
These permutations are exactly the isometries of the binary residue
tree: they preserve `v_2(i-j)` for distinct inputs.

The resulting kernels retain the usual zero-cut one-particle list

\[
 1,3,\ldots,2L-1
 \quad\text{with every entry repeated twice}.             \tag{2}
\]

This is a constructor and classification theorem. It supplies no
seed below the half-score threshold. The complete factoring score
remains `1/2`.

## 1. The complete binary power filtration

Work in the `L`-dimensional vector space
`R=F_2[T]/((T+1)^L)`. For a map `phi` from `{0,...,L-1}` to itself,
define the linear operator `U_phi` by `U_phi(T^i)=T^phi(i)`.
No multiplicativity is asserted.

For a bit mask `q`, put

\[
 F_q(T)=\sum_{i\supseteq q}T^i
       =T^q(T+1)^{L-1-q}.                              \tag{3}
\]

The polynomials `F_q` form a basis with distinct `T+1` orders.
Consequently `U_phi` preserves every subspace `(T+1)^k R` if and
only if

\[
 v_{T+1}\!\left(\sum_{i\supseteq q}T^{\phi(i)}\right)
       \ge L-1-q\qquad(0\le q<L).                    \tag{4}
\]

The order of the zero polynomial is infinity.

**Classification.** Condition (4) holds exactly for the maps (1).

To prove this, write
`g_k(i)=binom(phi(i),k) mod2`. Lucas's binary identity gives
`g_k(i)=prod_(j in k) phi_j(i)`. The coefficient of `(T+1)^k`
in the expression in (4) is

\[
 \sum_{i\supseteq q}g_k(i).                           \tag{5}
\]

Expand a Boolean function uniquely in square-free monomials:
`g_k(i)=sum_s c_s prod_(j in s) i_j`. The sum in (5) is then
`sum_(s superset (L-1-q)) c_s`. Descending Boolean inversion shows
that (5) vanishes whenever `k<L-1-q` if and only if `c_s=0`
for every integer mask `s>k`.

For `k=2^j`, this says that the algebraic normal form of `phi_j`
contains only masks at most `2^j`. Its only monomial containing
the bit `i_j` is therefore the single variable `i_j`; all remaining
monomials use lower bits. This proves necessity of (1).

Conversely, each factor `phi_j` in `g_k` has monomial masks at most
`2^j`. A product chooses one such mask per factor. The integer
value of their union is at most their sum, which is at most `k`.
Thus `g_k` has no monomial mask greater than `k`. This proves
sufficiency, including maps that are not injective.

When all `a_j=1`, successive recovery of input bits shows that
`phi` is a permutation. The inverse has the same form. Its linear
operator preserves the same flag, so the induced map on each
one-dimensional successive quotient is nonzero. Over `F_2` this
forces the sharper identity

\[
 v_{T+1}\!\left(\sum_{i\supseteq q}T^{\phi(i)}\right)
       =L-1-q.                                       \tag{6}
\]

If some `a_j=0`, the first such index makes the output prefix
non-surjective. The full map therefore cannot be a permutation.
This establishes the permutation count. For a permutation, the
first input bit at which two inputs differ remains their first
differing output bit, proving the tree-isometry description. The
converse follows by recovering the output bits successively from
the preserved residue partitions.

## 2. An actual crossing kernel for every classified map

Put `c=1+u`, `P_ij=[j subseteq i]`, and define

\[
 E_\phi(i)=2\phi(i)+1,\qquad E_0(i)=2i+1,
\]

\[
 A=\mathrm{diag}(c^{E_\phi(i)})P^{\mathsf T},\qquad
 B=P^{\mathsf T}\mathrm{diag}(c^{2L-E_\phi(i)})P,
\]

\[
 \mathcal D_\phi=
 \begin{pmatrix}A&u^{2L}P\\I&B\end{pmatrix}.          \tag{7}
\]

Every exponent is positive and smaller than `2L`. The same direct
algebra as in the weighted constructor gives

\[
 AB=c^{2L}P,\qquad BPA=c^{2L}I,\qquad
 \det\mathcal D_\phi=1.                              \tag{8}
\]

The endpoint at `u=1` is `[[0,P],[I,0]]`. Thus it supplies the same
crossing after independent matching-type restrictions on `(B,Y)`
and `(X,V)`. The fixed ordinary Pascal index gives the strict
endpoint filtration; monotonicity of `phi` is unnecessary.

For `q=i|j`, the lower-right entry is

\[
 B_{ij}=c\sum_{a\supseteq q}(c^2)^{L-1-\phi(a)}.       \tag{9}
\]

The complemented map `L-1-phi` still has the form (1), since
complementation adds one to each Boolean output coordinate. Apply
(4) with `T=c^2=1+u^2` to obtain

\[
 v_u(B_{ij})\ge2(L-1-(i|j)).                         \tag{10}
\]

This is equality for every entry when `phi` is a permutation.
In boundary order `(B,X,V,Y)`, use zero heights

\[
 h_0(i)=(0,E_0(i)-2L,0,E_0(i)).                     \tag{11}
\]

The internal `A` block is regular. The weighted valuation of an
internal `B` entry is at least

\[
  E_0(i)+E_0(j)-2L+2(L-1-(i|j))=2(i\mathbin{\&}j)\ge0.
                                                               \tag{12}
\]

Both internal Gaussian factors are therefore invertible over the
valuation ring. The remaining cross block is the same weighted
`[[0,I],[u^(2L)P,0]]` as in the interval constructor. Its complete
one-particle Smith list is (2), and the full Fock list consists of
all subset sums of those `2L` orders.

Conversely, regularity of `B` at (11), restricted to entries `B_i0`,
forces (10) for every mask. The classification in Section 1 is
therefore exact for these fixed internal regularity inequalities.
It does not classify every possible way to obtain the same Smith
list outside that regularity cone.

All entries of (7) have degree at most `2L`. The identities (8)
also give the same polynomial inverse and Clifford formulas as in
the weighted flag note. Those formulas remain valid for arbitrary
physical heights at either place.

## 3. Scope of the enlarged family

For `L=8`, take

\[
 \phi_0=i_0,\qquad \phi_1=i_1,\qquad
 \phi_2=i_0i_1.
\]

Then `phi=(0,1,2,7,0,1,2,7)` and
`E_phi=(1,3,5,15,1,3,5,15)`. This is a nonlinear, noninjective
exponent map. Nevertheless its entire zero profile is exactly (2).
It is outside the positive coordinate-weight formula.

The exponent list `E_phi` is generally **not** the infinity Smith
list. For example, at `L=2`, the permutation `phi=(1,0)` and the
candidate infinity heights
`(E_phi-L,L,E_phi-L,-L)` give actual one-particle orders
`(-3,-3,3,3)`, not `(-3,-1,1,3)`. The actual infinity matrix must
be evaluated. An arbitrary optimization must also retain the two
independent matching-type costs.

This distinction matters: preserving one full cut does not certify
a useful tradeoff at the other cut. Numerical searches in this
family are exploratory until an actual candidate has exact profile
and pressure certificates.

## 4. Reproducible exact controls

The [constructor](../tools/filtered_pascal_exponents.py) validates
the Boolean classification and builds the actual polynomial matrix.
The [checker](../tools/check_filtered_pascal_exponents.py) tests
the classification by complete small enumeration, verifies literal
upset-polynomial orders, and checks polynomial inverse identities,
Clifford orders, and original Fock profiles for finite controls.

Run:

```sh
python3 tools/check_filtered_pascal_exponents.py \
  --output evidence/filtered-pascal-exponents.json
```

The theorem and controls used substantial OpenAI assistant support.
No new Lean formalization or independent external review is claimed.
