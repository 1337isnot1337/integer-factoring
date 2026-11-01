# Explicit Pascal cut profiles for every exponent permutation

## Result and scope

Every permutation of a power-of-two interval gives an actual Pascal
crossing with an explicit complete zero spectrum. No linearity,
filtration preservation, or monotonicity assumption on the permutation
is necessary. A coefficientwise-product inequality supplies the
regularity that a projected-rank argument previously established for
binary linear maps.

A compatible cone of infinity heights also gives a complete explicit
infinity spectrum. Both profiles belong to the same polynomial kernel.
The two matching-pair physical energies are stated separately and must
both be charged. These formulas enlarge the constructive search space;
they do not supply a crossing below `1/sqrt(2)`. The factoring score
remains `1/2`.

## 1. A coefficientwise-product valuation inequality

Let `L=2^r`, let the coefficient field have characteristic two, and
consider polynomials of degree smaller than `L`. Write `v` for order
at `T=1`, and write `F circ G` for coefficientwise multiplication
in the monomial basis. Then

\[
 v(F\circ G)\ge v(F)+v(G)-(L-1).                       \tag{1}
\]

The inequality is immediate if either polynomial is zero. Otherwise
use the basis

\[
 H_q(T)=\sum_{i\supseteq q}T^i
       =T^q(T+1)^{L-1-q},\qquad 0\le q<L.             \tag{2}
\]

Its valuations are all distinct: `v(H_q)=L-1-q`. Consequently a
polynomial of valuation at least `a` is a linear combination of the
`H_q` with `q<=L-1-a`. Moreover,

\[
 H_i\circ H_j=H_{i\mathbin{|}j},\qquad i\mathbin{|}j\le i+j.
                                                               \tag{3}
\]

Expanding both factors in this basis proves (1), including arbitrary
field coefficients. Cancellation can only raise the resulting order.

## 2. Apply the inequality after an arbitrary permutation

Let `pi` be any permutation of `0,...,L-1`. Set

\[
 F_q(T)=\sum_{i\supseteq q}T^{\pi(i)},\qquad
 \nu(q)=v(F_q).                                       \tag{4}
\]

Every `F_q` is nonzero and has degree below `L`. Also
`nu(0)=L-1` and `nu(L-1)=0`. Because a permutation preserves
coefficientwise products of the indicator vectors in (2),

\[
 F_i\circ F_j=F_{i\mathbin{|}j}.
\]

Equation (1) gives the exact inequality needed below:

\[
 \boxed{\nu(i)+\nu(j)-\nu(i\mathbin{|}j)\le L-1.}      \tag{5}
\]

There is no claim that `nu` is monotone or supermodular. For example,
`pi=(0,1,2,3,4,6,5,7)` has
`nu=(7,1,1,4,3,1,2,0)`, which is not monotone under containment.
It still satisfies (5).

## 3. The same polynomial crossing and its zero cut

Let `P_ij=[j subset i]`, `c=1+u`, and put

\[
 E_\infty(i)=2\pi(i)+1,\qquad
 A=\operatorname{diag}(c^{E_\infty(i)})P^{\mathsf T},
\]

\[
 B=P^{\mathsf T}\operatorname{diag}(c^{2L-E_\infty(i)})P,
 \qquad
 \mathcal D=\begin{pmatrix}A&u^{2L}P\\I&B\end{pmatrix}.
                                                               \tag{6}
\]

These are the actual kernels from the
[weighted Pascal constructor](weighted-pascal-flag-crossings.md).
The identities `AB=c^(2L)P` and `BPA=c^(2L)I` imply determinant
one and its polynomial inverse. At `u=1`, the kernel is
`[[0,P],[I,0]]`, with the same full matching-type endpoint proof.

Define

\[
 E_0(i)=2L-1-2\nu(i).                                 \tag{7}
\]

These are positive odd integers. Since
`B_ij=c^(2L-1) F_(i|j)(c^(-2))` and
`c^(-2)-1=u^2/c^2`, its exact zero order is

\[
 v_u(B_{ij})=2\nu(i\mathbin{|}j).                     \tag{8}
\]

Choose physical zero heights in order `(B,X,V,Y)` as

\[
 h_0=(0,E_0-2L,0,E_0).                                \tag{9}
\]

The internal `A` factor is regular. The order of a weighted internal
`B` entry is

\[
 2\bigl(L-1-\nu(i)-\nu(j)+\nu(i\mathbin{|}j)\bigr)\ge0
                                                               \tag{10}
\]

by (5). Both internal Gaussian factors are therefore unimodular over
the zero valuation ring. The remaining crossblock is
`[[0,I],[u^(2L)P,0]]` with the column heights in (9).
Its two blocks are `diag(u^E_0)` and `P diag(u^E_0)`.
The constant invertibility of `P` gives the complete one-particle
profile: **two copies of every `E_0(i)`**. Its full Fock profile
consists of all subset sums. This step does not require monotonicity
of `E_0`.

For complementary masks, (5) gives
`E_0(i)+E_0(bar i)>=2L`. The changed zero energies therefore
change the physical charge as well. The
[binary linear formula](linear-pascal-matroid-energies.md)
is the special case `nu(i)=g(bar i)` obtained from minimum-weight
projected row bases.

## 4. An explicit compatible infinity cone

Choose any integer vector `beta` satisfying
`beta(i)<=beta(j)` whenever `i subset j`. Use infinity heights

\[
 h_\infty=(\beta-L,L,E_\infty-L,-L).                  \tag{11}
\]

Write `T=P^T`. The exact constant operations

\[
 \operatorname{diag}(I,T)\mathcal D\operatorname{diag}(T,P)
 =\begin{pmatrix}
    \operatorname{diag}(c^{E_\infty})&u^{2L}I\\
    I&\operatorname{diag}(c^{2L-E_\infty})
   \end{pmatrix}                                     \tag{12}
\]

remain unimodular after the heights (11). The factors on `X` and
`Y` have constant heights. The factor on `B` is regular precisely
under the stated containment inequalities for `beta`.

After interleaving modes, (12) is a direct sum of two-by-two blocks.
In the block indexed by `i`, the four entry valuations are two copies
each of `beta(i)-2L` and `E_infty(i)-2L`. Its unweighted
determinant is `c^(2L)+u^(2L)=1`. Its complete weighted infinity
orders are therefore

\[
 \boxed{\min\{\beta(i),E_\infty(i)\}-2L,\quad
        \max\{\beta(i),E_\infty(i)\}.}                \tag{13}
\]

The full infinity Fock profile consists of all subset sums of (13).
For example, `beta(i)=2i+1` is always admissible. A computed
`E_0` is an admissible choice of `beta` only when it has the required
containment monotonicity.

The physical energy at the `B,Y` pair and at the `X,V` pair is,
respectively,

\[
 G_{BY}(i)=E_0(i)+\beta(i)-2L,\qquad
 G_{XV}(i)=E_0(i)+E_\infty(i)-2L.                     \tag{14}
\]

These need not agree. Both type partition functions enter the
pressure. Formula (13) is an actual common-kernel spectrum; it is
not a favorable profile obtained by specializing another variable.

## 5. Effective construction and exact controls

The [constructor](../tools/permuted_pascal_kernels.py) evaluates (4)
by a Boolean zeta transform on literal binary polynomials. It then
constructs (6) and the infinity gauge (11). These one-particle data
and individual coefficient computations have polynomial cost in `L`,
hence singly exponential in `r=log_2 L`. A written full Fock table
still has exponentially many entries in `L`.

The [checker](../tools/check_permuted_pascal_kernels.py) verifies (1)
for every pair of binary polynomials through `L=8`, and (5) for all
40,320 permutations at `L=8`. It reconstructs actual zero and
infinity profiles, including nonmonotone examples, and compares
small cases against the original physical Fock tensor.

```sh
python3 tools/check_permuted_pascal_kernels.py \
  --output evidence/permuted-pascal-zero-spectra.json
```

The displayed spectra enlarge the candidate family. An improving
pressure certificate remains necessary before changing the factoring
score. No new Lean formalization is claimed.
