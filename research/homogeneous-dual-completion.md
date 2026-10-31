# An affine family from all homogeneous completions of a selected dual

## Result and scope

At weight 32, the characteristic-two quotient has 23 monomial states.
The reflected sector with four odd parts below 16 and total weight 32
has eight labels.  Its complete set of homogeneous dual choices has
120 independent scalar parameters: fifteen kernel directions for each
of the eight duals.

The resulting complete polynomial tensors span an **88-dimensional
affine family**.  Four exact dependencies reduce the fifteen kernel
directions to eleven distinct tensor directions per first label.  The
dimension statement is a finite binary linear-algebra certificate; it
holds over every coefficient field of characteristic two.

Every member has the literal crossing endpoint.  The family permits
mixing two-part and four-part operator parameters and thus goes beyond
changing the physical weights of a fixed monomial dual.  It does not
supply a strict pressure witness or improve the factoring exponent.

The construction uses the source modes from Anthropic's Appendix B,
the [modular quotient](modular-crossing.md), and the
[reflected neutral compiler](neutral-reflected-crossings.md).
These are ordinary proofs and exact finite checks, with no new Lean work.

## 1. All right inverses of a selected vacuum pairing

Let `V` be a finite homogeneous parameter space and let
`b_0,...,b_(D-1)` be the selected states.  Define

\[
 G(a)_j=[1]a[-h]b_j.
 \tag{1}
\]

Suppose `a_0,...,a_(D-1)` satisfy `G(a_i)=e_i` and are members of a
monomial basis of `V`.  For each remaining monomial `m`, put

\[
 k_m=m-\sum_jG(m)_j a_j.                                  \tag{2}
\]

The vectors `k_m` form a basis of `ker G`: their coefficients on the
remaining monomials form the identity matrix, and subtracting the
`a_j` components puts every vector in the kernel.  Thus every dual
family is uniquely

\[
 \widetilde a_x=a_x+\sum_m c_{x,m}k_m.                      \tag{3}
\]

The source crossing is linear in its first parameter.  Its value at
the crossing endpoint depends only on (1), so (3) preserves that
endpoint for all coefficients.  All parameters remain homogeneous of
weight `h`; their angle factors are unchanged.  No determinant inverse
is recomputed when the coefficients change.

## 2. The weight-32 specialization

Work in `W_2` and take

\[
 L=\{\mu:\text{four distinct odd parts below 16},\ |\mu|=32\}.
\]

There are eight such labels.  Set `b_mu=t^mu` and
`a_mu=t^(16-mu)`.  The reflected duality theorem supplies (1).
The full regular-partition basis of weight 32 has 23 elements,
so (2) gives fifteen independent kernel vectors and (3) has
`8*15=120` parameters.

Eight monomials in the full space have length two.  Their vacuum
pairing with every four-part state is zero: reaching the vacuum requires
one annihilation for each input part, and there are only two operator
factors.  These are particularly simple kernel directions.  The other
seven directions use the exact subtraction in (2).

The [constructor](../tools/homogeneous_dual_completion.py) lists the
labels, baseline reflected parameters, and all fifteen kernel vectors.
It keeps each completed tensor as a linear combination of exact
homogeneous source tensors.

## 3. A direct compiler for arbitrary homogeneous monomials

The reflected Pfaffian formula extends to an arbitrary regular monomial
parameter `a=t^alpha` of weight `h`, even if its number of parts differs
from the selected states.  Choose a power of two `Q` larger than every
part involved.  For every operator part `m` in `alpha`, use an `A` node
with index `Q-m` in the kernel table of
[neutral-reflected-crossings.md](neutral-reflected-crossings.md).
Use the actual parts of `v,y,u` in the other three groups.

The gauge product on the `A` nodes is `z^(-h)c^h`, and the product on
the `B` nodes is `c^h`.  Thus the total prefactor remains
`z^(-h)c^(2h)`, independently of `Q` and of the number of parts of `a`.
The same Clifford matching proof therefore gives the complete source
polynomial.  Changing `Q` changes the displayed intermediate kernel,
but not the resulting Laurent entry.

This step matters for the two-part directions: padding them with fictitious
operator factors would change the source tensor.  The compiler instead
uses their actual number of nodes.

The construction clears one common Laurent monomial, retains every
coefficient, and verifies a common exponent spacing from the complete
family.  For the recorded family the lowest exponent is `-64` and the
common spacing is 2.  The resulting seed polynomials have degree at
most 64.  The eight-letter baseline alone admits spacing 16; that
larger spacing must not be applied to the new directions.

## 4. The exact dimension of the tensor family

Fix a first logical label `x`.  A kernel vector gives a slab indexed
by `(y,u,v)`, with a polynomial at each of its `8^3` entries.  Pack all
coefficients into one binary vector and perform exact elimination.
The fifteen slabs have rank eleven.

In the constructor's zero-based kernel order, independent directions are

\[
 0,1,2,3,4,5,6,7,8,9,11.
\]

The four relations are

\[
 \begin{aligned}
 K_6+K_8+K_9+K_{10}&=0,\\
 K_6+K_8+K_{11}+K_{12}&=0,\\
 K_6+K_9+K_{11}+K_{13}&=0,\\
 K_8+K_9+K_{11}+K_{14}&=0.
 \end{aligned}                                           \tag{4}
\]

Each identity is checked coefficient by coefficient across the complete
slab.  A nonzero binary rank-eleven minor remains nonzero over every
characteristic-two extension field, and (4) gives the matching upper
bound there.

Different first labels occupy disjoint tensor slabs.  Therefore the
complete tensor direction space has dimension `8*11=88`.  This is a
dimension of a finite coefficient space, not a rank of either physical
flattening and not a factoring rate.

## Reproduction

```sh
python3 tools/check_homogeneous_dual_completion.py
```

The checker reconstructs the full affine data, verifies every vacuum
kernel, checks the two-particle directions, compares exact source
coefficients, proves the eleven-dimensional slab rank and all four
relations, and verifies literal endpoints after independent F4
specializations.  The retained
[coefficient data](../evidence/homogeneous-dual-completion-q16-h32.json)
are sufficient to regenerate any specialization without recomputing
the Pfaffians.

Any pressure search on this space must charge the full degree-64
support, both complete cut profiles, and the physical type weights.
Passing the endpoint is only an admissibility check.  No strict
stronger-than-half certificate is claimed.

This note, constructor, and checker used substantial OpenAI assistant support.
