# Nilpotent parameter jets of the Catalan crossing

The [full cup frame](catalan-cup-crossings.md) has a unitriangular decoder,
so its construction also works over a commutative coefficient ring with
nilpotents. This gives a uniform enlargement of its F2 coefficient space
without requiring a field-valued Gram inverse.

Fix an integer s>=1, work in R_s=F2[u]/(u^s), and set q=1+u. This q is a
unit, with inverse 1+u+...+u^(s-1). Apply the cup-frame, coordinate-decoder
and normalized braid construction over R_s, and write

```math
\mathcal G_n(t,u)=\sum_{j=0}^{s-1}u^j G_{n,j}(t).
```

These are Hasse coefficients: (1+u)^k contributes binom(k,j) modulo two.
Ordinary higher derivatives would lose this information in characteristic two.

## 1. Exact endpoints over the full truncated ring

The physical minor T(q) has diagonal one over every R_s. Its inverse exists
by triangular substitution, even though delta=q+q^-1 is not a unit.
The cup slide uses q^-1 but never delta^-1. Therefore the same ordinary
proof gives

```math
\mathcal G_n(0,u)=I,\qquad \mathcal G_n(1,u)=E.
```

In particular G_(n,0) is the original strict F2 crossing, while every
positive jet has both endpoints zero. A positive jet supplies additional
coefficient tensors whose total is zero, with the same four label alphabets.
No rank over a non-field is substituted for their actual F2 cut ranks.

All these coefficients can be exposed in one strict scalar family. Set
D=n^2+2 and use

```math
\widetilde G_n(t)=G_{n,0}(t)+
                  \sum_{j=1}^{s-1}t^{jD}G_{n,j}(t).   \tag{1}
```

The degree windows are disjoint; both endpoints are preserved; and the
degree is O(sn^2). If s is polynomial in n, the number of pieces is
polynomial. Cut reflections or other specified constant tensors may also
be appended by separate endpoint-vanishing windows, as in the earlier
coherent families. A favorable spanning/rank sum still has to be proved.

## 2. Series compilation without exponential ring tables

The [implementation](../tools/catalan_cup_jets.py) stores s binary coefficient
arrays. It does not allocate a multiplication table of size 2^s by 2^s.
For T(u)=sum T_j u^j, its inverse coefficients satisfy

```math
V_0=T_0^{-1},\qquad
V_j=V_0\sum_{i=1}^j T_i V_{j-i}.                    \tag{2}
```

All products in (2) are ordinary F2 matrix products, with the factor order
retained. Decoder tensor products and boundary maps use the same truncated
coefficient convolution. Selecting a common binary alphabet restricts the
full decoder rows; it does not silently recompute a different minor.

In the exact link-state representation, reconnection coefficients remain
one, while a closed loop has

```math
\delta(u)=(1+u)+(1+u)^{-1}=u^2+u^3+\cdots\pmod{u^s}.
```

Thus a local qI term takes the current u coefficient and its predecessor;
a loop takes the prefix sum starting two orders earlier. The compiler keeps
every intermediate jet. It does not discard mixed terms on the strength of
the final endpoint cancellation.

Finally put h=n^2/2 and gamma=(1+u)^(-h). Its coefficient is
binom(h+j-1,j) modulo two. If A is the decoded unnormalized braid,
the strict family is gamma^2 A+t(gamma+gamma^2)A. This accounts for both
normalization factors, including at the top retained jet.

All arrays are polynomial in s times those of the full cup compiler.
Construction and dense factors therefore cost 2^O(lambda)poly(s) at
n=lambda+O(log(lambda+2)), and (1) has polynomial degree for polynomial s.
This is an effective family construction, not a faster factoring theorem.

## 3. Evidence and present limit

The [checker](../tools/check_catalan_cup_jets.py) verifies inverse-series
identities, all jet endpoints, the baseline specialization and both actual
cut profiles. An independent scalar truncated-ring physical-word walk
compares complete sampled polynomial outputs, including the retained
order-eight coefficients. Its small multiplication tables are confined to
finite controls; the uniform compiler uses (2) and binary convolutions.
The [receipt](../evidence/catalan-cup-jets.json) pins the sources.

```sh
make verify-catalan-cup-jets
```

These jets substantially enlarge the small coefficient spaces. The recorded
exact gamma=7/5 product filter is exploratory: admitting a tensor requires
(r*c)^5<L^14. Its finite spanning outcomes are retained and do not establish
an all-length favorable decomposition. In particular the current five-label
jet-coefficient list misses the crossing after that filter, even through
order eight. Recombination inside the larger span and a uniform rank proof
remain open. The complete factoring score stays1/2; no new Lean is used.
