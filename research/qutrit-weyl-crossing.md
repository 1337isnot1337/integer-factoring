# A qutrit Weyl crossing and its constant coefficient

## Result and scope

An order-three Weyl construction gives an explicit binary matrix-valued
crossing, with physical alphabet three and auxiliary dimension three.
The phase cancellation is exact. However, its constant coefficient has
rank \(3^{2t-1}\) in both cuts at length \(t\). Thus the raw coefficient
decomposition costs at least \(d^2/3\), where \(d=3^t\).

This excludes an asymptotic saving from raw monomial coefficients of this
particular pencil. It does not exclude arbitrary recombinations of its
coefficients, different Weyl paths, or different boundaries. It is not a
factoring lower bound.

## 1. An exact binary tensor from noncommuting order-three matrices

Work first in \(F_4=F_2[\omega]/(\omega^2+\omega+1)\). Index matrices by
\(\mathbb Z/3\mathbb Z\), and define

\[
X_{r,s}=[r=s+1],\qquad Z_{r,s}=[r=s]\omega^s,\qquad
W_{a,b}=X^aZ^b .
\]

Then \(ZX=\omega XZ\). Physical indices are \(x,y,u,v\); the crossing is
\(E(x,y,u,v)=[x=v][y=u]\). Define an auxiliary matrix

\[
A(x,y,u,v)=\sum_{a,b=0}^{2}
 (W_{a,b})_{xy}(W_{a,b})_{uv} W_{a,b}.                       \tag{1}
\]

For one entry, the shifts force
\(a=x-y=u-v=r-s\), while the phase is
\(\omega^{b(y+v+s)}\). Since \(1+\omega+\omega^2=0\) and
\(1+1+1=1\) in characteristic two,

\[
\boxed{\quad
A(x,y,u,v)_{r,s}
=[x-y=u-v=r-s]\,[s=-y-v],
\quad}                                                     \tag{2}
\]

with all equalities in \(\mathbb Z/3\mathbb Z\). Thus \(A\) is defined
over \(F_2\); no extension-field descent is needed.

Set

\[
P(z)=A+z(EI_3+A),\qquad f=e=(1,0,0)^T.                     \tag{3}
\]

For length \(t\), concatenate physical label words and contract
\[
G_t(z)=f^TP(z;x_1,y_1,u_1,v_1)\cdots
           P(z;x_t,y_t,u_t,v_t)e.
\]

The local endpoint is \(P(1)=EI_3\), so \(G_t(1)=E^{\otimes t}\).
Consequently its \(t+1\) monomial coefficients are an exact crossing
decomposition on \(d=3^t\) labels.

## 2. The constant coefficient is a partial permutation in the row cut

Let the auxiliary path be \(h_0,\ldots,h_t\), with \(h_0=h_t=0\).
In the product of constant terms, (2) requires

\[
a_i=x_i-y_i=u_i-v_i=h_{i-1}-h_i,\qquad
h_i=-y_i-v_i.                                             \tag{4}
\]

Fix the row-cut label \((x,y)\). It determines \(a_i=x_i-y_i\) and
\[
h_i=-\sum_{j\le i}a_j.
\]
The final boundary is satisfied exactly when \(\sum_i a_i=0\).
If so, the unique opposite label is

\[
v_i=\sum_{j\le i}a_j-y_i,\qquad u_i=v_i+a_i.                \tag{5}
\]

Conversely, \((u,v)\) determines \(a_i=u_i-v_i\), then reconstructs
\((x,y)\) by the same equations. Therefore the row matrix is a partial
permutation, with one nonzero entry for each of the
\(3^{2t-1}\) pairs satisfying \(\sum_i(x_i-y_i)=0\). Hence

\[
\mathrm{rank}R([z^0]G_t)=3^{2t-1}.                   \tag{6}
\]

## 3. The column cut has the same exact rank

Fix the column-cut row label \((x,u)\). Substitute
\(y_i=x_i-a_i\), \(v_i=u_i-a_i\) into (4). Since \(3a_i=0\) in
the label group, the result is

\[
h_{i-1}=-x_i-u_i.                                         \tag{7}
\]

The initial boundary requires \(x_1+u_1=0\).
For every such label, (7) determines \(h_0,\ldots,h_{t-1}\);
the final boundary supplies \(h_t=0\). Then
\[
a_i=h_{i-1}-h_i,\qquad y_i=x_i-a_i,\qquad v_i=u_i-a_i
\]
gives the unique opposite label.

Conversely, \((y,v)\) determines \(h_i=-y_i-v_i\), and hence
the entire path and \((x,u)\). The column matrix is another partial
permutation, with \(3^{2t-1}\) nonzero rows. Thus

\[
\mathrm{rank}C([z^0]G_t)=3^{2t-1}.                   \tag{8}
\]

These are ranks over \(F_2\), and over every extension, because the matrices
are partial permutation matrices.

## 4. Consequence for this coefficient split

The constant coefficient alone contributes

\[
\sqrt{r([z^0]G_t)c([z^0]G_t)}=3^{2t-1}=d^2/3              \tag{9}
\]

to the total rank cost. Its logarithmic exponent tends to two, exceeding
the strict finite-seed exponent \(\sqrt2\). In the raw coefficient-width
interface, the constant piece forces asymptotic \(\rho\ge1\).

This does not transfer to arbitrary linear combinations of the coefficients:
such combinations can cancel the constant coefficient against other terms.
No optimality result for the full coefficient space is claimed.

## 5. Exact controls

The [standalone checker](../tools/check_qutrit_weyl_crossing.py) compares
all 729 entries of the literal \(F_4\) Weyl sum with (2), checks the endpoint,
and contracts the complete polynomial tensor at lengths one through three.
It independently builds the partial permutations from (5) and (7), compares
both complete matrices, and checks ranks with two pivot orders.

| Length | Labels | Complete coefficient rank pairs |
| --- | ---: | --- |
| 1 | 3 | \((3,3),(7,7)\) |
| 2 | 9 | \((27,27),(40,40),(46,57)\) |
| 3 | 27 | \((243,243),(391,391),(428,437),(310,495)\) |

~~~bash
python tools/check_qutrit_weyl_crossing.py --output /tmp/qutrit-weyl-crossing.json
~~~

The [recorded evidence](../evidence/qutrit-weyl-crossing.json) preserves
the checker hash and exact scope. The finite controls test the implementation;
the path bijections prove (6)--(9) at every length.

The construction and its analysis used substantial OpenAI assistant support.
No Lean formalization or new factoring exponent is claimed.
