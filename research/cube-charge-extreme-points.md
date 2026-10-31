# All extreme symbol charges are cube height functions

## Result and scope

The [nonlinear charge theorem](nonlinear-symbol-charge-lattices.md)
characterizes admissible paired potentials by Lip(f)+Lip(g)<=2. Its
extreme points have an exact description at every physical alphabet size:
put twice an integer cube height function on one paired species and zero
on the other. A height function changes by exactly one across every cube
edge and is anchored to zero at the empty occupation.

For the full sixteen-symbol alphabet, this gives 1,980 extreme paired
charges, of which 1,948 are nonlinear in the occupation bits. There is
also an explicit family with independently chosen signs on half of all
symbols. It contains 2^(q/2) vertices on each species and needs endpoint
poles of order at most three, independently of the interval length.
These are exact endpoint-preserving construction parameters. No negative
full-pressure witness or new factoring exponent is asserted.

## 1. A vertex theorem for connected bipartite graphs

Let G be a finite connected bipartite graph with at least two vertices,
and fix a root o. Its anchored unit edge-Lipschitz polytope is

\[
 B_G=\{f:f(o)=0,\ |f(u)-f(v)|\le1\text{ for every edge }uv\}.
 \tag{1}
\]

It is bounded, since |f(v)| is at most the graph distance from o, and
it has full dimension |V(G)|-1. For a feasible f, call an edge tight
when its absolute difference is one.

**Lemma.** A feasible f is a vertex if and only if its tight-edge graph
is connected.

If the tight graph is disconnected, choose a component not containing o.
Every edge from that component to its complement has positive slack.
Increase or decrease all values in the component by the same sufficiently
small positive amount. Both perturbed functions remain feasible, and
their midpoint is f, so f is not a vertex.

Conversely, a connected collection of tight edges determines all values
from f(o)=0 and the signed edge differences. If f is the midpoint of two
feasible functions, every tight inequality must remain an equality in
both functions. Those connected equalities force both functions to equal
f. This proves the lemma.

**Theorem.** The vertices of B_G are exactly the anchored integer
functions with absolute difference one on every graph edge.

At a vertex, connected tight paths from o give integral values. Along
such a path, the parity of the height is the parity of its length.
Bipartiteness makes that parity independent of the chosen path. The
endpoints of every graph edge therefore have opposite height parity.
Their integral difference is odd and has absolute value at most one,
so it has absolute value exactly one. Conversely, a function with all
edges tight has connected tight graph, and is a vertex by the lemma.

## 2. The complete paired-charge vertices

On the m-cube, the anchored edge-Lipschitz quantity is a norm. Write B
for its unit ball. The paired charge domain is

\[
 C=\{(f,g):\operatorname{Lip}(f)+\operatorname{Lip}(g)\le2\}
   =\operatorname{conv}\bigl((2B,0)\cup(0,2B)\bigr).
 \tag{2}
\]

To verify the equality, express a nonzero pair as a convex combination
of its two axis normalizations, with weights Lip(f)/2 and Lip(g)/2;
the remaining weight, if any, goes to the origin. The reverse inclusion
follows from convexity of the two norms.

Every vertex of C must therefore be an extreme point of one of the two
axis copies. Conversely, let h be a vertex of B, and choose a linear
functional that uniquely maximizes at h on B, normalized to have maximum
one there. At a point of C its value
on f is at most Lip(f), at most 2-Lip(g). Reaching its maximum two
forces g=0 and f=2h. Thus every (2h,0), and similarly every (0,2h), is
a vertex. This proves the claimed complete description.

For the reflection tensor, m=M/2 and q=2^m. Combining (2) with the exact
charge theorem gives every extreme potential for the complete alphabet,
without numerical linear programming or a list of tensor inequalities.
Convex combinations recover the full admissible domain. Interior charges
can still matter for a nonlinear objective; a search over vertices alone
does not certify the optimum pressure.

## 3. Exact enumeration through four cube dimensions

Visit cube vertices in increasing occupation weight. Every predecessor of
a new vertex has already been assigned. Its height must belong to the
intersection of the two-element sets consisting of each predecessor's
height minus one and plus one. Branch over that intersection.

Every edge joins adjacent occupation weights, so this procedure checks
every edge exactly when its upper endpoint is visited. It enumerates
every anchored height function exactly once. The resulting counts are:

| Cube dimension m | Alphabet q | Height vertices on one unit ball | Paired charge vertices | Additive paired vertices | Nonlinear paired vertices |
|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 2 | 4 | 4 | 0 |
| 2 | 4 | 6 | 12 | 8 | 4 |
| 3 | 8 | 38 | 76 | 16 | 60 |
| 4 | 16 | 990 | 1,980 | 32 | 1,948 |

An additive height function has an independently chosen coefficient plus
or minus one on each occupation bit, giving exactly 2^m unit-ball vertices.
The reflection interval theorem uses power-of-two M>=4; dimensions one
and three in the table are graph controls, not extra interval claims.

The checker also enumerates every integral point in the full distance
box through dimension three, including unsaturated and parity-mixed
functions, and independently selects those with connected tight graph.
Its vertex sets agree with the constructive enumeration. The evidence
retains the complete height catalogues and their hashes.

## 4. Exponentially many bounded-pole nonlinear choices

Assign height zero to every even-parity cube vertex. Independently assign
height plus or minus one to every odd-parity vertex. Every edge joins
opposite parities, so every edge difference has absolute value one.
All these functions are distinct anchored vertices. There are

\[
 2^{2^{m-1}}=2^{q/2}
 \tag{3}
\]

choices on one species. After multiplication by two, the potential has
values in {-2,0,2}; its oscillation is at most four. With the other
species' potential zero, every paired charge satisfies E<=4. The strict
entry condition v_chi(delta T)+E>0 therefore gives

\[
 v_\chi(\delta T)\ge-3.
 \tag{4}
\]

Thus chi cubed is always a sufficient common endpoint denominator for
this entire growing collection. The numerator degree need only increase
by that actual pole; effective degree accounting is unchanged. The
negative-charge boundary extension also preserves this pole bound,
because it changes only entries whose required order is nonnegative.

This construction separates the number and complexity of admissible
symbol charges from their required pole order. A particular uniform
choice still needs a specified computable sign rule; enumerating the
entire collection is not part of a claimed factoring algorithm.

## Reproduction

```sh
make verify-cube-charge-vertices
```

The constructor is [`cube_charge_vertices.py`](../tools/cube_charge_vertices.py),
and the complete finite record is
[`cube-charge-vertices.json`](../evidence/cube-charge-vertices.json).
This theorem describes the charge domain. It does not restrict the
independent physical-letter heights used in full-profile cost optimization.
