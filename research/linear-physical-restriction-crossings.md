# Constant linear physical restrictions of the elementary crossing

## Result and scope

Linear combinations of physical letters give a family outside the
[coordinate-code entropy theorem](arbitrary-binary-code-pressure-floor.md).
This note supplies an exact constructor and exhausts one complete finite
case: all constant rank-seven retractions of three binary elementary
sites over `F_2`.

There are 32,640 such projectors. They reduce to 60 symmetry orbits and
produce exactly two pairs of good-cut Smith profiles. Every member has

$$B(s)>686,\qquad
\rho>\frac{\log686}{4\log7}>\frac34
\quad(s\ge1). \tag{1}$$

The opposite cut assignment has `rho>=1`. These statements concern
constant embeddings and retractions without additional physical gauges.
They do not exclude larger ambient spaces, other dimensions, Laurent
physical frames, nonuniform gauges, or new scalar seeds. They are not
lower bounds for factoring.

## 1. A retraction preserves the exact endpoints

Let `V=(F_2^2)^(tensor t)`, with ambient physical dimension `n=2^t`.
Use the elementary tensor

$$T(z)=(1-z)C+zE,$$

where `C` is the straight tensor and `E` is the swap. Choose constant
maps `P:F_2^q->V` and `L:V->F_2^q` satisfying `LP=I`. Define

$$S(z)=(L\otimes L)T(z)^{\otimes t}(P\otimes P), \tag{2}$$

viewed in the row cut `(x,y)|(u,v)`. It has degree at most `t`, and

$$S(0)=I,\qquad S(1)=E_q. \tag{3}$$

The second identity follows because the swap commutes with equal
physical embeddings and retractions. Thus this is a valid polynomial
crossing for every supplied retraction. Its two literal cut matrices are
constructed by
[`linear_physical_restriction.py`](../tools/linear_physical_restriction.py).
The code checks `LP=I`, the full identity at zero, and every swap entry
at one.

An internal change of basis `P->PA`, `L->A^(-1)L` changes both cut matrices
only by constant invertible row and column factors. Their local Smith
orders are unchanged. Consequently, the profiles depend only on the
projector `F=PL`, rather than on its chosen image basis.

## 2. Every rank-seven projector is explicitly covered

For `t=3`, a rank-seven projector on `F_2^8` has a one-dimensional kernel
spanned by a unique nonzero vector `v`, and image `ker(f)` for a unique
nonzero covector `f`. The direct-sum condition is `f(v)=1`, and

$$F=I+vf. \tag{4}$$

Conversely every pair `(v,f)` satisfying this condition gives such a
projector. There are exactly `255*128=32640` pairs.

Choose a coordinate `j` with `f_j=1`. A concrete image basis and left
inverse are, for `i!=j`,

$$P_i=e_i+f_i e_j,\qquad L_i=e_i^{\mathsf T}+v_i f. \tag{5}$$

These satisfy `LP=I` and `PL=I+vf`. Thus no subspace, complementary
kernel, or left inverse is omitted by the finite enumeration.

The group `GL_2(F_2)^3` acts separately on the three binary sites, and
the site-permutation group `S_3` permutes those sites. Each local
`G` commutes with the elementary tensor in the form `G tensor G`,
because both the identity and the swap do. It follows that the induced
group on `V` preserves the profiles under

$$ (v,f)\longmapsto(Gv,fG^{-1}). \tag{6}$$

The checker uses local swaps and transvections, together with adjacent
site permutations. It verifies every generator on the complete ambient
polynomial tensor, then performs exact orbit traversal on all 32,640
pairs. There are 60 orbits. This symmetry reduction is exact and is not
a random sample of retractions.

## 3. Two complete profile pairs

Index the ambient coordinates by binary words `0,...,7` and define

$$Q(v)=v_0v_7+v_1v_6+v_2v_5+v_3v_4\in F_2.$$

This quadratic expression is invariant under the local generators and
site permutations; the same holds for the dual action. Direct exact
Smith elimination on all 60 orbit representatives gives the following
classification by `Q(v)Q(f)`.

| Type | `Q(v)Q(f)` | Projectors | Orbits | Row infinity order multiplicities, orders `0,...,5` |
| --- | --- | --- | --- | --- |
| A | 0 | 25,920 | 51 | `1,6,12,13,11,6` |
| B | 1 | 6,720 | 9 | `1,6,15,7,14,6` |

The row orders here belong to the degree-three reversed matrix. The
column zero orders are the same for both types: their multiplicities at
orders `0,...,4` are

$$1,6,15,18,9. \tag{7}$$

Each list has length 49. The implementation resolves all Smith pivots;
it does not treat an unresolved finite jet as a zero or infinite order.
Independent block-Toeplitz ranks verify both resulting profile pairs.
The finite classification and its orbit coverage constitute the exact
proof for this eight-dimensional case.

## 4. An all-scale rational certificate

Put `x=1/s` and define the integer counting polynomials

$$\begin{aligned}
R_A(x)&=1+6x+12x^2+13x^3+11x^4+6x^5,\\
R_B(x)&=1+6x+15x^2+7x^3+14x^4+6x^5,\\
K(x)&=1+6x+15x^2+18x^3+9x^4.
\end{aligned}$$

The pressure is `s^3 R_A(1/s)K(1/s)` or its type-B counterpart.
These are ordinary real counting polynomials; their coefficients are
not reduced modulo two. Since

$$R_B(x)-R_A(x)=3x^2(1-x)^2,$$

type A has the smaller pressure at every `s>=1`.

To certify (1), consider

$$\begin{aligned}
F(x)&=R_A(x)K(x)-686x^3\\
&=1+12x+63x^2-493x^3+386x^4+537x^5\\
&\qquad+543x^6+405x^7+207x^8+54x^9.
\end{aligned}$$

Write `F(x)=sum_j c_j x^j`. Its degree-48 Bernstein coefficients
have the same signs as the integers

$$d_k=\sum_{j\le k}c_j\binom{48-j}{k-j},\qquad0\le k\le48.$$

Every `d_k` is strictly positive, as recorded in the evidence. Hence
`F(x)>0` throughout `[0,1]`, proving `B(s)>686` at every finite scale
`s>=1`. No floating minimization enters this certificate.

For orientation only, numerical minimization of the exact type-A
polynomial gives pressure approximately `723.66049` and
`rho approximately 0.84591808`; type B is worse at every scale. These
display values are not used to establish the lower bound.

For the opposite orientation, the row zero orders are all zero by
(3). If the column matrix has size `q^2` and degree bound `t`, the sum
of its reversed Smith orders is `tq^2-deg(det)`, which is at most
`tq^2`. Arithmetic-geometric mean then gives column infinity pressure
at least `q^2`. The product is at least `q^4`, proving `rho>=1` for
every constant restriction (2), not only this finite case.

## 5. Replay

Run:

```sh
python3 tools/check_codimension_one_crossings.py
```

The [checker](../tools/check_codimension_one_crossings.py) records all
orbit representatives, exact profile counts, complete symmetry actions,
independent Toeplitz controls, and the positive Bernstein certificate.
The [evidence](../evidence/codimension-one-crossings.json) includes the
constructor and checker hashes.
