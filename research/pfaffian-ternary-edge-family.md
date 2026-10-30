# A three-parameter Pfaffian family on fixed-length codes

## Result and scope

The local tensor in
[complement-shell-crossings.md](complement-shell-crossings.md)
extends to exact crossings on arbitrary fixed-length occupation codes.
An energy restriction is no longer required. Six uniform edge
polynomials reduce to three products. The whole tensor is homogeneous
of degree \(L\) in those products, so its coefficient space has
dimension at most \(\binom{L+2}{2}\), regardless of the degrees of
the chosen univariate polynomials.

The resulting finite search checks every binary combination of the
independent coefficient tensors and both complete cut ranks. It also
exports the same binary basis for extension-field searches. Linear
injections and projections provide a further exact construction beyond
monomial label subsets.

The six finite binary spaces below fail the seed criterion after
arbitrary coefficient recombination. Integer shortest-path lower
bounds certify these exclusions. They do not exclude other codes,
larger domains, different choices at different sites, linear retracts,
or other fields.

The family extends the characteristic-two Fock tensor. The original
operators and profile derivation are attributed in the complement-shell
note. This extension, proof, and search code used OpenAI assistant
support. No new Lean proof or improved factoring exponent is asserted.

## 1. The six-edge family

Work over a characteristic-two field. Fix \(P\) complementary mode
pairs, containing \(N=2P\) physical modes. A word in
\(\{0,1,2,3\}^P\) records an occupied subset in each pair; its length
is the total number of occupied modes.

Choose \(d\) distinct length-\(L\) parameter labels and \(d\)
distinct length-\(L\) state labels. The collections may differ.
Use four blocks \(U,A,V,Y\) of occupied-mode vertices. The
\(A,V\) labels use the parameter collection, and \(U,Y\) use
the state collection.

For six commuting scalars or polynomials \(f,g,h,k,\ell,m\),
put these entries in each local alternating matrix:

| Blocks | Entry |
| --- | --- |
| \(U,A\) | \([u=a]f\) |
| \(A,Y\) | \([a\ne y]g\) |
| \(A,V\) | \([a=v]h\) |
| \(U,V\) | \([u\ne v]k\) |
| \(V,Y\) | \([v=y]\ell\) |
| \(U,Y\) | \([u=y]m\) |

Comparisons refer to the two mode bits within a pair. Within-block
entries are zero. Define \(\mathcal G(f,g,h,k,\ell,m)\) as the
product of the local Pfaffians, equivalently the Pfaffian of their
block-diagonal union.

The original complement-shell profile is the specialization

\[
(f,g,h,k,\ell,m)=(t+t^2,1+t,t,1,1,1),
\tag{1}
\]

whose product is \(t^L F(z)\), where \(t=z^M\).

## 2. Reduction to three products

In one perfect matching, let the six edge counts in table order be
\(a,b,c,d,e,f_0\). Each global vertex block has length \(L\), so

\[
a+b+c=L,\quad a+d+f_0=L,\quad
c+d+e=L,\quad b+e+f_0=L.
\]

These integer count equations imply

\[
a=e,\qquad b=d,\qquad c=f_0,\qquad a+b+c=L.
\tag{2}
\]

For example, adding the \(A,U\) counts and subtracting the \(V,Y\)
counts gives \(2(a-e)=0\) as an integer equation. This argument
does not divide by two in the coefficient field.

The polynomial weight of the matching is therefore
\((f\ell)^a(gk)^b(hm)^c\).

**Theorem 1.** There are fixed tensors \(T_{a,b,c}\), indexed by
\(a,b,c\ge0\), \(a+b+c=L\), such that

\[
\mathcal G(f,g,h,k,\ell,m)=G_L(f\ell,gk,hm),
\qquad
G_L(p,q,r)=\sum_{a+b+c=L}p^a q^b r^c T_{a,b,c}.
\tag{3}
\]

This all-size identity follows by grouping the exact perfect-matching
expansion using (2). Uniform edge choices across the sites are a
hypothesis. With different edge polynomials at different sites,
global count equality does not give (3).

At \(p=q=0,r=1\), every \(A\) vertex matches an identical \(V\)
mode, and every \(U\) vertex matches an identical \(Y\) mode.
Consequently

\[
G_L(0,0,1)=T_{0,0,L}=E_d,\qquad
E_d(x,y,u,v)=[x=v][y=u].
\tag{4}
\]

Any univariate polynomials with \(p(1)=q(1)=0,r(1)=1\) therefore
give a complete crossing decomposition by the coefficients of
\(G_L(p(t),q(t),r(t))\). This works over every characteristic-two
extension field.

The endpoint proof uses equal global lengths. The unrestricted
four-state local tensor need not itself be a crossing.

## 3. The complete uniform edge-polynomial space

Let

\[
\mathcal V=\operatorname{span}_{\mathbf F_2}
 \{T_{a,b,c}:a+b+c=L\}.
\tag{5}
\]

Its dimension is at most \(\binom{L+2}{2}\). Every coefficient of
every univariate specialization above lies in \(\mathcal V\), as
does every binary recombination of those coefficients. Searching
all decompositions of \(E_d\) with pieces in \(\mathcal V\)
therefore covers all uniform binary edge-polynomial choices,
including arbitrarily large polynomial degrees.

The checker constructs the tensors by setting

\[
(f,g,h,k,\ell,m)=(z,z^{L+1},1,1,1,1).
\tag{6}
\]

The coefficient at degree \(a+(L+1)b\) is
\(T_{a,b,L-a-b}\). Since \(a+b\le L\), no encoded exponents
collide. The constant coefficient is the crossing.

It constructs every coefficient in both cuts, then selects an
independent binary basis. The crossing is its first vector.
Binary linear independence persists after field extension because
a nonzero binary minor stays nonzero. Thus the same exported basis
describes the complete extended coefficient space.

## 4. Exact syndrome blocks

For a label \(A\), let \(\sigma(A)\in\mathbf F_2^P\) record
occupation parity at each complementary pair. A local Pfaffian with
an odd number of vertices vanishes, so any nonzero entry satisfies

\[
\sigma(x)+\sigma(y)+\sigma(u)+\sigma(v)=0.
\tag{7}
\]

Both cuts,

\[
(x,y)\mid(v,u),\qquad (x,u)\mid(v,y),
\]

split into direct sums grouped by the XOR of the two syndromes.
The rank of a complete piece is the sum of every block rank.
The checker does not omit blocks, sample minors, or substitute a
state-space dimension for a matrix rank.

It also groups entries with the same multiset of local four-state
patterns. This saves repeated polynomial evaluations while retaining
each matrix entry.

## 5. Exhaustive binary recombination

Let \(B_0,\ldots,B_{k-1}\) be the independent basis, with
\(B_0=E_d\). Each coefficient-space element is uniquely

\[
A_v=\sum_i v_i B_i,\qquad v\in\mathbf F_2^k.
\]

Write \(r_v,c_v\) for its two complete cut ranks. A decomposition
of \(E_d\) is a sequence of vectors whose XOR is \(e_0\);
its cost is \(\sum\sqrt{r_vc_v}\). This is a shortest-path
problem in the finite additive group \(\mathbf F_2^k\), with
edge weight \(\sqrt{r_vc_v}\) for step \(v\).

For a rigorous exclusion, replace every edge weight by

\[
w_v=\left\lfloor\sqrt{r_vc_v}\right\rfloor.
\tag{8}
\]

The integer shortest-path distance \(D\) is a lower bound for every
real decomposition cost. The checker certifies failure of the
strict seed target by the integer comparison

\[
D^{70}>d^{99},\qquad \sqrt2<99/70.
\tag{9}
\]

The displayed real costs describe explicit decompositions found by
the numerical search. The exclusions use (8)--(9), so numerical tie
handling between sums of radicals does not affect them.

A full \(P,L\) domain consists of all \(L\)-subsets of \(2P\)
modes and has \(d=\binom{2P}{L}\).

| Domain | \(d\) | \(\dim\mathcal V\) | Found cost | Integer lower bound \(D\) | \(d^{\sqrt2}\) |
| --- | ---: | ---: | ---: | ---: | ---: |
| \(M=16,L=4,j=32\) energy shell | 8 | 7 | 44.000 | 44 | 18.931 |
| \(M=32,L=3,j=47\) energy shell | 28 | 9 | 568.030 | 567 | 111.324 |
| Full \(P=2,L=2\) | 6 | 6 | 26.000 | 26 | 12.603 |
| Full \(P=3,L=2\) | 15 | 6 | 126.500 | 126 | 46.052 |
| Full \(P=3,L=3\) | 20 | 10 | 189.880 | 188 | 69.173 |
| Full \(P=4,L=2\) | 28 | 6 | 356.269 | 356 | 111.324 |

Every basis combination is checked in both cuts. The report includes
basis selection, complete rank lists, selected decompositions,
integer lower bounds, and threshold certificates.

These results go beyond checking the original monomial coefficients.
For example, the eight-label energy shell has original profile cost
\(87.984\). Grouping those coefficients reduces it to \(45\);
the broader ternary space contains a decomposition of cost \(44\).
Each remains above the required threshold.

## 6. Linear retracts give further exact constructions

Let the full fixed-length space have dimension \(D_0\).
Choose an injection \(I:\mathbf F_2^d\to\mathbf F_2^{D_0}\)
and projection \(P:\mathbf F_2^{D_0}\to\mathbf F_2^d\)
with \(PI=I_d\). Apply \(P\) to the \(x,u\) legs and \(I\)
to the \(v,y\) legs:

\[
\widetilde G_{x,y,u,v}
=\sum_{\alpha,\beta,\gamma,\delta}
P_{x,\alpha}I_{\beta,y}P_{u,\gamma}I_{\delta,v}
G_{\alpha,\beta,\gamma,\delta}.
\tag{10}
\]

At the crossing point, (4) gives

\[
\widetilde G(0,0,1)_{x,y,u,v}
=(PI)_{x,v}(PI)_{u,y}
=[x=v][u=y].
\tag{11}
\]

The parameter and state legs may use different retracts of the same
rank. Linear combinations of occupation labels are therefore allowed.
The coefficient-dimension bound survives because (10) is a fixed
linear map on tensors.

This is an exact construction interface, not an existence theorem
for a favorable retract. The finite exclusions above do not exclude
all retracts of the ambient spaces. An improved candidate must
supply its injection and projection and its complete two-cut ranks.

For fixed finite data, all coefficients and factorizations are
constructible by ordinary finite-field arithmetic. A growing family
must also establish its construction cost and full rank-sum bound.

## 7. Reproduction and basis export

Reproduce the binary certificates with:

    python3 tools/explore_complement_edge_pencils.py \
      --universal-span \
      --cases shell16:4:32,shell32:3:47,length2:2,length3:2,length3:3,length4:2 \
      --output evidence/pfaffian-ternary-F2-spans.json

Regenerate the four extension-field bases with:

    python3 tools/explore_complement_edge_pencils.py \
      --export-basis-dir evidence \
      --cases shell16:4:32,length2:2,length3:2,length3:3

The canonical files are:

- [pfaffian-basis-shell16_L4_j32.json](../evidence/pfaffian-basis-shell16_L4_j32.json)
- [pfaffian-basis-all_P2_L2.json](../evidence/pfaffian-basis-all_P2_L2.json)
- [pfaffian-basis-all_P3_L2.json](../evidence/pfaffian-basis-all_P3_L2.json)
- [pfaffian-basis-all_P3_L3.json](../evidence/pfaffian-basis-all_P3_L3.json)

Each file contains occupation labels, block sizes, independent
coefficient monomials, both cut matrices as hexadecimal binary row
masks, endpoint coordinates, and the exporter SHA-256 hash.
The first basis vector is the crossing.

The startup check compares all \(256\) original local entries with
the independently written complement-shell implementation. Each
code calculation checks the exact crossing on every retained
entry; (7) proves that entries outside the syndrome blocks vanish.

## 8. Complete extension-field screens over four elements

The same four exported coefficient spaces were independently screened over
\(\mathbb F_4=\mathbb F_2[t]/(t^2+t+1)\). Nonzero scalar multiples have
the same two cut ranks, so one representative from each projective class
suffices. A space of dimension \(k\) has \((4^k-1)/3\) such classes.

Choose an integer cost gate \(g\) satisfying \(g^{70}>d^{99}\).
Then \(g>d^{\sqrt2}\). A decomposition meeting the strict seed criterion
could use only pieces with \(r c<g^2\). The native screen enumerates every
projective vector and computes exact ranks for every vector below that gate.
It may stop elimination early on the others, only after their partial rank
product has already reached the gate. Thus rejected ranks are certified
lower bounds, not invented exact values.

| Domain | Dimension | Projective classes | Gate \(g\) | Classes below gate | Rank of their span | Crossing in that span? |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Eight-label energy shell | 7 | 5461 | 19 | 33 | 6 | No |
| Full \(P=2,L=2\) | 6 | 1365 | 13 | 29 | 5 | No |
| Full \(P=3,L=2\) | 6 | 1365 | 47 | 8 | 6 | Yes |
| Full \(P=3,L=3\) | 10 | 349525 | 70 | 40 | 8 | No |

The three failed span tests exclude every \(\mathbb F_4\) decomposition
immediately. In the 15-label case, include all three nonzero scalar
multiples of each affordable projective vector as steps in the additive
coefficient group. Integer lower bounds for the square-root edge costs give
a shortest-path lower bound of \(140.678999262545\) within that restricted
graph. Any omitted vector costs at least 47 by itself. Hence every full
decomposition costs at least 47, which already exceeds
\(15^{\sqrt2}\). The restricted graph distance is not asserted to be the
minimum cost over the entire coefficient space.

An independent Python implementation uses scalar field multiplication and
the opposite pivot direction. It rechecks every affordable vector and 40
other projective vectors in each space, as well as the endpoint ranks,
field-span certificates, and shortest-path bound. The all-vector coverage
comes from the complete native enumeration; the independent replay checks
its arithmetic and every vector used by the final certificates.

The records and certificates are in
[pfaffian-F4-span-audit.json](../evidence/pfaffian-F4-span-audit.json).
For example, reproduce the six-label case with:

    c++ -O3 -std=c++17 tools/screen_pfaffian_span_f4.cpp -o /tmp/pfaffian-screen
    python3 tools/export_pfaffian_span_native.py evidence/pfaffian-basis-all_P2_L2.json --output /tmp/pfaffian-input.txt
    /tmp/pfaffian-screen < /tmp/pfaffian-input.txt > /tmp/pfaffian-screen.jsonl
    python3 tools/check_pfaffian_span_f4.py --basis evidence/pfaffian-basis-all_P2_L2.json --screen /tmp/pfaffian-screen.jsonl --output /tmp/pfaffian-audit.json

The native program also accepts `--field2` for complete binary screens.
The corresponding capped-screen analyzer is
[analyze_pfaffian_binary_screen.py](../tools/analyze_pfaffian_binary_screen.py).
It reports a prospective successful seed separately from an independently
audited tensor certificate. These four-field exclusions concern the listed
finite coefficient spaces; they do not exclude other fields, larger codes,
or linear retracts.
