# Polynomial crossings on perfect matchings

## Result and scope

For every `l>=2`, there is a strict polynomial crossing of degree `l`
on the `(2l-1)!!` perfect matchings of `2l` vertices. Its entries have
an explicit common-edge formula. For `l>=3`, the tensor is an exact
reflected sector of the characteristic-two source construction.

Both cuts have uniform block decompositions. The row blocks satisfy
the reciprocal identity `G_r(w)G_r(w^-1)=I`. The opposite-cut profiles
are given by a labeled-cycle recurrence, without a dense tensor.

The full matching alphabet has a **bare-pressure floor `rho>=3/4`**
at every size. It therefore cannot improve the square-root exponent
through this unweighted criterion. This does not exclude physical
restrictions, coherent maps, independent physical heights, or different
dual parameters. The construction and its obstruction are ordinary
proofs; no Lean formalization or improved factoring bound is claimed.

The matching association scheme provides a standard setting for these
matrices; see Srinivasan, [*The perfect matching association scheme*](https://alco.centre-mersenne.org/articles/10.5802/alco.104/).
The entry, reciprocal, and profile identities below are proved directly
over characteristic two and do not use characteristic-zero spectral
decompositions.

## 1. The complete polynomial tensor

Write `M_l` for the perfect matchings of the complete graph on `2l`
vertices, and put `D_l=|M_l|=(2l-1)!!`, including `D_0=1`.
For four matchings `x,y,u,v`, set the entry to zero unless

\[
 x\setminus y=v\setminus u,\qquad
 y\setminus x=u\setminus v.                              \tag{1}
\]

When (1) holds, put

\[
 r=|x\cap y|,\qquad f=|x\cap y\cap u\cap v|,
 \qquad T_l(w)_{xyuv}=w^f(1+w)^{r-f}.                      \tag{2}
\]

The nonshared edges must cross unchanged. On the vertices covered by
the shared edges, the common matching can be replaced by another common
matching. Formula (2) assigns the coefficient for that replacement.

At `w=1`, an entry can survive only when `f=r`. Its common matching is
then unchanged, so (1) forces `u=y` and `v=x`. Conversely that entry is
one. Thus

\[
 T_l(1)_{xyuv}=[x=v][y=u].                                 \tag{3}
\]

The degree is at most `l`. For `l=2`, this is the three-letter Brauer
pencil after replacing the parameter by its square. The larger matching
alphabets contain additional common-edge blocks.

## 2. A lacunary reflected source realization

Take `l>=3`, `Q=2^(2l+1)`, and vertex labels `1,...,2l`. To edge
`e={i,j}` associate the odd state part and reflected parameter part

\[
 b_e=Q-1-2^i-2^j,\qquad a_e=Q-b_e=1+2^i+2^j.              \tag{4}
\]

For a perfect matching `M`, use `b_M=prod_(e in M)t_(b_e)` and the
reflected dual `a_M=prod_(e in M)t_(a_e)`. Every selected state and
dual has, respectively, the same weight

\[
 H=(l-1)Q-l+2,\qquad W=Q+l-2,\qquad H+W=lQ.               \tag{5}
\]

Here `W<=H`. The reflected vacuum-pairing proof applies at mode `-H`
even when `W!=H`; it gives `[1]a_M[-H]b_N=[M=N]`.
The source angle prefactor for a parameter of weight `W` and physical
states of weight `H` is `z^(-W)(1+z)^(H+W)`. This is the product of
the gauges in the [neutral Wick kernel](neutral-reflected-crossings.md).

The part `b_e` has every binary digit below `Q` except the two digits
in `e`. The three kernel blocks `AB`, `AY`, and `BU` are therefore
diagonal on the edge index. Their nonzero entries, together with `YU`,
are

\[
 AB=YU=1,\qquad AY=z^{-Q}(1+z)^{a_e},\qquad BU=(1+z)^{b_e}.
 \tag{6}
\]

The `BY` block vanishes: two edges cover at most four of the `2l>=6`
vertex digits, so `(b_e-1)&b_f` has a nonzero digit. In a complete
Wick matching the number of `AU` edges equals the number of `BY` edges,
because each of the four groups contains `l` fields. Thus no `AU` edge
can contribute either.

If an edge occurs on all four physical boundaries, its two Wick
matchings give

\[
 1+z^{-Q}(1+z)^{a_e+b_e}=z^{-Q}.                           \tag{7}
\]

If `k` common edges are removed and replaced on the same vertex set,
the sums of their `a_e` parts agree for every matching of that set.
The product of the corresponding `AY` and `BU` edges is therefore

\[
 z^{-kQ}(1+z)^{kQ}=(1+z^{-Q})^k.                         \tag{8}
\]

Equations (6)--(8), with `w=z^-Q`, give exactly (1)--(2). All large
binary powers have disappeared from the final polynomial tensor.
The source realization has `log H=O(l)` while `log D_l=Theta(l log l)`;
this fact alone does not reduce either cut rank.

The restriction `l>=3` is necessary for this particular source collapse.
At `l=2`, complementary edges can have a nonzero `BY` contraction.
The checker retains a literal counterexample to extending this collapse
to `l=2`. Formula (2) itself remains a strict crossing there.

## 3. The row cut and an exact reciprocal identity

In the cut `(x,y)|(u,v)`, fix the two ordered nonshared matchings.
The remaining `2r` vertices support one block

\[
 G_r(w)_{P,Q}=w^{|P\cap Q|}(1+w)^{r-|P\cap Q|},\qquad
 P,Q\in M_r.                                             \tag{9}
\]

Its columns are indexed with the nonshared core reversed. Let `F_m`
be the number of ordered pairs of disjoint perfect matchings on `2m`
vertices. Inclusion and exclusion give

\[
 F_m=D_m\sum_{j=0}^m(-1)^j\binom mjD_{m-j}.               \tag{10}
\]

The multiplicity of the block `G_r` in the row cut is

\[
 \binom{2l}{2r}F_{l-r}.                                  \tag{11}
\]

### Reciprocal theorem

Over every field of characteristic two,

\[
 \boxed{G_r(w)G_r(w^{-1})=I.}                            \tag{12}
\]

**Proof.** Put `c=1+w` and fix row and column matchings `P,Q`. After
factoring out `c^(2r)w^-r`, the matrix product is the partition sum
over intermediate perfect matchings `L`, with edge weights `w/c`
on `P\Q`, `1/c` on `Q\P`, and `w/c^2` on `P cap Q`.
All other edge weights are one.

Expand each exceptional edge weight as `1+(weight-1)`. Every selected
partial matching has an odd number `D_(r-k)` of completions, and hence
has completion coefficient one in characteristic two. The sum becomes
the matching polynomial of the union multigraph `P union Q`.

A common edge is an isolated two-vertex component. Its contribution is
`1+(w/c^2-1)=w/c^2`. Every other component is an alternating even
cycle. Its alternating edge increments are `a=1/c` and `b=w/c`, so
`a+b=1`. The matching polynomial of a cycle of length `2k` is

\[
 \operatorname{tr}\left[\left(
 \begin{pmatrix}1&a\\1&0\end{pmatrix}
 \begin{pmatrix}1&b\\1&0\end{pmatrix}\right)^k\right].      \tag{13}
\]

The matrix inside (13) has trace `1+a+b=0`. By its two-dimensional
characteristic equation, every power has trace zero in characteristic
two. Thus the product entry vanishes if `P!=Q`. When `P=Q`, its `r`
isolated common edges contribute `(w/c^2)^r`, canceling the prefactor.
This proves (12) as a rational identity, and hence as a Laurent
polynomial identity. **End of proof.**

Every zero Smith order `a` of `G_r` therefore satisfies `0<=a<=r`,
and its infinity orders are the negatives of its zero orders. The upper
bound follows because every entry of the inverse has pole order at most
`r`. A further useful exact expansion is

\[
 G_r(w)=\sum_{j=0}^r(1+w)^{r-j}I_jI_j^{\mathsf T},        \tag{14}
\]

where `I_j` is the incidence matrix between perfect matchings and
their `j`-edge partial matchings. Its `(P,Q)` coefficient is
`sum_j binom(|P cap Q|,j)c^(r-j)=c^(r-f)(c+1)^f`.

## 4. The opposite cut has explicit cycle profiles

In the cut `(x,u)|(y,v)`, fix the union multigraph `x union u`. Write
`c` for its number of common edges and `k_1,...,k_s` for the half-lengths
of its nontrivial alternating cycles. A common edge is forced. On each
longer cycle the two colors can either be preserved or exchanged.
The block is

\[
 w^c\bigotimes_{i=1}^s
 \begin{pmatrix}(1+w)^{k_i}&1\\1&(1+w)^{k_i}\end{pmatrix}.
 \tag{15}
\]

At zero, each two-dimensional factor has orders
`0` and `2^(v_2(k_i)+1)`: one entry is a unit, and the determinant is
`(1+w)^(2k_i)+1`. The block orders are consequently

\[
 c+\sum_i\epsilon_i2^{v_2(k_i)+1},\qquad \epsilon_i\in\{0,1\}.
 \tag{16}
\]

At infinity the leading coefficient of each factor is the identity.
All orders of (15) are `-c-sum_i k_i=-l`. Thus the entire opposite
cut has `D_l^2` copies of order `-l` at infinity.

Let `C_l(t)` count the zero orders, so that the coefficient of `t^a`
is the multiplicity of order `a`. Labeled cycle counting gives

\[
 C_l(t)=(2l)![x^l]\exp\left(\frac{tx}{2}
 +\sum_{k\ge2}\frac{(1+t^{2^{v_2(k)+1}})x^k}{4k}\right).
 \tag{17}
\]

Equivalently, with `C_0=1`, all calculations use the integer recurrence

\[
 C_l(t)=(2l-1)tC_{l-1}(t)
 +\sum_{k=2}^l\frac{(2l-1)!}{2(2l-2k)!}
       (1+t^{2^{v_2(k)+1}})C_{l-k}(t).                    \tag{18}
\]

In particular `C_l(1)=D_l^2`, as required for a complete profile.

## 5. A uniform bare-pressure floor

Use `s>1` and the unweighted pressure convention

\[
 B(s)=\left(\sum_i s^{-\sigma_i}\right)
      \left(\sum_j s^{-\tau_j}\right),\qquad
 \rho=\frac{\log B(s)}{4\log D_l}.                        \tag{19}
\]

First use row infinity and opposite-cut zero. Every disjoint ordered
pair supplies a scalar `G_0=1` row block. Hence the first factor is
at least `F_l`. A union multigraph of disjoint matchings has at most
`floor(l/2)` cycles and therefore at most `2^floor(l/2)` ordered
colorings. Each such union gives at least one order-zero channel in
(16). The second factor is at least `F_l/2^floor(l/2)`. Thus

\[
 B(s)\ge F_l^2/2^{\lfloor l/2\rfloor}.                  \tag{20}
\]

Fixing one matching, a random other matching includes a specified edge
with probability `1/(2l-1)`. The union bound yields

\[
 F_l\ge D_l^2\frac{l-1}{2l-1}\ge D_l^2/3\quad(l>=2).    \tag{21}
\]

For `l>=4`, `D_l>=9*2^floor(l/2)`: check `D_4=105>=36` and use
`D_(l+1)/D_l=2l+1`. Equations (20)--(21) then imply `B(s)>=D_l^3`.
For `l=3`, `F_3=120` gives `B(s)>=120^2/2>15^3` directly.
For `l=2`, the complete product is

\[
 (8+s^2)\,3(1+s^{-2}+s^{-4}),                            \tag{22}
\]

whose constant terms already sum to `27=D_2^3`.

For the other orientation, the opposite infinity factor is `D_l^2 s^l`.
Every row-zero order is at most `l`, by (12), so the row-zero factor
is at least `D_l^2 s^-l`. This orientation has `B(s)>=D_l^4`.
Both orientations therefore give

\[
 \boxed{\rho\ge3/4\qquad(l>=2).}                         \tag{23}
\]

More strongly, (20)--(21) give
`rho>=1-(floor(l/2)log2+2log3)/(4log D_l)`, which tends to one.
The floor concerns the complete alphabet with no physical heights.
It does not cover arbitrary codes or weighted physical distributions.

## 6. Exact finite controls

The [constructor](../tools/perfect_matching_crossing.py) uses integer
polynomial bits and the explicit entry formula. The
[independent checker](../tools/check_perfect_matching_crossings.py)
compares every one of the `15^4=50,625` entries with the complete
neutral Wick compiler at `Q=128`, state weight 255, and dual weight 129.
It checks both cut profiles by binary block-Toeplitz ranks, independently
of the native weighted Smith oracle. It also checks (12) and (14)
coefficientwise through `r=4`.

The exact fifteen-letter profiles are:

| Cut and place | Order multiplicities |
| --- | --- |
| Row zero | `0:124, 1:96, 2:4, 3:1` |
| Row infinity | `-3:1, -2:4, -1:96, 0:124` |
| Opposite zero | `0:60, 1:45, 2:60, 3:15, 5:45` |
| Opposite infinity | `-3:225` |

```sh
python3 tools/check_perfect_matching_crossings.py
```

The [evidence](../evidence/perfect-matching-crossings.json) retains the
exact counts and profiles. Separate floating evaluations of the complete
profiles give bare minima about `0.876858` at three letters and
`0.963738` at fifteen letters. Those displays are not global weighted
optimization certificates; the all-scale unweighted floor is (23).

This construction, proof, and implementation used substantial OpenAI
assistant support.
