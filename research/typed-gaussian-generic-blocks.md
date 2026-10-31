# Generic mode blocks obstruct every fixed-occupation coefficient family

## Result and scope

Restrict powers of the four-letter Gaussian crossing to a common word
type. If the numbers of plus, minus, and doubly occupied letters stay
fixed while the word length grows, every decomposition using arbitrary
linear combinations of the polynomial coefficients has rank cost
\(\Omega(D^2)\), where \(D\) is the number of typed words.

The proof identifies many disjoint matrix blocks in both flattenings.
Every block contains the same pure monomial. A coefficient recombination
that supplies this monomial must be nonzero in every block, so at least
one piece of every crossing decomposition has large rank in both cuts.
The argument gives an explicit finite lower bound and does not require
the coefficient family to be linearly independent.

This result concerns the coefficient span of the specific Gaussian seed
from the [Schur digit reduction](schur-gaussian-digit-reduction.md),
after the displayed type restriction. It does not exclude other local
seeds, non-Gaussian changes, or additional operations that leave this
coefficient span. It establishes no improved factoring exponent.

The derivation and checks were developed with substantial OpenAI
assistant support.

## 1. Typed words and the central monomial

The four physical letters are
\[
0=\varnothing,\qquad 1=+,\qquad 2=-,\qquad 3=\{+,-\}.
\]
Let \(a,b,c,d\) be their multiplicities in a word of length

~~~math
m=a+b+c+d,\qquad s=b+c+d,\qquad
D=\frac{m!}{a!\,b!\,c!\,d!}.                            \tag{1}
~~~

Write \(\mathcal L\) for the set of these \(D\) words. The local
Gaussian kernel is

~~~math
K(w)=
\begin{pmatrix}
0&1+w&1+w&1\\
1+w&0&w^2&1+w\\
1+w&w^2&0&w^2(1+w)\\
1&1+w&w^2(1+w)&0
\end{pmatrix},                                         \tag{2}
~~~

in boundary order \(B,X,V,Y\). The physical labels
\((x,y,u,v)\) select boundary states
\((\mathrm{flip}(u),\mathrm{flip}(x),v,y)\);
the flip exchanges the two species. The selected determinant gives
the tensor entry.

Take \(m\) copies with one common variable and restrict each boundary
to \(\mathcal L\). Denote the resulting polynomial tensor by

~~~math
T_m(w)=\sum_j A_jw^j,\qquad
T_m(1)=E,\qquad
E(x,y,u,v)=\mathbf1_{x=v}\mathbf1_{y=u}.                  \tag{3}
~~~

The row flattening uses \((x,y)\mid(u,v)\); the column flattening
uses \((x,u)\mid(y,v)\).

Assume \(m\ge2s\). Two words are called disjoint here when their
nonvacuum supports are disjoint. For any such \(x,y\), the crossing
entry is

~~~math
T_m(w;x,y,y,x)=w^J,\qquad J=2(b+c)+4d.                   \tag{4}
~~~

Indeed, a single plus or minus occupation of \(x\) uses the
\(X,V\) edge \(w^2\), and a double occupation uses the two-by-two
\(X,V\) determinant \(w^4\). The \(B,Y\) contractions of \(y\)
have weight one. Disjoint supports prevent the two pairs from sharing
a mode, so the factors multiply to (4).

## 2. Disjoint row blocks

Restrict both row and column labels to pairs of disjoint words.
At a given mode, at most one boundary on either side of the chosen cut
is occupied.

The diagonal entries of (2) vanish. Thus a selected fermion cannot
contract within its own boundary. Disjoint row supports also prevent
an internal contraction between the two row boundaries, and the same
holds on the column side. Every contraction must cross the cut.
The per-mode species pattern therefore determines which restricted
rows and columns can meet.

### Row cut

For a row \((x,y)\), define three disjoint mode sets:

~~~math
\begin{aligned}
P&=\{\text{plus-only modes of }x\}
       \cup\{\text{minus-only modes of }y\},\\
N&=\{\text{minus-only modes of }x\}
       \cup\{\text{plus-only modes of }y\},\\
B&=\{\text{double modes of }x\}
       \cup\{\text{double modes of }y\}.
\end{aligned}                                           \tag{5}
~~~

They have sizes \(b+c,b+c,2d\). The compatible column signature
is the one obtained from \((u,v)=(y,x)\). Different signatures
cannot interact. Consequently the restricted flattening is block
diagonal, with

~~~math
N_R=\frac{(m)_{2s}}{(b+c)!^2(2d)!}                       \tag{6}
~~~

blocks. Here \((m)_k=m!/(m-k)!\). Each block has

~~~math
\binom{b+c}{b}^{\!2}\binom{2d}{d}
~~~

rows and the same number of columns.

### Column cut

For a row \((x,u)\), the corresponding sets contain the plus-only
positions of both words, the minus-only positions of both words, and
their double positions. Their sizes are \(2b,2c,2d\). The compatible
column signature is obtained from \((y,v)=(u,x)\). This gives

~~~math
N_C=\frac{(m)_{2s}}{(2b)!(2c)!(2d)!}                     \tag{7}
~~~

blocks, each of size

~~~math
\binom{2b}{b}\binom{2c}{c}\binom{2d}{d}.
~~~

These are support statements for the whole polynomial tensor. They
remain true for every linear combination of its coefficients.

## 3. Arbitrary recombinations retain a large piece

Let \(\mathcal V=\mathrm{span}\{A_j\}\), over
\(\mathbb F_2\) or a field extension. Equation (4) defines a
linear functional

~~~math
\lambda\!\left(\sum_j t_jA_j\right)=t_J.                 \tag{8}
~~~

This is well-defined even if the coefficient family is dependent:
evaluation at any one of the entries in (4) extracts \(t_J\).
Moreover, every generic block in both flattenings contains an entry
with exactly the polynomial \(w^J\).

If \(B\in\mathcal V\) and \(\lambda(B)\ne0\), every one of
the blocks in (6)--(7) is nonzero. Therefore

~~~math
\mathrm{rank}_R B\ge N_R,\qquad
\mathrm{rank}_C B\ge N_C.                         \tag{9}
~~~

Now let \(E=\sum_{\ell}B_\ell\), with every
\(B_\ell\in\mathcal V\). Since \(\lambda(E)=1\), some
\(\lambda(B_\ell)\) is nonzero. Thus every such decomposition obeys

~~~math
\boxed{\quad
\sum_\ell
 \sqrt{\mathrm{rank}_R B_\ell\,
       \mathrm{rank}_C B_\ell}
\ \ge\ \sqrt{N_RN_C}.
\quad}                                                  \tag{10}
~~~

No positivity or absence of cancellation among the pieces is assumed.
The proof uses only the one entry functional and the block supports.

## 4. Consequence for fixed occupations

Fix \(b,c,d\), with \(s>0\), and let \(a=m-s\) grow.
Then

~~~math
D=\frac{(m)_s}{b!\,c!\,d!},
\qquad
\sqrt{N_RN_C}
 =\frac{(m)_{2s}}
 {(b+c)!\,(2d)!\sqrt{(2b)!\,(2c)!}}.                    \tag{11}
~~~

In particular,

~~~math
\frac{\sqrt{N_RN_C}}{D^2}
\longrightarrow
\frac{(b!c!d!)^2}
 {(b+c)!\,(2d)!\sqrt{(2b)!\,(2c)!}}
>0.                                                    \tag{12}
~~~

The crossing cost is therefore \(\Omega(D^2)\). In the convention
\(\mathrm{cost}=D^{2\rho+o(1)}\), every construction using this
fixed-occupation coefficient span has \(\rho\ge1\).

### One plus and one minus

For \(b=c=1,d=0\), \(D=m(m-1)\) and
\(N_R=N_C=(m)_4/4\). If \(m\ge7\), then

~~~math
\frac{(m)_4}{4}\ge\frac{5}{42}D^2>D^{10/7}>D^{\sqrt2}.   \tag{13}
~~~

The first inequality follows because
\((m-2)(m-3)/(4m(m-1))\) increases for \(m\ge7\).
For the second, use \(D\ge42\) and the exact integer comparison
\(5^7>42^3\). Finally \(10/7>\sqrt2\).
Thus every coefficient-span recombination in this family already
misses the post-half threshold at each \(m\ge7\).

The native finite rank probes show additional structure, but exact
formulas for all individual coefficient ranks are unnecessary for
this conclusion.

For lengths \(m=3,4,5,6\), the checker enumerates every element of
the five-dimensional coefficient span over \(\mathbb F_2\). It then
assigns each element the integer weight
\(\lfloor\sqrt{\mathrm{rank}_R B\mathrm{rank}_C B}\rfloor\)
and computes the shortest path to the crossing in the XOR graph.
Every true decomposition cost is at least this integer path cost:

| Word length \(m\) | Typed dimension \(D\) | Certified cost lower bound |
|---:|---:|---:|
| 3 | 6 | 26 |
| 4 | 12 | 90 |
| 5 | 20 | 230 |
| 6 | 30 | 483 |

Each displayed lower bound exceeds \(D^{17/12}>D^{\sqrt2}\),
checked by integer powers. These finite enumerations concern
\(\mathbb F_2\) coefficient recombinations. The generic-block
argument in Sections 2--3 also applies over field extensions; no
extension-field enumeration is asserted for the four small cases.

## 5. Growing occupations remain a separate regime

The finite bound (10) still applies when \(b,c,d\) grow with \(m\),
provided \(2s\le m\). The positive constant in (12), however, then
depends on \(m\), so the fixed-occupation conclusion cannot be
transferred unchanged.

For example, suppose the type proportions tend to
\((p_0,p_1,p_2,p_3)\), with \(p_0\ge1/2\), and write \(H\)
for entropy with natural logarithms. The multinomial counts give

~~~math
\begin{aligned}
H_R&=H(2p_0-1,\ p_1+p_2,\ p_1+p_2,\ 2p_3),\\
H_C&=H(2p_0-1,\ 2p_1,\ 2p_2,\ 2p_3).
\end{aligned}                                           \tag{14}
~~~

If \(H(p_0,p_1,p_2,p_3)>0\), (10) yields the lower bound

~~~math
\rho\ge
\frac{H_R+H_C}{4H(p_0,p_1,p_2,p_3)}-o(1).               \tag{15}
~~~

This need not exceed the post-half threshold. In particular, the
generic-block argument leaves room for families with positive
occupation densities and for constructions outside the coefficient
span considered here.

## 6. Exact reproduction

Run

~~~bash
python tools/check_typed_gaussian_generic_blocks.py \
  --output evidence/typed-gaussian-generic-blocks.json
~~~

The checker validates the generic signatures, pure monomial entries,
finite block counts, and the stated one-plus/one-minus threshold. Its
five generic controls include unbalanced types and double occupations.
The small-span shortest paths use exact integer floor weights and
integer threshold comparisons. The fixed-occupation asymptotic
assertion follows from Sections 2--4, independently of those finite
enumerations.
