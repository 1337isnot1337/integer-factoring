# A fast row-cut action for the twelve-site coefficient algebra

## Result and scope

The naive word-pair expansion of the twelve-site cut module is
infeasible: it has `924^2=853776` word pairs and millions of local
entries. This note records a different construction of the **row**
action of one original-pencil monomial. It uses the signature blocks
from the coefficient-table constructor and a Krawtchouk/Walsh-Hadamard
convolution on each block. The action reproduces the dense row ranks at
lengths six and eight, and reproduces the published ten-site
idempotent-`0` module dimension `34` with an idempotent check.

At twelve sites the first central idempotent is `e_0=4 C_9`. The fast
action is idempotent on sampled vectors, its first sampled image has
ordinary `F5` matrix rank `88` and matrix trace `0`, and a randomized
image span gives a **lower bound** for the row rank of `e_0`. A
`140`-column diagnostic run reached rank `140`; the default replay in
the checker uses eight columns and certifies the corresponding smaller
lower bound. This is not an exact rank, not a column-cut result, and not
a change to the complete bound.

The earlier `apply2` and intermediate `/tmp/action_correct.py`
prototypes are superseded. In particular, the previously reported
rank lower bound `98` came from the intermediate action and is
withdrawn.

## 1. Block action

Let `q` be the number of paths, and let `F[p,w]` be the path-entry
matrix on balanced words. The coefficient tensor is a sum over
signatures `(a,d)`. For a signature with `s=|d|` and `t=s/2`, the pairs
are indexed by `t`-subsets `e` of `d`:

```math
x_e=a\cup e,\qquad y_e=a\cup(d\setminus e).
```

For an input path-pair matrix `V`, put

```math
W[x,y]=\sum_{p_4,p_3}F[p_3,x]V[p_4,p_3]F[p_4,y].       \tag{1}
```

The local monomial block depends only on `x_i\oplus y_h`. If
`b=weight(x_i\oplus y_h)=2r`, its degree-`j` local coefficient is

```math
K_r=\binom{2r}{j},\qquad j=index-s+2r,
```

with the sign `(-1)^j` and the convention that the coefficient is zero
unless `0<=j<=s-2r`. The row action is

```math
(LW)[e]=\sum_f K(e\oplus f)W[f],                        \tag{2}
```

the XOR convolution on the `s`-bit signature group. It is evaluated by
a Walsh-Hadamard transform, not by a pair-by-pair matrix. The output
pair contribution is

```math
\operatorname{out}[p_1,p_2]
 =\sum_e (LW)[e]F[p_1,x_e]F[p_2,y_e].                  \tag{3}
```

Only signatures with `s>=index` contribute. At `index=9` this leaves
the signature weights `10` and `12`, which is why the twelve-site
action is much smaller than the full word-pair expansion.

## 2. Controls

The checker compares the fast action with the dense row operator from
the original-pencil constructor. The rank controls are

```text
length 6, indices 0..6:  25,25,25,20,10,1,1
length 8, indices 0..8:  169,169,169,112,33,56,21,1,1
length 10, e_0:          34 (idempotent, and matrix rank of a generic image)
```

The ten-site vector and basis are read from the published algebra model.
The length-ten image is checked to satisfy `T^2=T`, and its ordinary
matrix rank is the published primary dimension `34`.

At twelve sites the same construction gives, for `e_0=4 C_9`,

```text
first image: ordinary F5 matrix rank 88, matrix trace 0
sampled image span: rank at least 8 (default replay)
diagnostic 140-column run: rank 140 (not part of the default replay)
```

The published row trace of `e_0` is zero, so an idempotent row operator
has rank divisible by five. The sampled lower bound is therefore
consistent with, but does not prove, an exact rank. The next target is
the **column** action: the present formula contracts the `(p3,p4)`
pair and does not compute the column pairing.

## 3. Consequences

- The row action of the twelve-site coefficient algebra is no longer
  blocked by the word-pair expansion.
- The first primary is not a small row component: its row rank is at
  least `140`, far above the `5` lower bound from a single nonzero
  entry.
- No column rank, exact primary rank, two-cut product, or factoring
  exponent is claimed.

## Reproduction

```sh
python3 tools/check_bounded_spin_twelve_row_action.py --samples 8
make verify-bounded-spin-twelve-row-action
```

The checker uses the published path-entry and coefficient-table
constructors, the ten-site algebra model, and the twelve-site idempotent
certificate. It does not build a dense twelve-site cut operator.
