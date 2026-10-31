# Finite polynomial representatives for every regular local row frame

## Result and scope

Fix sorted row targets and a constant flag over a finite field. Every
regular formal row frame with that constant flag has an equivalent
polynomial representative with explicit bounds on each lower entry.
The representative is unique once the constant frame is fixed.

This is a complete local classification, rather than a fixed-depth jet
approximation. It can substantially reduce the
[nonconstant-frame synthesis problem](nonconstant-binary-row-flags.md).
The two local frames in that problem are independent, so normalization
at one place imposes no extra behavior at the other. The source tensor
still belongs to the same complete global endpoint fiber.

The theorem does not assert that a search of one constant flag covers
the other flags. Nor does it preserve any independently imposed global
degree or other-place condition on the frame itself. No such conditions
are imposed on the independent local frames used here.

## 1. The exact stabilizer of the row conditions

Let O=k[[z]] for a field k, and let

```math
\tau_1\le\cdots\le\tau_d
```

be integer row targets. A matrix A satisfies the row conditions if each
entry of row i has valuation at least tau_i. Put D=diag(z^tau_i). The
regular changes preserving these conditions in both directions form

```math
H=\{P\in\mathrm{GL}_d(O):
       v(P_{ij})\ge\max(0,\tau_i-\tau_j)\}.              \tag{1}
```

Indeed, these are precisely the regular P for which D^-1 P D is
integral. Its determinant is a unit, so its inverse is integral too.
This proves that H is a group and that its elements preserve the row
conditions in both directions. Every regular upper triangular matrix
with unit diagonal entries belongs to H.

Two frames W and W' are equivalent when W'=P W for P in H. They give
exactly the same conditions on every input matrix A.

## 2. The identity constant matrix

**Theorem.** Suppose W is regular and W(0)=I. Its H-equivalence class
has a unique lower unipotent representative L such that

```math
L_{ij}(z)\in z k[z],\qquad
\deg L_{ij}<\tau_i-\tau_j\quad(i>j).                    \tag{2}
```

A gap of zero or one forces that lower entry to be zero.

**Existence.** Since W(0)=I, all its trailing principal minors are
units. Reverse Gaussian row elimination therefore gives W=U L_0,
where U is upper triangular with unit diagonal entries, U(0)=I, and
L_0 is lower unipotent with L_0(0)=I. The matrix U belongs to H.

For each row i, process the earlier columns j in descending order.
Write its current entry as the sum of its terms of degree below
tau_i-tau_j and its remaining tail b(z). Subtract b(z) times row j
from row i. This elementary matrix belongs to H by (1). It removes
the tail because the pivot L_jj is one, and changes only columns less
than j. Thus it cannot disturb a previously completed column.
If the target gap is zero, the entire entry is removed. The constant
matrix remains I throughout. The resulting L satisfies (2).

**Uniqueness.** If L'=P L with both representatives satisfying (2),
then P=L' L^-1 is lower unipotent and belongs to H. Fix a row i and
proceed through j<i in descending order. Once all P_iq for q>j have
vanished, matrix multiplication gives

```math
L'_{ij}-L_{ij}=P_{ij}.
```

The left side is a polynomial of degree below tau_i-tau_j; the right
side is divisible by z^(tau_i-tau_j), by (1). They are zero. This
proves P=I and L'=L.

Only coefficients of W below max(tau)-min(tau) can affect this result.
The reverse elimination and unit inverses preserve the z-adic
filtration, and every retained degree in (2) is smaller than that
cutoff. The zero-gap case is trivial. Thus the theorem gives an exact
finite algorithm on jets of a formal frame.

## 3. Constant flags and exact catalogue sizes

For a general invertible W(0), choose a constant row-flag representative
F and an upper triangular constant matrix P_0 with P_0 W(0)=F. Then
P_0 W F^-1 has constant matrix I. Section 2 gives a representative L F
of the original local row-condition class.

All complete flags suffice. Equal row targets permit additional
constant row operations, so complete flags can overcount in that case.
No completeness claim requires those duplicates to be distinct.

Over a field of order q, the identity-constant representatives number

```math
q^{\sum_{i>j}\max(\tau_i-\tau_j-1,0)}.                  \tag{3}
```

The coefficient variables in row i number
sum_{j<i} max(tau_i-tau_j-1,0). For example:

| Targets | Free coefficients by row | Complete identity-constant cell over F4 |
|---|---|---|
| `(5,8,10,13)` | `0,2,5,13` | `4^20` |
| `(6,9,9,11)` | `0,2,2,6` | `4^10` |
| `(6,10,10,11)` | `0,3,3,4` | `4^10` |

For the latter two targets, complete formal cells can be easier to
enumerate than a uniform depth-two catalogue. Some apparently free
depth-two entries are stabilizer operations, while a few higher powers
in other entries are essential.

## 4. Exact controls

The implementation normalizes arbitrary polynomial representatives
through the required precision. Independent dictionary arithmetic checks
100 random products of a canonical L with noncommuting stabilizer
operations. Literal matrix products check the defining row-condition
equivalence, including intentional violations of a target.

One complete control normalizes all 4,096 frames I+zN over F4, for
strictly lower triangular N and targets `(0,1,1,2)`. Exactly four
canonical representatives occur, each with multiplicity 1,024. Another
enumerates every one of the 256 representatives for `(0,1,2,3)` and
checks that normalization fixes it. These finite controls supplement
the field-independent proof above.

```sh
python3 tools/check_formal_row_flags.py --output /tmp/formal-row-flags-replay.json
```

The source is [formal_row_flags.py](../tools/formal_row_flags.py), with
[independent controls](../tools/check_formal_row_flags.py) and a
[saved record](../evidence/formal-row-flag-controls.json).
