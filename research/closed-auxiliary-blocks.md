# Closed auxiliary blocks with literal crossing endpoints

## Result and scope

A two-state polynomial matrix at each physical bit gives a scalar crossing
on every fixed block of binary words. The auxiliary state begins and ends
at the same fixed basis vector. All internal states are contracted before
computing either physical cut rank.

The construction admits arbitrary local off-diagonal tensors and an
arbitrary common physical word code. It supplies exact new crossing
families and a reproducible finite search. No favorable uniform pressure
bound or improved factoring exponent is established here.

This differs from the multiplicative
[ordered physical gauges](ordered-physical-gauge-crossings.md): the local
matrix below is an additive endpoint-zero perturbation of a scalar matrix.
Neither construction's finite negative controls exclude the other family.

## 1. The local matrix and fixed boundary

For a physical bit tuple `a=(x,y,u,v)`, define

\[
I_a=[x=u][y=v],\qquad E_a=[x=v][y=u],\qquad
P_a(z)=I_a+z(I_a+E_a)
\]

over a characteristic-two field. Let `B_a,C_a` be arbitrary fixed field
elements, and set

\[
K_a(z)=P_a(z)I_2+(1+z)
\begin{pmatrix}0&B_a\\C_a&0\end{pmatrix}.                 \tag{1}
\]

For a block of length `r`, use the scalar polynomial

\[
T_r(z;x,y,u,v)
 =e_0^T K_{a_1}(z)\cdots K_{a_r}(z)e_0,
\qquad e_0=(1,0)^T.                                    \tag{2}
\]

Each factor has degree at most one, so (2) has degree at most `r`.
At `z=1`, equation (1) gives `K_a(1)=E_a I_2`. Therefore

\[
T_r(1;x,y,u,v)
 =\prod_{j=1}^r[x_j=v_j][y_j=u_j]
 =[x=v][y=u].                                           \tag{3}
\]

The boundary is fixed independently of the physical inputs. No auxiliary
entry is selected after evaluating a rank or an endpoint.

Let `C` be any set of distinct binary words of length `r`, used on all
four physical legs. Restricting (2) to `C^4` keeps (3) literally. The
physical alphabet is `q=|C|`; internal auxiliary states do not contribute
to its entropy. Conversely, no cut-rank saving is inferred merely from
the auxiliary dimension. Both cut matrices are formed from the fully
contracted scalar tensor and their profiles are recomputed.

For example, the two-position block is

\[
P_{a_1}P_{a_2}+(1+z)^2B_{a_1}C_{a_2}.
\]

The three-position block includes all three ways to leave and return to
auxiliary state zero. Those terms can cancel over the coefficient field;
the constructor retains their literal coefficients.

## 2. Repeated blocks and the effective interface

Once `r`, the local tensors, the coefficient field, and the physical code
are fixed, (2) is an ordinary finite scalar seed. Its tensor powers use
a reset boundary between blocks. Their degree is at most `rm`, and they
have at most `rm+1` scalar coefficient pieces at block length `m`.

Exact construction of the fixed seed is finite. Constructing all entries
and rank factors of its `m`-fold tensor power has the same exponential
table cost as the existing
[block Smith crossing interface](block-smith-crossings.md). A proposed
improvement must still satisfy the complete two-cut pressure criterion,
the physical type cost, and the effective factoring transfer. There is no
claim here that fixed-block performance describes a growing auxiliary
construction without resets.

For binary `B,C`, transposing every local matrix exchanges the two masks.
The identity

\[
e_0^T(K_1\cdots K_r)e_0
 =e_0^T(K_r^T\cdots K_1^T)e_0
\]

then exchanges `B,C` and reverses each physical word. Thus unordered mask
pairs suffice when the physical alphabet is the full word cube. For a
restricted code, this reduction also reverses the code; it must not be
used without that relabeling or reversal invariance.

## 3. Exact finite controls

The [constructor](../tools/closed_auxiliary_blocks.py) uses binary
polynomial bitmasks. The independent
[checker](../tools/check_closed_auxiliary_blocks.py) sums every closed
auxiliary path with literal coefficient convolution, independently of
the constructor's two-coordinate recurrence. It checks:

- 9,329 complete path entries and their literal endpoints;
- 52 independent block-Toeplitz rank identities;
- ten complete native profile pairs against separate exact polynomial
  elimination; and
- 8,192 entries of the mask-exchange and word-reversal identity.

The retained cases include full alphabets at lengths two and three and a
five-word code at length four. The
[control record](../evidence/closed-auxiliary-block-controls.json) stores
all finite profiles and source hashes.

```sh
python3 tools/check_closed_auxiliary_blocks.py > /tmp/closed-auxiliary-controls.json
cmp /tmp/closed-auxiliary-controls.json evidence/closed-auxiliary-block-controls.json
```

The separate [search driver](../tools/search_closed_auxiliary_blocks.py)
enumerates a stated finite mask catalogue, computes complete unweighted
profiles, and explores physical heights and types numerically. Its
evidence distinguishes a completed coefficient catalogue from an
unfinished weighted schedule. Floating minima are exploratory; they are
not global lower bounds or rational pressure certificates.

No Lean source is added or changed.
