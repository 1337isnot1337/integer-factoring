# Audit of assignment-dependent global cube pairing

## Result

The replacement in `e233eac` does **not** supply the missing transfer.
Its fixed common-map identity is valid, but choosing that map after the
complete piece assignment is a different identity. The proposed rule
fails on coefficients of the cited Gaussian word construction: 72 of
256 fixed parent assignments disagree with the endpoint.

The [original replacement proposal](two-axis-global-pairing-proposal-rejected.md)
and [original receipt](../evidence/two-axis-global-pairing-proposal-rejected.json)
are preserved as rejected historical work. The current
[checker](../tools/check_two_axis_global_pairing.py) and
[receipt](../evidence/two-axis-global-pairing-repair.json) retain all valid
fixed-map and arithmetic controls and add the actual-family failure.
The complete factoring score remains **1/2 because this proof step fails**;
this is a substantive gap, not merely a wait for independent acceptance.

The Gaussian construction, its profiles and its three-axis transfer are
unchanged. No factoring lower bound, new Lean or independent mathematical
review is claimed. This ordinary source audit used OpenAI assistant support.

## 1. The quantifier change inside the expansion

For each fixed pairing c, linearity and endpoint symmetry give

~~~
sum_I product_s P_c(T_(I_s)) = product_s E_4.            (1)
~~~

Here I is the complete piece assignment to the squares, and P_c maps the
eight physical cube legs to the four fixed parent superlabels. Applying
one **fixed** map to every summand is valid; the positive controls check
all three such maps.

The proposed algorithm first inspects I, chooses the hidden axis c(I)
maximizing the total width, then uses

~~~
sum_I product_s P_(c(I))(T_(I_s)).                      (2)
~~~

Equation (1) does not assert that (2) equals the endpoint. Its c depends
on the summand. Calling this map common across the squares does not make
it common across the expansion. With just one square, a complete piece
assignment is one selected coefficient tensor: the distinction is already
present at B=1.

The source multiaxis transfer's assignment-dependent axis choice does not
have this problem. It changes which factorization of the **same fixed
summand** is used to evaluate that summand; it does not permute the
summand's arguments. The proposed repair changes those arguments.

## 2. Exact counterexample from the actual Gaussian coefficients

Use the existing homogeneous Gaussian seed G from the
[cube construction](cube-homogeneous-gaussian.md). Take word length m=5,
weight one and the two-letter physical alphabet

~~~
Q={00001,00010},       lambda=1,       lambda'=2.
~~~

This is a permitted common constant-weight subalphabet. Identify its two
letters with bits zero and one. For an eight-letter assignment encoded by
mask in {0,...,255}, the actual word entry polynomial is

~~~
F_mask(x,y,z)=G_mask(x,y,z)*G_(mask xor255)(x,y,z).
~~~

The other three sites are empty and contribute one. Its 22 nonzero
coefficient tensors sum exactly to the eight-leg endpoint. The checker
acquires these polynomials from the original Pfaffian implementation,
without importing a precomputed list of counterexample values.

For every coefficient tensor, compute its exact F_2 rank at each of the
three cuts and set r_a=ceil(log_2 rank_a). These are valid exact finite
channel widths. Apply the proposal's rule: choose the first axis maximizing
r_a as the hidden axis and use the other two as the parent axes. Among the
22 pieces the chosen hidden-axis counts are (10,7,5).

For each one of the three pairings held **fixed** before summing all
pieces, the result agrees with E_4 on every one of the 256 parent
assignments. With the assignment-dependent choices, 72 assignments
instead disagree. The first is mask3:

~~~
X=(1,1),   Y=U=V=(0,0),
E_4=0,     adaptive coefficient sum=1.
~~~

Thus the adaptive rule is not another implementation of the fixed tensor
identity, even on coefficients of the actual inherited family. Arbitrary
parent tables can select such an entry, so a false local identity cannot
be hidden by appealing to a scalar contraction.

This finite control targets the claimed exact identity step. It is not an
asymptotic pressure lower bound for every possible regrouping, and does not
claim that these small exact ranks by themselves certify the source's
asymptotic pressure. That separate pressure question cannot repair a
linearity step that has not been proved for the adaptive operation.

## 3. What survives and what a replacement must supply

The previous controls remain true:

- Fixed common permutations preserve the endpoint and the summand sum.
- Removing a largest width total leaves at most two thirds of the total.
- The smaller of two parent totals is at most their average.
- The parent-cut incidence and separate numerical width caps are correct.
- A total-width assignment can satisfy rank feasibility even when a
  uniform per-piece fixed-cut bound would violate endpoint rank.

These are necessary ingredients, not the missing equality between (1)
and (2). The checker now explicitly reports
factoring_exponent_improved=false and rejects promotion of these positive
controls to a complete transfer.

A valid construction must keep a fixed parent embedding with the required
width bound, or prove a new exact adaptive reassignment identity together
with compatible arbitrary vertex-table wiring. Merely applying a common
map separately within each summand of the outer expansion is insufficient.
The original [per-cube audit](two-axis-cube-regrouping.md) and its rejection
controls remain valid and preserved.

Run make verify-two-axis-global-pairing. Its success means the audit,
including the actual coefficient-family counterexample, has passed.
