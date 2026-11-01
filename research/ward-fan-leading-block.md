# An exact vanishing block improves the first fan upper profile

For the explicit \(P=L=10,k=3\) duplicated fan with 84 selected
labels, one \(128\times128\) block vanishes on the entire boundary
\([y^3]G(x,y)\) in both rank cuts. The five coefficient caps at
\((a,b)=(3,3),\ldots,(7,3)\) improve from 7,056 to 6,928 in
each cut. Combining these exact restrictions with the previously
certified complete profile reduces its integer rank-sum upper
bound from
\[
 192,340\quad\hbox{to}\quad191,700.
\]
This is a verified finite upper-bound improvement for this code
and output map. The resulting cost remains far above the seed
criterion and does not improve the factoring score.

## The exact restricted block

The code, preimages, and coordinate outputs are the sparse
construction in [ward-fan-polynomial-queries.md](ward-fan-polynomial-queries.md).
The five type columns are ordered
\[
 e_1,\ e_1+e_2,\ e_2,\ e_2+e_3,\ e_3,
\]
with each column repeated twice. Every code label has a type
syndrome \(\sigma\in\mathbf F_2^5\), recording the groups whose
two physical sites are singly occupied.

Take the pair syndrome \(\tau=7\), meaning the first three type
bits are one. The ordered logical pair set is
\[
 \mathcal P_7=\{(a,b):\sigma_a+\sigma_b=7\}.
\]
It has exactly 128 elements. In the row cut, rows are indexed by
\((x,y)\in\mathcal P_7\), columns by \((v,u)\in\mathcal P_7\), and
the entry is \(G_{xyuv}\). In the column cut, rows are
\((x,u)\in\mathcal P_7\), columns are \((v,y)\in\mathcal P_7\).
The evidence stores the full ordered pair list and all 84 input
vectors, homogeneous preimages, and coordinate outputs.

Write
\[
 G(x,y)=y^3\Gamma(x,y),
 \qquad p=1+x,\quad q=1+y,\quad r=1.
\]
The Ward quotient is evaluated using the exact identity
\[
 Q\mathscr M Q
 =y^3\,W(y)\mathscr M Q,
\]
with the sparse homogeneous preimage on one state leg and the
input image on the other. The
[two-state physical oracle](ward-fan-polynomial-queries.md)
evaluates each required physical polynomial without expanding
the form subsets.

At \(y=0\), the local parameters are \(q=r=1,\eta=0\). The
quotient is a univariate polynomial of degree at most seven.
The computation checks every one of the 16,384 entries in
each cut and obtains
\[
 \Gamma_{\mathcal P_7,\mathcal P_7}(x,0)=0.
 \tag{1}
\]
This is an exact polynomial identity for the specified finite
block, not a sampled rank at a parameter value.

The binary joint-order theorem gives order at least six for
\(G\). Thus the only potentially nonzero boundary coefficients
are \([x^a y^3]G\), \(3\le a\le7\), and all five are zero on
these blocks. Each corresponding actual block rank is zero in
both cuts.

## The complete coefficient bounds

Each logical cut has dimension \(84^2=7056\) and is block
diagonal in pair syndrome. The sum of the dimensions of the
other syndrome blocks is \(7056-128=6928\). Hence (1) gives
\[
 \operatorname{rank}G[a,3]\le6928
 \qquad (3\le a\le7)
\]
in each cut, using only the computed block and the dimension
bounds on all uncomputed blocks.

The previous complete profile used the cap 7,056 in both cuts
for each of these five coefficients. All remaining bounds are
retained. The integer rank-sum upper bound therefore decreases
by exactly \(5\cdot128=640\). The full updated coefficient
profile is included in the evidence. The preceding rectangular
bounds allow arbitrary homogeneous syndrome-preserving output
maps, so they apply to the coordinate outputs explicitly
specified here.

## Independent controls

The checker first constructs the input vectors by creation
forms. It independently applies the literal product
\(Q_1Q_2Q_3\), including every local diagonal and flip term, to
each sparse preimage and checks that the result is the stated
input vector. It also checks all homogeneous grades, the exact
left inverse \(\Pi I=I\), and the isotropic input Gram matrix.

For 32 selected logical entries, it separately contracts the
original Pfaffian kernel with both full input vectors. This
calculation uses the truncated ring
\(\mathbf F_2[x,y]/(x^8,y^4)\), not the Ward factorization or
the two-state recurrence. The result agrees with the entire
boundary polynomial \(y^3\Gamma(x,0)\).

There are also 32 independent crossing-endpoint controls at
\((p,q,r)=(0,0,1)\). Sixteen entries equal one, and sixteen
equal zero. Both the Ward evaluator and the direct original
kernel contraction give these values. These positive controls
verify that the zero boundary block is not caused by an
identically zero evaluator or a failed input/output contraction.

The block computation makes 729,088 physical query requests.
Of these, 728,448 have a zero collapsed local pair. Exact caching
reduces the remaining work to six distinct physical polynomial
evaluations. No full \(84^2\times84^2\) coefficient matrix is
assembled.

## Reproduction

Run

~~~sh
python3 tools/check_ward_fan_leading_block.py --syndrome 7 --output /tmp/ward-fan-leading-block.json
~~~

The [checker](../tools/check_ward_fan_leading_block.py) generates
both complete restricted polynomial blocks, checks their
coefficient ranks by two elimination orders, and updates the
frozen complete upper profile. The
[evidence](../evidence/ward-fan-leading-block.json) includes all
source hashes and the hash of the prior upper profile.
