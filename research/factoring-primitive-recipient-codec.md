# Arbitrary primitive recipient sets have a sharp full-range threshold

## Result and scope

The Fourier-coset assumption in
[the joint-channel theorem](factoring-network-structure.md) is unnecessary
for recipient indices of full multiplicative order. For the source
parameters, **every** set of `k<h/4` odd recipient indices has the full
joint message range. A deterministic polynomial-time algorithm constructs
bounded input digits for any prescribed messages. At `k>=h/4`, the range
cannot be full, by cardinality alone.

This is an exact local theorem for the arithmetic network, not a new
factoring bound. It closes the possibility that a different selection of
a small number of odd recipients has a smaller joint alphabet. It does
not determine the exact rank above the threshold, and it does not apply
to a whole network after target-product and carry constraints are imposed.

The source definitions and their pinned revision are those in
[factoring-network-structure.md](factoring-network-structure.md). Write

```math
h=2^s\ge8,\quad L=2h^2,\quad B=2^L,\quad
M=B^4+1,\quad A=2h,\quad\zeta=2^{4L/h}\pmod M.
\tag{1}
```

For a recipient set `J`, ignore the known invertible gate-dependent
coordinate multipliers and put

```math
F_J(u)=\left(\sum_{i=0}^{h-1}u_i\zeta^{ij}\right)_{j\in J},
\qquad u\in\{0,\ldots,B-1\}^h.
\tag{2}
```

Theorem. If `J` consists of `k` distinct odd indices in `{0,...,A-1}`,
then (2) is onto `(Z/MZ)^k` if and only if `k<h/4`. For `k<h/4`, its
equality matrix has exact rank `M^k` over every coefficient field, and
a preimage of any output can be computed in time polynomial in `h,L`.
The empty recipient set is included trivially.

## 1. The exact kernel lattice

Let `1<=k<h/4`, and define

```math
\Lambda_J=\left\{v\in\mathbb Z^h:
     \sum_i v_i\zeta^{ij}=0\pmod M\text{ for every }j\in J\right\}.
\tag{3}
```

The first `k` columns of the evaluation matrix form a Vandermonde matrix
on the distinct points `zeta^j`. Every difference of two such points is
a unit modulo `M`, by the composite-ring interpolation lemma in the
companion note. Thus evaluation from `Z^h` onto `(Z/MZ)^k` is surjective,
and its kernel has index and determinant

```math
[\mathbb Z^h:\Lambda_J]=\det\Lambda_J=M^k.       \tag{4}
```

This is a statement about unrestricted integer coefficient vectors.
The bounded-digit surjectivity in the theorem requires the next steps.

There is an explicit integer basis with `O(L)`-bit entries. Use `M e_i`
for the first `k` basis vectors. For each remaining coordinate, set it
to one and the other free coordinates to zero; solve the first `k`
coordinates modulo `M` so that the evaluation vanishes. The resulting
basis matrix is block triangular, with determinant `M^k`. Lagrange
interpolation or unit-pivot Vandermonde elimination constructs it in
polynomial bit time.

## 2. Negacyclic symmetry turns one short vector into a full set

Let

```math
S(v_0,\ldots,v_{h-1})=(-v_{h-1},v_0,\ldots,v_{h-2}).
\tag{5}
```

For odd `j`, `(zeta^j)^h=-1`, so evaluating `Sv` multiplies the
evaluation of `v` by `zeta^j`. Hence `S Lambda_J = Lambda_J`.

For any nonzero vector `v`, the vectors

```math
v,Sv,\ldots,S^{h-1}v                              \tag{6}
```

are linearly independent over `Q`. To see this, identify vectors with
polynomials of degree below `h` modulo `X^h+1`. For power-of-two `h`,
that polynomial is irreducible: after replacing `X` by `X+1`, every
intermediate coefficient is even and the constant coefficient is two,
so Eisenstein's criterion at two applies. A dependence in (6) would say
that `X^h+1` divides the product of two nonzero polynomials of degree
below `h`, one of which is the polynomial representing `v`. Irreducibility
would force one factor to be divisible by `X^h+1`, a contradiction.

Let `W` be the integer matrix with the vectors in (6) as its columns.
Every row of `W` contains the absolute values of the coordinates of `v`
once each. Therefore

```math
\sum_a|W_{ia}|=\|v\|_1\quad\text{for every row }i.       \tag{7}
```

The columns need not be a basis of the entire kernel lattice. Their
integer span is a full-rank sublattice of it, and that is sufficient for
the rounding algorithm below.

## 3. A short vector is computable within the digit budget

Apply the deterministic LLL algorithm with parameter `3/4` to the
integer basis from Section 1. The algorithm and its polynomial bit-cost
analysis are from A. K. Lenstra, H. W. Lenstra Jr., and L. Lovasz,
*Factoring polynomials with rational coefficients*, Math. Ann. 261
(1982), 515--534, Section 1
([primary paper](https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/Lovasz/LovaszLenstraLenstrafactor.pdf),
[DOI](https://doi.org/10.1007/BF01457454)).

For completeness, the norm bound needed here follows directly from the
LLL conditions. If `b_i*` are the orthogonalized basis vectors, size
reduction and the Lovasz condition imply
`||b_(i+1)*||^2 >= ||b_i*||^2/2`. Multiplying the resulting inequalities
for the first basis vector and using the determinant gives

```math
\|v\|_2\le2^{(h-1)/4}(\det\Lambda_J)^{1/h},
\qquad
\|v\|_1\le\sqrt h\,2^{(h-1)/4}M^{k/h}.
\tag{8}
```

Here `v` is the first reduced basis vector and is nonzero. Since
`k<=h/4-1`, `log2 M<4L+1`, and `log2 h<=h`, the base-two logarithm of
the right side in (8) is at most

```math
\begin{aligned}
\tfrac12\log_2h+\tfrac{h-1}4
 +(\tfrac14-\tfrac1h)(4L+1)
&\le L-\frac{4L}{h}+\frac{3h}{4}\\
&=L-\frac{29h}{4}<L-1.
\end{aligned}
\tag{9}
```

Thus `||v||_1 < B/2`. In particular `||v||_1<=B-1`.
This large margin comes from the source's actual relation `L=2h^2`;
it is not a generic consequence of lattice determinant alone.

All input basis entries have `O(L)` bits. The LLL computation has
polynomial bit cost in `h,L`. Its vector has polynomially many bits,
as do `W`, its determinant, and the exact rational entries of `W^-1`.
The uniform algorithm never uses a nonconstructive shortest-vector
oracle.

## 4. Construct a preimage by rounding

Given the requested residues `w`, solve the first `k` coefficient
coordinates by the same invertible Vandermonde matrix, set the others
to zero, and choose canonical integer representatives. This gives a
polynomial-size integer vector `x` with `F_J(x)=w`.

Let `c=((B-1)/2,...,(B-1)/2)` and compute exactly

```math
z=W^{-1}(x-c),\qquad n_a=\lfloor z_a+1/2\rfloor,
\qquad u=x-Wn.                              \tag{10}
```

All entries of `u` are integers. Since each column of `W` belongs to
the kernel lattice, `F_J(u)=w`. Equations (7)--(8) give

```math
|u_i-c_i|
 =\left|\sum_aW_{ia}(z_a-n_a)\right|
 \le\tfrac12\|v\|_1\le(B-1)/2.
\tag{11}
```

Thus every coordinate lies in the required integer interval
`0<=u_i<B`. Exact rational matrix inversion, nearest-integer rounding,
and matrix multiplication all take polynomial bit time. Ties may use
the fixed floor rule in (10).

Restoring the original factors `16^(r_gate*j)` only changes the requested
residues before this construction. The resulting digits are settings of
the original input gate, rather than values in a larger relaxed domain.

## 5. Sharpness and independent-bundle consequence

For `k>=h/4`, the number of possible source digit vectors is `B^h`,
whereas

```math
M^k=(B^4+1)^k>B^{4k}\ge B^h.                \tag{12}
```

Surjectivity is impossible. This proves the exact threshold stated at
the start. It does **not** prove an exact formula for the smaller image
when `k>=h/4`.

As one consequence, partition the whole recipient layer into arbitrary
bundles, each containing fewer than `h/4` recipients. Let `o_b` be the
number of odd recipients in bundle `b`. Projecting its message image
onto those coordinates gives a surjection onto `R^(o_b)`. Its equality
matrix therefore has rank at least `M^(o_b)`. Since there are exactly
`h` odd recipients, any separately factored bundle interface has total
logarithmic channel cost at least

```math
\sum_b\log_2\mathrm{rank}(E_b)
 \ge h\log_2M>4hL.                           \tag{13}
```

This conclusion allows adaptive or irregular recipient sets; it does
not require Fourier cosets. The whole unsplit fan-out still has rank
`B^h` and needs only `hL` bits. Equation (13) concerns the declared
independent bundle factors. It neither sums ranks of arbitrary correlated
factors nor applies after later arithmetic constraints are contracted.

## 6. Sets with a common two-adic valuation

There is also a useful sufficient condition outside the primitive class.
Suppose every recipient is `q` times an odd integer, with `q` a power of
two dividing `h`. Put `r=h/q`. Folding input positions modulo `r` gives

```math
a_\ell=\sum_{t=0}^{q-1}(-1)^t u_{\ell+rt},\qquad 0\le\ell<r.
\tag{14}
```

Each coordinate independently fills an integer interval of exactly

```math
H=q(B-1)+1
```

values. Translating those intervals only translates each output by a
known residue. The reduced root `zeta^q` has `r`th power minus one, and
the reduced recipient indices are odd. Thus the same kernel-lattice
argument, now in dimension `r`, applies.

If a nonzero vector in the reduced kernel lattice is found with
`l1` norm at most `H-1`, the rounding procedure gives a complete
preimage algorithm for the reduced digit box. A sufficient condition
guaranteed by LLL is

```math
\sqrt r\,2^{(r-1)/4}M^{k/r}\le H-1.          \tag{15}
```

Every reduced digit can then be expressed as a sum of `q` independently
bounded signed original digits by a greedy interval decomposition.
This restores an actual source input vector in polynomial bit time.
In particular, `k<r/4` always suffices for the source parameters when
that range is nonempty: the estimate in (9) becomes only stronger because
`L=2h^2>=2r^2` and `H>=B`.

For `q>1`, (15) is a sufficient condition, not a sharp threshold. The
larger folded alphabet can make the boundary `k=r/4` surjective; the
exact Fourier-coset theorem already exhibits that phenomenon. For a
specific set, the actual norm of a computed vector is a stronger exact
certificate than the general upper bound (15).

## 7. Verification and scope of the implementation

Run

```sh
python3 tools/check_factoring_primitive_recipients.py \
  --output evidence/factoring-primitive-recipients.json
```

The implementation constructs the modular kernel basis, uses exact
rational LLL arithmetic, verifies its reduction conditions, checks the
kernel congruences and nonzero negacyclic determinant, and reconstructs
prescribed messages using exact rational rounding. It records the
actual short-vector norms and checks their digit budget. Bounded finite
controls compare independent direct evaluations with the decoded
source digits. The theorem for every odd recipient set follows from
Sections 1--5; finite examples do not replace that proof.

The new arithmetic statement leaves the factoring score at `1/2`.
Its constructive content is the exact full-range threshold and its
efficient source-digit preimage algorithm. Any further compression
would have to exploit larger joint bundles or the contracted
multiplication/carry constraints, rather than selecting different small
primitive recipient sets.

This analysis uses the source transform of Rudich and Rousseau and the
LLL algorithm of Lenstra, Lenstra, and Lovasz. The new application and
checker were developed with OpenAI assistant support. No Lean source
is changed.
