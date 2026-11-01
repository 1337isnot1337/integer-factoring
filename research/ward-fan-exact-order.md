# The 84-label fan has exact Ward vanishing order four

The [exhaustive lowest-boundary calculation](ward-fan-boundary.md)
proved that the specified \(P=L=10,k=3\), 84-label fan tensor
is divisible by \((q+r)^4\). An exact coefficient witness now
shows that the fourth power is its **exact** common order.
This statement concerns the explicit input/output map in the
frozen fan packages.

Write
\[
 G(1+x,1+y,1)=y^3\Gamma(x,y).
\]
At the logical index tuple
\[
 (x_{\rm label},y_{\rm label},u_{\rm label},v_{\rm label})
 =(71,79,79,71),
\]
the next quotient layer is
\[
 [y]\Gamma(x,y)=x^2.
 \tag{1}
\]
The tuples \((69,77,77,69)\) and \((71,77,77,71)\) give the
same polynomial. Label indices are zero-based and refer to
the exact basis stored in
[ward-fan-leading-block.json](../evidence/ward-fan-leading-block.json).
In all three cases both state inputs are the second
unicyclic residual vector.

Thus
\[
 [x^2y^4]G=1
\]
at a displayed entry. The tensor is not divisible by
\((q+r)^5\). Combined with the preceding full-boundary
identity, this proves exact common order four over every
characteristic-two extension field.

For the first displayed entry, homogeneity also gives the
coordinate-free leading statement
\[
 G_{\rm entry}(p,q,r)
 \equiv (q+r)^4r^4(p+r)^2
       \pmod{(q+r)^5}.
 \tag{2}
\]
Both sides have homogeneous degree ten.

## Exact thin-jet calculation

The evaluator keeps \(\Gamma\) modulo \(y^2\), with total
degree at most seven in \(x,y\). Multiplication takes place
in a commutative characteristic-two quotient ring. Its
binary encoding separates consecutive \(y\)-layers by
sixteen bit positions; surviving \(x\)-degrees are at most
seven, so multiplying two terms cannot alias an \(x\)
degree into a \(y\)-degree.

The [two-state fan oracle](ward-fan-polynomial-queries.md)
is valid over this quotient ring and uses no division.
Zero queries are pruned only by a valid local lower bound
on the total \(y\)-valuation. No assumption about a
state-sector support rule is used.

The bounded probe checked 2,048 logical entries, favoring
the two unicyclic residual input directions. It found the
three exact entries above. Sampling is used only to locate
the witness. Each displayed polynomial is then computed
exactly and independently verified.

In particular, the checker separately contracts the
original Pfaffian kernel with both full input vectors
in \(\mathbf F_2[x,y]/(x^8,y^5)\). It compares that result
with \(y^3\Gamma\). This calculation uses neither the Ward
preimages nor the two-state transfer.

The ring implementation also passes 120 comparisons with
direct squarefree form-variable multiplication. Eight
positive crossing-jet controls expand around
\(p=0,q=y,r=1\), where \(\eta=1+y\) is a unit. Their
constant crossing values are one, and their full
first-order jets agree with direct original-kernel
contraction.

The probe required 92,712 physical query requests and
145 distinct nontrivial physical jets. It took less
than one second in the recorded environment. It does
not assert the complete rank of the next coefficient
boundary.

## Reproduction

Run

~~~sh
python3 tools/check_ward_fan_jet_probe.py --samples 2048 --output /tmp/ward-fan-jet-probe.json
~~~

The [jet evaluator](../tools/ward_fan_jet_queries.py),
[probe checker](../tools/check_ward_fan_jet_probe.py), and
[evidence](../evidence/ward-fan-jet-probe.json) record the
exact coefficient witnesses and source hashes. The next
complete-boundary calculation, if performed, is a
separate rank computation.
