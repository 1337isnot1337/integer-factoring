# No naive large-prime lift of the bounded-spin coefficient table

## Result

The rank=trace identity for an idempotent holds in every characteristic.
If the characteristic-five coefficient table were the reduction of an
integral table, rebuilding it over a prime `p > 89^2` would turn each
central-idempotent trace into its exact integer rank. This checker tests
that shortcut at length ten, where the eight-site control is cheap.

Using the same entry-functional construction and a **full-rank witness**
at every prime:

```text
p=5      witness rank 8   commutative yes   unital yes
p=7      witness rank 8   commutative no    unital yes
p=11     witness rank 8   commutative no    unital yes
p=10007  witness rank 8   commutative no    unital yes
```

Thus the larger-prime tables are not the reduction of the characteristic-
five table in the same coefficient basis. The characteristic-five table
is not the reduction of an integral table over these primes. The
characteristic-five dimension collapse is essential to the ten-site
coefficient algebra.

## Consequence

- Trace forms over `F5` give cut-module ranks only modulo five. They cannot
  be lifted to exact integer ranks by changing the prime.
- Exact ranks require a direct `F5` rank computation on the cut modules, a
  characteristic-zero model, or an independent rank certificate. The
  word-pair expansion obstruction rules out the naive direct action; a
  transfer-matrix or quotient method is still required.
- The same test is a guard against treating the `F5` idempotent
  decomposition as an integral decomposition. The central-idempotent
  block structure is exact over `F5`; it does not by itself determine
  integer cut ranks.

No factoring exponent or score changes. The complete bound remains `1/2`.

## Reproduction

```sh
python3 tools/check_bounded_spin_characteristic_five_lift.py \
  --output evidence/bounded-spin-characteristic-five-lift.json
make verify-bounded-spin-characteristic-five-lift
```

The checker rebuilds the ten-site entry-functional table at four primes
with full-rank witnesses and checks commutativity and unitality. It is a
negative diagnostic, not a rank bound.
