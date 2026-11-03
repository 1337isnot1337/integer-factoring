# The identity minor is a hash join, not a Cartesian product

## Result and scope

The source's finite identity-minor diagnostic in
[the half-top structure note](half-top-structure.md) enumerates every
ordered pair of messages and compares them. The same sum is an equality
join. For `s=1` the source check performs `16` edges times `256^2 =
1,048,576` pair comparisons. A hash table computes the identical integer
answer in `8,192` message-plus-output operations, a factor of `128`
fewer. For `s=0` the factor is `2`.

This is an exact finite source-definition diagnostic. It shows that the
Cartesian enumeration in this minor is avoidable. It does **not** show
that the full edge tensors are functional, does not compress the
evaluator pad `V=P^2 lambda`, and does not change the factoring exponent.
The complete score remains `1/2`.

## Mechanism

For one edge the source builds the message list `m_0,...,m_{2^d-1}` and
counts

```math
\sum_{v=0}^{2^d-1}\sum_{w=0}^{2^d-1}[m_v=m_w].
```

Let `c_m` be the multiplicity of message `m`. The same sum is
`\sum_m c_m^2`. A hash table builds the multiplicities in one pass over
the messages and evaluates the squares in one pass over the distinct
messages. The work is proportional to messages plus output pairs, not to
the square of the message count. On the identity minor the messages are
distinct, so the output is one pair per message and the join is linear.

This is the local-solver principle named as an unexcluded route in
[the half-top structure note](half-top-structure.md): count compatible
settings without visiting their whole Cartesian product. It is distinct
from the source's deterministic-message identity, which is what makes
the join exact here.

## What it does not settle

The half-top lower bound has two terms: the local enumeration
`H_Q` and the pad term `V/r` with `V=P^2 lambda`. Replacing a pair scan
by a hash join addresses only the first term, and only when the actual
edge relation is functional or sparse. The pad term is unchanged. A
complete improvement still needs a joint boundary compression that
reduces the pad or a changed source family.

The next bounded diagnostic is to extract the actual small edge tensor
from the source definition and measure its message-multiplicity profile.
If a constant fraction of messages is unique, the same join is
output-sensitive on that edge. If the profile is diffuse, the Cartesian
lower bound survives for that tensor. The identity minor alone is not
evidence either way for the full tensor.

## Reproduction

```sh
python3 tools/check_hash_join_identity_minor.py \
  --output evidence/hash-join-identity-minor.json
make verify-hash-join-identity-minor
```

The checker reproduces the source's `s=0,1` identity minor exactly,
compares brute-force and hash answers, and records the operation counts.
It uses the same composite modulus `4294967297=641*6700417` at `s=1`.
