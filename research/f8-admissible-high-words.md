# The q96 F8 map has only 45 admissible single-word extensions

## Result

Enumerate all single-word vectors `{w:1}` that can be added to both the
embedding and retraction of the q96 F8 hyperoval map while preserving
`LP=I`, deletion closure, and the graded-letter axioms. Exactly 45 words
qualify. Their heights are

```text
26 (6 words), 27 (7), 28 (5), 29 (4), 30 (4), 31 (2), 32 (4),
33 (1), 34 (3), 36 (3), 38 (2), 40 (2), 42 (1), 44 (1).
```

There is no admissible single-word extension of height at most 25, and
no admissible single word of height 35, 37, 39, 41, or 43. The q99
extension is a **joint** addition of three words: `160` and `192` are
individually admissible (heights `26` and `28`), while `224` (height
`39`) becomes admissible only after the two pairs are added, because
its deletions need them. Thus `224` is not in the 45-word single-word
inventory.

## Low-support saturation

Every word of height at most 25 occurs in the support of at least one q96
letter. There is no word outside the old supports at low height, so the
single-word census also rules out a pair of free low words: both members
of a support-disjoint pair would have to be free, and none are. A new
low coordinate would have to use words already appearing in old
supports, with coefficients arranged to stay orthogonal; that is a
change of the letter supports, not an addition of a free word.

## Interpretation

The q99 height-39 extension is not one arbitrary choice from a large
space; its single-word inventory is finite and small. The low core has
no hidden single-word extension. The 45 words all add capacity above the
height-25 boundary, and the complete q99 profile already shows that a
three-word joint extension (two of them from this inventory) gains only
`1.27e-5`. Even if every admissible single word were added, the number
of added letters would be 45; the low-core plateau is not removed by
this inventory. This is a structural
limitation of the single-word high extension, not a proof that no
multi-word addition or different low block can cross.

The next low-block search must therefore change the letter supports
themselves, for example through a different rank-three image or a
different protected structure; adding free low or high words is not the
route.

## Reproduction

```sh
make verify-f8-admissible-high-words
```

The checker rebuilds the 45-word census and compares it with the
[recorded list](../evidence/f8-admissible-high-words.json), including
the base-profile hash and the q99 words. It does not screen multi-word
additions.
