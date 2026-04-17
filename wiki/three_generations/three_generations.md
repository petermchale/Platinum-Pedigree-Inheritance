# Extending the same pass to a third generation

This page is part of the [wiki](../index.md) and extends the
[nuclear-family walkthrough](../nuclear_family/nuclear_family.md) by
marrying **Kid3** — the paternal recombinant from the nuclear-family
page — to a fresh outside-marriage founder and adding two
grandchildren. Two points fall out of this setup:

1. `gtg-ped-map` handles G2→G3 with the **same single loop** it used
   for G1→G2, without ever constructing a joint inheritance vector
   across all founders.
2. A crossover in Kid3's own G2→G3 meiosis can layer on top of the
   ancestral G1-meiosis crossover she already carries, producing what
   the grandchild's letter trace records as a **recombinant of a
   recombinant**. Seeing this requires two grandchildren rather than
   one, because the diagnostic signature is that the ancestral
   transition appears in *both* G3 sibs while the per-meiosis one
   appears in only *one*.

All line numbers refer to commit `583ef82`. As in the other
walkthrough pages, each function link is followed by its call site in
the driver — `main()` in
[`map_builder.rs`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L989) — so you can step through the
driver source in parallel with this walkthrough.

The toy simulation adds three things on top of the nuclear-family
example:

- **Kid3** (already carrying a G1-ancestral A→B transition on her
  paternal homolog and C throughout on her maternal homolog, from the
  nuclear-family pass) marries **Spouse**, a fresh founder whose two
  homologs are labelled **E** and **F** by
  [`Iht::new`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/iht.rs#L172) (driver calls at
  [`map_builder.rs:1059`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L1059) for the master template
  and [`map_builder.rs:1111`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L1111) for each VCF site).
- **GK1** inherits Spouse's E homolog plus a Kid3 gamete produced by a
  crossover: Kid3's paternal homolog (letters A,A,A,A,B,B) for sites
  0–5, then Kid3's maternal homolog (letter C,C) for sites 6–7.
- **GK2** inherits Spouse's F homolog plus Kid3's paternal homolog
  unrecombined (letters A,A,A,A,B,B,B,B).

The whole simulation spans 8 VCF sites. Everything below
is reproducible by running

```
python wiki/generate_wiki.py --page three_generations
```

which regenerates both the figure PNGs referenced here and this
markdown file itself.

## 1. Pedigree with outside marriage

![Figure 1 — Three-generation pedigree with outside marriage](fig1.png)

Kid3 arrives at this pass already carrying per-site letter labels from
the nuclear-family page: paternal slot **A** for sites 0–3 and **B**
for sites 4–7 (an ancestral G1-Dad crossover), maternal slot **C**
throughout. She is not a founder of the three-generation pedigree, but
from the perspective of the G2→G3 sub-problem she plays exactly the
role Dad and Mom played in G1→G2: her two homologs are already tagged
with letters, and those letters are what `gtg-ped-map` will propagate
to GK1 and GK2.

Spouse, on the other hand, *is* a founder relative to this pedigree
branch, so [`Iht::new`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/iht.rs#L172) (called from the driver
at [`map_builder.rs:1059`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L1059) and
[`map_builder.rs:1111`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L1111)) hands him the next
fresh letter pair `(E, F)`. Nothing about Spouse depends on the G1 pass.

## 2. Unphased VCF rows for the G2→G3 pass

![Figure 2 — Unphased VCF rows for the G2->G3 pass](fig2.png)

These are the only genotype rows the tool sees for the new nuclear
unit. Kid3's row is read with Kid3 in the **parent** role. That is the
key structural point: `gtg-ped-map` does not treat G1 and G2
individuals differently; it just iterates over (parent, spouse, child)
triples in ancestor-first depth order given by
`family.get_individual_depths()` (see
[`ped.rs`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/ped.rs)).

## 3. Recursive informative-site deduction

![Figure 3 — Recursive informative-site deduction (G2 -> G3)](fig3.png)

The exact same function,
[`track_alleles_through_pedigree`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L295) (driver call at
[`map_builder.rs:1116`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L1116)), that handled G1→G2 now
handles G2→G3. For each (parent, spouse) pair it calls
[`unique_allele`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L243) (invoked from inside the walk
at [`map_builder.rs:315`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L315)) to find alleles that
one partner carries and the other does not:

- **Spouse-informative** (Spouse het × Kid3 hom): the unique paternal
  allele tags whichever Spouse homolog (`E` or `F`) each grandchild
  inherited. In this simulation these are sites `[2, 3]`.
- **Kid3-informative** (Kid3 het × Spouse hom): symmetric, but here
  the letter tagged on each grandchild's maternal slot is whichever of
  Kid3's *per-site* letters — `A`, `B`, or `C` — sits on the homolog
  carrying Kid3's unique allele at that site. These are sites
  `[0, 1, 4, 5, 6, 7]`.

When the walk reaches the Kid3–Spouse pair,
[`get_iht_markers`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L274) (called from inside the walk
at [`map_builder.rs:328`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L328)) reads Kid3's
already-assigned letters directly — those labels were written during
the earlier G1→G2 iteration of the same loop. That is what makes the
algorithm look recursive across generations even though it is a single
ancestor-first pass: by the time the loop reaches a G2 parent, her
letter labels are already finalized and they serve as the "founder
labels" for the G2→G3 sub-problem, even when those labels vary from
site to site. No joint inheritance vector over all founders
`{A, B, C, D, E, F}` is ever constructed; each grandkid ends with
exactly two letters per site — one per homolog — identical in shape
to the output of the nuclear-family pass.

Figure 3 keeps the same layout as the nuclear-family analogue: the two
indicator rows (`*` marks informative sites, `_` non-informative ones)
sit above the grandkid rows, each grandkid's paternal row (`p`) is
placed directly above its maternal row (`m`), and every column is
aligned so you can read each letter assignment straight up to the
indicator that produced it.

## 4. Ancestral vs per-meiosis crossovers

![Figure 4 — Ancestral vs per-meiosis crossovers after block collapse](fig4.png)

The same block-collapse, gap-fill and flip routines invoked for G1→G2 —
[`collapse_identical_iht`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L385) (driver call at
[`map_builder.rs:1191`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L1191)),
[`fill_missing_values`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L617) (driver call at
[`map_builder.rs:1200`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L1200)),
[`fill_missing_values_by_neighbor`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L540) (driver call
at [`map_builder.rs:1201`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L1201)), and
[`perform_flips_in_place`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L702) (driver calls at
[`map_builder.rs:1135`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L1135),
[`map_builder.rs:1193`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L1193), and
[`map_builder.rs:1203`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L1203)) — run on the G3 trace
without modification. Two transitions fall out, and the contrast
between them is the point of this page:

- **Ancestral A→B at sites 3/4** appears on *both* grandchildren's
  maternal rows. This is the G1-Dad crossover Kid3 already carried on
  her paternal homolog — it propagates to every descendant that
  inherits the affected segment, which is both GK1 and GK2 here. A
  single meiotic event (in G1) produces two rows in
  `{prefix}.recombinants.txt`.
- **Per-meiosis B→C at sites 5/6** appears on GK1's maternal row only.
  This is the new crossover in Kid3's G2→G3 meiosis, and because it
  happened in the gamete that made GK1, it is absent from GK2's trace.
  One meiotic event (in G2), one row.

Both transitions are emitted by
[`summarize_child_changes`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L673) (driver call at
[`map_builder.rs:1228`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/583ef829a14f276f3a6cb1e2f34195d45397904f/code/rust/src/bin/map_builder.rs#L1228)). The raw output
therefore contains three rows — two for the ancestral crossover, one
for the per-meiosis crossover — even though the underlying biology is
only two distinct meiotic events. Collapsing the shared ancestral
entries is the downstream reconciliation step that
[`methods.md §4.5`](../methods.md) describes; the current
implementation deliberately exposes the raw transitions and leaves the
reconciliation to the analyst.

## 5. Truth versus deduced

![Figure 5 — Truth vs deduced founder labels](fig5.png)

For every grandkid the deduced paternal and maternal label streams
match the ground truth at every site (0 mismatches out of
32 label slots), including both the ancestral A→B and the
per-meiosis B→C in GK1. The block map stored for G3 uses the letters
`{A, B, C, E, F}` — all three of Kid3's per-site letters reach G3,
because GK1's gamete from Kid3 crossed over between homologs. As in
the nuclear-family case, the block map contains only founder letters;
the 0/1 allele sequence of each haplotype is reconstructed downstream
by `gtg-concordance`, covered in the
[concordance walkthrough](../concordance/concordance.md).
