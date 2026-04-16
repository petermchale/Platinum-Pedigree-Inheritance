# Structural haplotype mapping in a nuclear family

This page is part of the [wiki](../index.md) and walks through
`gtg-ped-map`'s structural labelling algorithm on the simplest possible
pedigree: a two-generation nuclear family with two founders (dad and
mom) and three children. It complements the full
[`methods.md`](../methods.md) write-up by zooming in on the per-site
mechanics and pinning each panel to the exact Rust code that implements
it. All line numbers refer to commit `7a91d01`. Each function link is
followed by its call site in the driver — `main()` in
[`map_builder.rs`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L989) for `gtg-ped-map`, and `main()` in
[`gtg_concordance.rs`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/gtg_concordance.rs#L315) for `gtg-concordance` —
so you can step through the driver source in parallel with this
walkthrough.

The toy simulation hard-codes four founder haplotypes over
8 sites and three children whose transmissions are known a
priori. Everything below is reproducible by running

```
python wiki/generate_wiki.py --page nuclear_family
```

which regenerates both the panel PNGs referenced here and this markdown
file itself.

## 1. Ground truth

![Figure 1 — Ground-truth founder haplotypes](fig1.png)

Dad carries two haplotypes, arbitrarily labelled **A** and **B**; mom
carries **C** and **D**. These letter labels are assigned at startup by
[`Iht::new`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/iht.rs#L172) (driver calls at
[`map_builder.rs:1059`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L1059) for the master template
and [`map_builder.rs:1111`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L1111) for each VCF site),
which hands each founder a fresh pair of capital letters — `(A,B)`,
`(C,D)`, `(E,F)`, … — *without* associating any allele sequence with
them. The letters are pure structural placeholders whose only job is
to be carried from founders to descendants.

In this simulation:

- **Kid1** inherits `(A, C)` with no recombination.
- **Kid2** inherits `(B, D)` with no recombination.
- **Kid3** inherits `(A|B, C)` — the paternal slot crosses over between
  sites 3 and 4, so Kid3 carries dad's `A` haplotype on sites 0–3 and
  dad's `B` haplotype on sites 4–7.

The goal of `gtg-ped-map` is to recover exactly these letter
transmissions from the jointly-called VCF alone, without ever looking
at the underlying 0/1 allele sequence.

## 2. Unphased VCF input

![Figure 2 — Unphased VCF view](fig2.png)

This is the only genotype information `gtg-ped-map` sees (plus the PED
file that declares who is whose parent). Two observations matter:

- **Haplotypes cannot be distinguished from genotypes alone.** A `0/1`
  call for dad does not reveal which of his two homologs carries the
  `1`, so a single-individual view has no way to label `A` vs `B`.
- **Patterns across the family resolve the ambiguity.** If dad is `0/1`
  while mom is `0/0`, then any child that also carries a `1` must have
  inherited dad's `1`-carrying homolog — precisely the logic of the
  informative-site test in the next section.

Only biallelic SNVs enter the map; indels are filtered at read time via
[`is_indel`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L501), invoked from the VCF-reading loop at
[`map_builder.rs:164`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L164) inside `parse_vcf` (the
driver calls `parse_vcf` at [`map_builder.rs:1092`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L1092)).

## 3. Informative-site detection and letter deduction

![Figure 3 — Informative-site deduction (paternal and maternal)](fig3.png)

For each VCF record,
[`track_alleles_through_pedigree`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L295) (driver call at
[`map_builder.rs:1116`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L1116)) walks the pedigree in
ancestor-first depth order and, for every `(parent, spouse)` pair,
calls [`unique_allele`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L243) (from inside the walk at
[`map_builder.rs:315`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L315)) to ask whether the parent
carries an allele that the spouse does not. Two cases can arise:

- **Dad-informative** (dad het × mom hom): dad's unique allele tags
  whichever paternal homolog (`A` or `B`) each child inherited. In
  this simulation these are sites `[0, 1, 4, 5]`.
- **Mom-informative** (mom het × dad hom): symmetric, tagging `C` or
  `D`. These are sites `[2, 3, 6, 7]`.

When a child carries the parent's unique allele, the child's paternal
(or maternal) slot is filled with the parent's letter; otherwise the
slot is filled with the *other* letter of that parent's pair. Because
the depth-ordered walk always processes a parent before its children,
[`get_iht_markers`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L274) (called from inside the walk
at [`map_builder.rs:328`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L328)) reads the parent's
already-assigned letters when propagating to the next generation,
which is what makes the method look "recursive" across generations
while being expressed as a single loop.

Non-informative sites (both parents het, or both hom for the same
allele) contribute nothing at this stage and are rendered as `.` in
Figure 3. The two indicator rows (`*` marks informative sites, `_` marks
non-informative ones) sit directly above the kid rows, with every
column aligned, so you can read each letter assignment straight up to
the indicator that produced it. Each kid's paternal row (`p`) sits
directly above its maternal row (`m`). Kid3's `p` row already exhibits
the A→B recombination at sites 3/4, even though no block collapse has
happened yet.

## 4. Block collapse and noise filtering

![Figure 4 — Collapsed blocks with recombination](fig4.png)

Several Rust routines clean the per-site letter trace up before it is
written to disk:

1. [`backfill_sibs`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L804) (driver call at
   [`map_builder.rs:1122`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L1122)) uses the fact that
   siblings must together carry both founder homologs. If exactly one
   child is tagged at a site, the others can be inferred by elimination.
   In this toy simulation every informative site already tags all three
   kids, so backfill is a no-op here, but on real data it is essential
   in noisy regions.
2. [`collapse_identical_iht`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L385) (driver call at
   [`map_builder.rs:1191`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L1191)) merges adjacent sites
   with compatible letter assignments into blocks, while
   [`fill_missing_values`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L617) (driver call at
   [`map_builder.rs:1200`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L1200)) and
   [`fill_missing_values_by_neighbor`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L540) (driver
   call at [`map_builder.rs:1201`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L1201)) fill the `.`
   gaps visible in Figure 3 from flanking blocks.
3. [`count_matching_neighbors`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L935) (driver call at
   [`map_builder.rs:1172`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L1172)) and
   [`mask_child_alleles`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L970) (driver call at
   [`map_builder.rs:1187`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L1187)) identify isolated
   runs shorter than `--run` (default 10 markers) and mask them back
   to `?` as likely sequencing noise, so that collapse does not invent
   spurious recombinations.
4. [`perform_flips_in_place`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L702) enforces consistent
   founder-letter orientation across blocks, since the two letters in
   a founder's pair are interchangeable within any single block. The
   driver calls it three times — before and after block collapse,
   and again after gap fill — at
   [`map_builder.rs:1135`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L1135),
   [`map_builder.rs:1193`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L1193), and
   [`map_builder.rs:1203`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L1203).

After these steps, each kid's paternal and maternal slots are fully
resolved — shown on separate rows per kid in Figure 4. Kid3's
highlighted A→B transition on the paternal row is emitted to
`{prefix}.recombinants.txt` by
[`summarize_child_changes`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L673) (driver call at
[`map_builder.rs:1228`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/map_builder.rs#L1228)).

## 5. Truth versus deduced

![Figure 5 — Truth vs deduced founder labels](fig5.png)

For every kid the deduced paternal and maternal label streams match the
ground truth at every site (0 mismatches out of
48 label slots), including Kid3's recombination. The full
output of `gtg-ped-map` for this chromosome is the set of blocks shown
above plus the `recombinants.txt` entry for Kid3's switch — and
critically, **nothing else**. The block map stores only founder
letters; it does *not* store the 0/1 allele sequence of any haplotype.

Reconstructing which allele each letter represents at every VCF site
is the job of `gtg-concordance`, which will have its own wiki page
once migrated. For every block, `gtg-concordance` enumerates the
`2^F` founder-phase orientations produced by
[`Iht::founder_phase_orientations`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/iht.rs#L492) (driver call
at [`gtg_concordance.rs:256`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/gtg_concordance.rs#L256)), maps letters
to VCF alleles via [`assign_genotypes`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/iht.rs#L442) (driver
call at [`gtg_concordance.rs:267`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/7a91d0149a5c070255f938ebd9d0d8a91f459cf1/code/rust/src/bin/gtg_concordance.rs#L267)), and
picks the orientation that minimises mismatches against the observed
genotypes. The split of responsibilities is deliberate and strict:
`gtg-ped-map` writes only letters and only at informative sites, while
`gtg-concordance` is the sole place where letter→allele correspondence
is computed and written out.
