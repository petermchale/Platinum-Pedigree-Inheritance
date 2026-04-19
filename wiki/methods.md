# Methods — pedigree haplotype mapping with gtg-ped-map and gtg-concordance

The pipeline consists of two Rust binaries — `gtg-ped-map`
([`map_builder.rs`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L989)) and `gtg-concordance`
([`gtg_concordance.rs`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/gtg_concordance.rs#L315)) — that share a common
inheritance-tracking library (`code/rust/src/iht.rs`) and are driven by a
standard PED file plus a jointly-called VCF. The two binaries have
deliberately separate responsibilities: `gtg-ped-map` labels haplotypes
with structural letters at informative sites only, and
`gtg-concordance` is the sole place where letter→allele correspondence is
resolved and phased genotypes are written. Everything in between — the
per-site letter propagation, the block-collapse and flip machinery, and
the orientation search — is described below.

All line numbers refer to commit `9d28575`; every function link is
paired with its driver call site in `main()` so the manuscript can be
read in parallel with the source. Worked toy examples for each major
component live in the companion walkthrough pages:
[nuclear family](nuclear_family/nuclear_family.md),
[three-generation pedigree](three_generations/three_generations.md),
and [gtg-concordance](concordance/concordance.md).

## 1. Inputs

The tools consume two files:

- A **PED file** describing the pedigree, parsed by
  [`code/rust/src/ped.rs`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/ped.rs). Each row declares an individual's
  family, parents, sex, and affection status; founders are rows whose
  parent columns are `0`.
- A **jointly-called, tabix-indexed VCF** — typically DeepVariant on
  HiFi — with genotypes for every individual named in the PED. The VCF
  must be jointly called so that every record has a defined GT for
  every sample; per-sample VCFs merged after the fact will not work
  because the informative-site test below requires the parents' and
  children's genotypes at the same coordinate in a single record.

Only biallelic SNVs enter map construction. Indels are filtered at read
time by [`is_indel`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L501) (invoked from the VCF-reading loop
at [`map_builder.rs:164`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L164) inside `parse_vcf`; the driver
calls `parse_vcf` at [`map_builder.rs:1092`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1092)).
Multi-allelic records are likewise dropped at parse time — the letter
machinery is built around a two-homolog, two-allele model and does not
generalise to three or more alleles per site.

## 2. gtg-ped-map: structural haplotype labeling

### 2.1 Founder letter assignment

At startup the driver calls [`Iht::new`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/iht.rs#L172) once to build a
master inheritance template
([`map_builder.rs:1059`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1059)) that fixes the founder
alphabet for the whole run. For each founder in depth order it allocates
a fresh pair of capital letters — `(A,B)` for founder 0, `(C,D)` for
founder 1, `(E,F)` for founder 2, and so on — and assigns one letter to
the paternal homolog and one to the maternal homolog. Non-founders are
initialised with `?` slots that will be filled during the per-site walk.
The master template is never mutated: only its
[`legend()`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/iht.rs#L330) is read, to print the column header
(`Dad:A|B Mom:C|D Kid1:?|? …`) at the top of the IHT and marker
output files. Per-site bookkeeping happens on a separate object — the
driver re-invokes `Iht::new` per VCF record at
[`map_builder.rs:1111`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1111) to allocate a fresh `local_iht`
that [`track_alleles_through_pedigree`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L295) then *mutates* in
place to record which founder letter each child inherited at that
site. Two reasons a per-site copy is needed rather than reusing the
master:

1. **Each site needs its own mutable IHT vector** — that vector *is*
   the per-site output, so it cannot be shared across sites.
2. **Zygosity is per-chromosome.** The master is hard-coded to
   `ChromType::Autosome` (it only feeds the header), whereas
   `local_iht` is built with the chromosome's actual `zygosity`
   (autosome vs. chrX, decided at
   [`map_builder.rs:1086`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1086)), which changes how letters
   are laid out for males on chrX.

The two objects share a founder-letter convention only because both
flow through `Iht::new` with the same `family.founders()` /
`family.offspring()` ordering.

**No allele sequence is associated with these letters at this stage.**
That deliberate decoupling is what makes the rest of the pipeline
tractable: letter propagation becomes a pure combinatorial problem on
the pedigree graph, with the VCF acting only as an oracle for
"informative" versus "non-informative" at each coordinate. Allele
sequences re-enter the picture only in §3, inside `gtg-concordance`.

### 2.2 Informative-site detection

For every VCF record, [`track_alleles_through_pedigree`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L295)
(driver call at [`map_builder.rs:1116`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1116)) walks
individuals in ancestor-first depth order, using the depth ordering
returned by `family.get_individual_depths()`. For every
`(parent, spouse)` pair in the walk it calls
[`unique_allele`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L243) (from inside the walk at
[`map_builder.rs:315`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L315)) to test whether the parent
carries a VCF allele that the spouse does not. This is the
**informative-site** test, and it has three outcomes:

- **Dad-informative**: dad is heterozygous and mom is homozygous for the
  other allele. Dad's unique allele tags whichever paternal homolog
  (`A` or `B`) each child inherited; a child carrying dad's unique
  allele is labelled `A` on the paternal slot, a child without it is
  labelled `B` by elimination.
- **Mom-informative**: symmetric, tagging `C` or `D` on each child's
  maternal slot.
- **Non-informative**: both parents heterozygous, or both homozygous
  for the same allele. `unique_allele` returns `None` and the site
  contributes nothing to map construction. These sites are not
  discarded from the VCF — they are simply invisible to `gtg-ped-map`
  and re-enter the picture at the `gtg-concordance` stage.

The informative-site filter is not a quality filter; it is a logical
one. A non-informative site cannot, on its own, distinguish which
parental homolog a child inherited, regardless of sequencing quality.

### 2.3 Letter propagation across generations

When a child carries the parent's unique allele at an informative site,
the child's paternal (dad-informative) or maternal (mom-informative)
slot is filled with the parent's letter label.
[`get_iht_markers`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L274) (invoked from inside the walk at
[`map_builder.rs:328`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L328)) is the routine that performs the
lookup. Because the walk always processes a parent before its children,
by the time `get_iht_markers` reads a parent's letters those letters
have already been finalised.

This is what makes the method look "recursive" across generations while
remaining a single loop: a G2 individual's labels are fixed by the
founders-to-G2 iteration, and the G2→G3 iteration simply treats the G2
labels as input. Three-generation pedigrees therefore require no extra
code over two-generation ones — the same routine handles both. The
three-generation walkthrough page makes this explicit with a worked
example that includes an outside-marriage G2 founder whose letters
(`E`,`F`) are assigned by the same `Iht::new` call that issued the
original `(A,B)`/`(C,D)` pairs.

Each individual therefore ends up with exactly two letter labels — one
per homolog — drawn from the founder alphabet. The pipeline **never
constructs a joint inheritance vector across all founders**. The
inheritance state is strictly nuclear-family-local and is carried upward
one generation at a time, which is why memory and runtime scale linearly
in pedigree size rather than exponentially.

### 2.4 Sibship backfilling

Real VCFs have missing genotypes and depth-filtered sites. When a
founder has multiple children and only one child's haplotype is tagged
at a given site, [`backfill_sibs`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L804) (driver call at
[`map_builder.rs:1122`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1122)) infers the other founder
allele on the remaining children: a sibling who does not carry the
tagged allele must have inherited the other founder homolog, so the
missing `?` can be replaced with the other letter of the parent's pair.
On the toy nuclear-family simulation this is a no-op because every
informative site already tags all three children, but on real data it
is essential for keeping the map dense through noisy regions.

### 2.5 Block collapse and phase flipping

After the per-site letter trace is complete, adjacent sites with
compatible letter assignments are merged into contiguous blocks by
[`collapse_identical_iht`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L385) (driver call at
[`map_builder.rs:1191`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1191)). A block is a maximal run of
sites over which every individual carries the same letter pair; a
recombination event in any individual forces a block boundary.

Because the two letters in a founder's pair `(A,B)` are interchangeable
within any single block — nothing in the combinatorial machinery of §2.2
distinguishes "A-then-B" from "B-then-A" — neighbouring blocks can
disagree on which letter is "first" even when the underlying biology is
continuous. [`perform_flips_in_place`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L702) resolves this by
walking the block list left-to-right and, for every founder, swapping
that founder's letter pair inside a block when doing so reduces the
per-individual letter mismatch with the previous block. The driver
invokes it three times —
[`map_builder.rs:1135`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1135) (before block collapse),
[`map_builder.rs:1193`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1193) (after block collapse), and
[`map_builder.rs:1203`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1203) (after gap fill) — so that each
subsequent step sees consistently-oriented labels.

Gaps left by `?` slots (depth-filtered sites, backfilled misses) are
then filled by [`fill_missing_values`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L617) (driver call at
[`map_builder.rs:1200`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1200)) and
[`fill_missing_values_by_neighbor`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L540) (driver call at
[`map_builder.rs:1201`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1201)), which propagate letters from
flanking blocks into the unassigned interior under the assumption that
intra-block recombination is rare.

### 2.6 Short-run masking (noise suppression)

A single sequencing error in one child at one site produces an
`A → B → A` flicker in that child's per-site letter trace: one isolated
site disagreeing with otherwise-consistent flanking sites. If passed
through `collapse_identical_iht` unmodified, each flicker would open and
close a one-site block and surface as two adjacent recombination
breakpoints in `recombinants.txt` — dramatically inflating the apparent
crossover rate.

[`count_matching_neighbors`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L935) (driver call at
[`map_builder.rs:1172`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1172)) scans each individual's
per-site letter sequence and flags runs shorter than `--run` (default
10 markers) whose flanking letters agree.
[`mask_child_alleles`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L970) (driver call at
[`map_builder.rs:1187`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1187)) then sets those slots back to
`?` before block collapse, so that the masked flicker neither creates
spurious blocks nor contaminates adjacent ones. The default threshold is
a trade-off: higher `--run` values suppress more noise but also mask
genuine short gene-conversion tracts.

### 2.7 Recombination reporting

After block collapse and flip normalisation,
[`summarize_child_changes`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L673) (driver call at
[`map_builder.rs:1228`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1228)) walks each individual's block
sequence and emits every letter transition — `A → B` on the paternal
slot, `C → D` on the maternal slot — to `{prefix}.recombinants.txt`.

These are **candidate** transitions, not meiotic recombinations.
Ancestral recombinations that occurred in an earlier generation
propagate to every descendant who inherits the affected segment, so the
same physical crossover coordinate can appear in multiple rows of
`recombinants.txt`. As the repository's HAPLOTYPING.md note makes
explicit, reducing candidate transitions to unique meiotic events
requires a downstream reconciliation step that the current
implementation does not perform.

### 2.8 Outputs

`gtg-ped-map` writes three files per run, all keyed by the
`--prefix` argument:

- `{prefix}.iht.txt` — the final block map: contiguous physical
  intervals with a letter pair for every individual and a marker count.
  This is the sole input that `gtg-concordance` consumes.
- `{prefix}.markers.txt` — raw per-site letter assignments before
  collapse, useful for breakpoint refinement and for auditing decisions
  made by `mask_child_alleles`.
- `{prefix}.recombinants.txt` — the candidate letter transitions
  described in §2.7.

Nothing here contains 0/1 allele sequences — the block map is a pure
structural labelling. That is not an oversight: it is the precondition
that lets `gtg-concordance` do its job cleanly.

## 3. gtg-concordance: allele phasing and concordance QC

The driver [`main()`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/gtg_concordance.rs#L315) reads the block map with
[`parse_ihtv2_file`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/iht.rs#L606) (driver call at
[`gtg_concordance.rs:405`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/gtg_concordance.rs#L405)) and, for each block, fetches
every overlapping VCF record — including the non-informative sites that
`gtg-ped-map` skipped — via the per-site phasing loop at
[`gtg_concordance.rs:437`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/gtg_concordance.rs#L437). Each record is then phased
against the block's letter labels by the routine described below.
Records that fail depth or quality filters are routed directly to
`{prefix}.fail.vcf` at [`gtg_concordance.rs:444`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/gtg_concordance.rs#L444)
(low-quality) and [`gtg_concordance.rs:450`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/gtg_concordance.rs#L450) (no-call),
and never enter the orientation search.

### 3.1 Orientation search

For every record that clears the input filters,
[`find_best_phase_orientation`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/gtg_concordance.rs#L252) (driver call at
[`gtg_concordance.rs:454`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/gtg_concordance.rs#L454)) enumerates the `2^F`
founder-phase orientations produced by
[`Iht::founder_phase_orientations`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/iht.rs#L492) (invoked inside
`find_best_phase_orientation` at
[`gtg_concordance.rs:256`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/gtg_concordance.rs#L256)). An orientation is a choice,
for each of the `F` founders independently, of which of the founder's
two sorted VCF alleles is tagged with its first letter and which with
its second. The `2^F` enumeration is exhaustive: every possible
letter→allele correspondence consistent with the block's structural
labels is tried.

Under a given orientation, [`Iht::assign_genotypes`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/iht.rs#L442)
produces the **expected** genotype for every individual: each founder
contributes `letter1 → vcf_allele[0]`, `letter2 → vcf_allele[1]`, and
each descendant's two letters are looked up in that map to give an
expected `{a1, a2}` pair. The expected genotypes are compared to the
observed (unphased) ones by
[`compare_genotype_maps`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/gtg_concordance.rs#L213) (driver call at
[`gtg_concordance.rs:268`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/gtg_concordance.rs#L268)), which counts the number
of samples whose observed pair disagrees with the expected pair under
unphased set-equality. The orientation with the fewest mismatches wins.

### 3.2 Pass, fail, and the "impossible genotype" rule

Two outcomes are possible:

- **Zero mismatches**: the winning orientation explains every sample's
  genotype exactly. The driver uses the winning letter→allele map —
  re-applied via [`Iht::assign_genotypes`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/iht.rs#L442) on the
  passing branch at [`gtg_concordance.rs:514`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/gtg_concordance.rs#L514) — to
  rewrite the record as phased (`0|1`, paternal first) and write it to
  `{prefix}.pass.vcf` at
  [`gtg_concordance.rs:534`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/gtg_concordance.rs#L534).
- **Non-zero mismatches under every orientation**: no assignment of
  alleles to founder letters explains the whole family. The site is
  **impossible** under the current block labels. The driver takes the
  best-available orientation — `Iht::assign_genotypes` is still called
  on the failing branch at
  [`gtg_concordance.rs:487`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/gtg_concordance.rs#L487) to identify the offending
  samples — and writes the unphased record to `{prefix}.fail.vcf` at
  [`gtg_concordance.rs:507`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/gtg_concordance.rs#L507). The offending sample
  names are appended to `{prefix}.failed_sites.txt`, and per-block
  counts of sites where exactly one sample is the culprit ("singleton"
  failures) are tallied in `{prefix}.filtering_stats.txt`.

The impossible-genotype rule is the mechanism by which
`gtg-concordance` refuses to produce phased calls that the block's
structural labels cannot justify. It simultaneously catches sequencing
errors, Mendelian violations, and residual block-labelling mistakes
upstream — all three manifest as "no orientation works".

The division of labour is strict. `gtg-ped-map` stores only letters and
only at informative sites; `gtg-concordance` is the sole place where
letter→allele correspondence is computed and written out. Any
downstream tool that needs phased genotypes must consume
`{prefix}.pass.vcf`, never the `.iht.txt` structural map.

## 4. Real-world issues

### 4.1 Sequencing errors

Short-read and HiFi sequencing errors flip individual genotypes
stochastically. They enter the pipeline in two different places and are
handled by two different mechanisms.

Inside `gtg-ped-map`, a single erroneous site in one child presents as
an `A → B → A` flicker in that child's pre-collapse letter trace. The
[`count_matching_neighbors`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L935) /
[`mask_child_alleles`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L970) pair (driver calls at
[`map_builder.rs:1172`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1172) and
[`map_builder.rs:1187`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1187)) recognises runs shorter than
`--run` (default 10 markers) and masks them before
`collapse_identical_iht` runs. Without this filter, each isolated
error would surface as two adjacent recombination breakpoints.

Inside `gtg-concordance`, sequencing errors show up as sites where no
orientation satisfies the whole family. These are written to
`{prefix}.fail.vcf` by the impossible-genotype rule described in
§3.2. Per-block singleton-failure tallies in
`{prefix}.filtering_stats.txt` distinguish a single noisy sample
(random errors) from a block-wide labelling mistake (systematic
failure across many samples), because the two fail the orientation
search with very different footprints.

### 4.2 Depth and quality filtering

`extract_depth_statistics` computes per-sample mean and standard
deviation over each chromosome.
[`depth_filters`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L217) drops any site where a sample is below
`--depth` (default 10) or more than one standard deviation from its
per-sample mean. The upper bound matters as much as the lower one: an
inflated depth is the signature of a collapsed duplication or a
copy-number-variable region where the two-homolog model breaks down.
`--qual` (default 20) filters out low-confidence calls at the
VCF-record level. Both binaries apply these filters independently on
the same input VCF so that a site dropped in `gtg-ped-map` is also
dropped in `gtg-concordance`.

### 4.3 Founder-phase instability across blocks

Because `(A,B)` and `(B,A)` are indistinguishable within a single block
(§2.5), adjacent blocks can disagree on letter order even when the
underlying biology is continuous. [`perform_flips_in_place`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L702)
resolves this by pair-swapping founder letters inside a block to
minimise mismatches with the previous block. The driver calls it three
times — [`map_builder.rs:1135`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1135),
[`map_builder.rs:1193`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1193), and
[`map_builder.rs:1203`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1203) — once before collapse, once
after, and once after gap fill, so that every step operates on
consistently-oriented blocks.

### 4.4 Missing data in sibships

Multi-child sibships provide redundant information that is essential in
noisy regions. [`backfill_sibs`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L804) (driver call at
[`map_builder.rs:1122`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1122)) uses the fact that if one
child in a sibship carries founder allele X, the remaining children
must each carry either X or the other founder allele, and a majority
vote across the sibship can recover a usable label for a child whose
own site-level signal is missing. Singleton sibships (one child only)
cannot benefit from this and remain more vulnerable to dropout.

### 4.5 Ancestral vs. per-meiosis recombinations

[`summarize_child_changes`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L673) (driver call at
[`map_builder.rs:1228`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L1228)) reports every letter
transition it observes, but the same ancestral crossover in a G2
individual propagates to all of that individual's G3 descendants and
appears at the same coordinate in each of their rows. The
[three-generation walkthrough](three_generations/three_generations.md)
shows this explicitly: its
[Figure 4](three_generations/three_generations.md#4-ancestral-vs-per-meiosis-crossovers)
contrasts two transitions in the same trace — an ancestral A→B at
sites 3/4 that appears on *both* grandchildren's rows (one meiotic
event, two `recombinants.txt` entries) and a per-meiosis B→C at sites
5/6 that appears on one grandchild's row only (one meiotic event, one
entry). Reducing these shared transitions to a unique meiotic-event
count is not a well-defined operation at the block level — it
requires cross-sample breakpoint alignment and is left to downstream
analysis. The current `recombinants.txt` output is deliberately raw
so that downstream tools can choose their own reconciliation policy.

### 4.6 Chromosome X

Male hemizygosity is treated as a special case throughout
[`track_alleles_through_pedigree`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L346) and
[`Iht::new`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/iht.rs#L181). Males receive a single non-dot letter on
the X slot, derived from mom; dads cannot transmit X to sons, so the
paternal-slot label on a male child's X is left undefined; and
[`perform_flips_in_place`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9d28575a182a58891c5e446c12b8dc2c73a414c9/code/rust/src/bin/map_builder.rs#L702) skips X entirely to avoid
meaningless swaps on the hemizygous slot. The result is that the same
block-map data structure serves both autosomes and X without a second
code path, at the cost of a few explicit sex checks in the inner loops.

## 5. Reproducing the toy figures

Every figure in the companion walkthroughs was produced by
[`wiki/generate_wiki.py`](generate_wiki.py), which contains fully
deterministic, hard-coded simulations of the method on (a) a nuclear
family with one paternal recombinant child, (b) a three-generation
pedigree with an outside marriage, and (c) non-informative site phasing
with an injected sequencing error. Running the script with no arguments

```
python wiki/generate_wiki.py
```

regenerates every PNG, every walkthrough markdown file, and this
Methods document, pinning every permalink above to the current `HEAD`
commit SHA. Individual pages can be regenerated in isolation via
`python wiki/generate_wiki.py --page <name>`.
