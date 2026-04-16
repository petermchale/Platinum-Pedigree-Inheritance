# Methods — pedigree haplotype mapping with gtg-ped-map and gtg-concordance

The pipeline consists of two Rust binaries that share a common
inheritance-tracking library (`code/rust/src/iht.rs`) and are driven by a
standard PED file and a jointly-called VCF. All line numbers below refer to
commit `2acfb8b`; permalinks are given in square brackets.

## 1. Inputs

The tools consume:

- A **PED file** describing the pedigree
  ([`code/rust/src/ped.rs`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/ped.rs)).
- A **jointly-called, tabix-indexed VCF** — typically DeepVariant on HiFi —
  containing genotypes for every individual in the PED.

Only biallelic SNVs are used for map construction; indels are filtered at
read time (`is_indel` [`map_builder.rs:501`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L501)).

## 2. gtg-ped-map: structural haplotype labeling

### 2.1 Founder letter assignment

At startup, the pipeline calls `Iht::new`
([`iht.rs:172`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/iht.rs#L172)), which assigns a pair of capital letters
(`A`,`B`), (`C`,`D`), (`E`,`F`), … to the two homologs of each founder. **No
allele sequence is associated with these letters at this stage.** The
letters are pure structural placeholders whose job is to be tracked from
founders to descendants.

### 2.2 Informative-site detection

For every VCF record in a chromosome, `track_alleles_through_pedigree`
([`map_builder.rs:295`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L295)) walks individuals in
ancestor-first depth order (`family.get_individual_depths()`). For each
(parent, spouse) pair in the walk, `unique_allele`
([`map_builder.rs:243`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L243)) checks whether the parent has an
allele that the spouse does not. This is the **informative-site** test:

- **Dad-informative**: dad is heterozygous and mom is homozygous for the
  other allele → the unique allele tags whichever paternal homolog each
  child inherited.
- **Mom-informative**: symmetric.

Sites where both parents are heterozygous, or where both are homozygous for
the same allele, yield no unique allele and are skipped at this stage.

### 2.3 Letter propagation across generations

When a child carries the parent's unique allele at an informative site, the
child's paternal (for dad-informative) or maternal (for mom-informative)
slot is filled with the parent's letter label. Because the depth-ordered
walk always processes a parent before its children,
`get_iht_markers` ([`map_builder.rs:274`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L274)) simply reads
the parent's already-assigned letters when propagating to the next
generation. This is why the method looks "recursive" across generations
even though it is expressed as a single loop: a G2 individual's labels are
already finalized by the time the loop reaches G3, so the same routine
reuses them as the "founder labels" for the G2→G3 sub-problem (Component 2).

Each individual therefore ends up with exactly two letter labels, one per
homolog, drawn from the founder alphabet. The pipeline **never constructs
a joint inheritance vector across all founders**; the inheritance state is
strictly nuclear-family-local and carried upward one generation at a time.

### 2.4 Sibship backfilling

When a founder has multiple children and only one child's haplotype is
tagged at a given site, `backfill_sibs`
([`map_builder.rs:804`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L804)) infers the other founder
allele on the remaining children — a sibling who does not carry the tagged
allele must have inherited the other founder homolog. This reduces the
number of `?` slots before block collapse.

### 2.5 Block collapse and phase flipping

Adjacent sites with compatible letter assignments are merged by
`collapse_identical_iht`
([`map_builder.rs:385`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L385)). Because each founder's two
letters ('A','B') could appear in either order in a block, neighbouring
blocks may disagree on which letter is "first". `perform_flips_in_place`
([`map_builder.rs:702`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L702)) walks the block list and swaps
founder letter pairs where doing so reduces mismatches with the previous
block, enforcing a consistent phase across the chromosome. Missing values
are then filled by `fill_missing_values` and
`fill_missing_values_by_neighbor` ([`map_builder.rs:617`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L617),
[`:540`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L540)) using flanking blocks.

### 2.6 Short-run masking (noise suppression)

Isolated marker assignments that are flanked on both sides by a different
letter are the haplotype-map signature of a single discordant site, which
is almost always a sequencing error rather than a true recombination.
`count_matching_neighbors`
([`map_builder.rs:935`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L935)) scans each individual's
per-site letters and identifies runs shorter than `--run` (default 10
markers). `mask_child_alleles`
([`map_builder.rs:970`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L970)) then sets the offending slots
back to `?` so they do not survive into the final block map.

### 2.7 Recombination reporting

After block collapse, `summarize_child_changes`
([`map_builder.rs:673`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L673)) writes per-individual letter
transitions (e.g. `A → B` on the paternal slot) to
`{prefix}.recombinants.txt`. As the HAPLOTYPING.md note makes explicit,
these are candidate transitions: ancestral recombinations in an earlier
generation appear in all descendant children and must be reconciled
downstream.

### 2.8 Outputs

`gtg-ped-map` writes three files per run:

- `{prefix}.iht.txt` — the final block map: contigous blocks with a
  letter pair for every individual and a marker count.
- `{prefix}.markers.txt` — raw per-site letter assignments before
  collapse, useful for breakpoint refinement.
- `{prefix}.recombinants.txt` — putative switch points.

## 3. gtg-concordance: allele phasing and concordance QC

`gtg-concordance` consumes the `.iht.txt` block map and the original VCF.
For each block it fetches every overlapping VCF record — including the
non-informative sites that `gtg-ped-map` discarded — and attempts to phase
each one against the block's letter labels.

### 3.1 Orientation search

For every record, `find_best_phase_orientation`
([`gtg_concordance.rs:252`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/gtg_concordance.rs#L252)) enumerates the 2^F
founder phase orientations produced by `Iht::founder_phase_orientations`
([`iht.rs:492`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/iht.rs#L492)). Each orientation swaps the order of
one or more founders' letter pairs, which is equivalent to swapping which
of that founder's two VCF alleles is paired with which letter.

`assign_genotypes` ([`iht.rs:442`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/iht.rs#L442)) then produces the
**expected** genotype for every individual under that orientation: for each
founder it maps `letter1 → vcf_allele[0]`, `letter2 → vcf_allele[1]`, and
each child's two letters are looked up in that map to give an expected
allele pair. The expected genotypes are compared to the observed
(unphased) genotypes by `compare_genotype_maps`
([`gtg_concordance.rs:213`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/gtg_concordance.rs#L213)), and the orientation with
the fewest mismatches wins.

### 3.2 Passing, failing, phased output

If the winning orientation yields **zero** mismatches, the site is
consistent with the haplotype map. `gtg-concordance` rewrites the
genotypes as phased (`0|1`, with the paternal allele first) using the
winning letter→allele map and writes the record to `{prefix}.pass.vcf`
([`gtg_concordance.rs:509`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/gtg_concordance.rs#L509)).

If **no** orientation yields zero mismatches, the site is "impossible"
under the current block labels and is routed to `{prefix}.fail.vcf`
([`gtg_concordance.rs:507`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/gtg_concordance.rs#L507)). The list of offending
samples per site is appended to `{prefix}.failed_sites.txt`, and
singleton failures (one sample mismatching in isolation) are tallied per
block in `{prefix}.filtering_stats.txt`.

The pipeline's division of labour is deliberate and strict: `gtg-ped-map`
stores **only letters** and **only at informative sites**, while
`gtg-concordance` is the sole place where letter→allele correspondence is
computed and written out. Any downstream tool that needs phased genotypes
must consume `{prefix}.pass.vcf`, not the `.iht.txt` structural map.

## 4. Real-world issues

### 4.1 Sequencing errors

Short-read or HiFi sequencing errors flip individual genotypes stochastically.
In the pre-collapse letter trace, a single erroneous site in one kid looks
like an `A → B → A` flicker surrounded by otherwise consistent blocks. The
`count_matching_neighbors` / `mask_child_alleles` pair
([`map_builder.rs:935`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L935),
[`:970`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L970)) recognises runs shorter than `--run` (default
10 markers) and drops them before `collapse_identical_iht` runs. Without
this filter, each isolated sequencing error would present as a pair of
adjacent recombination breakpoints in `recombinants.txt`, dramatically
inflating the apparent crossover rate.

At the concordance stage, sequencing errors show up as sites where no
orientation satisfies the full family and are written to `fail.vcf`
(Component 3, Panel D). Singleton failures are surfaced separately because
a single sample consistently failing across many sites in one block is
diagnostic of a block-labeling error rather than random noise.

### 4.2 Depth and quality filtering

`extract_depth_statistics` computes per-sample mean/std dev over each
chromosome, and `depth_filters`
([`map_builder.rs:217`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L217)) drops any site where a sample is
below `--depth` (default 10) or more than one standard deviation from its
per-sample mean — the upper bound catches copy-number-inflated regions.
`--qual` (default 20) filters out low-confidence calls at the VCF-record
level. Both binaries apply these filters independently.

### 4.3 Founder-phase instability across blocks

Because the two letters in a founder's pair are interchangeable within any
single block, the labels on adjacent blocks can disagree by a "flip" even
when the underlying biology is continuous. `perform_flips_in_place`
([`map_builder.rs:702`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L702)) resolves this by comparing each
block to its left neighbour and swapping founder letter pairs to minimise
mismatches. This is called twice — once before and once after block
collapse — so that block merging sees consistently-oriented labels.

### 4.4 Missing data in sibships

Multi-child sibships provide redundant information that is essential in
noisy regions: `backfill_sibs`
([`map_builder.rs:804`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L804)) uses the fact that if one child
in a sibship carries founder allele X, the other children must carry
either X or the other founder allele, and the majority vote across the
sibship can recover a usable label for a child whose own site-level signal
is missing.

### 4.5 Ancestral vs. per-meiosis recombinations

`summarize_child_changes` reports every letter transition, but the same
ancestral crossover in a G2 individual will propagate to all of that
individual's G3 descendants, producing apparent transitions in each of
them at the same coordinate. As HAPLOTYPING.md notes, downstream filtering
is required to collapse these shared transitions into a unique meiotic
event count. The current implementation exposes the raw transitions and
leaves the reconciliation to the analyst.

### 4.6 Chromosome X

Male hemizygosity on chromosome X is handled as a special case throughout
`track_alleles_through_pedigree` ([`map_builder.rs:346-373`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/bin/map_builder.rs#L346))
and `Iht::new`
([`iht.rs:181-185`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/2acfb8b0e9dadbcd707e9adbf1a546ef91ff145e/code/rust/src/iht.rs#L181)). Males receive a single non-dot
letter on the X slot derived from mom; dads cannot transmit X to sons;
the flip routine skips X entirely to avoid meaningless swaps on the
hemizygous slot.

## 5. Reproducing the toy figures

The three multi-panel figures in this manuscript were produced by
`wiki/generate_wiki.py`, which contains fully
deterministic, hard-coded simulations of the method on (a) a nuclear
family with one paternal recombinant child, (b) a three-generation
pedigree with an outside marriage, and (c) non-informative site phasing
with an injected sequencing error. Running the script with no arguments
regenerates all three PNGs and this Methods document alongside it in
`wiki/`.
