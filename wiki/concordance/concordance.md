# Phasing alleles consistently with the haplotype map

This page is part of the [wiki](../index.md) and picks up where the
[nuclear-family walkthrough](../nuclear_family/nuclear_family.md) left
off. `gtg-ped-map` emits only founder letters, and only at informative
sites; it never reconstructs the 0/1 allele sequence of any haplotype.
That job belongs to `gtg-concordance`, which re-reads the VCF for each
IHT block and phases **every** variant using the block's letter map.
All line numbers refer to commit `9a53789`. As in the other
walkthrough pages, each function link is followed by its call site in
the driver — `main()` in
[`gtg_concordance.rs`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L315) — so you can step through
the driver source in parallel with this walkthrough.

The toy simulation reuses the [left half of the nuclear-family block](https://github.com/petermchale/Platinum-Pedigree-Inheritance/blob/main/wiki/nuclear_family/fig4_2.png)
(Kid1=(A,C), Kid2=(B,D), Kid3=(A,C)) and adds two sites where both
parents are heterozygous. These are NON-informative for `gtg-ped-map`
because [`unique_allele`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/map_builder.rs#L243) returns `None` at each
of them, but `gtg-concordance` still has to phase them. One of the two
sites carries an injected sequencing error so both the clean-pass
(`pass.vcf`) and error-quarantine (`fail.vcf`) code paths are
exercised. Everything below is reproducible by running

```
python wiki/generate_wiki.py --page concordance
```

which regenerates both the figure PNGs referenced here and this
markdown file itself.

## 1. Ingesting inheritance blocks and their variants

The driver reads the `{prefix}.iht.txt` file produced by
`gtg-ped-map` using
[`parse_ihtv2_file`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/iht.rs#L606) (driver call at
[`gtg_concordance.rs:405`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L405)). The parser walks
the file line by line: the first three whitespace-separated fields of
each row are the `chrom`, `start`, and `end` of the block's BED
interval, the trailing columns carry the per-block marker count, and
the middle columns — one per individual, in the header order — hold
the two letter labels separated by `|` or `/`. Using `founder_count`
(taken from the pedigree) as the split point, the first chunk of
individuals is stored in the block's `founders` map and the rest in
its `children` map, so every block ends up as an `IhtVec { bed,
iht: Iht { founders, children }, count }` entry. Inside this block
the letters are constant; the concrete labels used by the toy
simulation appear in the "Block letter labels" panel of Figure 2.

The driver then iterates over the returned `Vec<IhtVec>` one block at
a time. For each block it converts the BED coordinates to the
`(chrom_id, start, end)` triple expected by `rust-htslib` and calls
`reader.fetch(...)` at
[`gtg_concordance.rs:427`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L427) to position the VCF
reader on the first record that falls inside the block — silently
skipping the block if the fetch fails (e.g. the contig is absent from
the VCF index). Every record returned by the resulting
`reader.records()` iterator is fed through
[`parse_vcf_record`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L116) at
[`gtg_concordance.rs:440`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L440), which pulls each
sample's depth and (unphased, index-sorted) allele vector into a
`HashMap<String, (depth, Vec<GenotypeAllele>)>`. Low-quality records
(`record.qual() < args.qual`) and records with missing alleles are
short-circuited straight to `{prefix}.fail.vcf` at
[`gtg_concordance.rs:444`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L444) and
[`gtg_concordance.rs:450`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L450); everything else is
handed to the per-site phasing step —
[`find_best_phase_orientation`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L252) at
[`gtg_concordance.rs:454`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L454) — including variants
that `gtg-ped-map` could not use because neither parent has a unique
allele.

## 2. Non-informative sites inside a block

![Figure 1 — Non-informative sites inside a block](fig1.png)

The two sites in Figure 1 are both homozygous-absent for informative
patterns: dad is `0/1` and so is mom. `gtg-ped-map`'s
[`unique_allele`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/map_builder.rs#L243) test therefore returns `None`
at both sites, so neither contributes to block construction. They
still enter `gtg-concordance`'s phasing loop; the per-site machinery
described below is what turns their unphased genotypes into either a
phased `pass.vcf` record or a quarantined `fail.vcf` record.

Site `N2` carries an **injected sequencing error**: Kid1 is reported
as `1/1` even though the simulation truth is `0/1`. This is the case
Figure 3 covers, where exhaustive enumeration finds no consistent phase.

## 3. Deducing variant phase by exhaustive enumeration at a clean site

![Figure 2 — Deducing variant phase by exhaustive enumeration at site N1 (clean pass)](fig2.png)

At every record inside a block,
[`find_best_phase_orientation`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L252) (driver call at
[`gtg_concordance.rs:454`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L454)) enumerates the
`2^F=4` letter-allele mappings produced by
[`Iht::founder_phase_orientations`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/iht.rs#L492) (invoked
inside `find_best_phase_orientation` at
[`gtg_concordance.rs:256`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L256)). Each mapping
chooses which of dad's two sorted VCF alleles is tagged `A` vs `B`
and which of mom's two is tagged `C` vs `D` — recall that the
founders' letters are *unphased* (A/B, C/D), whereas each kid's
letter pair is *phased* (paternal letter | maternal letter) by
construction of the block. Under a given mapping,
[`Iht::assign_genotypes`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/iht.rs#L442) (driver call at
[`gtg_concordance.rs:487`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L487) on the failing branch
and [`gtg_concordance.rs:514`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L514) on the passing
branch) reads each kid's phased letter pair as a phased VCF genotype.
That genotype is then unphased (alleles sorted) and compared against
the observed unphased kid genotype by
[`compare_genotype_maps`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L213) (driver call at
[`gtg_concordance.rs:268`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L268)), which counts how
many kids disagree.

At site `N1`, exactly one of the four mappings explains every kid
simultaneously; the three others each force a mismatch somewhere.
Under the winning mapping, the kids' phased letter pairs immediately
give the phased `p|m` genotypes shown in the
"kids (expected, phased)" column.

## 4. Failing to deduce variant phase by exhaustive enumeration at an error site

![Figure 3 — Failing to deduce variant phase by exhaustive enumeration at site N2 (injected error)](fig3.png)

At site `N2` the injected error means **no** mapping produces zero
mismatches — the best any mapping can do is 1 sample(s)
disagreeing. `find_best_phase_orientation` therefore returns a
non-empty mismatch list, the driver writes the record to
`{prefix}.fail.vcf` at
[`gtg_concordance.rs:507`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L507) (alongside the
low-quality and no-call records that were already routed to fail at
[`gtg_concordance.rs:444`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L444) and
[`gtg_concordance.rs:450`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L450)), and the offending
sample names are appended to `{prefix}.failed_sites.txt`. If
exactly one sample is the culprit across the whole block, the failure
is counted as a "singleton" — a strong signal of a sequencing error
in that one sample rather than a systemic block-labelling problem.

This is the mechanism by which `gtg-concordance` filters sequencing
errors, Mendelian violations, and residual block-labelling mistakes
without ever producing a phased call that the block's structural
labels cannot justify.

## 5. Truth versus deduced phased genotypes

![Figure 4 — Truth vs deduced phased genotypes](fig4.png)

At the PASS site (`N1`), every kid's deduced paternal and maternal
phased alleles match the simulation truth exactly
(0 mismatches out of 6 phased-allele slots). At
the FAIL site (`N2`), no phased output is emitted — the record lands
in `fail.vcf` untouched and the truth row in Figure 4 is shown only to
document what `gtg-concordance` declined to commit to.

This closes the pipeline. `gtg-ped-map` (`map_builder.rs`) is a pure
structural-labelling tool that operates on informative sites only and
writes founder letters per individual per block.
`gtg-concordance` (`gtg_concordance.rs`) is a pure phasing/QC tool
that uses those structural labels to assign alleles at every site in
the block: the passing records are phased via
[`Iht::assign_genotypes`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/iht.rs#L442) and emitted to
`{prefix}.pass.vcf` at
[`gtg_concordance.rs:534`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/9a5378948944b796ffdd933454d93bac18bebb24/code/rust/src/bin/gtg_concordance.rs#L534), and the failing
records are quarantined to `{prefix}.fail.vcf` as described above.
The split is the answer to "does `gtg-ped-map` reconstruct the 0/1
allele sequence of each haplotype?": no, and deliberately — that is
handled exclusively by `gtg-concordance`.
