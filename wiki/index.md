# Platinum Pedigree Inheritance — wiki

This wiki accompanies the `gtg-ped-map` (`map_builder.rs`) and
`gtg-concordance` (`gtg_concordance.rs`) Rust binaries in this
repository. They run sequentially: `gtg-ped-map` is invoked first on
a PED file plus a jointly-called VCF and emits a per-block
founder-letter map (`{prefix}.iht.txt`); `gtg-concordance` is then
invoked on that map together with the same VCF, and emits two output
VCFs — a phased `{prefix}.pass.vcf` for records the block labels
fully explain, and an unphased `{prefix}.fail.vcf` for records they
cannot. The three walkthrough pages below isolate the ideas that make
that pipeline work, on toy pedigrees small enough that every
intermediate table can be read off by hand.

The structure is inspired by Andrej Karpathy's "LLM Wiki" pattern —
see his [tweet](https://x.com/karpathy/status/2040572272944324650) and
accompanying [gist](https://gist.github.com/karpathy/442a6bf555914893e9891c11519de94f)
— in which a repository ships a thin, hand-curated catalog page that
links out to self-contained topic pages, each bundled with its own
assets and regenerable from source. Here the "source" is
[`wiki/generate_wiki.py`](generate_wiki.py) plus the pinned Rust
commit. Every Rust permalink in the walkthroughs is pinned to the SHA
captured when the wiki was last regenerated; running

```
python wiki/generate_wiki.py
```

refreshes every PNG, every markdown page, and every permalink.
Individual pages can be regenerated in isolation via
`python wiki/generate_wiki.py --page <name>`.

## Walkthrough pages (recommended reading order)

1. [Structural haplotype mapping in a nuclear family](nuclear_family/nuclear_family.md)
   — the longest and densest walkthrough, and the one in which
   `gtg-ped-map`'s algorithm is derived end-to-end on a two-founder,
   three-child pedigree with a single paternal recombinant. The page
   works sequentially through ground-truth founder haplotypes, the
   unphased VCF input, informative-site detection with founder-letter
   tagging and within-block haplotype inference, the expansion of
   per-site labels into linkage blocks by minimising recombinants, an
   equivalent pairwise-comparison view of the same block collapse, and
   the short-run masking and gap-fill passes that suppress genotyping
   noise. Read this first; the other two pages assume its algorithmic
   vocabulary.

2. [Structural haplotype mapping across three generations](three_generations/three_generations.md)
   — extends the nuclear-family pedigree by marrying Kid3 (the
   paternal recombinant from that page) to an outside-marriage
   founder **Spouse**, adding two grandchildren **GK1** and **GK2**.
   The focus is deliberately narrow: the page does not re-derive the
   per-site mechanics, which generalise verbatim, but addresses two
   things that only become visible once a third generation is in the
   pedigree. First, `gtg-ped-map` does not bolt a separate G2→G3 pass
   onto the G1→G2 pass — letter assignment is a single ancestor-first
   walk over (parent, spouse, child) triples in depth order, and the
   page works through the depth assignments, the triples this
   pedigree visits, and the relation between this construction and
   the Lander–Green inheritance vector. Second, a grandchild's letter
   trace can carry two qualitatively different kinds of crossover —
   an **ancestral** one inherited unchanged from one of Kid3's
   homologs and therefore shared across every descendant of the
   affected segment, versus a **de novo** one introduced in Kid3's
   meiosis and therefore visible in only one grandchild — and the
   page spells out the two-grandchild diagnostic that tells them
   apart.

3. [Phasing alleles consistently with the haplotype map](concordance/concordance.md)
   — closes the pipeline. `gtg-ped-map` emits only founder letters at
   informative sites and never reconstructs the 0/1 allele sequence
   of any haplotype; that job belongs to `gtg-concordance`, which
   re-reads the VCF for every IHT block and phases **every** variant
   using the block's **letter map** — the pair of founder letters
   each individual carries inside that block, as recorded in
   `{prefix}.iht.txt` (e.g. `A|C` for a child who inherits the `A`
   homolog paternally and the `C` homolog maternally). The
   walkthrough covers how
   non-informative sites inside a block — sites where both parents
   are heterozygous and `gtg-ped-map` declined to touch them — are
   still phaseable via the `2^F` founder-phase orientation search,
   which enumerates every letter→allele correspondence consistent
   with the block's structural labels and keeps the one that
   minimises per-sample mismatches. `gtg-concordance`'s main loop
   then splits on that mismatch count — the two **branches** (the
   `if`/`else` arms of an `if mismatches_nonempty { ... } else { ...
   }` check) are the "passing branch" for zero mismatches and the
   "failing branch" for one or more. A clean site exercises the
   passing branch, which rewrites the record as `p|m` and emits it to
   `{prefix}.pass.vcf`; a sister site with an injected sequencing
   error exercises the failing branch, along with the per-block
   singleton tally written to `{prefix}.filtering_stats.txt` and the
   debug-level per-record genotype dump that localises the offending
   sample. An
   appendix walks `Iht::assign_genotypes` — the pure function the
   orientation search calls once per orientation — line by line.

## Conventions

- **Filenames describe content, not sequence.** Reading order is
  defined above; subdirectories use short descriptive names
  (`nuclear_family/`) instead of numeric prefixes.
- **Figures live next to the page that references them.** Each
  walkthrough subdirectory holds the markdown plus its PNG figures,
  so a page is self-contained and can be moved or renamed as a unit.
- **Rust permalinks are pinned to a commit SHA.** Regenerating the
  wiki refreshes every permalink to the current `HEAD`.
