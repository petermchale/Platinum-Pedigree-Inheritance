# Platinum Pedigree Inheritance — wiki

This wiki accompanies the `gtg-ped-map` and `gtg-concordance` Rust
binaries in this repository. It is organised as a small catalog of
self-contained pages that walk through the pipeline in increasing
order of complexity, plus a full Methods write-up for the manuscript.

The structure is inspired by Andrej Karpathy's "LLM Wiki" pattern —
see his [tweet](https://x.com/karpathy/status/2040572272944324650) and
accompanying [gist](https://gist.github.com/karpathy/442a6bf555914893e9891c11519de94f)
— in which a repository ships a thin, hand-curated catalog page that
links out to self-contained topic pages, each bundled with its own
assets and regenerable from source. Here the "source" is the following Python
script plus the pinned Rust commit:
[`wiki/generate_wiki.py`](generate_wiki.py).

Every page that references Rust source uses a permalink pinned to the
commit SHA captured when the wiki was last regenerated. The whole wiki
is reproducible by running

```
python wiki/generate_wiki.py
```

## Walkthrough pages (recommended reading order)

1. [Nuclear family — structural haplotype mapping](nuclear_family/nuclear_family.md)
   — the simplest case: two founders and three children, one paternal
   recombinant. Introduces founder-letter labelling, informative-site
   detection, block collapse, and recombination reporting.
2. [Three-generation pedigree with an outside marriage](three_generations/three_generations.md)
   — extends the nuclear family by marrying Kid2 to a fresh founder
   (Spouse) and adding a G3 sibship. Shows how the same single loop
   handles G2→G3 via ancestor-first depth ordering, without ever
   constructing a joint inheritance vector over all founders.
3. *(coming soon)* **gtg-concordance phasing at non-informative
   sites** — closes the pipeline by mapping founder letters back to
   VCF alleles, with the "impossible genotype" rule that routes
   sequencing errors to `fail.vcf`.

## Reference pages

- [Methods](methods.md) — the manuscript-style write-up with the
  real-world caveats (depth filtering, phase instability, sibship
  backfilling, chromosome X) that the walkthrough pages deliberately
  skip.

## Conventions

- **Filenames describe content, not sequence.** Reading order is
  defined above; subdirectories use short descriptive names
  (`nuclear_family/`) instead of numeric prefixes.
- **Panels live next to the page that references them.** Each
  walkthrough subdirectory holds the markdown plus its PNG panels, so
  a page is self-contained and can be moved or renamed as a unit.
- **Rust permalinks are pinned to a commit SHA.** Regenerating the
  wiki refreshes every permalink to the current `HEAD`.
