# Pedagogical Visualization of gtg-ped-map and gtg-concordance

## Context

The user is preparing a manuscript (Bioinformatics-style journal) about the haplotype-mapping pipeline in this repository and wants pedagogical ASCII / matplotlib figures that walk a reader through the method step-by-step. Two Rust binaries implement the method:

- `gtg-ped-map` (`code/rust/src/bin/map_builder.rs`) — identifies informative sites, assigns founder haplotype labels (A,B,C,D,…) to each kid's two haplotypes, merges adjacent consistent blocks, and writes a structural inheritance map (`.iht.txt`). Only founder **labels**, not allele sequences, appear in this map.
- `gtg-concordance` (`code/rust/src/bin/gtg_concordance.rs`) — reads the structural map and, at *every* VCF site in each block (including non-informative sites), searches the 2^F founder phase orientations to find the one that explains the observed genotypes. The best orientation provides both (a) a failure flag when no orientation explains an offspring genotype and (b) a phased genotype when one does.

The goal is a Python script (pure simulation — no Rust binary calls) producing three PNG figures for Adobe Illustrator plus a Methods section describing the real-world complications.

## Key answers the simulation must bake in

1. **Does gtg-ped-map reconstruct 0/1 allele sequences per haplotype?** No. `map_builder.rs` only assigns founder letters (A/B/C/D/…) at informative sites and merges blocks. Allele-to-label correspondence at non-informative sites is exclusively the job of `gtg-concordance` (`iht.rs:442 assign_genotypes` + `gtg_concordance.rs:252 find_best_phase_orientation`).
2. **Multi-generation tracking — recursive or single-step inheritance vector?** Recursive / generational. `track_alleles_through_pedigree` (`map_builder.rs:295`) walks individuals in ancestor-first depth order (`family.get_individual_depths()`), and `get_iht_markers` (`map_builder.rs:274`) reads the parent's already-assigned label — which for G2 parents was set during the previous pass of the same function on an earlier generation. Each child carries only two haplotype characters; there is no joint inheritance vector over all founders.
3. **Impossible-genotype logic.** `find_best_phase_orientation` enumerates `2^F` founder orientations (`iht.rs:492`), calls `assign_genotypes` (`iht.rs:442`) to propagate founder alleles to children via the block's letter labels, and compares to the observed genotypes. If no orientation yields zero mismatches, the site fails (sequencing error, Mendelian violation, or block mislabeling). Otherwise the winning orientation provides the phased child genotypes.

## Deliverables

A single Python script `code/python/visualize_haplotype_method.py` (new file) that, when run, writes three PNGs into a `figures/` subdirectory. No existing Python plotting code is present in `code/` — this will be a standalone educational script using `matplotlib` with a monospace figure style so the panels read as ASCII diagrams but are vector-quality.

I will build the script **one component at a time**, stopping to show the user each PNG + script diff before moving to the next component.

### Component 1 — Nuclear family (dad, mom, 3 kids)
**File:** `figures/component1_nuclear_family.png`

Multi-panel figure illustrating `gtg-ped-map` on a single nuclear family.

- Panel A: True haplotypes. Dad = (A,B), Mom = (C,D). Hard-coded SNV alleles at ~12 sites (0/1) for A, B, C, D. Kid1 inherits (A,C), Kid2 inherits (B,D), Kid3 inherits (A-prefix+B-suffix, C) — i.e. paternal recombinant.
- Panel B: Unphased VCF view. Each row an individual; each column a site; each cell an unphased genotype (0/0, 0/1, 1/1).
- Panel C: Dad-informative sites highlighted (dad het, mom hom). At each such site, identify which kids carry the dad-unique allele → label paternal slot with A or B (mirrors `unique_allele` at `map_builder.rs:243` and `track_alleles_through_pedigree` at `map_builder.rs:295`).
- Panel D: Same for mom-informative sites → maternal slot labeled C or D.
- Panel E: Combined per-kid letter assignments at each informative site. Show run-length collapse into blocks (mirrors `collapse_identical_iht` at `map_builder.rs:385`) and mark the Kid3 A→B switch as a putative recombination.
- Panel F: Truth-vs-deduced comparison, showing the deduced blocks reproduce the true transmitted blocks; include a caption explicitly stating that the inheritance map carries **letters, not 0/1 sequences**.

Each panel has a plain-English caption and inline comments in the script pointing to source lines (`map_builder.rs:243`, `:274`, `:295`, `:385`, plus GitHub permalinks using the current HEAD SHA).

### Component 2 — Three-generation pedigree with outside marriage
**File:** `figures/component2_three_generations.png`

- Panel A: G1 founders Dad (A,B), Mom (C,D) transmit to G2 child Kid2 (B,D). Kid2 marries an outside founder Spouse (E,F). They have two G3 grandkids with known transmitted combinations (e.g. Grandkid1 = (B,E), Grandkid2 = (recombinant B/A from Kid2's gametogenesis of her own haplotypes? — actually Kid2 only has B+D, so a G3 recombinant would be B-prefix + D-suffix from Kid2. Use that.)
- Panel B: Show `gtg-ped-map`'s depth-ordered walk. Annotate that the function `track_alleles_through_pedigree` processes G1→G2 first (labeling Kid2 as B,D), then G2→G3 using Kid2's now-resolved labels as "parent labels" — so the same single-site procedure recurses across generations.
- Panel C: Show the resulting letter map for grandkids drawn from {A,B,C,D,E,F}. Emphasize that each grandkid's slot still holds just two letters; there is no inheritance vector over all six founders.
- Panel D: Caption contrasting this with a Lander–Green style joint inheritance vector; reference `map_builder.rs:295-383` and the depth-first walk `family.get_individual_depths()`.

### Component 3 — gtg-concordance phasing at non-informative sites
**File:** `figures/component3_concordance.png`

Uses the Component 1 nuclear family. Picks sites where both parents are het (non-informative — dad and mom both 0/1) and shows `gtg-concordance`'s logic.

- Panel A: Pull the block's letter map (A,B | C,D | kids' letters) from Component 1. State: "the block says Kid3's paternal slot is A for the left half and B for the right half".
- Panel B: Enumerate the 2^4 founder-phase orientations (`iht.rs:492 founder_phase_orientations`). For each orientation, evaluate the expected kid genotypes with `assign_genotypes` (`iht.rs:442`) — show a small side-by-side "seen vs expected" table à la `format_genotype_maps` (`gtg_concordance.rs:151`).
- Panel C: Highlight the winning orientation (zero mismatches). Show the phased kid genotypes this yields (matches Table S5 of Eberle 2017 — parents 0/1 × 0/1 yields four distinct phased outcomes in kids depending on which founder allele each haplotype letter maps to).
- Panel D: Now inject a simulated sequencing error flipping one kid's genotype at one site. Show that *no* orientation yields zero mismatches → the site is written to `.fail.vcf` (`gtg_concordance.rs:509`). Caption this as the "impossible genotype" rule.
- Panel E: Truth-vs-deduced phased genotypes confirming the phasing deductions match the simulated truth.

Captions cite `gtg_concordance.rs:252 find_best_phase_orientation`, `iht.rs:442 assign_genotypes`, `iht.rs:492 founder_phase_orientations`, and `gtg_concordance.rs:437-538` (the passing/failing split).

### Component 4 — Methods section (markdown block inside the script)
The script will, at the end, emit `figures/methods.md` containing a plain-English Methods section that:
- Describes informative-site detection, letter-label propagation, block collapse, phase flipping, missing-value fill, and short-run masking (`map_builder.rs:385,540,617,702,970`).
- Describes `gtg-concordance`'s founder-orientation search and allele-to-label mapping.
- Discusses real-world issues: sequencing errors mimicking recombinations (addressed by `count_matching_neighbors` / `mask_child_alleles` at `map_builder.rs:935,970`), depth-based filtering (`map_builder.rs:217`), per-block phase instability handled by `perform_flips_in_place` (`map_builder.rs:702`), allele backfilling in multi-child sibships (`backfill_sibs` at `map_builder.rs:804`), and Mendelian-violation / low-quality filtering in concordance (`gtg_concordance.rs:442-452`).
- States the boundary of responsibility: structural letter map = `gtg-ped-map`; allele-sequence phasing at every site = `gtg-concordance`.

## Critical files to reference from the script

| File | Lines | Role |
|------|-------|------|
| `code/rust/src/bin/map_builder.rs` | 243–268 | `unique_allele` — informative-site detection |
| `code/rust/src/bin/map_builder.rs` | 274–293 | `get_iht_markers` / `find_valid_char` — parent-label lookup |
| `code/rust/src/bin/map_builder.rs` | 295–383 | `track_alleles_through_pedigree` — depth-first label propagation |
| `code/rust/src/bin/map_builder.rs` | 385–434 | `collapse_identical_iht` — block merging |
| `code/rust/src/bin/map_builder.rs` | 702–788 | `perform_flips_in_place` — founder-phase consistency across blocks |
| `code/rust/src/bin/map_builder.rs` | 804–914 | `backfill_sibs` — inferring the other founder allele in sibships |
| `code/rust/src/bin/map_builder.rs` | 935–987 | `count_matching_neighbors` / `mask_child_alleles` — sequencing-error / short-run masking |
| `code/rust/src/bin/gtg_concordance.rs` | 252–305 | `find_best_phase_orientation` |
| `code/rust/src/bin/gtg_concordance.rs` | 437–538 | pass/fail split and phased VCF emission |
| `code/rust/src/iht.rs` | 172–198 | `Iht::new` — founder-letter assignment scheme (A,B / C,D / …) |
| `code/rust/src/iht.rs` | 442–489 | `assign_genotypes` — founder→child allele propagation |
| `code/rust/src/iht.rs` | 492–524 | `founder_phase_orientations` — 2^F enumeration |

GitHub permalinks will use the form `https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/blob/<sha>/code/rust/src/...#L<line>` with `<sha>` captured from `git rev-parse HEAD` at script-run time.

## Implementation approach

- One Python file, `code/python/visualize_haplotype_method.py`.
- Stdlib + `matplotlib`; no numpy needed beyond what matplotlib pulls in.
- Each component = one top-level function (`component_1_nuclear_family()`, `component_2_three_generations()`, `component_3_concordance()`, `emit_methods_section()`).
- A small `draw_text_grid()` helper renders monospace-font text tables as matplotlib axes so panels look like ASCII art but export cleanly as vector-friendly PNG.
- Simulation data is fully deterministic and hard-coded (no RNG seeds) so the paper figures are reproducible.
- After each component is implemented, **stop and show the user** the generated PNG and the function body, before starting the next.

## Verification

After each component:
1. Run `python code/python/visualize_haplotype_method.py --component N`.
2. Open the PNG and verify visually that panels match the described pedigree.
3. For Component 3, print (to stdout alongside the figure) the deduced vs. truth phased genotypes, confirming exact match at non-error sites and an "impossible" flag at the error site.

After all components:
1. `python code/python/visualize_haplotype_method.py` regenerates all four artifacts.
2. (Optional sanity check the user may want to run) Build the Rust binaries, craft a tiny VCF+PED that matches the simulated Component 1 scenario, and confirm that `gtg-ped-map` + `gtg-concordance` produce letter blocks and phased genotypes consistent with the simulation. This is outside the script but listed here so the user knows it is possible.

## Resolved choices (confirmed with user)

- **Script location:** `code/python/visualize_haplotype_method.py` (new directory).
- **Figure style:** Monospace matplotlib text panels — each panel is a matplotlib axes filled with `pyplot.text` in a monospace font so the content reads as ASCII art but exports as a high-DPI vector-friendly PNG ready for Adobe Illustrator.
- **Table S5 PDF:** Not fetched. Component 3 is built directly from `iht.rs:442` and `gtg_concordance.rs:252`; the Methods section cites Eberle 2017 Table S5 textually.
