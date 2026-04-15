"""
Pedagogical visualizations of the haplotype-mapping method implemented in
gtg-ped-map (code/rust/src/bin/map_builder.rs) and gtg-concordance
(code/rust/src/bin/gtg_concordance.rs).

Each component produces a multi-panel PNG for inclusion in a Bioinformatics
journal manuscript. Panels use monospace pyplot text rendering so the content
reads as ASCII art but exports as high-DPI vector-friendly PNG suitable for
Adobe Illustrator.

Run:
    python wiki/generate_wiki.py --page nuclear_family
    python wiki/generate_wiki.py                     # all

Source-code cross-references in the panel captions use GitHub permalinks
pinned to the current HEAD SHA (captured at script start).
"""

from __future__ import annotations

import argparse
import subprocess
import textwrap
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

REPO = "Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance"
FIG_DIR = Path(__file__).resolve().parent


def head_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip()
    except Exception:
        return "main"


def permalink(path: str, line: int, sha: str) -> str:
    return f"https://github.com/{REPO}/blob/{sha}/{path}#L{line}"


SHA = head_sha()


# ---------------------------------------------------------------------------
# Panel helper — renders a monospace text block as a matplotlib axes.
# ---------------------------------------------------------------------------

def text_panel(
    ax: plt.Axes, # type: ignore
    title: str,
    body_lines: List[str],
    caption: str,
    highlight_spans: List[Tuple[int, int, int, str]] | None = None,
) -> None:
    """
    Draw a monospace-text panel.

    highlight_spans: list of (row_index, col_start, col_end, color) to shade
                     individual character spans inside the body text.
    """
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    ax.text(
        0.0, 0.98, title,
        ha="left", va="top",
        fontsize=11, fontweight="bold",
        transform=ax.transAxes,
    )

    # Body rendered line-by-line so we can measure row positions.
    n_lines = len(body_lines)
    top = 0.92
    bottom_body = 0.22
    if n_lines > 0:
        line_h = (top - bottom_body) / max(n_lines, 1)
    else:
        line_h = 0.0

    # We render each char at a fixed column pitch so highlight rectangles line
    # up with the monospace glyphs.
    char_w = 0.0115
    font_size = 9
    for i, line in enumerate(body_lines):
        y = top - i * line_h
        ax.text(
            0.02, y, line,
            ha="left", va="top",
            fontsize=font_size, family="monospace",
            transform=ax.transAxes,
        )

    if highlight_spans:
        for row, c0, c1, color in highlight_spans:
            y = top - row * line_h - line_h * 0.82
            x = 0.02 + c0 * char_w
            w = (c1 - c0) * char_w
            ax.add_patch(
                Rectangle(
                    (x, y), w, line_h * 0.82,
                    transform=ax.transAxes,
                    facecolor=color, edgecolor="none",
                    alpha=0.35, zorder=0,
                )
            )

    wrapped = "\n".join(textwrap.wrap(caption, width=95))
    ax.text(
        0.0, 0.18, wrapped,
        ha="left", va="top",
        fontsize=8, style="italic",
        transform=ax.transAxes,
    )


def _render_panel_image(
    body_lines: List[str],
    out_path: Path,
    highlight_spans: List[Tuple[int, int, int, str]] | None = None,
) -> None:
    """Render a monospace text block as a standalone panel PNG.

    Unlike `text_panel`, this helper emits no title and no caption — the
    figure is sized to fit exactly the body lines so it can be embedded
    inline in a markdown narrative, where the surrounding prose does the
    explanatory work.
    """
    n_lines = max(len(body_lines), 1)
    fig_w = 12.0
    fig_h = 0.30 * n_lines + 0.4
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.96, bottom=0.04)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    top = 0.96
    line_h = 0.92 / n_lines
    font_size = 12

    for i, line in enumerate(body_lines):
        y = top - i * line_h
        ax.text(
            0.01, y, line,
            ha="left", va="top",
            fontsize=font_size, family="monospace",
            transform=ax.transAxes,
        )

    if highlight_spans:
        axes_width_in = fig_w * 0.96
        char_width_in = (font_size * 0.6) / 72.0
        char_w = char_width_in / axes_width_in
        for row, c0, c1, color in highlight_spans:
            y = top - row * line_h - line_h * 0.9
            x = 0.01 + c0 * char_w
            w = (c1 - c0) * char_w
            ax.add_patch(
                Rectangle(
                    (x, y), w, line_h * 0.9,
                    transform=ax.transAxes,
                    facecolor=color, edgecolor="none",
                    alpha=0.35, zorder=0,
                )
            )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Component 1 — Nuclear family (dad, mom, 3 kids)
# ---------------------------------------------------------------------------
#
# Rust cross-references:
#   unique_allele (informative-site detection)
#     code/rust/src/bin/map_builder.rs:243-268
#   track_alleles_through_pedigree (label propagation)
#     code/rust/src/bin/map_builder.rs:295-383
#   collapse_identical_iht (block merging)
#     code/rust/src/bin/map_builder.rs:385-434
#   Iht::new (founder letter assignment A,B / C,D / …)
#     code/rust/src/iht.rs:172-198
#
# The Rust code only ever writes founder letters ('A','B','C','D') to each
# child's two haplotype slots. It does NOT reconstruct 0/1 allele sequences —
# that task is performed later by gtg-concordance (see Component 3).

NUM_SITES = 12

# True founder haplotypes as 0/1 strings over NUM_SITES.
# Sites 0-2 and 6-8 are dad-informative (dad het, mom hom).
# Sites 3-5 and 9-11 are mom-informative (mom het, dad hom).
HAP_A = "010001101010"  # dad's first haplotype
HAP_B = "101001010010"  # dad's second haplotype
HAP_C = "000101111101"  # mom's first haplotype
HAP_D = "000010111010"  # mom's second haplotype

# True transmission pattern:
#   Kid1 inherits (A, C)        — no recombinations
#   Kid2 inherits (B, D)        — no recombinations
#   Kid3 inherits (A|B, C)      — paternal recombination between sites 5 and 6
KID1_LABELS = [("A", "C")] * NUM_SITES
KID2_LABELS = [("B", "D")] * NUM_SITES
KID3_LABELS = [("A", "C")] * 6 + [("B", "C")] * 6


def _hap_to_alleles(hap: str) -> List[int]:
    return [int(c) for c in hap]


def _child_genotype(
    dad_label: str, mom_label: str, site: int
) -> Tuple[int, int]:
    paternal = _hap_to_alleles(
        {"A": HAP_A, "B": HAP_B}[dad_label]
    )[site]
    maternal = _hap_to_alleles(
        {"C": HAP_C, "D": HAP_D}[mom_label]
    )[site]
    return paternal, maternal


def _unphased_gt(a: int, b: int) -> str:
    if a == b:
        return f"{a}/{a}"
    return "0/1"


def _build_simulation() -> Dict:
    dad_hap_a = _hap_to_alleles(HAP_A)
    dad_hap_b = _hap_to_alleles(HAP_B)
    mom_hap_c = _hap_to_alleles(HAP_C)
    mom_hap_d = _hap_to_alleles(HAP_D)

    def gts(labels):
        return [_child_genotype(la, lb, i) for i, (la, lb) in enumerate(labels)]

    kid1_phased = gts(KID1_LABELS)
    kid2_phased = gts(KID2_LABELS)
    kid3_phased = gts(KID3_LABELS)

    def unphased_row(hap1, hap2):
        return [_unphased_gt(hap1[i], hap2[i]) for i in range(NUM_SITES)]

    return {
        "dad_unphased": unphased_row(dad_hap_a, dad_hap_b),
        "mom_unphased": unphased_row(mom_hap_c, mom_hap_d),
        "kid1_unphased": [_unphased_gt(a, b) for a, b in kid1_phased],
        "kid2_unphased": [_unphased_gt(a, b) for a, b in kid2_phased],
        "kid3_unphased": [_unphased_gt(a, b) for a, b in kid3_phased],
        "dad_a": dad_hap_a,
        "dad_b": dad_hap_b,
        "mom_c": mom_hap_c,
        "mom_d": mom_hap_d,
        "kid1_phased": kid1_phased,
        "kid2_phased": kid2_phased,
        "kid3_phased": kid3_phased,
        "kid_labels": {
            "Kid1": KID1_LABELS,
            "Kid2": KID2_LABELS,
            "Kid3": KID3_LABELS,
        },
    }


def _informative_sites_dad(sim: Dict) -> List[int]:
    """Sites where dad is het and mom is homozygous — the allele unique to
    dad labels whichever paternal haplotype (A or B) each kid carries.
    Mirrors unique_allele at map_builder.rs:243."""
    sites = []
    for i in range(NUM_SITES):
        dad_het = sim["dad_a"][i] != sim["dad_b"][i]
        mom_hom = sim["mom_c"][i] == sim["mom_d"][i]
        if dad_het and mom_hom:
            sites.append(i)
    return sites


def _informative_sites_mom(sim: Dict) -> List[int]:
    sites = []
    for i in range(NUM_SITES):
        mom_het = sim["mom_c"][i] != sim["mom_d"][i]
        dad_hom = sim["dad_a"][i] == sim["dad_b"][i]
        if mom_het and dad_hom:
            sites.append(i)
    return sites


def _dad_label_for_kid_at_site(sim: Dict, kid: str, site: int) -> str:
    """At a dad-informative site, deduce whether the kid carries A or B."""
    # dad's unique allele is the allele in dad that isn't in mom.
    dad_alleles = {sim["dad_a"][site], sim["dad_b"][site]}
    mom_alleles = {sim["mom_c"][site], sim["mom_d"][site]}
    unique = list(dad_alleles - mom_alleles)
    if not unique:
        return "?"
    unique = unique[0]
    # Which dad haplotype (A or B) carries this unique allele?
    if sim["dad_a"][site] == unique:
        carrier_label = "A"
    else:
        carrier_label = "B"
    # Does the kid carry the unique allele?
    kid_gt_key = {"Kid1": "kid1_unphased", "Kid2": "kid2_unphased", "Kid3": "kid3_unphased"}[kid]
    gt = sim[kid_gt_key][site]
    kid_alleles = {int(gt[0]), int(gt[2])}
    if unique in kid_alleles:
        return carrier_label
    # Kid carries the *other* dad haplotype.
    return "B" if carrier_label == "A" else "A"


def _mom_label_for_kid_at_site(sim: Dict, kid: str, site: int) -> str:
    mom_alleles = {sim["mom_c"][site], sim["mom_d"][site]}
    dad_alleles = {sim["dad_a"][site], sim["dad_b"][site]}
    unique = list(mom_alleles - dad_alleles)
    if not unique:
        return "?"
    unique = unique[0]
    if sim["mom_c"][site] == unique:
        carrier_label = "C"
    else:
        carrier_label = "D"
    kid_gt_key = {"Kid1": "kid1_unphased", "Kid2": "kid2_unphased", "Kid3": "kid3_unphased"}[kid]
    gt = sim[kid_gt_key][site]
    kid_alleles = {int(gt[0]), int(gt[2])}
    if unique in kid_alleles:
        return carrier_label
    return "D" if carrier_label == "C" else "C"


def component_1_nuclear_family(out_dir: Path) -> None:
    """Render the nuclear-family component as a dedicated wiki page.

    Panels and prose are written to `out_dir/nuclear_family/`; the
    top-level wiki index at `out_dir/index.md` references the resulting
    `nuclear_family.md` via a stable subdirectory name (no numeric
    prefix — reading order lives in the index, not in filenames).
    """
    sim = _build_simulation()
    dad_info = _informative_sites_dad(sim)
    mom_info = _informative_sites_mom(sim)

    kids = ["Kid1", "Kid2", "Kid3"]
    paternal_deduced: Dict[str, List[str]] = {k: ["?"] * NUM_SITES for k in kids}
    maternal_deduced: Dict[str, List[str]] = {k: ["?"] * NUM_SITES for k in kids}
    for s in dad_info:
        for k in kids:
            paternal_deduced[k][s] = _dad_label_for_kid_at_site(sim, k, s)
    for s in mom_info:
        for k in kids:
            maternal_deduced[k][s] = _mom_label_for_kid_at_site(sim, k, s)

    def carry_forward(labels: List[str]) -> List[str]:
        out = list(labels)
        last = "?"
        for i in range(NUM_SITES):
            if out[i] != "?":
                last = out[i]
            else:
                out[i] = last
        last = "?"
        for i in range(NUM_SITES - 1, -1, -1):
            if out[i] != "?":
                last = out[i]
            else:
                out[i] = last
        return out

    paternal_blocks = {k: carry_forward(paternal_deduced[k]) for k in kids}
    maternal_blocks = {k: carry_forward(maternal_deduced[k]) for k in kids}

    nf_dir = out_dir / "nuclear_family"
    nf_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Panel A — ground-truth founder haplotypes (no title, no site row).
    # ------------------------------------------------------------------
    body_a = [
        "Dad A:  " + " ".join(str(x) for x in sim["dad_a"]),
        "Dad B:  " + " ".join(str(x) for x in sim["dad_b"]),
        "Mom C:  " + " ".join(str(x) for x in sim["mom_c"]),
        "Mom D:  " + " ".join(str(x) for x in sim["mom_d"]),
        "",
        "True transmission:",
        "  Kid1 <- (A,C)       paternal = A, maternal = C  (no recomb)",
        "  Kid2 <- (B,D)       paternal = B, maternal = D  (no recomb)",
        "  Kid3 <- (A|B,C)     paternal recomb: A at sites 0-5, B at 6-11",
    ]
    _render_panel_image(body_a, nf_dir / "panel_a.png")

    # ------------------------------------------------------------------
    # Panel B — unphased VCF view (no title, no site row).
    # ------------------------------------------------------------------
    def _fmt_gt(g: str) -> str:
        return g[0] + g[2] if g[0] == g[2] else "01"

    body_b = [
        "Dad :   " + " ".join(_fmt_gt(g) for g in sim["dad_unphased"]),
        "Mom :   " + " ".join(_fmt_gt(g) for g in sim["mom_unphased"]),
        "Kid1:   " + " ".join(_fmt_gt(g) for g in sim["kid1_unphased"]),
        "Kid2:   " + " ".join(_fmt_gt(g) for g in sim["kid2_unphased"]),
        "Kid3:   " + " ".join(_fmt_gt(g) for g in sim["kid3_unphased"]),
        "",
        "(0/0 rendered as '00', 0/1 as '01', 1/1 as '11')",
    ]
    _render_panel_image(body_b, nf_dir / "panel_b.png")

    # ------------------------------------------------------------------
    # Panel C — merged paternal + maternal informative-site deduction.
    # Rows are explicitly labelled "Kid<n> p" / "Kid<n> m".
    # ------------------------------------------------------------------
    def mark_sites(site_list: List[int]) -> str:
        marks = ["_"] * NUM_SITES
        for s in site_list:
            marks[s] = "*"
        return "           " + " ".join(marks)

    body_c = [
        mark_sites(dad_info) + "   <- dad-informative (dad het x mom hom)",
        mark_sites(mom_info) + "   <- mom-informative (mom het x dad hom)",
        "",
        "Deduced founder-letter labels at informative sites only",
        "('.' = site not informative for that slot):",
        "",
    ]
    for k in kids:
        pat_row = f"  {k} p:    " + " ".join(
            paternal_deduced[k][i] if i in dad_info else "."
            for i in range(NUM_SITES)
        )
        mat_row = f"  {k} m:    " + " ".join(
            maternal_deduced[k][i] if i in mom_info else "."
            for i in range(NUM_SITES)
        )
        body_c.append(pat_row)
        body_c.append(mat_row)
    _render_panel_image(body_c, nf_dir / "panel_c.png")

    # ------------------------------------------------------------------
    # Panel D — collapsed blocks with paternal / maternal on separate
    # rows, highlighting Kid3's A -> B transition.
    # ------------------------------------------------------------------
    body_d = ["Deduced labels after block collapse and gap-fill:", ""]
    for k in kids:
        pat_row = f"  {k} p:    " + " ".join(
            paternal_blocks[k][i] for i in range(NUM_SITES)
        )
        mat_row = f"  {k} m:    " + " ".join(
            maternal_blocks[k][i] for i in range(NUM_SITES)
        )
        body_d.append(pat_row)
        body_d.append(mat_row)
    body_d += [
        "",
        "Kid3's paternal row switches A -> B between sites 5 and 6;",
        "this transition is written to {prefix}.recombinants.txt.",
    ]
    _render_panel_image(body_d, nf_dir / "panel_d.png")

    # ------------------------------------------------------------------
    # Panel E — truth vs deduced, paternal and maternal on separate rows.
    # ------------------------------------------------------------------
    body_e = ["Truth (T) vs deduced (D) founder-letter labels:", ""]
    mismatches = 0
    for k in kids:
        truth = sim["kid_labels"][k]
        pat_truth = [truth[i][0] for i in range(NUM_SITES)]
        mat_truth = [truth[i][1] for i in range(NUM_SITES)]
        pat_dedu = [paternal_blocks[k][i] for i in range(NUM_SITES)]
        mat_dedu = [maternal_blocks[k][i] for i in range(NUM_SITES)]
        body_e.append(f"  {k} p  T:  " + " ".join(pat_truth))
        body_e.append(f"  {k} p  D:  " + " ".join(pat_dedu))
        body_e.append(f"  {k} m  T:  " + " ".join(mat_truth))
        body_e.append(f"  {k} m  D:  " + " ".join(mat_dedu))
        body_e.append("")
        for i in range(NUM_SITES):
            if pat_truth[i] != pat_dedu[i]:
                mismatches += 1
            if mat_truth[i] != mat_dedu[i]:
                mismatches += 1
    total_slots = NUM_SITES * len(kids) * 2
    body_e.append(f"Total label mismatches: {mismatches} / {total_slots}")
    _render_panel_image(body_e, nf_dir / "panel_e.png")

    # ------------------------------------------------------------------
    # Markdown narrative.
    # ------------------------------------------------------------------
    _emit_component1_markdown(
        nf_dir / "nuclear_family.md",
        dad_info=dad_info,
        mom_info=mom_info,
        mismatches=mismatches,
        total_slots=total_slots,
    )

    print(f"[component 1] Wrote panel PNGs + markdown to {nf_dir}")
    print(f"[component 1] Dad-informative sites: {dad_info}")
    print(f"[component 1] Mom-informative sites: {mom_info}")
    print(f"[component 1] Label mismatches vs truth: {mismatches}")
    print(
        f"[component 1] Permalink example: "
        f"{permalink('code/rust/src/bin/map_builder.rs', 295, SHA)}"
    )


def _emit_component1_markdown(
    out_path: Path,
    *,
    dad_info: List[int],
    mom_info: List[int],
    mismatches: int,
    total_slots: int,
) -> None:
    """Write the Component-1 narrative that interleaves panel PNGs with
    prose explanations and Rust-source permalinks.
    """
    def link(path: str, line: int) -> str:
        return permalink(path, line, SHA)

    map_rs = "code/rust/src/bin/map_builder.rs"
    iht_rs = "code/rust/src/iht.rs"

    content = f"""\
# Structural haplotype mapping in a nuclear family

This page is part of the [wiki](../index.md) and walks through
`gtg-ped-map`'s structural labelling algorithm on the simplest possible
pedigree: a two-generation nuclear family with two founders (dad and
mom) and three children. It complements the full
[`methods.md`](../methods.md) write-up by zooming in on the per-site
mechanics and pinning each panel to the exact Rust code that implements
it. All line numbers refer to commit `{SHA[:7]}`.

The toy simulation hard-codes four founder haplotypes over
{NUM_SITES} sites and three children whose transmissions are known a
priori. Everything below is reproducible by running

```
python wiki/generate_wiki.py --page nuclear_family
```

which regenerates both the panel PNGs referenced here and this markdown
file itself.

## 1. Ground truth

![Panel A — Ground-truth founder haplotypes](panel_a.png)

Dad carries two haplotypes, arbitrarily labelled **A** and **B**; mom
carries **C** and **D**. These letter labels are assigned at startup by
[`Iht::new`]({link(iht_rs, 172)}), which hands each founder a fresh
pair of capital letters — `(A,B)`, `(C,D)`, `(E,F)`, … — *without*
associating any allele sequence with them. The letters are pure
structural placeholders whose only job is to be carried from founders
to descendants.

In this simulation:

- **Kid1** inherits `(A, C)` with no recombination.
- **Kid2** inherits `(B, D)` with no recombination.
- **Kid3** inherits `(A|B, C)` — the paternal slot crosses over between
  sites 5 and 6, so Kid3 carries dad's `A` haplotype on sites 0–5 and
  dad's `B` haplotype on sites 6–11.

The goal of `gtg-ped-map` is to recover exactly these letter
transmissions from the jointly-called VCF alone, without ever looking
at the underlying 0/1 allele sequence.

## 2. Unphased VCF input

![Panel B — Unphased VCF view](panel_b.png)

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
[`is_indel`]({link(map_rs, 501)}).

## 3. Informative-site detection and letter deduction

![Panel C — Informative-site deduction (paternal and maternal)](panel_c.png)

For each VCF record,
[`track_alleles_through_pedigree`]({link(map_rs, 295)}) walks the
pedigree in ancestor-first depth order and, for every `(parent, spouse)`
pair, calls [`unique_allele`]({link(map_rs, 243)}) to ask whether the
parent carries an allele that the spouse does not. Two cases can arise:

- **Dad-informative** (dad het × mom hom): dad's unique allele tags
  whichever paternal homolog (`A` or `B`) each child inherited. In
  this simulation these are sites `{dad_info}`.
- **Mom-informative** (mom het × dad hom): symmetric, tagging `C` or
  `D`. These are sites `{mom_info}`.

When a child carries the parent's unique allele, the child's paternal
(or maternal) slot is filled with the parent's letter; otherwise the
slot is filled with the *other* letter of that parent's pair. Because
the depth-ordered walk always processes a parent before its children,
[`get_iht_markers`]({link(map_rs, 274)}) reads the parent's already-
assigned letters when propagating to the next generation, which is
what makes the method look "recursive" across generations while being
expressed as a single loop.

Non-informative sites (both parents het, or both hom for the same
allele) contribute nothing at this stage and are rendered as `.` in the
panel. Each kid's paternal row (`p`) sits directly above its maternal
row (`m`), so the two deduction streams line up site-by-site. Kid3's
`p` row already exhibits the A→B recombination at sites 5/6, even
though no block collapse has happened yet.

## 4. Block collapse and noise filtering

![Panel D — Collapsed blocks with recombination](panel_d.png)

Several Rust routines clean the per-site letter trace up before it is
written to disk:

1. [`backfill_sibs`]({link(map_rs, 804)}) uses the fact that siblings
   must together carry both founder homologs. If exactly one child is
   tagged at a site, the others can be inferred by elimination. In this
   toy simulation every informative site already tags all three kids,
   so backfill is a no-op here, but on real data it is essential in
   noisy regions.
2. [`collapse_identical_iht`]({link(map_rs, 385)}) merges adjacent
   sites with compatible letter assignments into blocks, while
   [`fill_missing_values`]({link(map_rs, 617)}) and
   [`fill_missing_values_by_neighbor`]({link(map_rs, 540)}) fill the
   `.` gaps visible in Panel C from flanking blocks.
3. [`count_matching_neighbors`]({link(map_rs, 935)}) and
   [`mask_child_alleles`]({link(map_rs, 970)}) identify isolated runs
   shorter than `--run` (default 10 markers) and mask them back to `?`
   as likely sequencing noise, so that collapse does not invent
   spurious recombinations.
4. [`perform_flips_in_place`]({link(map_rs, 702)}) enforces consistent
   founder-letter orientation across blocks, since the two letters in a
   founder's pair are interchangeable within any single block.

After these steps, each kid's paternal and maternal slots are fully
resolved — shown on separate rows per kid in the panel above. Kid3's
highlighted A→B transition on the paternal row is emitted to
`{{prefix}}.recombinants.txt` by
[`summarize_child_changes`]({link(map_rs, 673)}).

## 5. Truth versus deduced

![Panel E — Truth vs deduced founder labels](panel_e.png)

For every kid the deduced paternal and maternal label streams match the
ground truth at every site ({mismatches} mismatches out of
{total_slots} label slots), including Kid3's recombination. The full
output of `gtg-ped-map` for this chromosome is the set of blocks shown
above plus the `recombinants.txt` entry for Kid3's switch — and
critically, **nothing else**. The block map stores only founder
letters; it does *not* store the 0/1 allele sequence of any haplotype.

Reconstructing which allele each letter represents at every VCF site
is the job of `gtg-concordance`, which will have its own wiki page
once migrated. For every block, `gtg-concordance` enumerates the
`2^F` founder-phase orientations produced by
[`Iht::founder_phase_orientations`]({link(iht_rs, 492)}), maps letters
to VCF alleles via [`assign_genotypes`]({link(iht_rs, 442)}), and picks
the orientation that minimises mismatches against the observed
genotypes. The split of responsibilities is deliberate and strict:
`gtg-ped-map` writes only letters and only at informative sites, while
`gtg-concordance` is the sole place where letter→allele correspondence
is computed and written out.
"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(content)


# ---------------------------------------------------------------------------
# Component 2 — Three-generation pedigree with outside marriage
# ---------------------------------------------------------------------------
#
# Rust cross-references:
#   track_alleles_through_pedigree walks individuals in depth order
#     code/rust/src/bin/map_builder.rs:295-383
#   family.get_individual_depths() provides the ancestor-first ordering
#     (see code/rust/src/ped.rs for the depth computation)
#   Iht::new assigns founder letters in A,B / C,D / E,F / … pairs
#     code/rust/src/iht.rs:172-198
#
# Point of this component: the method is NOT a single-shot joint inheritance
# vector over all founders. It iterates over (parent-pair, child) nuclear
# units in depth order, so G1→G2 labels resolve first and the same routine
# then uses those already-assigned labels as "parent labels" for G2→G3.
# Each descendant still carries only two characters, one per homolog.

# Spouse (outside founder) haplotypes. Designed so every site is either
# Kid2-informative (Kid2 het × Spouse hom) or Spouse-informative
# (Kid2 hom × Spouse het) at the G2→G3 level.
HAP_E = "000100010101"
HAP_F = "010000000010"

# Grandkid transmissions (paternal slot, maternal slot).
# Kid2 is the mother; Spouse is the father.
#   Grandkid1 <- (E, B)           no recombination
#   Grandkid2 <- (F, D)           no recombination
#   Grandkid3 <- (E, B|D)         recombination in Kid2's gametogenesis:
#                                 B on sites 0-5, D on sites 6-11
GRANDKID1_LABELS = [("E", "B")] * NUM_SITES
GRANDKID2_LABELS = [("F", "D")] * NUM_SITES
GRANDKID3_LABELS = [("E", "B")] * 6 + [("E", "D")] * 6


def _letter_to_allele(letter: str, site: int) -> int:
    return {
        "A": _hap_to_alleles(HAP_A),
        "B": _hap_to_alleles(HAP_B),
        "C": _hap_to_alleles(HAP_C),
        "D": _hap_to_alleles(HAP_D),
        "E": _hap_to_alleles(HAP_E),
        "F": _hap_to_alleles(HAP_F),
    }[letter][site]


def _genotype_from_labels(labels: List[Tuple[str, str]]) -> List[str]:
    out = []
    for i, (pat, mat) in enumerate(labels):
        a = _letter_to_allele(pat, i)
        b = _letter_to_allele(mat, i)
        out.append(_unphased_gt(a, b))
    return out


def component_2_three_generations(out_path: Path) -> None:
    sim = _build_simulation()
    spouse_e = _hap_to_alleles(HAP_E)
    spouse_f = _hap_to_alleles(HAP_F)

    # Kid2's genotype row (same as Component 1).
    kid2_unphased = sim["kid2_unphased"]

    # Grandkid genotypes from true transmissions.
    gk1_unphased = _genotype_from_labels(GRANDKID1_LABELS)
    gk2_unphased = _genotype_from_labels(GRANDKID2_LABELS)
    gk3_unphased = _genotype_from_labels(GRANDKID3_LABELS)

    # --- G2 → G3 informative-site analysis ---
    # At a Kid2-informative site (Kid2 het × Spouse hom), Kid2's unique allele
    # tells us whether she passed her B or her D haplotype to a given grandkid.
    # That identification uses Kid2's *already-computed* labels (B,D) — so
    # the method is recursive across generations.
    kid2_info_sites = []
    spouse_info_sites = []
    for i in range(NUM_SITES):
        kid2_het = sim["kid2_phased"][i][0] != sim["kid2_phased"][i][1]
        spouse_het = spouse_e[i] != spouse_f[i]
        if kid2_het and not spouse_het:
            kid2_info_sites.append(i)
        elif spouse_het and not kid2_het:
            spouse_info_sites.append(i)

    grandkids = ["GK1", "GK2", "GK3"]
    gk_unphased = {"GK1": gk1_unphased, "GK2": gk2_unphased, "GK3": gk3_unphased}

    # Deduce maternal labels (B or D) from Kid2-informative sites.
    # At each Kid2-informative site, Kid2's unique allele corresponds to
    # whichever of Kid2's two haplotypes (B or D) carries that allele. We use
    # the *true* B/D haplotype strings because in Component 1 we verified
    # that the Rust pass assigns exactly those labels to Kid2 at every site.
    maternal_deduced_gk: Dict[str, List[str]] = {g: ["?"] * NUM_SITES for g in grandkids}
    for s in kid2_info_sites:
        kid2_alleles = {sim["kid2_phased"][s][0], sim["kid2_phased"][s][1]}
        spouse_alleles = {spouse_e[s], spouse_f[s]}
        unique = list(kid2_alleles - spouse_alleles)
        if not unique:
            continue
        unique = unique[0]
        b_allele = _letter_to_allele("B", s)
        carrier = "B" if b_allele == unique else "D"
        for g in grandkids:
            gt = gk_unphased[g][s]
            gk_alleles = {int(gt[0]), int(gt[2])}
            if unique in gk_alleles:
                maternal_deduced_gk[g][s] = carrier
            else:
                maternal_deduced_gk[g][s] = "D" if carrier == "B" else "B"

    # Deduce paternal labels (E or F) from Spouse-informative sites.
    paternal_deduced_gk: Dict[str, List[str]] = {g: ["?"] * NUM_SITES for g in grandkids}
    for s in spouse_info_sites:
        kid2_alleles = {sim["kid2_phased"][s][0], sim["kid2_phased"][s][1]}
        spouse_alleles = {spouse_e[s], spouse_f[s]}
        unique = list(spouse_alleles - kid2_alleles)
        if not unique:
            continue
        unique = unique[0]
        carrier = "E" if spouse_e[s] == unique else "F"
        for g in grandkids:
            gt = gk_unphased[g][s]
            gk_alleles = {int(gt[0]), int(gt[2])}
            if unique in gk_alleles:
                paternal_deduced_gk[g][s] = carrier
            else:
                paternal_deduced_gk[g][s] = "F" if carrier == "E" else "E"

    def carry_forward(labels: List[str]) -> List[str]:
        out = list(labels)
        last = "?"
        for i in range(NUM_SITES):
            if out[i] != "?":
                last = out[i]
            else:
                out[i] = last
        last = "?"
        for i in range(NUM_SITES - 1, -1, -1):
            if out[i] != "?":
                last = out[i]
            else:
                out[i] = last
        return out

    paternal_blocks_gk = {g: carry_forward(paternal_deduced_gk[g]) for g in grandkids}
    maternal_blocks_gk = {g: carry_forward(maternal_deduced_gk[g]) for g in grandkids}

    truth_labels = {
        "GK1": GRANDKID1_LABELS,
        "GK2": GRANDKID2_LABELS,
        "GK3": GRANDKID3_LABELS,
    }
    mismatches = 0
    for g in grandkids:
        for i in range(NUM_SITES):
            truth = truth_labels[g][i]
            dedu = (paternal_blocks_gk[g][i], maternal_blocks_gk[g][i])
            if truth != dedu:
                mismatches += 1

    # ------------------------------------------------------------------
    # Figure layout: 3 rows × 2 cols.
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(3, 2, figsize=(16, 14))
    plt.subplots_adjust(
        left=0.03, right=0.985, top=0.96, bottom=0.03,
        wspace=0.08, hspace=0.35,
    )

    site_header = "site:  " + " ".join(f"{i:>1}" for i in range(NUM_SITES))

    # --- Panel A: pedigree sketch ---
    body_a = [
        "          G1-Dad (A,B) --- G1-Mom (C,D)",
        "                     |",
        "          +----------+----------+",
        "          |          |          |",
        "        Kid1        Kid2       Kid3        <-- G2 (resolved first)",
        "       (A,C)       (B,D)     (A|B,C)",
        "                     |",
        "                Kid2 --- Spouse (E,F)       <-- outside marriage",
        "                     |",
        "          +----------+----------+",
        "          |          |          |",
        "         GK1        GK2        GK3         <-- G3 (resolved next)",
        "        (E,B)     (F,D)    (E, B|D recomb)",
        "",
        "At G2 level, gtg-ped-map has already labeled Kid2 as (B,D).",
        "At G3 level, Kid2 is treated as the 'mother' parent whose",
        "haplotype labels are B,D; Spouse is a fresh founder labeled E,F.",
    ]
    cap_a = (
        "Panel A. Three-generation pedigree. G1 founders Dad and Mom carry "
        "the (A,B) and (C,D) letter pairs assigned by Iht::new "
        "(iht.rs:172). In Component 1 we showed that Kid2 receives labels "
        "(B,D) via track_alleles_through_pedigree (map_builder.rs:295). "
        "Kid2 then marries Spouse, an outside founder who is brand-new to "
        "the pedigree and therefore receives the next fresh letter pair, "
        "(E,F). The same labeling routine is then applied to each G3 kid."
    )
    text_panel(axes[0, 0], "A. Pedigree with outside marriage", body_a, cap_a)

    # --- Panel B: genotype table for Kid2, Spouse, and grandkids ---
    def row(label, gts):
        return label + " ".join(
            g[0] + g[2] if g[0] == g[2] else "01" for g in gts
        )
    body_b = [
        site_header,
        row("Kid2  : ", kid2_unphased),
        row("Spouse: ", [_unphased_gt(spouse_e[i], spouse_f[i]) for i in range(NUM_SITES)]),
        row("GK1   : ", gk1_unphased),
        row("GK2   : ", gk2_unphased),
        row("GK3   : ", gk3_unphased),
        "",
        "Kid2's true haplotypes (already labelled by the G1->G2 pass):",
        "  B     : " + " ".join(str(x) for x in _hap_to_alleles(HAP_B)),
        "  D     : " + " ".join(str(x) for x in _hap_to_alleles(HAP_D)),
        "Spouse's haplotypes (labelled E,F by Iht::new on this fresh founder):",
        "  E     : " + " ".join(str(x) for x in spouse_e),
        "  F     : " + " ".join(str(x) for x in spouse_f),
    ]
    cap_b = (
        "Panel B. Unphased VCF rows for the G2→G3 nuclear unit. Kid2's rows "
        "reuse the same data shown in Component 1 — the key point is that "
        "her letter labels (B,D) are already known when the Rust routine "
        "arrives at G3, because individuals are processed in ancestor-first "
        "depth order via family.get_individual_depths() in ped.rs."
    )
    text_panel(axes[0, 1], "B. Input rows for the G2→G3 pass", body_b, cap_b)

    # --- Panel C: informative-site split at G2→G3 ---
    def mark_sites(site_list: List[int]) -> str:
        marks = ["_"] * NUM_SITES
        for s in site_list:
            marks[s] = "*"
        return "        " + " ".join(marks)

    body_c = [
        site_header,
        mark_sites(kid2_info_sites) + "   <- Kid2-informative (Kid2 het x Spouse hom)",
        mark_sites(spouse_info_sites) + "   <- Spouse-informative (Spouse het x Kid2 hom)",
        "",
        "At each Kid2-informative site the unique allele tags whether",
        "the grandkid inherited Kid2's B or D homolog.",
        "At each Spouse-informative site it tags E or F.",
        "",
        "Deduced maternal labels (B or D) at Kid2-informative sites only:",
    ]
    for g in grandkids:
        body_c.append(
            "  " + g + ":  " + " ".join(
                maternal_deduced_gk[g][i] if i in kid2_info_sites else "."
                for i in range(NUM_SITES)
            )
        )
    body_c += ["", "Deduced paternal labels (E or F) at Spouse-informative sites only:"]
    for g in grandkids:
        body_c.append(
            "  " + g + ":  " + " ".join(
                paternal_deduced_gk[g][i] if i in spouse_info_sites else "."
                for i in range(NUM_SITES)
            )
        )
    cap_c = (
        "Panel C. Informative-site split at the G2→G3 level. The exact same "
        "unique_allele routine (map_builder.rs:243) runs, but the 'parent' "
        "labels it propagates are Kid2's already-resolved B/D labels plus "
        "the freshly-created E/F labels on Spouse. This is what makes the "
        "approach recursive without ever constructing a joint inheritance "
        "vector over all founders."
    )
    text_panel(
        axes[1, 0],
        "C. Recursive informative-site deduction",
        body_c, cap_c,
    )

    # --- Panel D: collapsed grandkid block map ---
    highlight_spans = []
    body_d = [
        site_header,
        "",
        "Deduced (paternal | maternal) after block-fill:",
    ]
    for g in grandkids:
        row_str = "  " + g + ":  " + " ".join(
            paternal_blocks_gk[g][i] + maternal_blocks_gk[g][i] for i in range(NUM_SITES)
        )
        body_d.append(row_str)
    body_d += [
        "",
        "GK3's maternal slot switches B -> D at the site 5/6 boundary.",
        "This reflects a crossover in Kid2's gametogenesis and is",
        "emitted to {prefix}.recombinants.txt.",
    ]
    # Highlight the switch in the GK3 row.
    gk3_row = 4  # after header + blank + 'Deduced...' + GK1 + GK2
    prefix_len = len("  GK3:  ")
    col_site5 = prefix_len + 5 * 3
    col_site6 = prefix_len + 6 * 3
    highlight_spans.append((gk3_row, col_site5, col_site5 + 2, "#ffcc66"))
    highlight_spans.append((gk3_row, col_site6, col_site6 + 2, "#ff9966"))

    cap_d = (
        "Panel D. Per-grandkid haplotype block map drawn from letters "
        "{B,D,E,F}. Note what is NOT here: there is no slot for A or C, "
        "because neither reached G3 (Kid2 carried only B and D). Each "
        "grandkid's record still holds exactly two characters — one per "
        "homolog — identical in shape to Component 1's output. The Rust "
        "code never stores a joint 6-founder inheritance vector."
    )
    text_panel(
        axes[1, 1],
        "D. Grandkid letter map + G2→G3 recombination",
        body_d, cap_d, highlight_spans=highlight_spans,
    )

    # --- Panel E: truth-vs-deduced ---
    body_e = [
        "Truth vs deduced labels (per grandkid, per site):",
        site_header,
    ]
    for g in grandkids:
        truth_row = "  " + g + "T: " + " ".join(a + b for a, b in truth_labels[g])
        dedu_row = "  " + g + "D: " + " ".join(
            paternal_blocks_gk[g][i] + maternal_blocks_gk[g][i] for i in range(NUM_SITES)
        )
        body_e.append(truth_row)
        body_e.append(dedu_row)
    body_e += [
        "",
        f"Total label mismatches: {mismatches} / {NUM_SITES * len(grandkids)}",
    ]
    cap_e = (
        "Panel E. Deduced letters match the simulated truth for all three "
        "grandkids, including the mid-chromosome B→D recombination in GK3."
    )
    text_panel(axes[2, 0], "E. Truth vs deduced", body_e, cap_e)

    # --- Panel F: recursion vs single-shot inheritance vector ---
    body_f = [
        "Recursive pass (gtg-ped-map):",
        "  for ind in family.get_individual_depths():  # ancestor-first",
        "      for spouse, child in nuclear_unit(ind):",
        "          marker = unique_allele(parent, spouse)",
        "          label  = get_iht_markers(parent)  # read parent's letters",
        "          assign(child, label)              # one slot per homolog",
        "",
        "  => Each individual ends with a 2-character label, drawn from",
        "     {A,B,C,D,E,F}. The same function handles G1->G2 and G2->G3.",
        "",
        "Contrast with a Lander-Green-style joint inheritance vector,",
        "which would encode, per site, a bit for every founder-gamete",
        "combination across the whole pedigree jointly.  gtg-ped-map does",
        "NOT do this - it is strictly nuclear-family-local, iterated",
        "in depth order.",
    ]
    cap_f = (
        "Panel F. Why the G3 pass looks 'recursive' without any explicit "
        "recursion: track_alleles_through_pedigree (map_builder.rs:295) "
        "iterates over individuals in ancestor-first depth order. When it "
        "reaches a G2 parent, that parent already carries resolved letter "
        "labels from the earlier iteration, and they serve as the 'founder "
        "labels' for the G2→G3 sub-problem. No joint optimisation over all "
        "founders is ever performed."
    )
    text_panel(axes[2, 1], "F. Recursion vs joint inheritance vector", body_f, cap_f)

    fig.suptitle(
        "Component 2 — Depth-ordered haplotype mapping into a third generation "
        "(gtg-ped-map)",
        fontsize=13, fontweight="bold",
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)

    print(f"[component 2] Wrote {out_path}")
    print(f"[component 2] Kid2-informative sites: {kid2_info_sites}")
    print(f"[component 2] Spouse-informative sites: {spouse_info_sites}")
    print(f"[component 2] Grandkid label mismatches vs truth: {mismatches}")


# ---------------------------------------------------------------------------
# Component 3 — gtg-concordance phasing at non-informative sites
# ---------------------------------------------------------------------------
#
# Rust cross-references:
#   find_best_phase_orientation (enumerates 2^F founder orientations,
#   picks the one with fewest observed/expected mismatches)
#     code/rust/src/bin/gtg_concordance.rs:252-305
#   assign_genotypes (letter -> VCF allele propagation)
#     code/rust/src/iht.rs:442-489
#   founder_phase_orientations (2^F enumeration)
#     code/rust/src/iht.rs:492-524
#   pass / fail split and phased VCF emission
#     code/rust/src/bin/gtg_concordance.rs:437-538
#
# Setup: reuse the nuclear family from Component 1 and focus on the LEFT
# block of the chromosome where kid labels are Kid1=(A,C), Kid2=(B,D),
# Kid3=(A,C). Add three EXTRA sites inside that block where both parents
# are heterozygous (0/1 x 0/1). Such sites are NON-informative for
# gtg-ped-map (neither parent has a unique allele), so they contribute
# nothing to building the block, but gtg-concordance still has to phase
# them using the block's letter labels.

# Non-informative site definitions (dad and mom both 0/1 at every site).
# Each entry gives the true founder haplotype alleles (A,B,C,D) at that site.
# Kid3 is taken from the LEFT block (Kid3 = (A,C) before the recombination).
NON_INFORMATIVE_SITES = [
    # (label, A_allele, B_allele, C_allele, D_allele)
    ("N1", 0, 1, 0, 1),   # A=0,B=1,C=0,D=1 -> Kid1 0/0, Kid2 1/1, Kid3 0/0
    ("N2", 1, 0, 1, 0),   # A=1,B=0,C=1,D=0 -> Kid1 1/1, Kid2 0/0, Kid3 1/1
    ("N3", 0, 1, 1, 0),   # A=0,B=1,C=1,D=0 -> Kid1 0/1, Kid2 0/1, Kid3 0/1
]

# Block letter labels (from Component 1, left-half of the chromosome).
BLOCK_LABELS = {
    "Kid1": ("A", "C"),
    "Kid2": ("B", "D"),
    "Kid3": ("A", "C"),
}


def _site_genotypes(a: int, b: int, c: int, d: int) -> Dict[str, Tuple[int, int]]:
    """Compute the TRUE (paternal, maternal) phased allele pair for each kid
    from the founder haplotype alleles at one site."""
    return {
        "Dad": (a, b),
        "Mom": (c, d),
        "Kid1": (a, c),
        "Kid2": (b, d),
        "Kid3": (a, c),
    }


def _sorted_unphased(pair: Tuple[int, int]) -> Tuple[int, int]:
    return (min(pair), max(pair))


def _expected_under_orientation(
    dad_letters: Tuple[str, str],
    mom_letters: Tuple[str, str],
    dad_vcf: Tuple[int, int],
    mom_vcf: Tuple[int, int],
    block_labels: Dict[str, Tuple[str, str]],
) -> Dict[str, Tuple[int, int]]:
    """Mirrors iht.rs:442 assign_genotypes.

    `dad_letters` and `mom_letters` are the founder letter pair under the
    current orientation (e.g. ('A','B') or ('B','A')). The VCF allele pairs
    are stored sorted (min, max). letter[0] gets mapped to vcf[0], letter[1]
    gets mapped to vcf[1]. Kids' letters are then looked up in this map.
    """
    letter_to_allele: Dict[str, int] = {
        dad_letters[0]: dad_vcf[0],
        dad_letters[1]: dad_vcf[1],
        mom_letters[0]: mom_vcf[0],
        mom_letters[1]: mom_vcf[1],
    }
    expected: Dict[str, Tuple[int, int]] = {
        "Dad": (letter_to_allele[dad_letters[0]], letter_to_allele[dad_letters[1]]),
        "Mom": (letter_to_allele[mom_letters[0]], letter_to_allele[mom_letters[1]]),
    }
    for kid, (lp, lm) in block_labels.items():
        expected[kid] = (letter_to_allele[lp], letter_to_allele[lm])
    return expected


# All 2^2 = 4 orientations of dad/mom letter pairs.
ORIENTATIONS: List[Tuple[Tuple[str, str], Tuple[str, str]]] = [
    (("A", "B"), ("C", "D")),
    (("B", "A"), ("C", "D")),
    (("A", "B"), ("D", "C")),
    (("B", "A"), ("D", "C")),
]


def _count_mismatches(
    observed: Dict[str, Tuple[int, int]],
    expected: Dict[str, Tuple[int, int]],
) -> int:
    n = 0
    for key in observed:
        if _sorted_unphased(observed[key]) != _sorted_unphased(expected[key]):
            n += 1
    return n


def _fmt_gt(pair: Tuple[int, int]) -> str:
    s = _sorted_unphased(pair)
    return f"{s[0]}/{s[1]}"


def _fmt_phased(pair: Tuple[int, int]) -> str:
    return f"{pair[0]}|{pair[1]}"


def component_3_concordance(out_path: Path) -> None:
    # -------------- Per-site orientation analysis --------------
    # For each site, compute truth, observed (maybe with injected error),
    # run the 4 orientations, and record mismatches + winner.

    # Inject a sequencing error: flip Kid1's observed genotype at site N3
    # from the truth (0/1) to 1/1.
    error_injection = {"N3": {"Kid1": (1, 1)}}

    per_site = []
    for label, a, b, c, d in NON_INFORMATIVE_SITES:
        truth = _site_genotypes(a, b, c, d)

        observed = {k: _sorted_unphased(v) for k, v in truth.items()}
        if label in error_injection:
            for k, v in error_injection[label].items():
                observed[k] = _sorted_unphased(v)

        dad_vcf = _sorted_unphased(observed["Dad"])
        mom_vcf = _sorted_unphased(observed["Mom"])

        orient_results = []
        for dad_letters, mom_letters in ORIENTATIONS:
            expected = _expected_under_orientation(
                dad_letters, mom_letters, dad_vcf, mom_vcf, BLOCK_LABELS
            )
            n_mis = _count_mismatches(observed, expected)
            orient_results.append((dad_letters, mom_letters, expected, n_mis))

        best = min(orient_results, key=lambda r: r[3])
        per_site.append({
            "label": label,
            "truth": truth,
            "observed": observed,
            "orient_results": orient_results,
            "best": best,
            "error": label in error_injection,
        })

    # ------------------------------------------------------------------
    # Figure layout: 3 rows x 2 cols = 6 panels.
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(3, 2, figsize=(16, 14))
    plt.subplots_adjust(
        left=0.03, right=0.985, top=0.96, bottom=0.03,
        wspace=0.08, hspace=0.35,
    )

    # --- Panel A: block letter labels (reminder from Component 1) ---
    body_a = [
        "Block letter labels in the left half of the chromosome",
        "(from gtg-ped-map output, Component 1):",
        "",
        "  Dad  : (A, B)    <- fresh founder letters",
        "  Mom  : (C, D)    <- fresh founder letters",
        "  Kid1 : (A, C)",
        "  Kid2 : (B, D)",
        "  Kid3 : (A, C)    <- left-half, before the paternal recombination",
        "",
        "These letters are constant inside one IHT block. gtg-concordance",
        "runs over EVERY VCF site inside the block, including sites that",
        "gtg-ped-map could not use to build the block (because neither",
        "parent has a unique allele).",
    ]
    cap_a = (
        "Panel A. Input to gtg-concordance for this block: the letter map "
        "written by gtg-ped-map into {prefix}.iht.txt and parsed by "
        "parse_ihtv2_file (iht.rs:606). gtg-concordance re-reads the VCF "
        "for the same region and attempts to phase EVERY variant, not only "
        "the informative ones that built the block."
    )
    text_panel(axes[0, 0], "A. Block letter map (input to gtg-concordance)", body_a, cap_a)

    # --- Panel B: the three non-informative sites ---
    body_b = [
        "Three additional sites INSIDE the block where both parents are het:",
        "",
        "             Dad   Mom   Kid1  Kid2  Kid3    <- observed genotypes",
    ]
    for s in per_site:
        obs = s["observed"]
        tag = "  (ERROR)" if s["error"] else ""
        body_b.append(
            f"  {s['label']} :     "
            f"{_fmt_gt(obs['Dad'])}   {_fmt_gt(obs['Mom'])}   "
            f"{_fmt_gt(obs['Kid1'])}   {_fmt_gt(obs['Kid2'])}   "
            f"{_fmt_gt(obs['Kid3'])}{tag}"
        )
    body_b += [
        "",
        "At every site here, dad and mom are both 0/1:",
        "  unique_allele(dad, mom) = {} and unique_allele(mom, dad) = {}",
        "  => NON-informative for gtg-ped-map (nothing added to the block)",
        "  => still phaseable by gtg-concordance via the block letter map",
        "",
        "Site N3 carries an INJECTED sequencing error: Kid1 observed as",
        "1/1 when the simulation truth is 0/1.",
    ]
    cap_b = (
        "Panel B. Non-informative sites — the cases Component 1 could not "
        "touch. Because dad and mom are both heterozygous, unique_allele "
        "(map_builder.rs:243) returns None and gtg-ped-map emits nothing "
        "at these sites. gtg-concordance, however, visits every VCF record "
        "in the block (gtg_concordance.rs:437)."
    )
    text_panel(axes[0, 1], "B. Non-informative sites inside the block", body_b, cap_b)

    # --- Panel C: orientation enumeration at site N1 ---
    s1 = per_site[0]
    body_c = [
        f"Site {s1['label']}:  Dad=0/1  Mom=0/1  Kid1=0/0  Kid2=1/1  Kid3=0/0",
        "",
        "founder_phase_orientations() (iht.rs:492) enumerates 2^F=4",
        "orientations.  For each one, assign_genotypes() (iht.rs:442)",
        "maps letter -> VCF allele and propagates to the kids:",
        "",
        "   orient      letter->allele map     expected kids      #mis",
        "   --------    ---------------------  -----------------  ----",
    ]
    for i, (dl, ml, exp, nmis) in enumerate(s1["orient_results"]):
        map_str = (
            f"{dl[0]}->{_sorted_unphased(s1['observed']['Dad'])[0]}, "
            f"{dl[1]}->{_sorted_unphased(s1['observed']['Dad'])[1]}, "
            f"{ml[0]}->{_sorted_unphased(s1['observed']['Mom'])[0]}, "
            f"{ml[1]}->{_sorted_unphased(s1['observed']['Mom'])[1]}"
        )
        kids = (
            f"K1={_fmt_gt(exp['Kid1'])} "
            f"K2={_fmt_gt(exp['Kid2'])} "
            f"K3={_fmt_gt(exp['Kid3'])}"
        )
        star = "  <- winner" if nmis == 0 else ""
        body_c.append(
            f"   ({dl[0]},{dl[1]}),({ml[0]},{ml[1]})  {map_str}  {kids}    {nmis}{star}"
        )

    # Phased output for the winning orientation.
    best_exp = s1["best"][2]
    body_c += [
        "",
        "Phased output under the winning orientation:",
        f"  Kid1 : {_fmt_phased(best_exp['Kid1'])}   (paternal=A, maternal=C)",
        f"  Kid2 : {_fmt_phased(best_exp['Kid2'])}   (paternal=B, maternal=D)",
        f"  Kid3 : {_fmt_phased(best_exp['Kid3'])}   (paternal=A, maternal=C)",
    ]
    cap_c = (
        "Panel C. Orientation search at non-informative site N1. Only "
        "one of the four 2^F orientations yields zero mismatches; that "
        "orientation fixes the letter→allele map at this site and the "
        "block letter labels then immediately give the kids' phased "
        "genotypes. Compare this panel to format_genotype_maps "
        "(gtg_concordance.rs:151) — it prints essentially the same table."
    )
    text_panel(axes[1, 0], "C. Four orientations at site N1", body_c, cap_c)

    # --- Panel D: orientation search at the ERROR site ---
    s3 = per_site[2]  # N3
    body_d = [
        f"Site {s3['label']} (ERROR): "
        f"Dad={_fmt_gt(s3['observed']['Dad'])} "
        f"Mom={_fmt_gt(s3['observed']['Mom'])} "
        f"Kid1={_fmt_gt(s3['observed']['Kid1'])} (err)",
        f"                 Kid2={_fmt_gt(s3['observed']['Kid2'])}  "
        f"Kid3={_fmt_gt(s3['observed']['Kid3'])}",
        "",
        "   orient                expected kids       #mis",
        "   --------              ------------------  ----",
    ]
    for dl, ml, exp, nmis in s3["orient_results"]:
        kids = (
            f"K1={_fmt_gt(exp['Kid1'])} "
            f"K2={_fmt_gt(exp['Kid2'])} "
            f"K3={_fmt_gt(exp['Kid3'])}"
        )
        body_d.append(
            f"   ({dl[0]},{dl[1]}),({ml[0]},{ml[1]})         {kids}      {nmis}"
        )
    min_mis = min(r[3] for r in s3["orient_results"])
    body_d += [
        "",
        f"Minimum mismatch across all 4 orientations: {min_mis}",
        "",
        "Since no orientation yields 0 mismatches, the site is written",
        "to {prefix}.fail.vcf by gtg_concordance.rs:507 and the",
        "offending sample(s) are logged to {prefix}.failed_sites.txt.",
        "",
        "If exactly ONE sample is the culprit across the whole block, the",
        "failure is counted as a 'singleton' — a strong signal of a",
        "sequencing error in that one sample at that one site.",
    ]
    cap_d = (
        "Panel D. The 'impossible genotype' rule. At site N3 we flipped "
        "Kid1's observed genotype 0/1 → 1/1. No orientation can explain "
        "this under the block's letter labels, so find_best_phase_orientation "
        "(gtg_concordance.rs:252) returns a non-empty mismatch list and the "
        "site is routed to fail.vcf. This is how gtg-concordance filters "
        "sequencing errors, Mendelian violations, and residual block-"
        "labelling mistakes."
    )
    text_panel(axes[1, 1], "D. Sequencing error → impossible genotype", body_d, cap_d)

    # --- Panel E: truth-vs-deduced phased genotypes for all 3 sites ---
    body_e = [
        "Phased output vs simulated truth for all three non-informative",
        "sites (kids only):",
        "",
        "             TRUTH                       DEDUCED (gtg-concordance)",
    ]
    all_correct = True
    for s in per_site:
        truth = s["truth"]
        if s["best"][3] == 0:
            dedu = s["best"][2]
            status = "PASS"
        else:
            dedu = None
            status = "FAIL (->fail.vcf)"
        truth_str = (
            f"K1={_fmt_phased(truth['Kid1'])} "
            f"K2={_fmt_phased(truth['Kid2'])} "
            f"K3={_fmt_phased(truth['Kid3'])}"
        )
        if dedu is None:
            dedu_str = "(no phasing emitted)"
        else:
            dedu_str = (
                f"K1={_fmt_phased(dedu['Kid1'])} "
                f"K2={_fmt_phased(dedu['Kid2'])} "
                f"K3={_fmt_phased(dedu['Kid3'])}"
            )
            for k in ["Kid1", "Kid2", "Kid3"]:
                if _sorted_unphased(truth[k]) != _sorted_unphased(dedu[k]):
                    all_correct = False
        body_e.append(
            f"  {s['label']} :     {truth_str}    {dedu_str}   [{status}]"
        )
    body_e += [
        "",
        f"All PASS sites phased correctly: {all_correct}",
        "The FAIL site (N3) is NOT phased — it is written to",
        "{prefix}.fail.vcf alongside low-quality and no-call records.",
    ]
    cap_e = (
        "Panel E. Deduced phased genotypes at the PASS sites match the "
        "simulated truth exactly; the error site is quarantined. This "
        "closes the pipeline: gtg-ped-map builds a structural letter map "
        "from informative sites only, and gtg-concordance uses that map "
        "to phase every other site in each block — falling back to a "
        "fail-vcf when no orientation is consistent."
    )
    text_panel(axes[2, 0], "E. Truth vs deduced phased genotypes", body_e, cap_e)

    # --- Panel F: Boundary statement between the two tools ---
    body_f = [
        "Responsibility split:",
        "",
        "  gtg-ped-map (map_builder.rs)",
        "    - uses ONLY informative sites",
        "    - writes LETTERS per individual per block",
        "    - no 0/1 allele sequences leave this tool",
        "",
        "  gtg-concordance (gtg_concordance.rs)",
        "    - visits EVERY VCF record inside each block",
        "    - enumerates 2^F founder orientations",
        "    - picks the one with zero mismatches",
        "    - emits phased genotypes (pass.vcf) or routes to fail.vcf",
        "",
        "This is the answer to 'does gtg-ped-map reconstruct the 0/1",
        "allele sequence of each haplotype?': no, that is handled",
        "exclusively by gtg-concordance via assign_genotypes.",
    ]
    cap_f = (
        "Panel F. Responsibility split between the two binaries. The "
        "separation is clean: gtg-ped-map is a pure structural-labeling "
        "tool operating on informative sites, and gtg-concordance is a "
        "pure phasing/QC tool that uses those structural labels to assign "
        "alleles at every site."
    )
    text_panel(axes[2, 1], "F. Tool boundary", body_f, cap_f)

    fig.suptitle(
        "Component 3 — Phasing non-informative sites with gtg-concordance",
        fontsize=13, fontweight="bold",
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)

    print(f"[component 3] Wrote {out_path}")
    for s in per_site:
        print(
            f"[component 3] site {s['label']}: "
            f"best-orientation mismatches={s['best'][3]} "
            f"error={'yes' if s['error'] else 'no'}"
        )


# ---------------------------------------------------------------------------
# Component 4 — Methods section (markdown)
# ---------------------------------------------------------------------------
#
# Plain-English description of the Rust pipeline plus the real-world noise
# sources it has to cope with. Written to figures/methods.md as a drop-in
# draft for the Bioinformatics manuscript. Source-code line references use
# GitHub permalinks pinned to the current HEAD SHA (captured at script start).

METHODS_TEMPLATE = """\
# Methods — pedigree haplotype mapping with gtg-ped-map and gtg-concordance

The pipeline consists of two Rust binaries that share a common
inheritance-tracking library (`code/rust/src/iht.rs`) and are driven by a
standard PED file and a jointly-called VCF. All line numbers below refer to
commit `{sha_short}`; permalinks are given in square brackets.

## 1. Inputs

The tools consume:

- A **PED file** describing the pedigree
  ([`code/rust/src/ped.rs`]({ped_rs})).
- A **jointly-called, tabix-indexed VCF** — typically DeepVariant on HiFi —
  containing genotypes for every individual in the PED.

Only biallelic SNVs are used for map construction; indels are filtered at
read time (`is_indel` [`map_builder.rs:501`]({map_rs_501})).

## 2. gtg-ped-map: structural haplotype labeling

### 2.1 Founder letter assignment

At startup, the pipeline calls `Iht::new`
([`iht.rs:172`]({iht_rs_172})), which assigns a pair of capital letters
(`A`,`B`), (`C`,`D`), (`E`,`F`), … to the two homologs of each founder. **No
allele sequence is associated with these letters at this stage.** The
letters are pure structural placeholders whose job is to be tracked from
founders to descendants.

### 2.2 Informative-site detection

For every VCF record in a chromosome, `track_alleles_through_pedigree`
([`map_builder.rs:295`]({map_rs_295})) walks individuals in
ancestor-first depth order (`family.get_individual_depths()`). For each
(parent, spouse) pair in the walk, `unique_allele`
([`map_builder.rs:243`]({map_rs_243})) checks whether the parent has an
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
`get_iht_markers` ([`map_builder.rs:274`]({map_rs_274})) simply reads
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
([`map_builder.rs:804`]({map_rs_804})) infers the other founder
allele on the remaining children — a sibling who does not carry the tagged
allele must have inherited the other founder homolog. This reduces the
number of `?` slots before block collapse.

### 2.5 Block collapse and phase flipping

Adjacent sites with compatible letter assignments are merged by
`collapse_identical_iht`
([`map_builder.rs:385`]({map_rs_385})). Because each founder's two
letters ('A','B') could appear in either order in a block, neighbouring
blocks may disagree on which letter is "first". `perform_flips_in_place`
([`map_builder.rs:702`]({map_rs_702})) walks the block list and swaps
founder letter pairs where doing so reduces mismatches with the previous
block, enforcing a consistent phase across the chromosome. Missing values
are then filled by `fill_missing_values` and
`fill_missing_values_by_neighbor` ([`map_builder.rs:617`]({map_rs_617}),
[`:540`]({map_rs_540})) using flanking blocks.

### 2.6 Short-run masking (noise suppression)

Isolated marker assignments that are flanked on both sides by a different
letter are the haplotype-map signature of a single discordant site, which
is almost always a sequencing error rather than a true recombination.
`count_matching_neighbors`
([`map_builder.rs:935`]({map_rs_935})) scans each individual's
per-site letters and identifies runs shorter than `--run` (default 10
markers). `mask_child_alleles`
([`map_builder.rs:970`]({map_rs_970})) then sets the offending slots
back to `?` so they do not survive into the final block map.

### 2.7 Recombination reporting

After block collapse, `summarize_child_changes`
([`map_builder.rs:673`]({map_rs_673})) writes per-individual letter
transitions (e.g. `A → B` on the paternal slot) to
`{{prefix}}.recombinants.txt`. As the HAPLOTYPING.md note makes explicit,
these are candidate transitions: ancestral recombinations in an earlier
generation appear in all descendant children and must be reconciled
downstream.

### 2.8 Outputs

`gtg-ped-map` writes three files per run:

- `{{prefix}}.iht.txt` — the final block map: contigous blocks with a
  letter pair for every individual and a marker count.
- `{{prefix}}.markers.txt` — raw per-site letter assignments before
  collapse, useful for breakpoint refinement.
- `{{prefix}}.recombinants.txt` — putative switch points.

## 3. gtg-concordance: allele phasing and concordance QC

`gtg-concordance` consumes the `.iht.txt` block map and the original VCF.
For each block it fetches every overlapping VCF record — including the
non-informative sites that `gtg-ped-map` discarded — and attempts to phase
each one against the block's letter labels.

### 3.1 Orientation search

For every record, `find_best_phase_orientation`
([`gtg_concordance.rs:252`]({conc_rs_252})) enumerates the 2^F
founder phase orientations produced by `Iht::founder_phase_orientations`
([`iht.rs:492`]({iht_rs_492})). Each orientation swaps the order of
one or more founders' letter pairs, which is equivalent to swapping which
of that founder's two VCF alleles is paired with which letter.

`assign_genotypes` ([`iht.rs:442`]({iht_rs_442})) then produces the
**expected** genotype for every individual under that orientation: for each
founder it maps `letter1 → vcf_allele[0]`, `letter2 → vcf_allele[1]`, and
each child's two letters are looked up in that map to give an expected
allele pair. The expected genotypes are compared to the observed
(unphased) genotypes by `compare_genotype_maps`
([`gtg_concordance.rs:213`]({conc_rs_213})), and the orientation with
the fewest mismatches wins.

### 3.2 Passing, failing, phased output

If the winning orientation yields **zero** mismatches, the site is
consistent with the haplotype map. `gtg-concordance` rewrites the
genotypes as phased (`0|1`, with the paternal allele first) using the
winning letter→allele map and writes the record to `{{prefix}}.pass.vcf`
([`gtg_concordance.rs:509`]({conc_rs_509})).

If **no** orientation yields zero mismatches, the site is "impossible"
under the current block labels and is routed to `{{prefix}}.fail.vcf`
([`gtg_concordance.rs:507`]({conc_rs_507})). The list of offending
samples per site is appended to `{{prefix}}.failed_sites.txt`, and
singleton failures (one sample mismatching in isolation) are tallied per
block in `{{prefix}}.filtering_stats.txt`.

The pipeline's division of labour is deliberate and strict: `gtg-ped-map`
stores **only letters** and **only at informative sites**, while
`gtg-concordance` is the sole place where letter→allele correspondence is
computed and written out. Any downstream tool that needs phased genotypes
must consume `{{prefix}}.pass.vcf`, not the `.iht.txt` structural map.

## 4. Real-world issues

### 4.1 Sequencing errors

Short-read or HiFi sequencing errors flip individual genotypes stochastically.
In the pre-collapse letter trace, a single erroneous site in one kid looks
like an `A → B → A` flicker surrounded by otherwise consistent blocks. The
`count_matching_neighbors` / `mask_child_alleles` pair
([`map_builder.rs:935`]({map_rs_935}),
[`:970`]({map_rs_970})) recognises runs shorter than `--run` (default
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
([`map_builder.rs:217`]({map_rs_217})) drops any site where a sample is
below `--depth` (default 10) or more than one standard deviation from its
per-sample mean — the upper bound catches copy-number-inflated regions.
`--qual` (default 20) filters out low-confidence calls at the VCF-record
level. Both binaries apply these filters independently.

### 4.3 Founder-phase instability across blocks

Because the two letters in a founder's pair are interchangeable within any
single block, the labels on adjacent blocks can disagree by a "flip" even
when the underlying biology is continuous. `perform_flips_in_place`
([`map_builder.rs:702`]({map_rs_702})) resolves this by comparing each
block to its left neighbour and swapping founder letter pairs to minimise
mismatches. This is called twice — once before and once after block
collapse — so that block merging sees consistently-oriented labels.

### 4.4 Missing data in sibships

Multi-child sibships provide redundant information that is essential in
noisy regions: `backfill_sibs`
([`map_builder.rs:804`]({map_rs_804})) uses the fact that if one child
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
`track_alleles_through_pedigree` ([`map_builder.rs:346-373`]({map_rs_346}))
and `Iht::new`
([`iht.rs:181-185`]({iht_rs_181})). Males receive a single non-dot
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
"""


def emit_methods_section(out_path: Path) -> None:
    def link(path: str, line: int) -> str:
        return permalink(path, line, SHA)

    content = METHODS_TEMPLATE.format(
        sha_short=SHA[:7],
        ped_rs=f"https://github.com/{REPO}/blob/{SHA}/code/rust/src/ped.rs",
        map_rs_217=link("code/rust/src/bin/map_builder.rs", 217),
        map_rs_243=link("code/rust/src/bin/map_builder.rs", 243),
        map_rs_274=link("code/rust/src/bin/map_builder.rs", 274),
        map_rs_295=link("code/rust/src/bin/map_builder.rs", 295),
        map_rs_346=link("code/rust/src/bin/map_builder.rs", 346),
        map_rs_385=link("code/rust/src/bin/map_builder.rs", 385),
        map_rs_501=link("code/rust/src/bin/map_builder.rs", 501),
        map_rs_540=link("code/rust/src/bin/map_builder.rs", 540),
        map_rs_617=link("code/rust/src/bin/map_builder.rs", 617),
        map_rs_673=link("code/rust/src/bin/map_builder.rs", 673),
        map_rs_702=link("code/rust/src/bin/map_builder.rs", 702),
        map_rs_804=link("code/rust/src/bin/map_builder.rs", 804),
        map_rs_935=link("code/rust/src/bin/map_builder.rs", 935),
        map_rs_970=link("code/rust/src/bin/map_builder.rs", 970),
        conc_rs_213=link("code/rust/src/bin/gtg_concordance.rs", 213),
        conc_rs_252=link("code/rust/src/bin/gtg_concordance.rs", 252),
        conc_rs_507=link("code/rust/src/bin/gtg_concordance.rs", 507),
        conc_rs_509=link("code/rust/src/bin/gtg_concordance.rs", 509),
        iht_rs_172=link("code/rust/src/iht.rs", 172),
        iht_rs_181=link("code/rust/src/iht.rs", 181),
        iht_rs_442=link("code/rust/src/iht.rs", 442),
        iht_rs_492=link("code/rust/src/iht.rs", 492),
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(content)
    print(f"[component 4] Wrote {out_path}")


# ---------------------------------------------------------------------------
# Wiki index (LLM-wiki style catalog, Karpathy pattern)
# ---------------------------------------------------------------------------
#
# The wiki's top-level `index.md` is a thin catalog: one entry per page with
# a one-line summary, grouped by category, and a recommended reading order.
# Reading order lives here rather than in filenames — subdirectory names are
# descriptive ("nuclear_family") rather than numeric, so adding or reordering
# pages does not rename files.

def emit_wiki_index(out_path: Path) -> None:
    content = """\
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
2. *(coming soon)* **Three-generation pedigree with an outside
   marriage** — shows how the same routine recurses across generations
   via ancestor-first depth ordering, without ever constructing a
   joint inheritance vector over all founders.
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
"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(content)
    print(f"[wiki] Wrote {out_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

PAGE_CHOICES = ["nuclear_family", "three_generations", "concordance", "methods"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--page", choices=PAGE_CHOICES,
        help="Render only one wiki page (default: all).",
    )
    parser.add_argument(
        "--outdir", type=Path, default=FIG_DIR,
        help="Directory to write the wiki into.",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    pages = {
        "nuclear_family": lambda: component_1_nuclear_family(args.outdir),
        "three_generations": lambda: component_2_three_generations(
            args.outdir / "component2_three_generations.png"
        ),
        "concordance": lambda: component_3_concordance(
            args.outdir / "component3_concordance.png"
        ),
        "methods": lambda: emit_methods_section(
            args.outdir / "methods.md"
        ),
    }

    if args.page:
        pages[args.page]()
    else:
        for fn in pages.values():
            fn()

    # The wiki index is the catalog, so it is always refreshed.
    emit_wiki_index(args.outdir / "index.md")


if __name__ == "__main__":
    main()
