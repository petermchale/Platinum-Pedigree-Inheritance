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
from typing import Dict, List, Optional, Tuple

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

NUM_SITES = 9

# Parent haplotypes are labelled with Greek letters (alpha/beta for dad,
# gamma/delta for mom) to make clear that these are *physical homolog
# identities*, NOT the per-site Latin labels (A/B/C/D) that gtg-ped-map
# assigns later. The whole point of the wiki page is that A/B are
# per-site, per-block algorithm labels with no fixed correspondence to
# any specific physical homolog.
#
# Sites 0,1,4,5,8: dad-informative (dad het, mom hom)
# Sites 2,3,6,7:   mom-informative (mom het, dad hom)
# Site 8 is dad-informative AND has Kid1's genotype missing in the VCF;
# this is the case backfill_sibs handles by sibling elimination.
HAP_DAD_ALPHA = "100110010"
HAP_DAD_BETA  = "010101011"
HAP_MOM_GAMMA = "101110010"
HAP_MOM_DELTA = "100010100"

# True transmission pattern (Greek = physical homolog identity):
#   Kid1: (alpha, gamma) at every site            — no recombinations
#   Kid2: (beta,  delta) at every site            — no recombinations
#   Kid3: (alpha, gamma) at sites 0-3,            — paternal recombination
#         (beta,  gamma) at sites 4-8               between sites 3 and 4
KID1_LABELS = [("α", "γ")] * NUM_SITES
KID2_LABELS = [("β", "δ")] * NUM_SITES
KID3_LABELS = [("α", "γ")] * 4 + [("β", "γ")] * 5

# Sites with missing genotypes in the simulated VCF (rendered as ./.).
# At site 8, Kid1's genotype is missing — track_alleles_through_pedigree
# cannot tag Kid1, but backfill_sibs recovers its letter from siblings.
KID_MISSING = {"Kid1": {8}, "Kid2": set(), "Kid3": set()}


def _hap_to_alleles(hap: str) -> List[int]:
    return [int(c) for c in hap]


_GREEK_TO_HAP = {
    "α": HAP_DAD_ALPHA,
    "β": HAP_DAD_BETA,
    "γ": HAP_MOM_GAMMA,
    "δ": HAP_MOM_DELTA,
}


def _child_genotype(
    dad_label: str, mom_label: str, site: int
) -> Tuple[int, int]:
    paternal = _hap_to_alleles(_GREEK_TO_HAP[dad_label])[site]
    maternal = _hap_to_alleles(_GREEK_TO_HAP[mom_label])[site]
    return paternal, maternal


def _unphased_gt(a: int, b: int) -> str:
    if a == b:
        return f"{a}/{a}"
    return "0/1"


def _build_simulation() -> Dict:
    dad_alpha = _hap_to_alleles(HAP_DAD_ALPHA)
    dad_beta = _hap_to_alleles(HAP_DAD_BETA)
    mom_gamma = _hap_to_alleles(HAP_MOM_GAMMA)
    mom_delta = _hap_to_alleles(HAP_MOM_DELTA)

    def gts(labels):
        return [_child_genotype(la, lb, i) for i, (la, lb) in enumerate(labels)]

    kid1_phased = gts(KID1_LABELS)
    kid2_phased = gts(KID2_LABELS)
    kid3_phased = gts(KID3_LABELS)

    def unphased_row(hap1, hap2):
        return [_unphased_gt(hap1[i], hap2[i]) for i in range(NUM_SITES)]

    def kid_unphased(kid_name, phased):
        out = []
        for i, (a, b) in enumerate(phased):
            if i in KID_MISSING[kid_name]:
                out.append("./.")
            else:
                out.append(_unphased_gt(a, b))
        return out

    return {
        "dad_unphased": unphased_row(dad_alpha, dad_beta),
        "mom_unphased": unphased_row(mom_gamma, mom_delta),
        "kid1_unphased": kid_unphased("Kid1", kid1_phased),
        "kid2_unphased": kid_unphased("Kid2", kid2_phased),
        "kid3_unphased": kid_unphased("Kid3", kid3_phased),
        "dad_alpha": dad_alpha,
        "dad_beta": dad_beta,
        "mom_gamma": mom_gamma,
        "mom_delta": mom_delta,
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
    dad labels whichever paternal homolog each kid carries.
    Mirrors unique_allele at map_builder.rs:243."""
    sites = []
    for i in range(NUM_SITES):
        dad_het = sim["dad_alpha"][i] != sim["dad_beta"][i]
        mom_hom = sim["mom_gamma"][i] == sim["mom_delta"][i]
        if dad_het and mom_hom:
            sites.append(i)
    return sites


def _informative_sites_mom(sim: Dict) -> List[int]:
    sites = []
    for i in range(NUM_SITES):
        mom_het = sim["mom_gamma"][i] != sim["mom_delta"][i]
        dad_hom = sim["dad_alpha"][i] == sim["dad_beta"][i]
        if mom_het and dad_hom:
            sites.append(i)
    return sites


def _per_site_parent_labels(
    sim: Dict,
    info_sites: List[int],
    parent_hap_a_key: str,
    parent_hap_b_key: str,
    other_hap_a_key: str,
    other_hap_b_key: str,
    kid_unphased_key: Dict[str, str],
    letter_first: str,
    letter_second: str,
    slot: str,  # "p" or "m"
) -> Dict[str, List[str]]:
    """Faithfully simulate, per informative site, what
    track_alleles_through_pedigree (carriers tagged with letter_first)
    plus backfill_sibs (non-carriers get letter_second; swap-by-majority
    if needed) would write to each kid's slot.

    Returns {kid: [label_or_'?'] * NUM_SITES} where '?' means the site is
    not informative for this slot.
    """
    kids = ["Kid1", "Kid2", "Kid3"]
    out: Dict[str, List[str]] = {k: ["?"] * NUM_SITES for k in kids}

    for s in info_sites:
        parent_alleles = {sim[parent_hap_a_key][s], sim[parent_hap_b_key][s]}
        other_alleles = {sim[other_hap_a_key][s], sim[other_hap_b_key][s]}
        unique_set = parent_alleles - other_alleles
        if not unique_set:
            continue
        unique = next(iter(unique_set))

        # Step 1 (track_alleles_through_pedigree): tag carriers with
        # letter_first. Missing-genotype kids stay untagged.
        per_kid: Dict[str, str] = {}
        for k in kids:
            gt = sim[kid_unphased_key[k]][s]
            if gt == "./.":
                per_kid[k] = "?"
                continue
            kid_alleles = {int(gt[0]), int(gt[2])}
            per_kid[k] = letter_first if unique in kid_alleles else "?"

        # Step 2 (backfill_sibs): if at least one carrier was tagged,
        # assign letter_second to every other (non-missing) sibling.
        any_carrier = any(per_kid[k] == letter_first for k in kids)
        if any_carrier:
            # Backfill: untagged sibs (whether non-carrier or missing)
            # get letter_second, deducing the missing kid by elimination.
            for k in kids:
                if per_kid[k] == "?":
                    per_kid[k] = letter_second

        # Step 3 (backfill swap-by-majority, map_builder.rs:881-905):
        # if letter_first appears fewer times than letter_second across
        # children, swap them so the majority class carries letter_first.
        first_count = sum(1 for k in kids if per_kid[k] == letter_first)
        second_count = sum(1 for k in kids if per_kid[k] == letter_second)
        if first_count < second_count:
            for k in kids:
                if per_kid[k] == letter_first:
                    per_kid[k] = letter_second
                elif per_kid[k] == letter_second:
                    per_kid[k] = letter_first

        for k in kids:
            out[k][s] = per_kid[k]

    return out


def _flip_blocks(
    per_site: Dict[str, List[str]],
    info_sites: List[int],
    letters: Tuple[str, str],
) -> Dict[str, List[str]]:
    """Mimic perform_flips_in_place: walk informative sites left to right
    and, whenever the next site's labels would agree better with the
    running labels under an A<->B swap, flip them. Then propagate each
    site's labels forward and backward to fill non-informative gaps
    (mimicking gap-fill / collapse).
    """
    kids = list(per_site.keys())
    out: Dict[str, List[str]] = {k: ["?"] * NUM_SITES for k in kids}
    a, b = letters
    prev: Optional[Dict[str, str]] = None

    for s in info_sites:
        cur = {k: per_site[k][s] for k in kids}
        if prev is not None:
            flipped = {
                k: (b if cur[k] == a else a if cur[k] == b else cur[k])
                for k in kids
            }
            same = sum(1 for k in kids if cur[k] == prev[k])
            flip = sum(1 for k in kids if flipped[k] == prev[k])
            if flip > same:
                cur = flipped
        for k in kids:
            out[k][s] = cur[k]
        prev = cur

    for k in kids:
        last = "?"
        for i in range(NUM_SITES):
            if out[k][i] != "?":
                last = out[k][i]
            else:
                out[k][i] = last
        last = "?"
        for i in range(NUM_SITES - 1, -1, -1):
            if out[k][i] != "?":
                last = out[k][i]
            else:
                out[k][i] = last

    return out


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
    kid_unphased_key = {
        "Kid1": "kid1_unphased",
        "Kid2": "kid2_unphased",
        "Kid3": "kid3_unphased",
    }

    paternal_deduced = _per_site_parent_labels(
        sim, dad_info,
        parent_hap_a_key="dad_alpha", parent_hap_b_key="dad_beta",
        other_hap_a_key="mom_gamma", other_hap_b_key="mom_delta",
        kid_unphased_key=kid_unphased_key,
        letter_first="A", letter_second="B", slot="p",
    )
    maternal_deduced = _per_site_parent_labels(
        sim, mom_info,
        parent_hap_a_key="mom_gamma", parent_hap_b_key="mom_delta",
        other_hap_a_key="dad_alpha", other_hap_b_key="dad_beta",
        kid_unphased_key=kid_unphased_key,
        letter_first="C", letter_second="D", slot="m",
    )

    paternal_blocks = _flip_blocks(paternal_deduced, dad_info, ("A", "B"))
    maternal_blocks = _flip_blocks(maternal_deduced, mom_info, ("C", "D"))

    nf_dir = out_dir / "nuclear_family"
    nf_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Panel A — ground-truth founder haplotypes (no title, no site row).
    # ------------------------------------------------------------------
    def _kid_pat_row(kid_phased):
        return " ".join(str(a) for a, _ in kid_phased)

    def _kid_mat_row(kid_phased):
        return " ".join(str(b) for _, b in kid_phased)

    body_a = [
        "Figure 1 — Ground-truth founder haplotypes",
        "",
        "Dad  α:  " + " ".join(str(x) for x in sim["dad_alpha"]),
        "Dad  β:  " + " ".join(str(x) for x in sim["dad_beta"]),
        "Mom  γ:  " + " ".join(str(x) for x in sim["mom_gamma"]),
        "Mom  δ:  " + " ".join(str(x) for x in sim["mom_delta"]),
        "",
        "Kid1 p:  " + _kid_pat_row(sim["kid1_phased"]),
        "Kid1 m:  " + _kid_mat_row(sim["kid1_phased"]),
        "Kid2 p:  " + _kid_pat_row(sim["kid2_phased"]),
        "Kid2 m:  " + _kid_mat_row(sim["kid2_phased"]),
        "Kid3 p:  " + _kid_pat_row(sim["kid3_phased"]),
        "Kid3 m:  " + _kid_mat_row(sim["kid3_phased"]),
        "",
        "True transmission (Greek = physical homolog identity):",
        "  Kid1 <- (α, γ)        no recombination",
        "  Kid2 <- (β, δ)        no recombination",
        "  Kid3 <- (α|β, γ)      paternal recomb: α at sites 0-3, β at 4-8",
    ]
    _render_panel_image(body_a, nf_dir / "fig1.png")

    # ------------------------------------------------------------------
    # Panel B — unphased VCF view (no title, no site row).
    # ------------------------------------------------------------------
    def _fmt_gt(g: str) -> str:
        if g == "./.":
            return ".."
        return g[0] + g[2] if g[0] == g[2] else "01"

    body_b = [
        "Figure 2 — Unphased VCF view",
        "",
        "Dad :   " + " ".join(_fmt_gt(g) for g in sim["dad_unphased"]),
        "Mom :   " + " ".join(_fmt_gt(g) for g in sim["mom_unphased"]),
        "Kid1:   " + " ".join(_fmt_gt(g) for g in sim["kid1_unphased"]),
        "Kid2:   " + " ".join(_fmt_gt(g) for g in sim["kid2_unphased"]),
        "Kid3:   " + " ".join(_fmt_gt(g) for g in sim["kid3_unphased"]),
        "",
        "(0/0 rendered as '00', 0/1 as '01', 1/1 as '11', ./. as '..')",
    ]
    _render_panel_image(body_b, nf_dir / "fig2.png")

    # ------------------------------------------------------------------
    # Panel C — merged paternal + maternal informative-site deduction.
    # Rows are explicitly labelled "Kid<n> p" / "Kid<n> m".
    # ------------------------------------------------------------------
    # Align the "* / _" indicator rows with the kid-label columns.
    # The kid-row prefix "  Kid1 p:    " is 13 characters wide, so pad the
    # indicator rows with 13 leading spaces so the marks sit directly above
    # the founder-letter columns of the rows they explain.
    row_prefix_width = len("  Kid1 p:    ")

    def mark_sites(site_list: List[int]) -> str:
        marks = ["_"] * NUM_SITES
        for s in site_list:
            marks[s] = "*"
        return (" " * row_prefix_width) + " ".join(marks)

    body_c = [
        "Figure 3 — Informative-site deduction (paternal and maternal)",
        "",
        "Deduced founder-letter labels at informative sites only",
        "('.' = site not informative for that slot):",
        "",
        mark_sites(dad_info) + "   <- dad-informative (dad het x mom hom)",
        mark_sites(mom_info) + "   <- mom-informative (mom het x dad hom)",
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
    _render_panel_image(body_c, nf_dir / "fig3.png")

    # ------------------------------------------------------------------
    # Panel D — collapsed blocks with paternal / maternal on separate
    # rows, highlighting Kid3's A -> B transition.
    # ------------------------------------------------------------------
    body_d = [
        "Figure 4 — Collapsed blocks with recombination",
        "",
        "Deduced labels after perform_flips_in_place + collapse + gap-fill:",
        "",
    ]
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
        "Kid3's paternal row switches A -> B between sites 3 and 4;",
        "this transition is written to {prefix}.recombinants.txt.",
    ]
    _render_panel_image(body_d, nf_dir / "fig4.png")

    # ------------------------------------------------------------------
    # Panel E — truth vs deduced, paternal and maternal on separate rows.
    # ------------------------------------------------------------------
    def partition_match(truth_at_site: List[str], dedu_at_site: List[str]) -> bool:
        """Two labelings agree as partitions iff for every kid pair, truth
        and deduced agree on whether the two kids share a homolog."""
        n = len(truth_at_site)
        for i in range(n):
            for j in range(i + 1, n):
                t_eq = truth_at_site[i] == truth_at_site[j]
                d_eq = dedu_at_site[i] == dedu_at_site[j]
                if t_eq != d_eq:
                    return False
        return True

    body_e = [
        "Figure 5 — Truth vs deduced founder labels",
        "",
        "Truth uses Greek (physical homologs); deduced uses Latin",
        "(per-block algorithm letters). They match if they induce the",
        "same partition of kids per site.",
        "",
    ]
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

    # Partition mismatches: per site, per slot, across all 3 kids.
    for i in range(NUM_SITES):
        pat_truth_site = [sim["kid_labels"][k][i][0] for k in kids]
        pat_dedu_site = [paternal_blocks[k][i] for k in kids]
        mat_truth_site = [sim["kid_labels"][k][i][1] for k in kids]
        mat_dedu_site = [maternal_blocks[k][i] for k in kids]
        if not partition_match(pat_truth_site, pat_dedu_site):
            mismatches += 1
        if not partition_match(mat_truth_site, mat_dedu_site):
            mismatches += 1
    total_partitions = NUM_SITES * 2  # paternal + maternal slot per site
    body_e.append(
        f"Partition mismatches (paternal+maternal per site): "
        f"{mismatches} / {total_partitions}"
    )
    _render_panel_image(body_e, nf_dir / "fig5.png")

    # ------------------------------------------------------------------
    # Markdown narrative.
    # ------------------------------------------------------------------
    _emit_component1_markdown(
        nf_dir / "nuclear_family.md",
        dad_info=dad_info,
        mom_info=mom_info,
        mismatches=mismatches,
        total_partitions=total_partitions,
    )

    print(f"[component 1] Wrote panel PNGs + markdown to {nf_dir}")
    print(f"[component 1] Dad-informative sites: {dad_info}")
    print(f"[component 1] Mom-informative sites: {mom_info}")
    print(f"[component 1] Partition mismatches vs truth: {mismatches}")
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
    total_partitions: int,
) -> None:
    """Write the Component-1 narrative that interleaves panel PNGs with
    prose explanations and Rust-source permalinks.
    """
    def link(path: str, line: int) -> str:
        return permalink(path, line, SHA)

    map_rs = "code/rust/src/bin/map_builder.rs"
    iht_rs = "code/rust/src/iht.rs"
    concord_rs = "code/rust/src/bin/gtg_concordance.rs"

    content = f"""\
# Structural haplotype mapping in a nuclear family

This page is part of the [wiki](../index.md) and walks through
`gtg-ped-map`'s structural labelling algorithm on the simplest possible
pedigree: a two-generation nuclear family with two founders (dad and
mom) and three children. It complements the full
[`methods.md`](../methods.md) write-up by zooming in on the per-site
mechanics and pinning each figure to the exact Rust code that implements
it. All line numbers refer to commit `{SHA[:7]}`. Each function link is
followed by its call site in the driver — `main()` in
[`map_builder.rs`]({link(map_rs, 989)}) for `gtg-ped-map`, and `main()` in
[`gtg_concordance.rs`]({link(concord_rs, 315)}) for `gtg-concordance` —
so you can step through the driver source in parallel with this
walkthrough.

The toy simulation hard-codes four founder haplotypes over
{NUM_SITES} sites and three children whose transmissions are known a
priori. Everything below is reproducible by running

```
python wiki/generate_wiki.py --page nuclear_family
```

which regenerates both the figure PNGs referenced here and this markdown
file itself.

## 1. Ground truth

![Figure 1 — Ground-truth founder haplotypes](fig1.png)

Each column of Figure 1 corresponds to one biallelic SNV, and each
`0`/`1` entry is the allele carried by that homolog at that site:
`0` is the reference (REF) allele as recorded in the VCF, `1` is the
alternate (ALT) allele. (Indels and multi-allelic sites are filtered
out before this stage; see §2.)

Dad carries two physical homologs, named **α** and **β** here purely as
labels for the figure; mom carries **γ** and **δ**. We use Greek
letters at this stage to emphasise that these names refer to specific
*physical homologs* in the founders' cells. The Latin labels (`A`,
`B`, `C`, `D`) that `gtg-ped-map` eventually writes are something
different: they are *per-site, per-block* algorithm tags. The
phrase has two stages, both relevant downstream:

- **Per-site** refers to the raw per-VCF-record output described in
  §3: every site independently picks which of the parent's two
  letters goes to which group of kids (the grouping is the partition
  defined in §3 by the carrier test described there), so the same
  kid can be tagged `A` at one site and `B` at the next even though
  it inherited the same physical homolog. Figure 3 makes this
  visible in Kid2's paternal row.
- **Per-block** refers to what survives after the across-site
  reconciliation described in §4, which is what `gtg-ped-map`
  actually writes to disk: each contiguous block of sites that share
  the same partition gets one fixed, self-consistent labeling — but
  the block as a whole can still be flipped `A`↔`B` without losing
  any structural information, because the two letters in a founder's
  pair are interchangeable within any single block. That residual
  per-block freedom is what `gtg-concordance` resolves later by
  enumerating `2^n` founder-phase orientations (where `n` is the
  number of founders in the pedigree, i.e. one factor of 2 per
  founder for the independent A↔B / C↔D / … swap) and picking the
  one that best matches the observed alleles.

In neither stage are Latin letters pinned to a specific physical
homolog by `gtg-ped-map` itself, so it is a recurring source of
confusion to read `A` as a fixed name for dad's `α` homolog; it is
not.

In this simulation:

- **Kid1** inherits `(α, γ)` with no recombination.
- **Kid2** inherits `(β, δ)` with no recombination.
- **Kid3** inherits `(α|β, γ)` — the paternal slot crosses over between
  sites 3 and 4, so Kid3 carries dad's `α` homolog on sites 0–3 and
  dad's `β` homolog on sites 4–8.

At program startup, [`Iht::new`]({link(iht_rs, 172)}) (driver calls at
[`map_builder.rs:1059`]({link(map_rs, 1059)}) for the master template
and [`map_builder.rs:1111`]({link(map_rs, 1111)}) for each VCF site)
hands each founder a fresh pair of Latin letters — `(A,B)`, `(C,D)`,
`(E,F)`, … — *without* associating any allele or any physical homolog
with them. The letters are pure structural placeholders. The two
`Iht::new` call sites play different roles: the first builds a
**master template** that is never mutated — only its
[`legend()`]({link(iht_rs, 330)}) is read, to print the column header
(`Dad:A|B Mom:C|D Kid1:?|? …`) at the top of the output files. The
second allocates a fresh `local_iht` per VCF record that
[`track_alleles_through_pedigree`]({link(map_rs, 295)}) then *mutates*
in place to record which founder letter each child inherited at that
site. A per-site copy is needed rather than reusing the master because
(i) each site's IHT vector is itself an output, so it cannot be shared
across sites, and (ii) the master is hard-coded to
`ChromType::Autosome`, whereas `local_iht` uses the chromosome's
actual zygosity (autosome vs. chrX, decided at
[`map_builder.rs:1086`]({link(map_rs, 1086)})), which changes how
letters are laid out for males on chrX.

The goal of `gtg-ped-map` is to recover exactly the Greek-labelled
transmissions above — but expressed in Latin letters and only as
*partitions* of the children, not as physical-homolog identities —
from the jointly-called *unphased* VCF alone (see §2), without ever
looking at the underlying 0/1 allele sequence.

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
[`is_indel`]({link(map_rs, 501)}), invoked from the VCF-reading loop at
[`map_builder.rs:164`]({link(map_rs, 164)}) inside `parse_vcf` (the
driver calls `parse_vcf` at [`map_builder.rs:1092`]({link(map_rs, 1092)})).

## 3. Informative-site detection, letter tagging, and sibling backfill

![Figure 3 — Informative-site deduction (paternal and maternal)](fig3.png)

This section describes what `gtg-ped-map` does at *each VCF record
independently*. The two routines involved —
[`track_alleles_through_pedigree`]({link(map_rs, 295)}) (driver call at
[`map_builder.rs:1116`]({link(map_rs, 1116)})) and
[`backfill_sibs`]({link(map_rs, 804)}) (driver call at
[`map_builder.rs:1122`]({link(map_rs, 1122)})) — are called once per
site, and produce the per-site Latin labels rendered in Figure 3. No
across-site reasoning has happened yet at this stage.

**Step 1 — informative-site detection.**
[`track_alleles_through_pedigree`]({link(map_rs, 295)}) walks the
pedigree in ancestor-first depth order and, for every
`(parent, spouse)` pair, calls
[`unique_allele`]({link(map_rs, 243)}) (from inside the walk at
[`map_builder.rs:315`]({link(map_rs, 315)})) to ask whether the parent
carries an allele that the spouse does not. Two cases can arise:

- **Dad-informative** (dad het × mom hom): dad's unique allele tags
  whichever paternal homolog each child inherited. In this simulation
  these are sites `{dad_info}`.
- **Mom-informative** (mom het × dad hom): symmetric, tagging the
  child's maternal slot. These are sites `{mom_info}`.

**Step 2 — tag carriers with the first letter.** At an informative
site the parent has two distinct alleles, exactly one of which is
*unique* to that parent (i.e. absent from the spouse). The
operational test for each child is then a single allele lookup at
that site: does the child's genotype contain the parent's unique
allele, yes or no? Define a child to be a **carrier** (of the
parent's unique allele, at this site) iff the answer is yes. If the
child is a carrier, it must have inherited the parental homolog that
carries the unique allele (since the spouse could not have donated
it); if not, the child must have inherited the parent's *other*
homolog (the one carrying the allele common to both parents). So the
children are partitioned into two groups by the carrier test. The
two letters of the parent's pair are handed out one per group, but
`track_alleles_through_pedigree` only writes a letter to the carrier
group: it always picks the *first* letter of the parent's pair (`A`
for dad-informative sites, `C` for mom-informative sites,
[`map_builder.rs:333`]({link(map_rs, 333)})) and writes it to every
carrier; the non-carriers are left as `?` and resolved in Step 3.

This per-site choice of "first letter to carriers" is arbitrary in
two senses. First, the parent's `(A, B)` pair was created at startup
with no physical-homolog identity attached. Second, the *carrier
group* itself is defined by whichever physical homolog happens to
carry the unique allele at that particular site, and that can flip
between sites. So the same kid can be tagged `A` at one site and `B`
at the next while the underlying transmission is unchanged — these
are independent draws of an arbitrary label, not real switches. The
IHT therefore records the *partition* (which kids inherited the same
parental homolog) reliably, but identifying `A` with one specific
physical homolog is a per-site, per-block free choice that downstream
code (`perform_flips_in_place`, see §4, and ultimately
`gtg-concordance`'s `2^F`-orientation enumeration) is responsible for
reconciling.

**Step 3 — sibling backfill.**
[`backfill_sibs`]({link(map_rs, 804)}) is then called for the same
site. It exploits the fact that across siblings, both founder
homologs must be represented somewhere (when a parent has more than
one child). Concretely, it does three things per parent per site:

1. **Backfill non-carriers.** For every sibling left as `?` after
   Step 2, write the parent's *other* letter (`B` for dad, `D` for
   mom). This works whether the sibling was a known non-carrier or
   had a *missing genotype* — the latter is recovered purely by
   sibling elimination, since once carriers are identified and only
   one founder homolog remains, that homolog must have gone to the
   untagged child. Site 8 in Figure 3 shows this case: Kid1's
   genotype is `./.` in the VCF, so Step 2 cannot tag it, but Kid2
   and Kid3 are tagged carriers (`A`), so backfill assigns Kid1 the
   other letter (`B`).
2. **Swap by majority** ([`map_builder.rs:881`]({link(map_rs, 881)})).
   After backfilling, count how many sibling slots carry each of the
   two letters. If the letter assigned to carriers (`A` or `C`) ends
   up in the *minority*, swap the two letters across all siblings so
   the majority class always carries the first letter. This is a
   deterministic per-site convention so that two sites whose carrier
   groups happen to be *the same set of kids* — but where one site
   has a 2-carrier majority and the other a 1-carrier minority —
   nevertheless emerge with consistent labels, simplifying later
   block reconciliation.
3. **Skip families with one child** (the `children.len() > 1` guard
   at [`map_builder.rs:818`]({link(map_rs, 818)})), where there is no
   second sibling to provide the elimination signal.

Because the depth-ordered walk in Step 1 always processes a parent
before its children, [`get_iht_markers`]({link(map_rs, 274)}) (called
from inside the walk at [`map_builder.rs:328`]({link(map_rs, 328)}))
reads the parent's already-assigned letters when propagating to the
next generation, which is what makes the method look "recursive"
across generations while being expressed as a single loop.

Non-informative sites (both parents het, or both hom for the same
allele) contribute nothing at this stage and are rendered as `.` in
Figure 3. The two indicator rows (`*` marks informative sites, `_`
marks non-informative ones) sit directly above the kid rows, with
every column aligned, so you can read each letter assignment straight
up to the indicator that produced it. Each kid's paternal row (`p`)
sits directly above its maternal row (`m`). Notice in Figure 3 that
Kid2's paternal row is *not* a uniform run of one letter even though
Kid2 inherits the same physical homolog (`β`) at every dad-informative
site: the per-site `A`/`B` assignment shifts as the carrier group and
its majority shift across sites. That apparent letter switching is
the per-site arbitrariness described in Step 2 — it is not a real
recombination signal, and §4's flip pass is what reconciles it.

## 4. Block collapse, noise filtering, and flip reconciliation

![Figure 4 — Collapsed blocks with recombination](fig4.png)

The per-site labels in Figure 3 are correct as *partitions* of the
children but, as flagged in §3, the letter convention can flip from
site to site. Several Rust routines reconcile and clean the trace up
before it is written to disk:

1. [`collapse_identical_iht`]({link(map_rs, 385)}) (driver call at
   [`map_builder.rs:1191`]({link(map_rs, 1191)})) merges adjacent
   sites with compatible letter assignments into blocks, while
   [`fill_missing_values`]({link(map_rs, 617)}) (driver call at
   [`map_builder.rs:1200`]({link(map_rs, 1200)})) and
   [`fill_missing_values_by_neighbor`]({link(map_rs, 540)}) (driver
   call at [`map_builder.rs:1201`]({link(map_rs, 1201)})) fill the
   `.` gaps visible in Figure 3 from flanking blocks.
2. [`count_matching_neighbors`]({link(map_rs, 935)}) (driver call at
   [`map_builder.rs:1172`]({link(map_rs, 1172)})) and
   [`mask_child_alleles`]({link(map_rs, 970)}) (driver call at
   [`map_builder.rs:1187`]({link(map_rs, 1187)})) identify isolated
   runs shorter than `--run` (default 10 markers) and mask them back
   to `?` as likely sequencing noise, so that collapse does not
   invent spurious recombinations.
3. [`perform_flips_in_place`]({link(map_rs, 702)}) enforces a
   consistent founder-letter orientation across blocks, since the
   two letters in a founder's pair are interchangeable within any
   single block. This is the routine that finally pins each block's
   `A`/`B` (or `C`/`D`) convention so that consecutive blocks agree
   on every kid that did *not* recombine. The driver calls it three
   times — before and after block collapse, and again after gap
   fill — at
   [`map_builder.rs:1135`]({link(map_rs, 1135)}),
   [`map_builder.rs:1193`]({link(map_rs, 1193)}), and
   [`map_builder.rs:1203`]({link(map_rs, 1203)}).

After these steps, each kid's paternal and maternal slots are fully
resolved — shown on separate rows per kid in Figure 4. Within each
block all three kids' labels agree on a single partition; across
adjacent blocks, only the kid(s) that genuinely recombined change
letter. Kid3's highlighted `A`→`B` transition on the paternal row is
emitted to `{{prefix}}.recombinants.txt` by
[`summarize_child_changes`]({link(map_rs, 673)}) (driver call at
[`map_builder.rs:1228`]({link(map_rs, 1228)})).

## 5. Truth versus deduced

![Figure 5 — Truth vs deduced founder labels](fig5.png)

The truth row uses Greek (physical homologs); the deduced row uses
Latin (per-block algorithm letters). They cannot be compared
character-by-character because, as discussed in §3 and §4, the
algorithm guarantees the *partition* of children at each site — not
which physical homolog each Latin letter corresponds to. Two
labelings are therefore counted as agreeing at a site iff they induce
the same partition of the three children: for every kid pair `(i,j)`,
truth and deduced agree on whether kid *i* and kid *j* share a
homolog. By that criterion the deduced trace matches the ground truth
at every site ({mismatches} partition mismatches out of
{total_partitions} partition slots = paternal + maternal per site),
including Kid3's recombination at sites 3/4.

The full output of `gtg-ped-map` for this chromosome is the set of
blocks shown above plus the `recombinants.txt` entry for Kid3's
switch — and critically, **nothing else**. The block map stores only
founder letters; it does *not* store the 0/1 allele sequence of any
haplotype.

Reconstructing which allele each letter represents at every VCF site
is the job of `gtg-concordance`, which will have its own wiki page
once migrated. For every block, `gtg-concordance` enumerates the
`2^F` founder-phase orientations produced by
[`Iht::founder_phase_orientations`]({link(iht_rs, 492)}) (driver call
at [`gtg_concordance.rs:256`]({link(concord_rs, 256)})), maps letters
to VCF alleles via [`assign_genotypes`]({link(iht_rs, 442)}) (driver
call at [`gtg_concordance.rs:267`]({link(concord_rs, 267)})), and
picks the orientation that minimises mismatches against the observed
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

G3_NUM_SITES = 8

# Kid3 arrives at the G2->G3 pass carrying the SAME labels the
# nuclear-family walkthrough deduced for her: paternal slot A at
# sites 0-3 and B at sites 4-7 (the ancestral G1-Dad crossover),
# maternal slot C throughout. In this sub-problem Kid3 plays the
# parent role, and those per-site letters are what gtg-ped-map will
# propagate to the G3 grandchildren.
G3_KID3_PAT_LETTERS = ["A"] * 4 + ["B"] * 4
G3_KID3_MAT_LETTERS = ["C"] * G3_NUM_SITES

# Allele sequences on Kid3's two homologs. These are chosen fresh for
# this page so that Kid3-informative sites (Kid3 het x Spouse hom)
# flank BOTH crossover breakpoints in GK1's maternal-slot trace — the
# ancestral one at sites 3/4 and the new one at sites 5/6. They do
# not need to match the nuclear-family page's allele sequences
# because what propagates between passes is the letter labels, not
# the alleles.
G3_KID3_PAT_ALLELES = "01001010"
G3_KID3_MAT_ALLELES = "10000101"

# Spouse is an outside-marriage founder; Iht::new hands him the next
# fresh letter pair (E,F). Spouse-informative sites are placed at the
# two Kid3-hom sites (2 and 3) so that the grandchildren's paternal
# slot can also be tagged.
G3_SPOUSE_E_ALLELES = "00100000"
G3_SPOUSE_F_ALLELES = "00010000"

# Grandkid transmissions — paternal slot from Spouse, maternal slot
# from Kid3. GK1 illustrates a "recombinant of a recombinant":
#
#   GK1 — paternal = E throughout; maternal = Kid3's paternal homolog
#         at sites 0-5 (letters A,A,A,A,B,B) then Kid3's maternal
#         homolog at sites 6-7 (letters C,C). GK1's maternal-slot
#         letter trace therefore carries TWO transitions: the
#         ancestral A->B at sites 3/4 (inherited from Kid3's own
#         paternal recombinant) and a NEW B->C at sites 5/6
#         (introduced by a crossover in Kid3's G2->G3 meiosis).
#
#   GK2 — paternal = F throughout; maternal = Kid3's paternal
#         homolog throughout (letters A,A,A,A,B,B,B,B). No new
#         crossover, so GK2 carries only the ancestral A->B at
#         sites 3/4.
#
# Both grandchildren show the ancestral A->B at exactly the same
# coordinate, but only GK1 shows the per-meiosis B->C. This is the
# footprint methods.md §4.5 describes for ancestral vs per-meiosis
# recombinations: ancestral events appear in every descendant that
# inherits the affected segment; per-meiosis events appear in a
# single one.
G3_GK1_LABELS = (
    [("E", "A")] * 4 + [("E", "B")] * 2 + [("E", "C")] * 2
)
G3_GK2_LABELS = [("F", "A")] * 4 + [("F", "B")] * 4


def component_2_three_generations(out_dir: Path) -> None:
    """Render the three-generation component as a dedicated wiki page.

    Panels and prose are written to `out_dir/three_generations/`, mirroring
    the layout used for the nuclear-family page.
    """
    num_sites = G3_NUM_SITES
    kid3_pat_letters = G3_KID3_PAT_LETTERS
    kid3_mat_letters = G3_KID3_MAT_LETTERS
    kid3_pat_alleles = _hap_to_alleles(G3_KID3_PAT_ALLELES)
    kid3_mat_alleles = _hap_to_alleles(G3_KID3_MAT_ALLELES)
    spouse_e_alleles = _hap_to_alleles(G3_SPOUSE_E_ALLELES)
    spouse_f_alleles = _hap_to_alleles(G3_SPOUSE_F_ALLELES)

    grandkids = ["GK1", "GK2"]
    truth_labels: Dict[str, List[Tuple[str, str]]] = {
        "GK1": G3_GK1_LABELS,
        "GK2": G3_GK2_LABELS,
    }

    def allele_from_pat_letter(letter: str, site: int) -> int:
        return spouse_e_alleles[site] if letter == "E" else spouse_f_alleles[site]

    def allele_from_mat_letter(letter: str, site: int) -> int:
        if letter == kid3_pat_letters[site]:
            return kid3_pat_alleles[site]
        return kid3_mat_alleles[site]

    def gk_gt(g: str, site: int) -> str:
        pat_letter, mat_letter = truth_labels[g][site]
        return _unphased_gt(
            allele_from_pat_letter(pat_letter, site),
            allele_from_mat_letter(mat_letter, site),
        )

    kid3_unphased = [
        _unphased_gt(kid3_pat_alleles[i], kid3_mat_alleles[i])
        for i in range(num_sites)
    ]
    spouse_unphased = [
        _unphased_gt(spouse_e_alleles[i], spouse_f_alleles[i])
        for i in range(num_sites)
    ]
    gk_unphased: Dict[str, List[str]] = {
        g: [gk_gt(g, i) for i in range(num_sites)] for g in grandkids
    }

    kid3_info_sites: List[int] = []
    spouse_info_sites: List[int] = []
    for i in range(num_sites):
        kid3_het = kid3_pat_alleles[i] != kid3_mat_alleles[i]
        spouse_het = spouse_e_alleles[i] != spouse_f_alleles[i]
        if kid3_het and not spouse_het:
            kid3_info_sites.append(i)
        elif spouse_het and not kid3_het:
            spouse_info_sites.append(i)

    def deduce_maternal_at(site: int) -> Dict[str, str]:
        kid3_alleles = {kid3_pat_alleles[site], kid3_mat_alleles[site]}
        spouse_alleles = {spouse_e_alleles[site], spouse_f_alleles[site]}
        unique = list(kid3_alleles - spouse_alleles)
        if not unique:
            return {g: "?" for g in grandkids}
        u = unique[0]
        if kid3_pat_alleles[site] == u:
            carrier_letter = kid3_pat_letters[site]
            other_letter = kid3_mat_letters[site]
        else:
            carrier_letter = kid3_mat_letters[site]
            other_letter = kid3_pat_letters[site]
        out: Dict[str, str] = {}
        for g in grandkids:
            gt = gk_unphased[g][site]
            gk_alleles = {int(gt[0]), int(gt[2])}
            out[g] = carrier_letter if u in gk_alleles else other_letter
        return out

    def deduce_paternal_at(site: int) -> Dict[str, str]:
        kid3_alleles = {kid3_pat_alleles[site], kid3_mat_alleles[site]}
        spouse_alleles = {spouse_e_alleles[site], spouse_f_alleles[site]}
        unique = list(spouse_alleles - kid3_alleles)
        if not unique:
            return {g: "?" for g in grandkids}
        u = unique[0]
        carrier = "E" if spouse_e_alleles[site] == u else "F"
        other = "F" if carrier == "E" else "E"
        out: Dict[str, str] = {}
        for g in grandkids:
            gt = gk_unphased[g][site]
            gk_alleles = {int(gt[0]), int(gt[2])}
            out[g] = carrier if u in gk_alleles else other
        return out

    paternal_deduced: Dict[str, List[str]] = {g: ["?"] * num_sites for g in grandkids}
    maternal_deduced: Dict[str, List[str]] = {g: ["?"] * num_sites for g in grandkids}
    for s in kid3_info_sites:
        d = deduce_maternal_at(s)
        for g in grandkids:
            maternal_deduced[g][s] = d[g]
    for s in spouse_info_sites:
        d = deduce_paternal_at(s)
        for g in grandkids:
            paternal_deduced[g][s] = d[g]

    def carry_forward(labels: List[str]) -> List[str]:
        out = list(labels)
        last = "?"
        for i in range(G3_NUM_SITES):
            if out[i] != "?":
                last = out[i]
            else:
                out[i] = last
        last = "?"
        for i in range(G3_NUM_SITES - 1, -1, -1):
            if out[i] != "?":
                last = out[i]
            else:
                out[i] = last
        return out

    paternal_blocks = {g: carry_forward(paternal_deduced[g]) for g in grandkids}
    maternal_blocks = {g: carry_forward(maternal_deduced[g]) for g in grandkids}

    mismatches = 0
    for g in grandkids:
        for i in range(G3_NUM_SITES):
            pt, mt = truth_labels[g][i]
            if pt != paternal_blocks[g][i]:
                mismatches += 1
            if mt != maternal_blocks[g][i]:
                mismatches += 1
    total_slots = G3_NUM_SITES * len(grandkids) * 2

    tg_dir = out_dir / "three_generations"
    tg_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Figure 1 — pedigree with outside marriage + ground-truth haplotypes.
    # ------------------------------------------------------------------
    def _letters_row(letters: List[str]) -> str:
        return " ".join(letters)

    body_1 = [
        "Figure 1 — Three-generation pedigree with outside marriage",
        "",
        "        G1-Dad (A,B) --- G1-Mom (C,D)",
        "                      |",
        "           +----------+----------+",
        "           |          |          |",
        "         Kid1        Kid2       Kid3                <-- G2 (resolved first)",
        "        (A,C)       (B,D)     (A|B, C)",
        "                                 |",
        "                           Kid3 --- Spouse (E,F)    <-- outside marriage",
        "                                 |",
        "                           +-----+-----+",
        "                           |           |",
        "                          GK1         GK2           <-- G3 (resolved next)",
        "                       (E, A|B|C)  (F, A|B)",
        "",
        "Kid3 arrives at this pass with her per-site letter labels already",
        "fixed by the G1->G2 walk:",
        "  Kid3 paternal letters:  " + _letters_row(kid3_pat_letters),
        "  Kid3 maternal letters:  " + _letters_row(kid3_mat_letters),
        "",
        "Spouse is a fresh founder — Iht::new hands him letter pair (E,F):",
        "  Spouse E alleles:  " + " ".join(str(x) for x in spouse_e_alleles),
        "  Spouse F alleles:  " + " ".join(str(x) for x in spouse_f_alleles),
        "",
        "True G2->G3 transmissions:",
        "  GK1 <- (E, A|B|C)   ancestral A->B at 3/4 inherited from Kid3,",
        "                      plus a NEW B->C at 5/6 from a Kid3 meiosis crossover.",
        "  GK2 <- (F, A|B)     same ancestral A->B at 3/4, no new crossover.",
    ]
    _render_panel_image(body_1, tg_dir / "fig1.png")

    # ------------------------------------------------------------------
    # Figure 2 — unphased VCF rows for the G2->G3 nuclear unit.
    # ------------------------------------------------------------------
    def _fmt_gt(g: str) -> str:
        return g[0] + g[2] if g[0] == g[2] else "01"

    body_2 = [
        "Figure 2 — Unphased VCF rows for the G2->G3 pass",
        "",
        "Kid3  :   " + " ".join(_fmt_gt(g) for g in kid3_unphased),
        "Spouse:   " + " ".join(_fmt_gt(g) for g in spouse_unphased),
        "GK1   :   " + " ".join(_fmt_gt(g) for g in gk_unphased["GK1"]),
        "GK2   :   " + " ".join(_fmt_gt(g) for g in gk_unphased["GK2"]),
        "",
        "(0/0 rendered as '00', 0/1 as '01', 1/1 as '11')",
    ]
    _render_panel_image(body_2, tg_dir / "fig2.png")

    # ------------------------------------------------------------------
    # Figure 3 — recursive informative-site deduction.
    # Row layout mirrors nuclear_family/fig3: paternal-indicator row first,
    # then maternal-indicator row, then per-grandkid (p, m) rows with the
    # paternal row always directly above its matching maternal row.
    # ------------------------------------------------------------------
    row_prefix_width = len("  GK1 p:    ")

    def mark_sites(site_list: List[int]) -> str:
        marks = ["_"] * num_sites
        for s in site_list:
            marks[s] = "*"
        return (" " * row_prefix_width) + " ".join(marks)

    body_3 = [
        "Figure 3 — Recursive informative-site deduction (G2 -> G3)",
        "",
        "Deduced founder-letter labels at informative sites only",
        "('.' = site not informative for that slot):",
        "",
        mark_sites(spouse_info_sites) + "   <- Spouse-informative (Spouse het x Kid3 hom)",
        mark_sites(kid3_info_sites)   + "   <- Kid3-informative   (Kid3 het x Spouse hom)",
    ]
    for g in grandkids:
        pat_row = f"  {g} p:    " + " ".join(
            paternal_deduced[g][i] if i in spouse_info_sites else "."
            for i in range(num_sites)
        )
        mat_row = f"  {g} m:    " + " ".join(
            maternal_deduced[g][i] if i in kid3_info_sites else "."
            for i in range(num_sites)
        )
        body_3.append(pat_row)
        body_3.append(mat_row)
    _render_panel_image(body_3, tg_dir / "fig3.png")

    # ------------------------------------------------------------------
    # Figure 4 — collapsed blocks showing the ancestral crossover shared
    # between GK1 and GK2, plus the per-meiosis crossover unique to GK1.
    # ------------------------------------------------------------------
    body_4 = [
        "Figure 4 — Ancestral vs per-meiosis crossovers after block collapse",
        "",
        "Deduced labels after block collapse and gap-fill:",
        "",
    ]
    for g in grandkids:
        pat_row = f"  {g} p:    " + " ".join(
            paternal_blocks[g][i] for i in range(num_sites)
        )
        mat_row = f"  {g} m:    " + " ".join(
            maternal_blocks[g][i] for i in range(num_sites)
        )
        body_4.append(pat_row)
        body_4.append(mat_row)
    body_4 += [
        "",
        "Both grandkids' maternal rows switch A -> B between sites 3 and 4;",
        "this is the ancestral crossover from Kid3's own paternal homolog,",
        "inherited intact by every descendant that touches that segment.",
        "",
        "Only GK1's maternal row switches B -> C between sites 5 and 6;",
        "this is the NEW crossover in Kid3's G2->G3 meiosis — a recombinant",
        "of a recombinant — and it appears in one G3 sib only.",
        "",
        "Both transitions are emitted verbatim to {prefix}.recombinants.txt;",
        "collapsing shared ancestral events into a unique meiotic-event count",
        "is the downstream-reconciliation step methods.md section 4.5 describes.",
    ]
    _render_panel_image(body_4, tg_dir / "fig4.png")

    # ------------------------------------------------------------------
    # Figure 5 — truth vs deduced, paternal and maternal on separate rows.
    # ------------------------------------------------------------------
    body_5 = [
        "Figure 5 — Truth vs deduced founder labels",
        "",
        "Truth (T) vs deduced (D) founder-letter labels:",
        "",
    ]
    for g in grandkids:
        pat_truth = [truth_labels[g][i][0] for i in range(num_sites)]
        mat_truth = [truth_labels[g][i][1] for i in range(num_sites)]
        pat_dedu = [paternal_blocks[g][i] for i in range(num_sites)]
        mat_dedu = [maternal_blocks[g][i] for i in range(num_sites)]
        body_5.append(f"  {g} p  T:  " + " ".join(pat_truth))
        body_5.append(f"  {g} p  D:  " + " ".join(pat_dedu))
        body_5.append(f"  {g} m  T:  " + " ".join(mat_truth))
        body_5.append(f"  {g} m  D:  " + " ".join(mat_dedu))
        body_5.append("")
    body_5.append(f"Total label mismatches: {mismatches} / {total_slots}")
    _render_panel_image(body_5, tg_dir / "fig5.png")

    # ------------------------------------------------------------------
    # Markdown narrative.
    # ------------------------------------------------------------------
    _emit_component2_markdown(
        tg_dir / "three_generations.md",
        kid3_info_sites=kid3_info_sites,
        spouse_info_sites=spouse_info_sites,
        g3_num_sites=num_sites,
        mismatches=mismatches,
        total_slots=total_slots,
    )

    print(f"[component 2] Wrote panel PNGs + markdown to {tg_dir}")
    print(f"[component 2] Kid3-informative sites: {kid3_info_sites}")
    print(f"[component 2] Spouse-informative sites: {spouse_info_sites}")
    print(f"[component 2] Grandkid label mismatches vs truth: {mismatches}")


def _emit_component2_markdown(
    out_path: Path,
    *,
    kid3_info_sites: List[int],
    spouse_info_sites: List[int],
    g3_num_sites: int,
    mismatches: int,
    total_slots: int,
) -> None:
    """Write the Component-2 narrative that interleaves panel PNGs with
    prose explanations and Rust-source permalinks.
    """
    def link(path: str, line: int) -> str:
        return permalink(path, line, SHA)

    map_rs = "code/rust/src/bin/map_builder.rs"
    iht_rs = "code/rust/src/iht.rs"
    ped_rs = "code/rust/src/ped.rs"

    content = f"""\
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

All line numbers refer to commit `{SHA[:7]}`. As in the other
walkthrough pages, each function link is followed by its call site in
the driver — `main()` in
[`map_builder.rs`]({link(map_rs, 989)}) — so you can step through the
driver source in parallel with this walkthrough.

The toy simulation adds three things on top of the nuclear-family
example:

- **Kid3** (already carrying a G1-ancestral A→B transition on her
  paternal homolog and C throughout on her maternal homolog, from the
  nuclear-family pass) marries **Spouse**, a fresh founder whose two
  homologs are labelled **E** and **F** by
  [`Iht::new`]({link(iht_rs, 172)}) (driver calls at
  [`map_builder.rs:1059`]({link(map_rs, 1059)}) for the master template
  and [`map_builder.rs:1111`]({link(map_rs, 1111)}) for each VCF site).
- **GK1** inherits Spouse's E homolog plus a Kid3 gamete produced by a
  crossover: Kid3's paternal homolog (letters A,A,A,A,B,B) for sites
  0–5, then Kid3's maternal homolog (letter C,C) for sites 6–7.
- **GK2** inherits Spouse's F homolog plus Kid3's paternal homolog
  unrecombined (letters A,A,A,A,B,B,B,B).

The whole simulation spans {g3_num_sites} VCF sites. Everything below
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
branch, so [`Iht::new`]({link(iht_rs, 172)}) (called from the driver
at [`map_builder.rs:1059`]({link(map_rs, 1059)}) and
[`map_builder.rs:1111`]({link(map_rs, 1111)})) hands him the next
fresh letter pair `(E, F)`. Nothing about Spouse depends on the G1 pass.

## 2. Unphased VCF rows for the G2→G3 pass

![Figure 2 — Unphased VCF rows for the G2->G3 pass](fig2.png)

These are the only genotype rows the tool sees for the new nuclear
unit. Kid3's row is read with Kid3 in the **parent** role. That is the
key structural point: `gtg-ped-map` does not treat G1 and G2
individuals differently; it just iterates over (parent, spouse, child)
triples in ancestor-first depth order given by
`family.get_individual_depths()` (see
[`ped.rs`](https://github.com/{REPO}/blob/{SHA}/{ped_rs})).

## 3. Recursive informative-site deduction

![Figure 3 — Recursive informative-site deduction (G2 -> G3)](fig3.png)

The exact same function,
[`track_alleles_through_pedigree`]({link(map_rs, 295)}) (driver call at
[`map_builder.rs:1116`]({link(map_rs, 1116)})), that handled G1→G2 now
handles G2→G3. For each (parent, spouse) pair it calls
[`unique_allele`]({link(map_rs, 243)}) (invoked from inside the walk
at [`map_builder.rs:315`]({link(map_rs, 315)})) to find alleles that
one partner carries and the other does not:

- **Spouse-informative** (Spouse het × Kid3 hom): the unique paternal
  allele tags whichever Spouse homolog (`E` or `F`) each grandchild
  inherited. In this simulation these are sites `{spouse_info_sites}`.
- **Kid3-informative** (Kid3 het × Spouse hom): symmetric, but here
  the letter tagged on each grandchild's maternal slot is whichever of
  Kid3's *per-site* letters — `A`, `B`, or `C` — sits on the homolog
  carrying Kid3's unique allele at that site. These are sites
  `{kid3_info_sites}`.

When the walk reaches the Kid3–Spouse pair,
[`get_iht_markers`]({link(map_rs, 274)}) (called from inside the walk
at [`map_builder.rs:328`]({link(map_rs, 328)})) reads Kid3's
already-assigned letters directly — those labels were written during
the earlier G1→G2 iteration of the same loop. That is what makes the
algorithm look recursive across generations even though it is a single
ancestor-first pass: by the time the loop reaches a G2 parent, her
letter labels are already finalized and they serve as the "founder
labels" for the G2→G3 sub-problem, even when those labels vary from
site to site. No joint inheritance vector over all founders
`{{A, B, C, D, E, F}}` is ever constructed; each grandkid ends with
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
[`collapse_identical_iht`]({link(map_rs, 385)}) (driver call at
[`map_builder.rs:1191`]({link(map_rs, 1191)})),
[`fill_missing_values`]({link(map_rs, 617)}) (driver call at
[`map_builder.rs:1200`]({link(map_rs, 1200)})),
[`fill_missing_values_by_neighbor`]({link(map_rs, 540)}) (driver call
at [`map_builder.rs:1201`]({link(map_rs, 1201)})), and
[`perform_flips_in_place`]({link(map_rs, 702)}) (driver calls at
[`map_builder.rs:1135`]({link(map_rs, 1135)}),
[`map_builder.rs:1193`]({link(map_rs, 1193)}), and
[`map_builder.rs:1203`]({link(map_rs, 1203)})) — run on the G3 trace
without modification. Two transitions fall out, and the contrast
between them is the point of this page:

- **Ancestral A→B at sites 3/4** appears on *both* grandchildren's
  maternal rows. This is the G1-Dad crossover Kid3 already carried on
  her paternal homolog — it propagates to every descendant that
  inherits the affected segment, which is both GK1 and GK2 here. A
  single meiotic event (in G1) produces two rows in
  `{{prefix}}.recombinants.txt`.
- **Per-meiosis B→C at sites 5/6** appears on GK1's maternal row only.
  This is the new crossover in Kid3's G2→G3 meiosis, and because it
  happened in the gamete that made GK1, it is absent from GK2's trace.
  One meiotic event (in G2), one row.

Both transitions are emitted by
[`summarize_child_changes`]({link(map_rs, 673)}) (driver call at
[`map_builder.rs:1228`]({link(map_rs, 1228)})). The raw output
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
match the ground truth at every site ({mismatches} mismatches out of
{total_slots} label slots), including both the ancestral A→B and the
per-meiosis B→C in GK1. The block map stored for G3 uses the letters
`{{A, B, C, E, F}}` — all three of Kid3's per-site letters reach G3,
because GK1's gamete from Kid3 crossed over between homologs. As in
the nuclear-family case, the block map contains only founder letters;
the 0/1 allele sequence of each haplotype is reconstructed downstream
by `gtg-concordance`, covered in the
[concordance walkthrough](../concordance/concordance.md).
"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(content)


# ---------------------------------------------------------------------------
# Component 3 — gtg-concordance phasing at non-informative sites
# ---------------------------------------------------------------------------
#
# Rust cross-references (all permalinks assembled in `_emit_component3_markdown`):
#   main()                             gtg_concordance.rs:315
#   parse_ihtv2_file                   iht.rs:606           (driver: gtg_concordance.rs:405)
#   find_best_phase_orientation        gtg_concordance.rs:252 (driver: gtg_concordance.rs:454)
#   Iht::founder_phase_orientations    iht.rs:492           (driver: gtg_concordance.rs:256)
#   Iht::assign_genotypes              iht.rs:442           (driver: gtg_concordance.rs:487 / :514)
#   compare_genotype_maps              gtg_concordance.rs:213 (driver: gtg_concordance.rs:268)
#   {prefix}.fail.vcf writes            gtg_concordance.rs:444 / :450 / :507
#   {prefix}.pass.vcf write             gtg_concordance.rs:534
#
# Setup: reuse the nuclear family from Component 1 and focus on the LEFT
# block of the chromosome where kid labels are Kid1=(A,C), Kid2=(B,D),
# Kid3=(A,C). Add two EXTRA sites inside that block where both parents
# are heterozygous (0/1 x 0/1). Such sites are NON-informative for
# gtg-ped-map (neither parent has a unique allele), so they contribute
# nothing to building the block, but gtg-concordance still has to phase
# them using the block's letter labels.

# Non-informative site definitions (dad and mom both 0/1 at every site).
# Each entry gives the true founder haplotype alleles (A,B,C,D) at that site.
# Kid3 is taken from the LEFT block (Kid3 = (A,C) before the recombination).
NON_INFORMATIVE_SITES = [
    # (label, A_allele, B_allele, C_allele, D_allele)
    ("N1", 0, 1, 0, 1),   # clean pass: exactly one orientation consistent
    ("N2", 0, 1, 1, 0),   # with injected error: no orientation consistent
]

# Block letter labels (from Component 1, left-half of the chromosome).
BLOCK_LABELS = {
    "Kid1": ("A", "C"),
    "Kid2": ("B", "D"),
    "Kid3": ("A", "C"),
}

# Inject a sequencing error at site N2: Kid1's observed genotype is 1/1
# although the simulation truth is 0/1.
CONCORDANCE_ERROR_INJECTION: Dict[str, Dict[str, Tuple[int, int]]] = {
    "N2": {"Kid1": (1, 1)},
}


def _site_genotypes(a: int, b: int, c: int, d: int) -> Dict[str, Tuple[int, int]]:
    """Return TRUE (paternal, maternal) allele pairs at one site."""
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
    """Mirror of Iht::assign_genotypes (iht.rs:442).

    Under the orientation (dad_letters, mom_letters), dad's sorted VCF
    allele pair maps onto (dad_letters[0], dad_letters[1]) and similarly
    for mom. Kids' letters are then looked up in that letter->allele map.
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


# All 2^2 = 4 orientations of the two founder letter pairs.
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


def component_3_concordance(out_dir: Path) -> None:
    """Render the gtg-concordance component as a dedicated wiki page.

    Panels and prose are written to `out_dir/concordance/`, mirroring the
    layout used for the nuclear-family and three-generation pages.
    """
    # -------------- Per-site orientation analysis --------------
    per_site = []
    for label, a, b, c, d in NON_INFORMATIVE_SITES:
        truth = _site_genotypes(a, b, c, d)

        observed = {k: _sorted_unphased(v) for k, v in truth.items()}
        if label in CONCORDANCE_ERROR_INJECTION:
            for k, v in CONCORDANCE_ERROR_INJECTION[label].items():
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
            "error": label in CONCORDANCE_ERROR_INJECTION,
        })

    c_dir = out_dir / "concordance"
    c_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Figure 1 — block letter map (the gtg-ped-map output that
    # gtg-concordance reads in as its input).
    # ------------------------------------------------------------------
    body_1 = [
        "Figure 1 — Block letter map (input to gtg-concordance)",
        "",
        "Letter labels in the left half of the chromosome,",
        "as written by gtg-ped-map into {prefix}.iht.txt:",
        "",
        "Dad  p:  A    <- founder, paternal homolog",
        "Dad  m:  B    <- founder, maternal homolog",
        "Mom  p:  C    <- founder, paternal homolog",
        "Mom  m:  D    <- founder, maternal homolog",
        "Kid1 p:  A",
        "Kid1 m:  C",
        "Kid2 p:  B",
        "Kid2 m:  D",
        "Kid3 p:  A    <- left half, before the paternal recomb",
        "Kid3 m:  C",
        "",
        "These letters are constant inside one IHT block.",
        "gtg-concordance re-reads the VCF for the same region and",
        "attempts to phase EVERY variant in the block, not only the",
        "informative ones that built it.",
    ]
    _render_panel_image(body_1, c_dir / "fig1.png")

    # ------------------------------------------------------------------
    # Figure 2 — the two non-informative sites inside the block.
    # ------------------------------------------------------------------
    site_cols = ["Site", "Dad", "Mom", "Kid1", "Kid2", "Kid3"]
    site_widths = [8, 6, 6, 6, 6, 6]

    def _cols(items: List[str], widths: List[int], trailing: str = "") -> str:
        parts = [f"{item:<{w}}" for item, w in zip(items, widths)]
        return "".join(parts).rstrip() + trailing

    body_2 = [
        "Figure 2 — Non-informative sites inside the block",
        "",
        "Two additional sites INSIDE the block where both parents are het:",
        "",
        _cols(site_cols, site_widths, "   <- observed genotypes"),
    ]
    for s in per_site:
        obs = s["observed"]
        tag = "   (ERROR)" if s["error"] else ""
        body_2.append(
            _cols(
                [
                    s["label"],
                    _fmt_gt(obs["Dad"]),
                    _fmt_gt(obs["Mom"]),
                    _fmt_gt(obs["Kid1"]),
                    _fmt_gt(obs["Kid2"]),
                    _fmt_gt(obs["Kid3"]),
                ],
                site_widths,
                tag,
            )
        )
    body_2 += [
        "",
        "At every site here, dad and mom are both 0/1:",
        "unique_allele(dad, mom) = {} and unique_allele(mom, dad) = {}",
        "=> NON-informative for gtg-ped-map (nothing added to the block)",
        "=> still phaseable by gtg-concordance via the block letter map",
        "",
        "Site N2 carries an INJECTED sequencing error: Kid1 observed as",
        "1/1 when the simulation truth is 0/1.",
    ]
    _render_panel_image(body_2, c_dir / "fig2.png")

    # ------------------------------------------------------------------
    # Figure 3 — orientation enumeration at the clean site N1.
    # ------------------------------------------------------------------
    s1 = per_site[0]
    obs1 = s1["observed"]

    orient_cols_full = ["orient", "letter->allele map", "expected kids", "#mis"]
    orient_widths_full = [13, 24, 24, 5]

    body_3 = [
        f"Figure 3 — Four orientations at site {s1['label']} (clean pass)",
        "",
        f"Site {s1['label']}:  "
        f"Dad={_fmt_gt(obs1['Dad'])}  "
        f"Mom={_fmt_gt(obs1['Mom'])}  "
        f"Kid1={_fmt_gt(obs1['Kid1'])}  "
        f"Kid2={_fmt_gt(obs1['Kid2'])}  "
        f"Kid3={_fmt_gt(obs1['Kid3'])}",
        "",
        "For each 2^F=4 orientation, assign_genotypes maps letter -> VCF",
        "allele and propagates to the kids:",
        "",
        _cols(orient_cols_full, orient_widths_full),
        _cols(["-" * (w - 1) for w in orient_widths_full], orient_widths_full),
    ]
    for dl, ml, exp, nmis in s1["orient_results"]:
        dad_pair = _sorted_unphased(obs1["Dad"])
        mom_pair = _sorted_unphased(obs1["Mom"])
        orient = f"({dl[0]},{dl[1]}),({ml[0]},{ml[1]})"
        map_str = (
            f"{dl[0]}->{dad_pair[0]},{dl[1]}->{dad_pair[1]},"
            f"{ml[0]}->{mom_pair[0]},{ml[1]}->{mom_pair[1]}"
        )
        kids = (
            f"K1={_fmt_gt(exp['Kid1'])} "
            f"K2={_fmt_gt(exp['Kid2'])} "
            f"K3={_fmt_gt(exp['Kid3'])}"
        )
        star = "   <- winner" if nmis == 0 else ""
        body_3.append(
            _cols([orient, map_str, kids, str(nmis)], orient_widths_full, star)
        )

    best_exp = s1["best"][2]
    body_3 += [
        "",
        "Phased output under the winning orientation:",
        f"Kid1 p|m:  {_fmt_phased(best_exp['Kid1'])}   (paternal=A, maternal=C)",
        f"Kid2 p|m:  {_fmt_phased(best_exp['Kid2'])}   (paternal=B, maternal=D)",
        f"Kid3 p|m:  {_fmt_phased(best_exp['Kid3'])}   (paternal=A, maternal=C)",
    ]
    _render_panel_image(body_3, c_dir / "fig3.png")

    # ------------------------------------------------------------------
    # Figure 4 — orientation search at the error site N2.
    # ------------------------------------------------------------------
    s2 = per_site[1]
    obs2 = s2["observed"]

    orient_cols_slim = ["orient", "expected kids", "#mis"]
    orient_widths_slim = [13, 24, 5]

    body_4 = [
        f"Figure 4 — Four orientations at site {s2['label']} (injected error)",
        "",
        f"Site {s2['label']} (ERROR):  "
        f"Dad={_fmt_gt(obs2['Dad'])}  "
        f"Mom={_fmt_gt(obs2['Mom'])}  "
        f"Kid1={_fmt_gt(obs2['Kid1'])} (err)  "
        f"Kid2={_fmt_gt(obs2['Kid2'])}  "
        f"Kid3={_fmt_gt(obs2['Kid3'])}",
        "",
        _cols(orient_cols_slim, orient_widths_slim),
        _cols(["-" * (w - 1) for w in orient_widths_slim], orient_widths_slim),
    ]
    for dl, ml, exp, nmis in s2["orient_results"]:
        orient = f"({dl[0]},{dl[1]}),({ml[0]},{ml[1]})"
        kids = (
            f"K1={_fmt_gt(exp['Kid1'])} "
            f"K2={_fmt_gt(exp['Kid2'])} "
            f"K3={_fmt_gt(exp['Kid3'])}"
        )
        body_4.append(_cols([orient, kids, str(nmis)], orient_widths_slim))
    min_mis = min(r[3] for r in s2["orient_results"])
    body_4 += [
        "",
        f"Minimum mismatch across all 4 orientations: {min_mis}",
        "",
        "No orientation yields 0 mismatches, so the site is written to",
        "{prefix}.fail.vcf and the offending sample(s) are logged to",
        "{prefix}.failed_sites.txt. If exactly ONE sample is the culprit",
        "across the whole block, the failure is counted as a 'singleton',",
        "a strong signal of a sequencing error in that one sample.",
    ]
    _render_panel_image(body_4, c_dir / "fig4.png")

    # ------------------------------------------------------------------
    # Figure 5 — truth vs deduced phased genotypes, paternal and
    # maternal on separate rows per kid.
    # ------------------------------------------------------------------
    body_5 = [
        "Figure 5 — Truth vs deduced phased genotypes",
        "",
        "Truth (T) vs deduced (D) phased genotypes at the two sites:",
        "",
    ]
    all_pass = True
    for s in per_site:
        truth = s["truth"]
        if s["best"][3] == 0:
            dedu = s["best"][2]
            status = "PASS -> pass.vcf"
        else:
            dedu = None
            status = "FAIL -> fail.vcf"
            all_pass = False
        body_5.append(f"Site {s['label']}    [{status}]")
        for kid in ["Kid1", "Kid2", "Kid3"]:
            tp, tm = truth[kid]
            if dedu is None:
                dp_str = "."
                dm_str = "."
            else:
                dp_str = str(dedu[kid][0])
                dm_str = str(dedu[kid][1])
            body_5.append(f"{kid} p   T: {tp}   D: {dp_str}")
            body_5.append(f"{kid} m   T: {tm}   D: {dm_str}")
        body_5.append("")

    total_slots = len(NON_INFORMATIVE_SITES) * 3 * 2
    pass_slots = sum(
        3 * 2 for s in per_site if s["best"][3] == 0
    )
    mismatches = 0
    for s in per_site:
        if s["best"][3] != 0:
            continue
        dedu = s["best"][2]
        for kid in ["Kid1", "Kid2", "Kid3"]:
            tp, tm = s["truth"][kid]
            dp, dm = dedu[kid]
            if tp != dp:
                mismatches += 1
            if tm != dm:
                mismatches += 1
    body_5.append(
        f"At PASS sites: {mismatches} phased-allele mismatches out of "
        f"{pass_slots} slots."
    )
    body_5.append(
        "FAIL sites are NOT phased — they are written as-is to fail.vcf."
    )
    _render_panel_image(body_5, c_dir / "fig5.png")

    # ------------------------------------------------------------------
    # Markdown narrative.
    # ------------------------------------------------------------------
    _emit_component3_markdown(
        c_dir / "concordance.md",
        per_site=per_site,
        mismatches=mismatches,
        pass_slots=pass_slots,
        total_slots=total_slots,
    )

    print(f"[component 3] Wrote panel PNGs + markdown to {c_dir}")
    for s in per_site:
        print(
            f"[component 3] site {s['label']}: "
            f"best-orientation mismatches={s['best'][3]} "
            f"error={'yes' if s['error'] else 'no'}"
        )


def _emit_component3_markdown(
    out_path: Path,
    *,
    per_site: List[Dict],
    mismatches: int,
    pass_slots: int,
    total_slots: int,
) -> None:
    """Write the Component-3 narrative that interleaves panel PNGs with
    prose explanations and Rust-source permalinks.
    """
    def link(path: str, line: int) -> str:
        return permalink(path, line, SHA)

    conc_rs = "code/rust/src/bin/gtg_concordance.rs"
    iht_rs = "code/rust/src/iht.rs"
    map_rs = "code/rust/src/bin/map_builder.rs"

    n1 = per_site[0]
    n2 = per_site[1]
    n2_min_mis = min(r[3] for r in n2["orient_results"])

    content = f"""\
# Closing the loop with `gtg-concordance`

This page is part of the [wiki](../index.md) and picks up where the
[nuclear-family walkthrough](../nuclear_family/nuclear_family.md) left
off. `gtg-ped-map` emits only founder letters, and only at informative
sites; it never reconstructs the 0/1 allele sequence of any haplotype.
That job belongs to `gtg-concordance`, which re-reads the VCF for each
IHT block and phases **every** variant using the block's letter map.
All line numbers refer to commit `{SHA[:7]}`. As in the other
walkthrough pages, each function link is followed by its call site in
the driver — `main()` in
[`gtg_concordance.rs`]({link(conc_rs, 315)}) — so you can step through
the driver source in parallel with this walkthrough.

The toy simulation reuses the left half of the nuclear-family block
(Kid1=(A,C), Kid2=(B,D), Kid3=(A,C)) and adds two sites where both
parents are heterozygous. These are NON-informative for `gtg-ped-map`
because [`unique_allele`]({link(map_rs, 243)}) returns `None` at each
of them, but `gtg-concordance` still has to phase them. One of the two
sites carries an injected sequencing error so both the clean-pass
(`pass.vcf`) and error-quarantine (`fail.vcf`) code paths are
exercised. Everything below is reproducible by running

```
python wiki/generate_wiki.py --page concordance
```

which regenerates both the figure PNGs referenced here and this
markdown file itself.

## 1. Block letter map — the input to `gtg-concordance`

![Figure 1 — Block letter map (input to gtg-concordance)](fig1.png)

The driver reads the `{{prefix}}.iht.txt` file produced by
`gtg-ped-map` using
[`parse_ihtv2_file`]({link(iht_rs, 606)}) (driver call at
[`gtg_concordance.rs:405`]({link(conc_rs, 405)})). Each entry carries
one block's letter labels for every individual in the pedigree. Inside
this block the letters are constant; for the left half of the
nuclear-family chromosome they are the ones shown in Figure 1. Every
VCF record whose position falls inside the block is then handed to
the per-site phasing loop at
[`gtg_concordance.rs:437`]({link(conc_rs, 437)}) — including variants
that `gtg-ped-map` could not use because neither parent has a unique
allele.

## 2. Non-informative sites inside the block

![Figure 2 — Non-informative sites inside the block](fig2.png)

The two sites in Figure 2 are both homozygous-absent for informative
patterns: dad is `0/1` and so is mom. `gtg-ped-map`'s
[`unique_allele`]({link(map_rs, 243)}) test therefore returns `None`
at both sites, so neither contributes to block construction. They
still enter `gtg-concordance`'s phasing loop; the per-site machinery
described below is what turns their unphased genotypes into either a
phased `pass.vcf` record or a quarantined `fail.vcf` record.

Site `N2` carries an **injected sequencing error**: Kid1 is reported
as `1/1` even though the simulation truth is `0/1`. This is the case
that the "impossible genotype" rule in Figure 4 is designed to catch.

## 3. Orientation enumeration at a clean site

![Figure 3 — Four orientations at site N1 (clean pass)](fig3.png)

At every record inside a block,
[`find_best_phase_orientation`]({link(conc_rs, 252)}) (driver call at
[`gtg_concordance.rs:454`]({link(conc_rs, 454)})) enumerates the
`2^F=4` orientations produced by
[`Iht::founder_phase_orientations`]({link(iht_rs, 492)}) (invoked
inside `find_best_phase_orientation` at
[`gtg_concordance.rs:256`]({link(conc_rs, 256)})). Each orientation is
a choice of which of dad's two sorted VCF alleles is tagged `A` vs `B`
and which of mom's two is tagged `C` vs `D`. Under a given
orientation,
[`Iht::assign_genotypes`]({link(iht_rs, 442)}) (driver call at
[`gtg_concordance.rs:487`]({link(conc_rs, 487)}) on the failing branch
and [`gtg_concordance.rs:514`]({link(conc_rs, 514)}) on the passing
branch) turns each kid's letter pair into an expected genotype. A
straight equality check —
[`compare_genotype_maps`]({link(conc_rs, 213)}) (driver call at
[`gtg_concordance.rs:268`]({link(conc_rs, 268)})) — counts how many
samples disagree with the observation.

At site `N1`, exactly one of the four orientations explains every
sample simultaneously; the three others each force a mismatch
somewhere among the kids. The winning orientation fixes the
letter→allele map at this site, and the block's letter labels
immediately give the phased `p|m` genotypes shown beneath the
orientation table.

## 4. The "impossible genotype" rule

![Figure 4 — Four orientations at site N2 (injected error)](fig4.png)

At site `N2` the injected error means **no** orientation produces zero
mismatches — the best any orientation can do is {n2_min_mis} sample(s)
disagreeing. `find_best_phase_orientation` therefore returns a
non-empty mismatch list, the driver writes the record to
`{{prefix}}.fail.vcf` at
[`gtg_concordance.rs:507`]({link(conc_rs, 507)}) (alongside the
low-quality and no-call records that were already routed to fail at
[`gtg_concordance.rs:444`]({link(conc_rs, 444)}) and
[`gtg_concordance.rs:450`]({link(conc_rs, 450)})), and the offending
sample names are appended to `{{prefix}}.failed_sites.txt`. If
exactly one sample is the culprit across the whole block, the failure
is counted as a "singleton" — a strong signal of a sequencing error
in that one sample rather than a systemic block-labelling problem.

This is the mechanism by which `gtg-concordance` filters sequencing
errors, Mendelian violations, and residual block-labelling mistakes
without ever producing a phased call that the block's structural
labels cannot justify.

## 5. Truth versus deduced phased genotypes

![Figure 5 — Truth vs deduced phased genotypes](fig5.png)

At the PASS site (`N1`), every kid's deduced paternal and maternal
phased alleles match the simulation truth exactly
({mismatches} mismatches out of {pass_slots} phased-allele slots). At
the FAIL site (`N2`), no phased output is emitted — the record lands
in `fail.vcf` untouched and the truth row in Figure 5 is shown only to
document what `gtg-concordance` declined to commit to.

This closes the pipeline. `gtg-ped-map` (`map_builder.rs`) is a pure
structural-labelling tool that operates on informative sites only and
writes founder letters per individual per block.
`gtg-concordance` (`gtg_concordance.rs`) is a pure phasing/QC tool
that uses those structural labels to assign alleles at every site in
the block: the passing records are phased via
[`Iht::assign_genotypes`]({link(iht_rs, 442)}) and emitted to
`{{prefix}}.pass.vcf` at
[`gtg_concordance.rs:534`]({link(conc_rs, 534)}), and the failing
records are quarantined to `{{prefix}}.fail.vcf` as described above.
The split is the answer to "does `gtg-ped-map` reconstruct the 0/1
allele sequence of each haplotype?": no, and deliberately — that is
handled exclusively by `gtg-concordance`.
"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(content)
    print(f"[wiki] Wrote {out_path}")

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

The pipeline consists of two Rust binaries — `gtg-ped-map`
([`map_builder.rs`]({map_rs_main})) and `gtg-concordance`
([`gtg_concordance.rs`]({conc_rs_main})) — that share a common
inheritance-tracking library (`code/rust/src/iht.rs`) and are driven by a
standard PED file plus a jointly-called VCF. The two binaries have
deliberately separate responsibilities: `gtg-ped-map` labels haplotypes
with structural letters at informative sites only, and
`gtg-concordance` is the sole place where letter→allele correspondence is
resolved and phased genotypes are written. Everything in between — the
per-site letter propagation, the block-collapse and flip machinery, and
the orientation search — is described below.

All line numbers refer to commit `{sha_short}`; every function link is
paired with its driver call site in `main()` so the manuscript can be
read in parallel with the source. Worked toy examples for each major
component live in the companion walkthrough pages:
[nuclear family](nuclear_family/nuclear_family.md),
[three-generation pedigree](three_generations/three_generations.md),
and [gtg-concordance](concordance/concordance.md).

## 1. Inputs

The tools consume two files:

- A **PED file** describing the pedigree, parsed by
  [`code/rust/src/ped.rs`]({ped_rs}). Each row declares an individual's
  family, parents, sex, and affection status; founders are rows whose
  parent columns are `0`.
- A **jointly-called, tabix-indexed VCF** — typically DeepVariant on
  HiFi — with genotypes for every individual named in the PED. The VCF
  must be jointly called so that every record has a defined GT for
  every sample; per-sample VCFs merged after the fact will not work
  because the informative-site test below requires the parents' and
  children's genotypes at the same coordinate in a single record.

Only biallelic SNVs enter map construction. Indels are filtered at read
time by [`is_indel`]({map_rs_501}) (invoked from the VCF-reading loop
at [`map_builder.rs:164`]({map_rs_164}) inside `parse_vcf`; the driver
calls `parse_vcf` at [`map_builder.rs:1092`]({map_rs_1092})).
Multi-allelic records are likewise dropped at parse time — the letter
machinery is built around a two-homolog, two-allele model and does not
generalise to three or more alleles per site.

## 2. gtg-ped-map: structural haplotype labeling

### 2.1 Founder letter assignment

At startup the driver calls [`Iht::new`]({iht_rs_172}) once to build a
master inheritance template
([`map_builder.rs:1059`]({map_rs_1059})) that fixes the founder
alphabet for the whole run. For each founder in depth order it allocates
a fresh pair of capital letters — `(A,B)` for founder 0, `(C,D)` for
founder 1, `(E,F)` for founder 2, and so on — and assigns one letter to
the paternal homolog and one to the maternal homolog. Non-founders are
initialised with `?` slots that will be filled during the per-site walk.
The master template is never mutated: only its
[`legend()`]({iht_rs_330}) is read, to print the column header
(`Dad:A|B Mom:C|D Kid1:?|? …`) at the top of the IHT and marker
output files. Per-site bookkeeping happens on a separate object — the
driver re-invokes `Iht::new` per VCF record at
[`map_builder.rs:1111`]({map_rs_1111}) to allocate a fresh `local_iht`
that [`track_alleles_through_pedigree`]({map_rs_295}) then *mutates* in
place to record which founder letter each child inherited at that
site. Two reasons a per-site copy is needed rather than reusing the
master:

1. **Each site needs its own mutable IHT vector** — that vector *is*
   the per-site output, so it cannot be shared across sites.
2. **Zygosity is per-chromosome.** The master is hard-coded to
   `ChromType::Autosome` (it only feeds the header), whereas
   `local_iht` is built with the chromosome's actual `zygosity`
   (autosome vs. chrX, decided at
   [`map_builder.rs:1086`]({map_rs_1086})), which changes how letters
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

For every VCF record, [`track_alleles_through_pedigree`]({map_rs_295})
(driver call at [`map_builder.rs:1116`]({map_rs_1116})) walks
individuals in ancestor-first depth order, using the depth ordering
returned by `family.get_individual_depths()`. For every
`(parent, spouse)` pair in the walk it calls
[`unique_allele`]({map_rs_243}) (from inside the walk at
[`map_builder.rs:315`]({map_rs_315})) to test whether the parent
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
[`get_iht_markers`]({map_rs_274}) (invoked from inside the walk at
[`map_builder.rs:328`]({map_rs_328})) is the routine that performs the
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
at a given site, [`backfill_sibs`]({map_rs_804}) (driver call at
[`map_builder.rs:1122`]({map_rs_1122})) infers the other founder
allele on the remaining children: a sibling who does not carry the
tagged allele must have inherited the other founder homolog, so the
missing `?` can be replaced with the other letter of the parent's pair.
On the toy nuclear-family simulation this is a no-op because every
informative site already tags all three children, but on real data it
is essential for keeping the map dense through noisy regions.

### 2.5 Block collapse and phase flipping

After the per-site letter trace is complete, adjacent sites with
compatible letter assignments are merged into contiguous blocks by
[`collapse_identical_iht`]({map_rs_385}) (driver call at
[`map_builder.rs:1191`]({map_rs_1191})). A block is a maximal run of
sites over which every individual carries the same letter pair; a
recombination event in any individual forces a block boundary.

Because the two letters in a founder's pair `(A,B)` are interchangeable
within any single block — nothing in the combinatorial machinery of §2.2
distinguishes "A-then-B" from "B-then-A" — neighbouring blocks can
disagree on which letter is "first" even when the underlying biology is
continuous. [`perform_flips_in_place`]({map_rs_702}) resolves this by
walking the block list left-to-right and, for every founder, swapping
that founder's letter pair inside a block when doing so reduces the
per-individual letter mismatch with the previous block. The driver
invokes it three times —
[`map_builder.rs:1135`]({map_rs_1135}) (before block collapse),
[`map_builder.rs:1193`]({map_rs_1193}) (after block collapse), and
[`map_builder.rs:1203`]({map_rs_1203}) (after gap fill) — so that each
subsequent step sees consistently-oriented labels.

Gaps left by `?` slots (depth-filtered sites, backfilled misses) are
then filled by [`fill_missing_values`]({map_rs_617}) (driver call at
[`map_builder.rs:1200`]({map_rs_1200})) and
[`fill_missing_values_by_neighbor`]({map_rs_540}) (driver call at
[`map_builder.rs:1201`]({map_rs_1201})), which propagate letters from
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

[`count_matching_neighbors`]({map_rs_935}) (driver call at
[`map_builder.rs:1172`]({map_rs_1172})) scans each individual's
per-site letter sequence and flags runs shorter than `--run` (default
10 markers) whose flanking letters agree.
[`mask_child_alleles`]({map_rs_970}) (driver call at
[`map_builder.rs:1187`]({map_rs_1187})) then sets those slots back to
`?` before block collapse, so that the masked flicker neither creates
spurious blocks nor contaminates adjacent ones. The default threshold is
a trade-off: higher `--run` values suppress more noise but also mask
genuine short gene-conversion tracts.

### 2.7 Recombination reporting

After block collapse and flip normalisation,
[`summarize_child_changes`]({map_rs_673}) (driver call at
[`map_builder.rs:1228`]({map_rs_1228})) walks each individual's block
sequence and emits every letter transition — `A → B` on the paternal
slot, `C → D` on the maternal slot — to `{{prefix}}.recombinants.txt`.

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

- `{{prefix}}.iht.txt` — the final block map: contiguous physical
  intervals with a letter pair for every individual and a marker count.
  This is the sole input that `gtg-concordance` consumes.
- `{{prefix}}.markers.txt` — raw per-site letter assignments before
  collapse, useful for breakpoint refinement and for auditing decisions
  made by `mask_child_alleles`.
- `{{prefix}}.recombinants.txt` — the candidate letter transitions
  described in §2.7.

Nothing here contains 0/1 allele sequences — the block map is a pure
structural labelling. That is not an oversight: it is the precondition
that lets `gtg-concordance` do its job cleanly.

## 3. gtg-concordance: allele phasing and concordance QC

The driver [`main()`]({conc_rs_main}) reads the block map with
[`parse_ihtv2_file`]({iht_rs_606}) (driver call at
[`gtg_concordance.rs:405`]({conc_rs_405})) and, for each block, fetches
every overlapping VCF record — including the non-informative sites that
`gtg-ped-map` skipped — via the per-site phasing loop at
[`gtg_concordance.rs:437`]({conc_rs_437}). Each record is then phased
against the block's letter labels by the routine described below.
Records that fail depth or quality filters are routed directly to
`{{prefix}}.fail.vcf` at [`gtg_concordance.rs:444`]({conc_rs_444})
(low-quality) and [`gtg_concordance.rs:450`]({conc_rs_450}) (no-call),
and never enter the orientation search.

### 3.1 Orientation search

For every record that clears the input filters,
[`find_best_phase_orientation`]({conc_rs_252}) (driver call at
[`gtg_concordance.rs:454`]({conc_rs_454})) enumerates the `2^F`
founder-phase orientations produced by
[`Iht::founder_phase_orientations`]({iht_rs_492}) (invoked inside
`find_best_phase_orientation` at
[`gtg_concordance.rs:256`]({conc_rs_256})). An orientation is a choice,
for each of the `F` founders independently, of which of the founder's
two sorted VCF alleles is tagged with its first letter and which with
its second. The `2^F` enumeration is exhaustive: every possible
letter→allele correspondence consistent with the block's structural
labels is tried.

Under a given orientation, [`Iht::assign_genotypes`]({iht_rs_442})
produces the **expected** genotype for every individual: each founder
contributes `letter1 → vcf_allele[0]`, `letter2 → vcf_allele[1]`, and
each descendant's two letters are looked up in that map to give an
expected `{{a1, a2}}` pair. The expected genotypes are compared to the
observed (unphased) ones by
[`compare_genotype_maps`]({conc_rs_213}) (driver call at
[`gtg_concordance.rs:268`]({conc_rs_268})), which counts the number
of samples whose observed pair disagrees with the expected pair under
unphased set-equality. The orientation with the fewest mismatches wins.

### 3.2 Pass, fail, and the "impossible genotype" rule

Two outcomes are possible:

- **Zero mismatches**: the winning orientation explains every sample's
  genotype exactly. The driver uses the winning letter→allele map —
  re-applied via [`Iht::assign_genotypes`]({iht_rs_442}) on the
  passing branch at [`gtg_concordance.rs:514`]({conc_rs_514}) — to
  rewrite the record as phased (`0|1`, paternal first) and write it to
  `{{prefix}}.pass.vcf` at
  [`gtg_concordance.rs:534`]({conc_rs_534}).
- **Non-zero mismatches under every orientation**: no assignment of
  alleles to founder letters explains the whole family. The site is
  **impossible** under the current block labels. The driver takes the
  best-available orientation — `Iht::assign_genotypes` is still called
  on the failing branch at
  [`gtg_concordance.rs:487`]({conc_rs_487}) to identify the offending
  samples — and writes the unphased record to `{{prefix}}.fail.vcf` at
  [`gtg_concordance.rs:507`]({conc_rs_507}). The offending sample
  names are appended to `{{prefix}}.failed_sites.txt`, and per-block
  counts of sites where exactly one sample is the culprit ("singleton"
  failures) are tallied in `{{prefix}}.filtering_stats.txt`.

The impossible-genotype rule is the mechanism by which
`gtg-concordance` refuses to produce phased calls that the block's
structural labels cannot justify. It simultaneously catches sequencing
errors, Mendelian violations, and residual block-labelling mistakes
upstream — all three manifest as "no orientation works".

The division of labour is strict. `gtg-ped-map` stores only letters and
only at informative sites; `gtg-concordance` is the sole place where
letter→allele correspondence is computed and written out. Any
downstream tool that needs phased genotypes must consume
`{{prefix}}.pass.vcf`, never the `.iht.txt` structural map.

## 4. Real-world issues

### 4.1 Sequencing errors

Short-read and HiFi sequencing errors flip individual genotypes
stochastically. They enter the pipeline in two different places and are
handled by two different mechanisms.

Inside `gtg-ped-map`, a single erroneous site in one child presents as
an `A → B → A` flicker in that child's pre-collapse letter trace. The
[`count_matching_neighbors`]({map_rs_935}) /
[`mask_child_alleles`]({map_rs_970}) pair (driver calls at
[`map_builder.rs:1172`]({map_rs_1172}) and
[`map_builder.rs:1187`]({map_rs_1187})) recognises runs shorter than
`--run` (default 10 markers) and masks them before
`collapse_identical_iht` runs. Without this filter, each isolated
error would surface as two adjacent recombination breakpoints.

Inside `gtg-concordance`, sequencing errors show up as sites where no
orientation satisfies the whole family. These are written to
`{{prefix}}.fail.vcf` by the impossible-genotype rule described in
§3.2. Per-block singleton-failure tallies in
`{{prefix}}.filtering_stats.txt` distinguish a single noisy sample
(random errors) from a block-wide labelling mistake (systematic
failure across many samples), because the two fail the orientation
search with very different footprints.

### 4.2 Depth and quality filtering

`extract_depth_statistics` computes per-sample mean and standard
deviation over each chromosome.
[`depth_filters`]({map_rs_217}) drops any site where a sample is below
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
underlying biology is continuous. [`perform_flips_in_place`]({map_rs_702})
resolves this by pair-swapping founder letters inside a block to
minimise mismatches with the previous block. The driver calls it three
times — [`map_builder.rs:1135`]({map_rs_1135}),
[`map_builder.rs:1193`]({map_rs_1193}), and
[`map_builder.rs:1203`]({map_rs_1203}) — once before collapse, once
after, and once after gap fill, so that every step operates on
consistently-oriented blocks.

### 4.4 Missing data in sibships

Multi-child sibships provide redundant information that is essential in
noisy regions. [`backfill_sibs`]({map_rs_804}) (driver call at
[`map_builder.rs:1122`]({map_rs_1122})) uses the fact that if one
child in a sibship carries founder allele X, the remaining children
must each carry either X or the other founder allele, and a majority
vote across the sibship can recover a usable label for a child whose
own site-level signal is missing. Singleton sibships (one child only)
cannot benefit from this and remain more vulnerable to dropout.

### 4.5 Ancestral vs. per-meiosis recombinations

[`summarize_child_changes`]({map_rs_673}) (driver call at
[`map_builder.rs:1228`]({map_rs_1228})) reports every letter
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
[`track_alleles_through_pedigree`]({map_rs_346}) and
[`Iht::new`]({iht_rs_181}). Males receive a single non-dot letter on
the X slot, derived from mom; dads cannot transmit X to sons, so the
paternal-slot label on a male child's X is left undefined; and
[`perform_flips_in_place`]({map_rs_702}) skips X entirely to avoid
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
"""


def emit_methods_section(out_path: Path) -> None:
    def link(path: str, line: int) -> str:
        return permalink(path, line, SHA)

    content = METHODS_TEMPLATE.format(
        sha_short=SHA[:7],
        ped_rs=f"https://github.com/{REPO}/blob/{SHA}/code/rust/src/ped.rs",
        map_rs_main=link("code/rust/src/bin/map_builder.rs", 989),
        conc_rs_main=link("code/rust/src/bin/gtg_concordance.rs", 315),
        map_rs_164=link("code/rust/src/bin/map_builder.rs", 164),
        map_rs_217=link("code/rust/src/bin/map_builder.rs", 217),
        map_rs_243=link("code/rust/src/bin/map_builder.rs", 243),
        map_rs_274=link("code/rust/src/bin/map_builder.rs", 274),
        map_rs_295=link("code/rust/src/bin/map_builder.rs", 295),
        map_rs_315=link("code/rust/src/bin/map_builder.rs", 315),
        map_rs_328=link("code/rust/src/bin/map_builder.rs", 328),
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
        map_rs_1059=link("code/rust/src/bin/map_builder.rs", 1059),
        map_rs_1086=link("code/rust/src/bin/map_builder.rs", 1086),
        map_rs_1092=link("code/rust/src/bin/map_builder.rs", 1092),
        map_rs_1111=link("code/rust/src/bin/map_builder.rs", 1111),
        map_rs_1116=link("code/rust/src/bin/map_builder.rs", 1116),
        map_rs_1122=link("code/rust/src/bin/map_builder.rs", 1122),
        map_rs_1135=link("code/rust/src/bin/map_builder.rs", 1135),
        map_rs_1172=link("code/rust/src/bin/map_builder.rs", 1172),
        map_rs_1187=link("code/rust/src/bin/map_builder.rs", 1187),
        map_rs_1191=link("code/rust/src/bin/map_builder.rs", 1191),
        map_rs_1193=link("code/rust/src/bin/map_builder.rs", 1193),
        map_rs_1200=link("code/rust/src/bin/map_builder.rs", 1200),
        map_rs_1201=link("code/rust/src/bin/map_builder.rs", 1201),
        map_rs_1203=link("code/rust/src/bin/map_builder.rs", 1203),
        map_rs_1228=link("code/rust/src/bin/map_builder.rs", 1228),
        conc_rs_213=link("code/rust/src/bin/gtg_concordance.rs", 213),
        conc_rs_252=link("code/rust/src/bin/gtg_concordance.rs", 252),
        conc_rs_256=link("code/rust/src/bin/gtg_concordance.rs", 256),
        conc_rs_268=link("code/rust/src/bin/gtg_concordance.rs", 268),
        conc_rs_405=link("code/rust/src/bin/gtg_concordance.rs", 405),
        conc_rs_437=link("code/rust/src/bin/gtg_concordance.rs", 437),
        conc_rs_444=link("code/rust/src/bin/gtg_concordance.rs", 444),
        conc_rs_450=link("code/rust/src/bin/gtg_concordance.rs", 450),
        conc_rs_454=link("code/rust/src/bin/gtg_concordance.rs", 454),
        conc_rs_487=link("code/rust/src/bin/gtg_concordance.rs", 487),
        conc_rs_507=link("code/rust/src/bin/gtg_concordance.rs", 507),
        conc_rs_514=link("code/rust/src/bin/gtg_concordance.rs", 514),
        conc_rs_534=link("code/rust/src/bin/gtg_concordance.rs", 534),
        iht_rs_172=link("code/rust/src/iht.rs", 172),
        iht_rs_181=link("code/rust/src/iht.rs", 181),
        iht_rs_330=link("code/rust/src/iht.rs", 330),
        iht_rs_442=link("code/rust/src/iht.rs", 442),
        iht_rs_492=link("code/rust/src/iht.rs", 492),
        iht_rs_606=link("code/rust/src/iht.rs", 606),
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
2. [Three-generation pedigree with an outside marriage](three_generations/three_generations.md)
   — extends the nuclear family by marrying Kid2 to a fresh founder
   (Spouse) and adding a G3 sibship. Shows how the same single loop
   handles G2→G3 via ancestor-first depth ordering, without ever
   constructing a joint inheritance vector over all founders.
3. [Closing the loop with gtg-concordance](concordance/concordance.md)
   — closes the pipeline by mapping founder letters back to VCF
   alleles at the sites `gtg-ped-map` could not touch. Introduces the
   `2^F` founder-phase orientation search and the "impossible
   genotype" rule that routes sequencing errors to `fail.vcf`.

## Reference pages

- [Methods](methods.md) — the manuscript-style write-up with the
  real-world caveats (depth filtering, phase instability, sibship
  backfilling, chromosome X) that the walkthrough pages deliberately
  skip.

## Conventions

- **Filenames describe content, not sequence.** Reading order is
  defined above; subdirectories use short descriptive names
  (`nuclear_family/`) instead of numeric prefixes.
- **Figures live next to the page that references them.** Each
  walkthrough subdirectory holds the markdown plus its PNG figures, so
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
        "three_generations": lambda: component_2_three_generations(args.outdir),
        "concordance": lambda: component_3_concordance(args.outdir),
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
