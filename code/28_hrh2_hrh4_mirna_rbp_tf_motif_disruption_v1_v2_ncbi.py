#!/usr/bin/env python3

from pathlib import Path
import csv
import re
from collections import Counter, defaultdict

PROJECT = Path("/mnt/d/histamine_receptor_pharmacogenomics")

SEQ_TABLE = PROJECT / "06_regulatory_genomics/sequence_context_windows/HRH2_HRH4_regulatory_variant_sequence_windows_v1_v2_ncbi.tsv"

OUTDIR = PROJECT / "06_regulatory_genomics/motif_disruption"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_ALL = OUTDIR / "HRH2_HRH4_regulatory_motif_disruption_scan_v1_v2_ncbi.tsv"
OUT_CHANGED = OUTDIR / "HRH2_HRH4_regulatory_motif_disruption_changed_targets_v1_v2_ncbi.tsv"
OUT_SUMMARY = OUTDIR / "HRH2_HRH4_regulatory_motif_disruption_summary_v1_v2_ncbi.tsv"
OUT_TOP = OUTDIR / "HRH2_HRH4_top_motif_disruption_candidates_v1_v2_ncbi.tsv"

TABLE11 = PROJECT / "13_tables/manuscript_ready/Table_11_HRH2_HRH4_miRNA_RBP_TF_motif_disruption_candidates_v1_v2_ncbi.tsv"
RESULTS_MD = PROJECT / "14_manuscript/results_sections/results_hrh2_hrh4_motif_disruption_v1_v2_ncbi.md"

def read_tsv(path):
    with path.open("r", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

def safe_float(x):
    try:
        return float(x)
    except Exception:
        return 0.0

def parse_variant_key(vkey):
    parts = str(vkey).split(":", 3)
    if len(parts) != 4:
        return "", 0, "", ""
    return parts[0], int(parts[1]), parts[2].upper(), parts[3].upper()

def revcomp(seq):
    comp = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")
    return seq.translate(comp)[::-1].upper().replace("U", "T")

def rna_to_dna(seq):
    return seq.upper().replace("U", "T")

# Focused miRNA seed panel. Seeds are miRNA positions 2-8, represented in RNA/DNA alphabet.
# Motif searched in DNA windows is the reverse complement of the seed.
MIRNA_SEEDS = {
    "let-7_family": "GAGGTAG",
    "miR-155-5p": "TAATGCT",
    "miR-146a-5p": "GAGAACT",
    "miR-21-5p": "AGCTTAT",
    "miR-223-3p": "GTCAGTT",
    "miR-124-3p": "AAGGCAC",
    "miR-9-5p": "CTTTGGT",
    "miR-132-3p": "AACAGTC",
    "miR-181_family": "ACATTCA",
    "miR-29_family": "AGCACCA",
    "miR-16_family": "AGCAGCA",
    "miR-125_family": "CCCTGAG",
    "miR-34_family": "GGCAGTG",
    "miR-17_20_93_106_family": "AAAGTGC",
    "miR-200_family": "AATACTG",
    "miR-30_family": "GTAAACA",
    "miR-26_family": "TCAAGTA",
    "miR-27_family": "TCACAGT",
    "miR-221_222_family": "GCTACAT",
    "miR-142-3p": "GTAGTGT",
    "miR-150-5p": "CTCCCAA",
    "miR-126-3p": "CGTACCG",
    "miR-31-5p": "GGCAAGA",
    "miR-10_family": "ACCCTGT",
}

MOTIFS = []

# miRNA seed complement motifs.
for name, seed in MIRNA_SEEDS.items():
    seed_dna = rna_to_dna(seed)
    target = revcomp(seed_dna)
    MOTIFS.append({
        "motif_id": name,
        "motif_group": "miRNA_seed_7mer",
        "pattern": re.escape(target),
        "description": f"7mer reverse-complement seed site for {name}",
    })

# 3'UTR/RBP/cis-element motifs in DNA alphabet.
TRANSCRIPT_MOTIFS = {
    "AU_rich_element_ATTTA": "ATTTA",
    "AU_rich_element_TATTTAT": "TATTTAT",
    "AU_rich_element_TTATTTAT": "TTATTTAT",
    "polyA_signal_AATAAA": "AATAAA",
    "polyA_signal_ATTAAA": "ATTAAA",
    "polyA_signal_TATAAA": "TATAAA",
    "CPE_like_TTTTAT": "TTTTAT",
    "CPE_like_TTTTA": "TTTTA",
    "U_rich_TTTTT": "TTTTT",
    "GU_rich_GTGTG": "GTGTG",
    "GU_rich_TGTGT": "TGTGT",
    "C_rich_CCCCC": "CCCCC",
    "C_rich_CCTCC": "CCTCC",
    "C_rich_CCCTC": "CCCTC",
}

for name, pat in TRANSCRIPT_MOTIFS.items():
    MOTIFS.append({
        "motif_id": name,
        "motif_group": "UTR_RBP_cis_element",
        "pattern": re.escape(pat),
        "description": name,
    })

# Compact core TF/cCRE motif screen. These are consensus-core screens, not definitive TF binding calls.
TF_PATTERNS = {
    "AP1_core_TGACTCA": "TGACTCA",
    "AP1_core_TGAGTCA": "TGAGTCA",
    "CREB_core_TGACGTCA": "TGACGTCA",
    "Ebox_CACGTG": "CACGTG",
    "Ebox_CAGCTG": "CAGCTG",
    "SP1_GCbox_GGGCGG": "GGGCGG",
    "SP1_GCbox_CCGCCC": "CCGCCC",
    "STAT_like_TTCNNNGAA": "TTC[ACGT]{3}GAA",
    "IRF_like_GAAA_N2_GAAA": "GAAA[ACGT]{2}GAAA",
}

for name, pat in TF_PATTERNS.items():
    MOTIFS.append({
        "motif_id": name,
        "motif_group": "TF_cCRE_core_motif",
        "pattern": pat,
        "description": name,
    })

def find_motifs(seq, motif):
    hits = []
    if not seq:
        return hits
    pattern = re.compile(motif["pattern"])
    for m in pattern.finditer(seq.upper()):
        hits.append({
            "motif_id": motif["motif_id"],
            "motif_group": motif["motif_group"],
            "description": motif["description"],
            "start": m.start(),
            "end": m.end(),
            "match": m.group(0),
        })
    return hits

def hit_key(hit):
    return f"{hit['motif_group']}|{hit['motif_id']}|{hit['start']}|{hit['end']}|{hit['match']}"

def motif_summary(hits):
    counts = Counter(f"{h['motif_group']}:{h['motif_id']}" for h in hits)
    return ";".join(f"{k}:{v}" for k, v in sorted(counts.items()))

def motif_names(hits):
    return ";".join(sorted(set(f"{h['motif_group']}:{h['motif_id']}" for h in hits)))

def variant_span_in_window(row, ref=True):
    chrom, pos, ref_allele, alt_allele = parse_variant_key(row.get("variant_key", ""))
    win_start = int(row.get("ucsc_window_start", 0) or 0)
    local_start = pos - 1 - win_start
    allele = ref_allele if ref else alt_allele
    local_end = local_start + max(1, len(allele))
    return local_start, local_end

def overlaps(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start)) > 0

def changed_overlapping_variant(changes, row, ref=True):
    v_start, v_end = variant_span_in_window(row, ref=ref)
    return [h for h in changes if overlaps(h["start"], h["end"], v_start, v_end)]

def classify_direction(ref_hits, alt_hits):
    ref_keys = {hit_key(h): h for h in ref_hits}
    alt_keys = {hit_key(h): h for h in alt_hits}

    lost = [ref_keys[k] for k in sorted(set(ref_keys) - set(alt_keys))]
    gained = [alt_keys[k] for k in sorted(set(alt_keys) - set(ref_keys))]

    return lost, gained

def score_row(row, lost, gained, lost_var, gained_var):
    score = 0

    all_changed = lost + gained
    groups = Counter(h["motif_group"] for h in all_changed)

    if groups.get("miRNA_seed_7mer", 0):
        score += 4
    if groups.get("UTR_RBP_cis_element", 0):
        score += 3
    if groups.get("TF_cCRE_core_motif", 0):
        score += 2

    if lost_var or gained_var:
        score += 3

    labels = " ".join([row.get("labels", ""), row.get("source_layers", ""), row.get("evidence_layers", "")]).lower()
    if "3prime" in labels or "utr" in labels:
        score += 2
    if "ccre" in labels or row.get("ccre_class"):
        score += 2

    tier = row.get("regulatory_deepening_priority", "")
    if "Tier_1" in tier:
        score += 3
    elif "Tier_2" in tier:
        score += 2
    elif "Tier_3" in tier:
        score += 1

    pop_class = row.get("population_structure_class", "")
    if pop_class == "very_strong_population_structure":
        score += 3
    elif pop_class == "strong_population_structure":
        score += 2
    elif pop_class in {"moderate_population_structure", "private_or_near_private_superpopulation_signal"}:
        score += 1

    if score >= 14:
        priority = "MotifTier_1_high_priority"
    elif score >= 10:
        priority = "MotifTier_2_moderate_priority"
    elif score >= 6:
        priority = "MotifTier_3_context_priority"
    else:
        priority = "MotifTier_4_background"

    return score, priority

def main():
    rows = read_tsv(SEQ_TABLE)

    out_rows = []

    for row in rows:
        ref_seq = row.get("ref_sequence_window", "").upper()
        alt_seq = row.get("alt_sequence_window", "").upper()

        ref_hits = []
        alt_hits = []
        for motif in MOTIFS:
            ref_hits.extend(find_motifs(ref_seq, motif))
            alt_hits.extend(find_motifs(alt_seq, motif))

        lost, gained = classify_direction(ref_hits, alt_hits)
        lost_var = changed_overlapping_variant(lost, row, ref=True)
        gained_var = changed_overlapping_variant(gained, row, ref=False)

        score, priority = score_row(row, lost, gained, lost_var, gained_var)

        changed_groups = Counter(h["motif_group"] for h in lost + gained)

        out = {
            "gene": row.get("gene", ""),
            "variant_key": row.get("variant_key", ""),
            "gtex_variant_id_b38": row.get("gtex_variant_id_b38", ""),
            "regulatory_deepening_priority": row.get("regulatory_deepening_priority", ""),
            "regulatory_deepening_score": row.get("regulatory_deepening_score", ""),
            "motif_disruption_priority": priority,
            "motif_disruption_score": score,
            "has_any_motif_change": bool(lost or gained),
            "has_variant_overlapping_motif_change": bool(lost_var or gained_var),
            "changed_motif_group_counts": ";".join(f"{k}:{v}" for k, v in sorted(changed_groups.items())),
            "lost_motifs": motif_names(lost),
            "gained_motifs": motif_names(gained),
            "lost_variant_overlapping_motifs": motif_names(lost_var),
            "gained_variant_overlapping_motifs": motif_names(gained_var),
            "ref_motif_summary": motif_summary(ref_hits),
            "alt_motif_summary": motif_summary(alt_hits),
            "evidence_layers": row.get("evidence_layers", ""),
            "source_layers": row.get("source_layers", ""),
            "labels": row.get("labels", ""),
            "candidate_tier": row.get("candidate_tier", ""),
            "selected_hgvsc": row.get("selected_hgvsc", ""),
            "ccre_name": row.get("ccre_name", ""),
            "ccre_class": row.get("ccre_class", ""),
            "ccre_label": row.get("ccre_label", ""),
            "best_anchor": row.get("best_anchor", ""),
            "best_anchor_group": row.get("best_anchor_group", ""),
            "best_ld_r2": row.get("best_ld_r2", ""),
            "global_AF": row.get("global_AF", ""),
            "superpop_AF_range": row.get("superpop_AF_range", ""),
            "population_structure_class": row.get("population_structure_class", ""),
            "max_pairwise_fst_hudson": row.get("max_pairwise_fst_hudson", ""),
            "max_pairwise_fst_pair": row.get("max_pairwise_fst_pair", ""),
            "top_superpopulation": row.get("top_superpopulation", ""),
            "top_superpopulation_AF": row.get("top_superpopulation_AF", ""),
            "ucsc_window_chrom": row.get("ucsc_window_chrom", ""),
            "ucsc_window_start": row.get("ucsc_window_start", ""),
            "ucsc_window_end": row.get("ucsc_window_end", ""),
            "observed_ref_in_ucsc": row.get("observed_ref_in_ucsc", ""),
            "reference_match": row.get("reference_match", ""),
            "ref_sequence_window": row.get("ref_sequence_window", ""),
            "alt_sequence_window": row.get("alt_sequence_window", ""),
        }
        out_rows.append(out)

    out_rows = sorted(
        out_rows,
        key=lambda r: (
            int(r["motif_disruption_score"]),
            safe_float(r.get("max_pairwise_fst_hudson", "")),
            safe_float(r.get("superpop_AF_range", "")),
        ),
        reverse=True,
    )

    fields = [
        "gene", "variant_key", "gtex_variant_id_b38",
        "regulatory_deepening_priority", "regulatory_deepening_score",
        "motif_disruption_priority", "motif_disruption_score",
        "has_any_motif_change", "has_variant_overlapping_motif_change",
        "changed_motif_group_counts",
        "lost_motifs", "gained_motifs",
        "lost_variant_overlapping_motifs", "gained_variant_overlapping_motifs",
        "ref_motif_summary", "alt_motif_summary",
        "evidence_layers", "source_layers", "labels", "candidate_tier",
        "selected_hgvsc", "ccre_name", "ccre_class", "ccre_label",
        "best_anchor", "best_anchor_group", "best_ld_r2",
        "global_AF", "superpop_AF_range",
        "population_structure_class", "max_pairwise_fst_hudson", "max_pairwise_fst_pair",
        "top_superpopulation", "top_superpopulation_AF",
        "ucsc_window_chrom", "ucsc_window_start", "ucsc_window_end",
        "observed_ref_in_ucsc", "reference_match",
        "ref_sequence_window", "alt_sequence_window",
    ]

    changed = [r for r in out_rows if str(r["has_any_motif_change"]) == "True"]

    write_tsv(OUT_ALL, out_rows, fields)
    write_tsv(OUT_CHANGED, changed, fields)
    write_tsv(OUT_TOP, out_rows[:60], fields)
    write_tsv(TABLE11, out_rows[:40], fields)

    summary_rows = []

    def add_summary(group, subset):
        summary_rows.append({
            "group": group,
            "targets": len(subset),
            "targets_with_any_motif_change": sum(1 for r in subset if str(r["has_any_motif_change"]) == "True"),
            "targets_with_variant_overlapping_motif_change": sum(1 for r in subset if str(r["has_variant_overlapping_motif_change"]) == "True"),
            "motif_priority_counts": ";".join(f"{k}:{v}" for k, v in sorted(Counter(r["motif_disruption_priority"] for r in subset).items())),
            "changed_group_counts": ";".join(
                f"{k}:{v}" for k, v in sorted(Counter(
                    group_count.split(":")[0]
                    for r in subset
                    for group_count in str(r.get("changed_motif_group_counts", "")).split(";")
                    if group_count
                ).items())
            ),
        })

    add_summary("ALL", out_rows)
    for gene in sorted(set(r["gene"] for r in out_rows)):
        add_summary(f"gene={gene}", [r for r in out_rows if r["gene"] == gene])
    for tier in sorted(set(r["regulatory_deepening_priority"] for r in out_rows)):
        add_summary(f"regulatory_tier={tier}", [r for r in out_rows if r["regulatory_deepening_priority"] == tier])
    for mtier in sorted(set(r["motif_disruption_priority"] for r in out_rows)):
        add_summary(f"motif_tier={mtier}", [r for r in out_rows if r["motif_disruption_priority"] == mtier])

    write_tsv(
        OUT_SUMMARY,
        summary_rows,
        [
            "group", "targets", "targets_with_any_motif_change",
            "targets_with_variant_overlapping_motif_change",
            "motif_priority_counts", "changed_group_counts"
        ]
    )

    gene_counts = Counter(r["gene"] for r in out_rows)
    motif_tier_counts = Counter(r["motif_disruption_priority"] for r in out_rows)
    changed_count = len(changed)
    variant_overlap_count = sum(1 for r in out_rows if str(r["has_variant_overlapping_motif_change"]) == "True")

    changed_group_counts = Counter()
    for r in changed:
        for item in str(r.get("changed_motif_group_counts", "")).split(";"):
            if not item:
                continue
            group_name = item.split(":", 1)[0]
            changed_group_counts[group_name] += 1

    top_lines = []
    for r in out_rows[:15]:
        top_lines.append(
            f"- {r['gene']} {r['variant_key']}: {r['motif_disruption_priority']}; "
            f"score={r['motif_disruption_score']}; motif_changes={r['changed_motif_group_counts'] or 'none'}; "
            f"lost={r['lost_motifs'] or 'none'}; gained={r['gained_motifs'] or 'none'}."
        )

    md = f"""# HRH2/HRH4 sequence-level miRNA/RBP/TF motif-disruption screen

Validated REF/ALT sequence windows for HRH2/HRH4 regulatory-deepening targets were screened for sequence-level motif changes. The analysis used a focused miRNA seed-site panel, 3prime UTR/RBP cis-element motifs and compact TF/cCRE core motifs. This is a prioritization screen rather than a definitive binding assay; motif changes require validation using full transcript context, strand-aware regulatory annotation, expression data and experimental assays.

The scan evaluated {len(out_rows)} sequence-validated targets. Gene-level counts were: {', '.join(f'{k}: {v}' for k, v in sorted(gene_counts.items()))}. Overall, {changed_count} targets showed at least one REF-versus-ALT motif change, and {variant_overlap_count} targets showed a motif change overlapping the variant allele itself. Motif-disruption priority tiers were: {', '.join(f'{k}: {v}' for k, v in sorted(motif_tier_counts.items()))}. Changed motif groups were: {', '.join(f'{k}: {v}' for k, v in sorted(changed_group_counts.items()))}.

This module provides the first sequence-level regulatory-disruption layer for the HRH2/HRH4 regulatory axis. The most interpretable candidates are those where motif change co-occurs with exact UTR context, cCRE overlap, formal population structure or LD support. These candidates should be prioritized for database-based miRNA/eQTL/sQTL follow-up and experimental regulatory reporter testing.

## Top motif-disruption candidates

{chr(10).join(top_lines)}

## Output files

- Full motif scan: `06_regulatory_genomics/motif_disruption/HRH2_HRH4_regulatory_motif_disruption_scan_v1_v2_ncbi.tsv`
- Changed targets: `06_regulatory_genomics/motif_disruption/HRH2_HRH4_regulatory_motif_disruption_changed_targets_v1_v2_ncbi.tsv`
- Motif summary: `06_regulatory_genomics/motif_disruption/HRH2_HRH4_regulatory_motif_disruption_summary_v1_v2_ncbi.tsv`
- Manuscript-ready Table 11: `13_tables/manuscript_ready/Table_11_HRH2_HRH4_miRNA_RBP_TF_motif_disruption_candidates_v1_v2_ncbi.tsv`
"""
    RESULTS_MD.write_text(md, encoding="utf-8")

    print(f"Wrote full motif scan: {OUT_ALL}")
    print(f"Wrote changed motif targets: {OUT_CHANGED}")
    print(f"Wrote motif summary: {OUT_SUMMARY}")
    print(f"Wrote top motif candidates: {OUT_TOP}")
    print(f"Wrote Table 11: {TABLE11}")
    print(f"Wrote Results subsection: {RESULTS_MD}")
    print()
    print(f"Targets scanned: {len(out_rows)}")
    print(f"Targets with any motif change: {changed_count}")
    print(f"Targets with variant-overlapping motif change: {variant_overlap_count}")
    print("Gene counts:")
    for k, v in sorted(gene_counts.items()):
        print(f"  {k}: {v}")
    print("Motif tier counts:")
    for k, v in sorted(motif_tier_counts.items()):
        print(f"  {k}: {v}")
    print("Changed motif group counts:")
    for k, v in sorted(changed_group_counts.items()):
        print(f"  {k}: {v}")

if __name__ == "__main__":
    main()
