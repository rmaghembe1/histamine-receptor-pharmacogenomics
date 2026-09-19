#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import math
from collections import defaultdict, Counter

PROJECT = Path("/mnt/d/histamine_receptor_pharmacogenomics")

BASE_EVIDENCE = PROJECT / "11_results/monitoring_evaluation_v4_focused/hrh_focused_functional_candidate_evidence_matrix_v2.tsv"
PAIRWISE = PROJECT / "11_results/monitoring_evaluation_v4_focused/hrh_focused_candidate_pairwise_population_tests_v2.tsv"

OUTDIR = PROJECT / "11_results/monitoring_evaluation_v4_precision"
FIGDIR = PROJECT / "12_figures/monitoring_evaluation_v4_precision"

OUT_EVIDENCE = OUTDIR / "hrh_precision_functional_candidate_evidence_matrix_v1.tsv"
OUT_LAYER = OUTDIR / "hrh_precision_evidence_layer_summary_v1.tsv"
OUT_RECEPTOR = OUTDIR / "hrh_precision_receptor_axis_summary_v1.tsv"
OUT_VALIDATED = OUTDIR / "hrh_precision_validated_vs_exploratory_evidence_table_v1.tsv"
OUT_AUDIT = OUTDIR / "hrh_precision_motif_gtex_source_audit_v1.tsv"
OUT_TRUE_MOTIF = OUTDIR / "hrh_true_motif_disruption_candidates_v1.tsv"
OUT_TRUE_GTEX = OUTDIR / "hrh_true_gtex_eqtl_intersection_candidates_v1.tsv"

TABLE56 = PROJECT / "13_tables/manuscript_ready/Table_56_precision_evidence_layer_summary_v1.tsv"
TABLE57 = PROJECT / "13_tables/manuscript_ready/Table_57_precision_receptor_axis_summary_v1.tsv"
TABLE58 = PROJECT / "13_tables/manuscript_ready/Table_58_precision_validated_vs_exploratory_evidence_v1.tsv"
TABLE59 = PROJECT / "13_tables/manuscript_ready/Table_59_true_motif_and_gtex_regulatory_candidates_v1.tsv"

RESULTS_MD = PROJECT / "14_manuscript/results_sections/results_precision_population_functional_monitoring_v1.md"

CORR_ROOT = PROJECT / (
    "14_manuscript/PLOS_ONE_revision_PONE-D-26-38698_20260918/"
    "07_corrected_analysis/R2-04_motif_indel_aware_v1/"
    "02_corrected_precision"
)

OUTDIR = CORR_ROOT
FIGDIR = CORR_ROOT / "figures"

OUT_EVIDENCE = OUTDIR / "hrh_precision_functional_candidate_evidence_matrix_indel_aware_v1.tsv"
OUT_LAYER = OUTDIR / "hrh_precision_evidence_layer_summary_indel_aware_v1.tsv"
OUT_RECEPTOR = OUTDIR / "hrh_precision_receptor_axis_summary_indel_aware_v1.tsv"
OUT_VALIDATED = OUTDIR / "hrh_precision_validated_vs_exploratory_indel_aware_v1.tsv"
OUT_AUDIT = OUTDIR / "hrh_precision_explicit_source_audit_indel_aware_v1.tsv"
OUT_TRUE_MOTIF = OUTDIR / "hrh_true_motif_disruption_candidates_indel_aware_v1.tsv"
OUT_TRUE_GTEX = OUTDIR / "hrh_true_gtex_eqtl_intersection_candidates_indel_aware_v1.tsv"

TABLE56 = OUTDIR / "tables/Table_56_precision_evidence_layer_summary_indel_aware_v1.tsv"
TABLE57 = OUTDIR / "tables/Table_57_precision_receptor_axis_summary_indel_aware_v1.tsv"
TABLE58 = OUTDIR / "tables/Table_58_precision_validated_vs_exploratory_indel_aware_v1.tsv"
TABLE59 = OUTDIR / "tables/Table_59_true_motif_and_gtex_regulatory_candidates_indel_aware_v1.tsv"

RESULTS_MD = OUTDIR / "results_precision_population_functional_indel_aware_v1.md"


FIG_PNG = FIGDIR / "Figure_11_precision_population_functional_monitoring_v1_900dpi.png"
FIG_SVG = FIGDIR / "Figure_11_precision_population_functional_monitoring_v1.svg"

GENES = ["HRH1", "HRH2", "HRH3", "HRH4"]

EVIDENCE_CLASSES = [
    "protein",
    "splice",
    "biogenesis_fate",
    "motif_disruption_true",
    "motif_context_only",
    "ccre",
    "gtex_true_hit",
    "gtex_lookup_only",
    "ld_context",
    "docking_variant",
    "population",
    "regulatory_priority",
]

COLORS = {
    "protein": "#1b9e77",
    "splice": "#d95f02",
    "biogenesis_fate": "#7570b3",
    "motif_disruption_true": "#e7298a",
    "motif_context_only": "#f781bf",
    "ccre": "#66a61e",
    "gtex_true_hit": "#e6ab02",
    "gtex_lookup_only": "#ffff99",
    "ld_context": "#a6761d",
    "docking_variant": "#666666",
    "population": "#1f78b4",
    "regulatory_priority": "#b15928",
}

def read_tsv(path):
    with path.open("r", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

def norm(x):
    return str(x or "").strip()

def ffloat(x):
    try:
        return float(x)
    except Exception:
        return None

def fmt(x, digits=6):
    if x is None or x == "":
        return ""
    return f"{float(x):.{digits}f}"

def row_text(row):
    return " ".join(norm(v) for v in row.values())

def variant_key_from_text(text):
    s = norm(text)
    m = re.search(r"\b(?:chr)?([0-9XYM]+):([0-9]+):([ACGTN]+):([ACGTN]+)\b", s, flags=re.I)
    if not m:
        return ""
    chrom, pos, ref, alt = m.groups()
    return f"{chrom.replace('chr','')}:{pos}:{ref.upper()}:{alt.upper()}"

def get_gene(row):
    for v in row.values():
        vv = norm(v)
        if vv in set(GENES):
            return vv
    txt = row_text(row)
    for g in GENES:
        if g in txt:
            return g
    return ""

def boolish(x):
    return norm(x).lower() in {"true", "1", "yes", "y"}

def compact_note(row):
    keep = []
    for k, v in row.items():
        kk = k.lower()
        vv = norm(v)
        if not vv:
            continue
        if any(x in kk for x in [
            "gene", "variant", "motif", "changed", "lost", "gained", "tier",
            "score", "tissue", "pvalue", "p_value", "effect", "nes",
            "phenotype", "ccre", "hgvsc", "hgvsp"
        ]):
            keep.append(f"{k}={vv}")
    return " | ".join(keep[:16])

def source_kind(path):
    s = str(path.relative_to(PROJECT)).lower()
    if "motif_disruption_changed_targets" in s or "top_motif_disruption_candidates" in s or "table_11" in s:
        return "motif_candidate_table"
    if "motif_disruption_scan" in s:
        return "motif_scan"
    if "gtex_target_intersections" in s or "gtex_gene_level_eqtl_sqtl_intersections" in s or "table_12_hrh2_hrh4_gtex_gene_level" in s:
        return "gtex_true_intersection"
    if "gtex_eqtl_sqtl_lookup_hits" in s:
        return "gtex_true_hit_table"
    if "gtex_eqtl_sqtl_lookup_audit" in s:
        return "gtex_lookup_audit"
    if "gtex_gene_level_eqtl_sqtl_records" in s:
        return "gtex_gene_level_records"
    return "other"

def discover_sources():
    sources = [
        PROJECT / (
            "14_manuscript/PLOS_ONE_revision_PONE-D-26-38698_20260918/"
            "07_corrected_analysis/R2-04_motif_indel_aware_v1/"
            "01_corrected_motif/"
            "HRH2_HRH4_regulatory_motif_disruption_scan_indel_aware_v1.tsv"
        ),
        PROJECT / (
            "06_regulatory_genomics/gtex_eqtl_sqtl_lookup/"
            "HRH2_HRH4_GTEx_eQTL_sQTL_lookup_audit_v1_v2_ncbi.tsv"
        ),
        PROJECT / (
            "06_regulatory_genomics/gtex_eqtl_sqtl_lookup/"
            "HRH2_HRH4_GTEx_target_intersections_v1_v2_ncbi.tsv"
        ),
    ]

    missing = [str(x) for x in sources if not x.is_file()]
    if missing:
        raise SystemExit("Missing explicit source(s): " + "; ".join(missing))

    return sources


def extract_annotations():
    annotations = defaultdict(lambda: {
        "gene": "",
        "classes": set(),
        "sources": set(),
        "notes": [],
        "motif_true_notes": [],
        "gtex_true_notes": [],
    })

    true_motif_rows = []
    true_gtex_rows = []
    audit = []

    for path in discover_sources():
        kind = source_kind(path)

        try:
            rows = read_tsv(path)
        except Exception as e:
            audit.append({
                "source": str(path.relative_to(PROJECT)),
                "source_kind": kind,
                "status": f"read_failed:{type(e).__name__}",
                "rows": "",
                "keys_seen": "",
                "true_motif_keys": "",
                "motif_context_keys": "",
                "true_gtex_keys": "",
                "gtex_lookup_keys": "",
            })
            continue

        keys_seen = set()
        true_motif_keys = set()
        motif_context_keys = set()
        true_gtex_keys = set()
        gtex_lookup_keys = set()

        for row in rows:
            key = ""
            for val in row.values():
                key = variant_key_from_text(val)
                if key:
                    break
            if not key:
                continue

            keys_seen.add(key)
            gene = get_gene(row)
            txt = row_text(row).lower()
            note = compact_note(row)

            has_any_change = any([
                boolish(row.get("has_any_motif_change")),
                boolish(row.get("has_variant_overlapping_motif_change")),
                "lost_motifs=" in txt and "lost_motifs=none" not in txt,
                "gained_motifs=" in txt and "gained_motifs=none" not in txt,
                "changed_motif_group_counts=" in txt and "changed_motif_group_counts=none" not in txt,
            ])

            # True motif disruption only if the row explicitly has a motif change.
            if kind in {"motif_candidate_table", "motif_scan"}:
                if has_any_change:
                    annotations[key]["classes"].add("motif_disruption_true")
                    annotations[key]["classes"].add("regulatory_priority")
                    annotations[key]["motif_true_notes"].append(note)
                    true_motif_keys.add(key)
                    true_motif_rows.append({
                        "variant_key": key,
                        "gene": gene,
                        "source": str(path.relative_to(PROJECT)),
                        "motif_status": "true_motif_disruption",
                        "note": note,
                    })
                else:
                    # Keep as context only if table is a motif table but no actual change.
                    annotations[key]["classes"].add("motif_context_only")
                    annotations[key]["classes"].add("regulatory_priority")
                    motif_context_keys.add(key)

            # True GTEx only from intersection/hit tables or gene-level records with pValue and target_variant_key.
            if kind in {"gtex_true_intersection", "gtex_true_hit_table"}:
                annotations[key]["classes"].add("gtex_true_hit")
                annotations[key]["classes"].add("regulatory_priority")
                annotations[key]["gtex_true_notes"].append(note)
                true_gtex_keys.add(key)
                true_gtex_rows.append({
                    "variant_key": key,
                    "gene": gene,
                    "source": str(path.relative_to(PROJECT)),
                    "gtex_status": "true_gtex_intersection_or_hit",
                    "note": note,
                })
            elif kind in {"gtex_lookup_audit", "gtex_gene_level_records"}:
                # Lookup-only does not become expression-supported evidence.
                annotations[key]["classes"].add("gtex_lookup_only")
                gtex_lookup_keys.add(key)

            if gene and not annotations[key]["gene"]:
                annotations[key]["gene"] = gene
            annotations[key]["sources"].add(str(path.relative_to(PROJECT)))
            if note:
                annotations[key]["notes"].append(note)

        audit.append({
            "source": str(path.relative_to(PROJECT)),
            "source_kind": kind,
            "status": "read_ok",
            "rows": len(rows),
            "keys_seen": len(keys_seen),
            "true_motif_keys": len(true_motif_keys),
            "motif_context_keys": len(motif_context_keys),
            "true_gtex_keys": len(true_gtex_keys),
            "gtex_lookup_keys": len(gtex_lookup_keys),
        })

    return annotations, audit, true_motif_rows, true_gtex_rows

def significant_pairwise(pairwise):
    return [
        r for r in pairwise
        if ffloat(r.get("q_value")) is not None
        and ffloat(r.get("q_value")) < 0.05
        and ffloat(r.get("abs_delta_af")) is not None
        and ffloat(r.get("abs_delta_af")) >= 0.10
    ]

def best_pairwise(pairwise):
    best = {}
    for r in pairwise:
        key = r.get("variant_key", "")
        q = ffloat(r.get("q_value"))
        d = ffloat(r.get("abs_delta_af"))
        if not key or q is None or d is None:
            continue
        if key not in best or (q, -d) < (ffloat(best[key].get("q_value")), -ffloat(best[key].get("abs_delta_af"))):
            best[key] = r
    return best

def rebuild_evidence(base_rows, annotations):
    out_by_key = {}

    for r in base_rows:
        key = r["variant_key"]
        nr = {
            "variant_key": key,
            "gene": r.get("gene", ""),
            **{c: "False" for c in EVIDENCE_CLASSES},
            "evidence_sources": r.get("evidence_sources", ""),
            "inclusion_reasons": r.get("inclusion_reasons", ""),
            "evidence_notes": r.get("evidence_notes", ""),
        }

        # Convert old evidence classes into precision classes.
        for c in ["protein", "splice", "biogenesis_fate", "ccre", "ld_context", "docking_variant", "population", "regulatory_priority"]:
            if norm(r.get(c)).lower() == "true" or c in norm(r.get("evidence_classes")).split(";"):
                nr[c] = "True"

        out_by_key[key] = nr

    for key, ann in annotations.items():
        if key not in out_by_key:
            out_by_key[key] = {
                "variant_key": key,
                "gene": ann["gene"],
                **{c: "False" for c in EVIDENCE_CLASSES},
                "evidence_sources": "",
                "inclusion_reasons": "added_by_precision_motif_gtex_repair",
                "evidence_notes": "",
            }

        r = out_by_key[key]

        if ann["gene"] and not r.get("gene"):
            r["gene"] = ann["gene"]

        for c in ann["classes"]:
            if c in EVIDENCE_CLASSES:
                r[c] = "True"

        sources = set(x for x in norm(r.get("evidence_sources")).split(";") if x)
        sources.update(ann["sources"])
        r["evidence_sources"] = ";".join(sorted(sources))

        notes = []
        if norm(r.get("evidence_notes")):
            notes.append(norm(r.get("evidence_notes")))
        notes.extend(ann["motif_true_notes"][:3])
        notes.extend(ann["gtex_true_notes"][:3])
        r["evidence_notes"] = " || ".join(notes[:10])

    out = []
    for r in out_by_key.values():
        classes = [c for c in EVIDENCE_CLASSES if r.get(c) == "True"]
        r["evidence_classes"] = ";".join(classes)
        r["evidence_class_count"] = str(len(classes))
        out.append(r)

    return sorted(out, key=lambda r: (-int(r["evidence_class_count"]), r.get("gene", ""), r["variant_key"]))

def layer_summary(evidence, sig_keys):
    rows = []
    for c in EVIDENCE_CLASSES:
        subset = [r for r in evidence if r.get(c) == "True"]
        genes = Counter(r.get("gene", "") for r in subset)
        rows.append({
            "evidence_layer": c,
            "candidate_variants": len(subset),
            "population_significant_variants_fdr_lt_0_05_delta_ge_0_10": len([r for r in subset if r["variant_key"] in sig_keys]),
            "gene_counts": ";".join(f"{k}:{v}" for k, v in sorted(genes.items()) if k),
            "example_variants": ";".join(r["variant_key"] for r in subset[:12]),
        })
    return rows

def receptor_axis(evidence, pairwise, sig_rows):
    sig_by_gene = defaultdict(set)
    for r in sig_rows:
        sig_by_gene[r.get("gene", "")].add(r["variant_key"])

    out = []
    for gene in GENES:
        subset = [r for r in evidence if r.get("gene") == gene]
        gene_pairs = [r for r in pairwise if r.get("gene") == gene]
        counts = {c: sum(1 for r in subset if r.get(c) == "True") for c in EVIDENCE_CLASSES}
        dominant = sorted(counts.items(), key=lambda x: (-x[1], x[0]))[:6]
        max_delta = max([ffloat(r.get("abs_delta_af")) for r in gene_pairs if ffloat(r.get("abs_delta_af")) is not None] or [0])
        best_q = min([ffloat(r.get("q_value")) for r in gene_pairs if ffloat(r.get("q_value")) is not None] or [1])
        out.append({
            "gene": gene,
            "candidate_variants": len(subset),
            **counts,
            "population_significant_pairwise_tests": len([r for r in sig_rows if r.get("gene") == gene]),
            "population_significant_variants": len(sig_by_gene.get(gene, set())),
            "max_abs_delta_af": fmt(max_delta, 6),
            "best_q_value": fmt(best_q, 12),
            "dominant_evidence_axes": ";".join(f"{k}:{v}" for k, v in dominant if v > 0),
        })
    return out

def validated_table(evidence, pairwise, sig_keys):
    best = best_pairwise(pairwise)
    rows = []

    for e in evidence:
        key = e["variant_key"]
        b = best.get(key, {})
        sig = key in sig_keys

        has_true_gtex = e.get("gtex_true_hit") == "True"
        has_true_motif = e.get("motif_disruption_true") == "True"
        has_ccre = e.get("ccre") == "True"
        has_bio = e.get("biogenesis_fate") == "True"
        has_protein = e.get("protein") == "True"
        has_dock = e.get("docking_variant") == "True"
        has_splice = e.get("splice") == "True"

        functional_count = sum([has_true_gtex, has_true_motif, has_ccre, has_bio, has_protein, has_dock, has_splice])

        if sig and has_true_gtex and has_true_motif:
            grade = "A_population_divergent_expression_motif_supported"
            claim = "strong_population_regulatory_expression_candidate"
        elif sig and has_true_gtex:
            grade = "B_population_divergent_expression_supported"
            claim = "statistically_supported_population_expression"
        elif sig and has_true_motif and has_ccre:
            grade = "C_population_divergent_motif_cCRE_candidate"
            claim = "statistically_supported_population_regulatory_candidate"
        elif sig and functional_count >= 2:
            grade = "D_population_divergent_multilayer_candidate"
            claim = "statistically_supported_population_functional_candidate"
        elif has_true_gtex:
            grade = "E_expression_supported_candidate"
            claim = "expression_supported"
        elif has_true_motif:
            grade = "F_true_motif_disruption_candidate"
            claim = "motif_disruption_hypothesis"
        elif sig:
            grade = "G_population_divergent_candidate"
            claim = "statistically_supported_population_divergence"
        elif functional_count >= 2:
            grade = "H_multilayer_exploratory_candidate"
            claim = "multi_layer_hypothesis"
        else:
            grade = "I_single_layer_or_background_candidate"
            claim = "exploratory_or_background"

        rows.append({
            "variant_key": key,
            "gene": e.get("gene", ""),
            "evidence_classes": e.get("evidence_classes", ""),
            "evidence_class_count": e.get("evidence_class_count", ""),
            "best_population_comparison": b.get("comparison", ""),
            "best_abs_delta_af": b.get("abs_delta_af", ""),
            "best_q_value": b.get("q_value", ""),
            "population_significant_fdr_lt_0_05_delta_ge_0_10": str(sig),
            "interpretive_grade": grade,
            "manuscript_claim_strength": claim,
            "evidence_notes": e.get("evidence_notes", ""),
        })

    return sorted(rows, key=lambda r: (r["interpretive_grade"], -(ffloat(r.get("best_abs_delta_af")) or 0), r["gene"], r["variant_key"]))

def make_figure(layer_rows, receptor_rows, validated):
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as e:
        print(f"WARNING: figure skipped: {type(e).__name__}: {e}")
        return

    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10.5))

    layers = [r["evidence_layer"] for r in layer_rows]
    vals = [int(r["candidate_variants"]) for r in layer_rows]
    sig_vals = [int(r["population_significant_variants_fdr_lt_0_05_delta_ge_0_10"]) for r in layer_rows]
    cols = [COLORS.get(x, "#999999") for x in layers]

    ax = axes[0, 0]
    ax.bar(range(len(layers)), vals, color=cols)
    ax.set_xticks(range(len(layers)))
    ax.set_xticklabels([x.replace("_", "\n") for x in layers], fontsize=7)
    ax.set_ylabel("Candidate count")
    ax.set_title("A. Precision evidence-layer coverage")

    ax = axes[0, 1]
    ax.bar(range(len(layers)), sig_vals, color=cols)
    ax.set_xticks(range(len(layers)))
    ax.set_xticklabels([x.replace("_", "\n") for x in layers], fontsize=7)
    ax.set_ylabel("Population-divergent candidates")
    ax.set_title("B. Population-supported evidence layers")

    ax = axes[1, 0]
    genes = [r["gene"] for r in receptor_rows]
    selected_classes = [
        "protein", "splice", "biogenesis_fate", "motif_disruption_true",
        "ccre", "gtex_true_hit", "ld_context", "docking_variant"
    ]
    bottom = np.zeros(len(genes))
    for c in selected_classes:
        y = np.array([int(r.get(c, 0)) for r in receptor_rows])
        ax.bar(genes, y, bottom=bottom, color=COLORS.get(c, "#999999"), label=c.replace("_", " "))
        bottom += y
    ax.set_ylabel("Candidate count")
    ax.set_title("C. Receptor precision functional axes")
    ax.legend(fontsize=7, frameon=False, ncol=2)

    ax = axes[1, 1]
    grade_counts = Counter(r["interpretive_grade"] for r in validated)
    grades = sorted(grade_counts)
    y = [grade_counts[g] for g in grades]
    ax.barh(range(len(grades)), y, color="#6a51a3")
    ax.set_yticks(range(len(grades)))
    ax.set_yticklabels(grades, fontsize=7)
    ax.set_xlabel("Variant count")
    ax.set_title("D. Precision claim-strength grading")

    fig.suptitle("Precision monitoring of HRH population-functional evidence layers", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_PNG, dpi=900)
    fig.savefig(FIG_SVG)
    plt.close(fig)

def write_md(layer_rows, receptor_rows, validated, true_motif, true_gtex):
    grade_counts = Counter(r["interpretive_grade"] for r in validated)

    true_motif_row = next((r for r in layer_rows if r["evidence_layer"] == "motif_disruption_true"), {})
    context_row = next((r for r in layer_rows if r["evidence_layer"] == "motif_context_only"), {})
    true_gtex_row = next((r for r in layer_rows if r["evidence_layer"] == "gtex_true_hit"), {})
    lookup_row = next((r for r in layer_rows if r["evidence_layer"] == "gtex_lookup_only"), {})

    md = f"""# Precision monitoring of HRH population-functional evidence layers

A precision repair was performed after the broad motif/GTEx repair revealed over-integration of lookup-only GTEx records and motif-context-only records. The precision matrix separates true motif-disruption evidence from motif-context evidence, and true GTEx/eQTL intersections from lookup-only audit/provenance records.

## Corrected regulatory evidence layers

- True motif-disruption candidates: {true_motif_row.get('candidate_variants','0')} candidates; population-supported: {true_motif_row.get('population_significant_variants_fdr_lt_0_05_delta_ge_0_10','0')}; gene counts: {true_motif_row.get('gene_counts','')}.
- Motif-context-only candidates: {context_row.get('candidate_variants','0')} candidates; population-supported: {context_row.get('population_significant_variants_fdr_lt_0_05_delta_ge_0_10','0')}; gene counts: {context_row.get('gene_counts','')}.
- True GTEx/eQTL candidates: {true_gtex_row.get('candidate_variants','0')} candidates; population-supported: {true_gtex_row.get('population_significant_variants_fdr_lt_0_05_delta_ge_0_10','0')}; gene counts: {true_gtex_row.get('gene_counts','')}.
- GTEx lookup-only candidates: {lookup_row.get('candidate_variants','0')} candidates; these should not be described as expression-supported unless also present in true hit/intersection tables.

## Receptor-axis summary

""" + "\n".join(
        f"- {r['gene']}: {r['candidate_variants']} candidates; dominant axes {r['dominant_evidence_axes']}; "
        f"{r['population_significant_variants']} population-supported candidates."
        for r in receptor_rows
    ) + """

## Claim-strength summary

""" + "\n".join(f"- {k}: {v}" for k, v in sorted(grade_counts.items())) + f"""

## Manuscript interpretation

This precision matrix should be used for manuscript claims. Population divergence is statistically supported where pairwise superpopulation tests pass FDR q < 0.05 and absolute allele-frequency difference >= 0.10. True GTEx/eQTL candidates can be described as expression-supported. True motif-disruption candidates can be described as predicted regulatory-mechanism candidates. Motif-context-only and GTEx lookup-only records should remain in supplementary/provenance material and should not be over-interpreted as functional evidence.

## Output files

- Precision evidence matrix: `{OUT_EVIDENCE.relative_to(PROJECT)}`
- Precision layer summary: `{OUT_LAYER.relative_to(PROJECT)}`
- Precision receptor summary: `{OUT_RECEPTOR.relative_to(PROJECT)}`
- Precision validated table: `{OUT_VALIDATED.relative_to(PROJECT)}`
- True motif candidates: `{OUT_TRUE_MOTIF.relative_to(PROJECT)}`
- True GTEx/eQTL candidates: `{OUT_TRUE_GTEX.relative_to(PROJECT)}`
- Precision Figure 11: `{FIG_PNG.relative_to(PROJECT)}`
"""
    RESULTS_MD.write_text(md, encoding="utf-8")

def main():
    if not BASE_EVIDENCE.exists():
        raise SystemExit(f"Missing base evidence matrix: {BASE_EVIDENCE}")
    if not PAIRWISE.exists():
        raise SystemExit(f"Missing pairwise tests: {PAIRWISE}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)

    base = read_tsv(BASE_EVIDENCE)
    pairwise = read_tsv(PAIRWISE)

    annotations, audit, true_motif_rows, true_gtex_rows = extract_annotations()
    evidence = rebuild_evidence(base, annotations)

    # Corrected precision-rebuild invariants.
    gene_counts = Counter(r.get("gene", "") for r in evidence)

    assert len(evidence) == 2011, (
        f"Expected 2011 precision rows, observed {len(evidence)}"
    )
    assert gene_counts["HRH1"] == 678
    assert gene_counts["HRH2"] == 400
    assert gene_counts["HRH3"] == 409
    assert gene_counts["HRH4"] == 524

    true_keys = {r["variant_key"] for r in true_motif_rows}
    true_gene_counts = Counter(
        r.get("gene", "") for r in true_motif_rows
    )

    assert len(true_keys) == 14
    assert true_gene_counts["HRH2"] == 8
    assert true_gene_counts["HRH4"] == 6

    by_key = {r["variant_key"]: r for r in evidence}

    h2 = by_key["5:175709501:TCTGCAGCTGCGTGC:T"]
    h4a = by_key["18:24479080:GC:G"]
    h4b = by_key["18:24456680:A:AT"]

    assert h2["motif_disruption_true"] == "True"

    assert h4a["motif_disruption_true"] == "False"
    assert h4a["motif_context_only"] == "True"
    assert h4a["gtex_lookup_only"] == "True"

    assert h4b["motif_disruption_true"] == "False"
    assert h4b["motif_context_only"] == "True"
    assert h4b["gtex_true_hit"] == "True"
    assert h4b["ld_context"] == "True"
    assert h4b["population"] == "True"


    sig_rows = significant_pairwise(pairwise)
    sig_keys = {r["variant_key"] for r in sig_rows}

    layers = layer_summary(evidence, sig_keys)
    receptors = receptor_axis(evidence, pairwise, sig_rows)
    validated = validated_table(evidence, pairwise, sig_keys)

    evidence_fields = [
        "variant_key", "gene", "evidence_class_count", "evidence_classes",
        *EVIDENCE_CLASSES, "evidence_sources", "inclusion_reasons", "evidence_notes",
    ]
    write_tsv(OUT_EVIDENCE, evidence, evidence_fields)

    layer_fields = [
        "evidence_layer", "candidate_variants",
        "population_significant_variants_fdr_lt_0_05_delta_ge_0_10",
        "gene_counts", "example_variants",
    ]
    write_tsv(OUT_LAYER, layers, layer_fields)
    write_tsv(TABLE56, layers, layer_fields)

    receptor_fields = [
        "gene", "candidate_variants", *EVIDENCE_CLASSES,
        "population_significant_pairwise_tests", "population_significant_variants",
        "max_abs_delta_af", "best_q_value", "dominant_evidence_axes",
    ]
    write_tsv(OUT_RECEPTOR, receptors, receptor_fields)
    write_tsv(TABLE57, receptors, receptor_fields)

    validated_fields = [
        "variant_key", "gene", "evidence_classes", "evidence_class_count",
        "best_population_comparison", "best_abs_delta_af", "best_q_value",
        "population_significant_fdr_lt_0_05_delta_ge_0_10",
        "interpretive_grade", "manuscript_claim_strength", "evidence_notes",
    ]
    write_tsv(OUT_VALIDATED, validated, validated_fields)
    write_tsv(TABLE58, validated[:220], validated_fields)

    write_tsv(OUT_AUDIT, audit, [
        "source", "source_kind", "status", "rows", "keys_seen",
        "true_motif_keys", "motif_context_keys", "true_gtex_keys", "gtex_lookup_keys",
    ])
    write_tsv(OUT_TRUE_MOTIF, true_motif_rows, ["variant_key", "gene", "source", "motif_status", "note"])
    write_tsv(OUT_TRUE_GTEX, true_gtex_rows, ["variant_key", "gene", "source", "gtex_status", "note"])

    combined = []
    for r in true_motif_rows:
        combined.append({
            "evidence_type": "true_motif_disruption",
            "variant_key": r["variant_key"],
            "gene": r["gene"],
            "source": r["source"],
            "note": r["note"],
        })
    for r in true_gtex_rows:
        combined.append({
            "evidence_type": "true_gtex_eqtl",
            "variant_key": r["variant_key"],
            "gene": r["gene"],
            "source": r["source"],
            "note": r["note"],
        })
    write_tsv(TABLE59, combined, ["evidence_type", "variant_key", "gene", "source", "note"])

    make_figure(layers, receptors, validated)
    write_md(layers, receptors, validated, true_motif_rows, true_gtex_rows)

    t_motif = next((r for r in layers if r["evidence_layer"] == "motif_disruption_true"), {})
    m_context = next((r for r in layers if r["evidence_layer"] == "motif_context_only"), {})
    t_gtex = next((r for r in layers if r["evidence_layer"] == "gtex_true_hit"), {})
    g_lookup = next((r for r in layers if r["evidence_layer"] == "gtex_lookup_only"), {})

    print("58C_STATUS=PASS")
    print(f"PRECISION_TOTAL={len(evidence)}")
    print(
        "PRECISION_BY_GENE="
        f"HRH1:{gene_counts['HRH1']};"
        f"HRH2:{gene_counts['HRH2']};"
        f"HRH3:{gene_counts['HRH3']};"
        f"HRH4:{gene_counts['HRH4']}"
    )
    print(f"TRUE_MOTIF_TOTAL={len(true_keys)}")
    print(
        "TRUE_MOTIF_BY_GENE="
        f"HRH2:{true_gene_counts['HRH2']};"
        f"HRH4:{true_gene_counts['HRH4']}"
    )

    print(f"Wrote precision evidence matrix: {OUT_EVIDENCE}")
    print(f"Wrote precision layer summary: {OUT_LAYER}")
    print(f"Wrote precision receptor summary: {OUT_RECEPTOR}")
    print(f"Wrote precision validated table: {OUT_VALIDATED}")
    print(f"Wrote true motif candidates: {OUT_TRUE_MOTIF}")
    print(f"Wrote true GTEx candidates: {OUT_TRUE_GTEX}")
    print(f"Wrote source audit: {OUT_AUDIT}")
    print("Wrote Tables 56-59")
    print(f"Wrote precision Figure 11: {FIG_PNG}, {FIG_SVG}")
    print(f"Wrote Results subsection: {RESULTS_MD}")
    print()
    print(f"Base evidence rows: {len(base)}")
    print(f"Precision evidence rows: {len(evidence)}")
    print(f"True motif-disruption candidates: {t_motif.get('candidate_variants', '0')}")
    print(f"Motif-context-only candidates: {m_context.get('candidate_variants', '0')}")
    print(f"True GTEx/eQTL candidates: {t_gtex.get('candidate_variants', '0')}")
    print(f"GTEx lookup-only candidates: {g_lookup.get('candidate_variants', '0')}")
    print("Precision receptor-axis summary:")
    for r in receptors:
        print(
            f"  {r['gene']}: candidates={r['candidate_variants']}; "
            f"true_motif={r['motif_disruption_true']}; true_gtex={r['gtex_true_hit']}; "
            f"sig_variants={r['population_significant_variants']}; axes={r['dominant_evidence_axes']}"
        )

if __name__ == "__main__":
    main()
