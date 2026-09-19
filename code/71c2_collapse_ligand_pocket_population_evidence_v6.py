#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import datetime
from collections import defaultdict

PROJECT = Path("/mnt/d/histamine_receptor_pharmacogenomics")
OUTDIR = PROJECT / "11_results/monitoring_evaluation_v6_global_functional_diversity"

IN_MATRIX = OUTDIR / "HRH_ligand_pocket_population_diversity_matrix_v6.tsv"
IN_PAIRWISE_SIG = OUTDIR / "HRH_ligand_pocket_population_significant_pairwise_tests_v6.tsv"

OUT_COLLAPSED = OUTDIR / "HRH_ligand_pocket_population_collapsed_variant_ligand_evidence_v6.tsv"
OUT_REDOCK_PANEL = OUTDIR / "HRH_ligand_pocket_population_redocking_candidate_panel_v6.tsv"
OUT_SUMMARY = OUTDIR / "HRH_ligand_pocket_population_collapsed_summary_v6.tsv"
OUT_REPORT = OUTDIR / "HRH_ligand_pocket_population_collapsed_report_v6.md"
OUT_CHECKLIST = OUTDIR / "HRH_ligand_pocket_population_collapsed_checklist_v6.tsv"

CONTACT_RANK = {
    "direct_or_exact_contact": 5,
    "near_or_pocket_proximal": 4,
    "structural_context": 3,
    "ligand_or_contact_context": 2,
    "annotation_inferred_structural_context": 1,
    "unclassified_contact_source": 0,
}

def read_tsv(path):
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

def as_float(x, default=0.0):
    try:
        if x in {None, ""}:
            return default
        return float(x)
    except Exception:
        return default

def as_int(x, default=0):
    try:
        if x in {None, ""}:
            return default
        return int(float(x))
    except Exception:
        return default

def source_reliability(row):
    src = row.get("source_path", "").lower()
    txt = row.get("source_row_excerpt", "").lower()
    score = 0

    if "contact_perturbation_proxy" in src:
        score += 40
    if "mutant_redocking_candidate_panel" in src:
        score += 35
    if "final_interpretive" in src or "top_structural_candidates" in src:
        score += 30
    if "variant_docking_integration" in src:
        score += 25
    if "manuscript_ready" in src:
        score += 15
    if "population_matrix_annotation" in src:
        score -= 20

    # Penalize unresolved/unmapped rows when choosing best mechanistic evidence.
    if "unresolved" in txt or "unmapped" in txt or "not resolved" in txt:
        score -= 35
    if "distal-from-docking-pocket" in txt or "distal from ligand" in txt:
        score -= 15
    if "near-pocket" in txt or "pocket-proximal" in txt:
        score += 20
    if "exact_contact" in txt or "exact representative contact" in txt or "direct contact" in txt:
        score += 40
    if "structurally-near" in txt:
        score += 10

    return score

def extract_protein_change(text):
    if not text:
        return ""
    m = re.search(r"\b(p\.[A-Za-z]{3}\d+[A-Za-z]{3})\b", text)
    if m:
        return m.group(1)
    m = re.search(r"\bENSP[^\s|:]*:?(p\.[A-Za-z]{3}\d+[A-Za-z]{3})\b", text)
    if m:
        return m.group(1)
    return ""

def classify_actionability(best_class, best_excerpt, sig_count, abs_delta):
    text = (best_excerpt or "").lower()
    if best_class == "direct_or_exact_contact":
        return "A_direct_contact_population_context" if sig_count > 0 else "A_direct_contact_low_population_signal"
    if best_class == "near_or_pocket_proximal":
        return "B_near_pocket_population_supported" if sig_count > 0 else "B_near_pocket_low_population_signal"
    if best_class == "structural_context":
        if sig_count > 0:
            return "C_structural_context_population_supported"
        return "C_structural_context_low_population_signal"
    if "unresolved" in text or "unmapped" in text or "distal" in text:
        return "D_distal_or_unresolved_do_not_redock_first"
    if sig_count > 0 and abs_delta >= 0.10:
        return "D_population_signal_but_weak_structural_contact"
    return "E_background_or_annotation_only"

def choose_best(rows):
    def key(r):
        return (
            CONTACT_RANK.get(r.get("contact_class", ""), 0),
            source_reliability(r),
            as_float(r.get("contact_priority_score")),
            as_int(r.get("significant_pairwise_count_fdr_0_05_delta_0_10")),
            as_float(r.get("top_pairwise_abs_delta_af")),
        )
    return sorted(rows, key=key, reverse=True)[0]

def main():
    if not IN_MATRIX.exists():
        raise SystemExit(f"Missing input matrix: {IN_MATRIX}")

    rows = read_tsv(IN_MATRIX)
    sig_rows = read_tsv(IN_PAIRWISE_SIG) if IN_PAIRWISE_SIG.exists() else []

    grouped = defaultdict(list)
    for r in rows:
        key = (r.get("variant_key", ""), r.get("gene", ""), r.get("ligand", ""))
        if key[0] and key[1] and key[2]:
            grouped[key].append(r)

    sig_by_variant = defaultdict(list)
    for r in sig_rows:
        sig_by_variant[r.get("variant_key", "")].append(r)

    collapsed = []
    for (variant_key, gene, ligand), group in grouped.items():
        best = choose_best(group)
        classes = sorted(set(r.get("contact_class", "") for r in group))
        sources = sorted(set(r.get("source_path", "") for r in group))
        protein_changes = sorted(set(extract_protein_change(r.get("source_row_excerpt", "")) for r in group if extract_protein_change(r.get("source_row_excerpt", ""))))

        sig_count = as_int(best.get("significant_pairwise_count_fdr_0_05_delta_0_10"))
        abs_delta = as_float(best.get("top_pairwise_abs_delta_af"))
        best_class = best.get("contact_class", "")
        action_class = classify_actionability(best_class, best.get("source_row_excerpt", ""), sig_count, abs_delta)

        sig_details = []
        for s in sig_by_variant.get(variant_key, []):
            sig_details.append(
                f"{s.get('comparison')}:{s.get('pop1_af')}vs{s.get('pop2_af')};delta={s.get('abs_delta_af')};q={s.get('fdr_q')}"
            )

        collapsed.append({
            "variant_key": variant_key,
            "gene": gene,
            "ligand": ligand,
            "protein_change_detected": ";".join(protein_changes),
            "best_contact_class": best_class,
            "best_contact_rank": str(CONTACT_RANK.get(best_class, 0)),
            "best_contact_priority_score": best.get("contact_priority_score", ""),
            "best_source_path": best.get("source_path", ""),
            "evidence_row_count": str(len(group)),
            "evidence_contact_classes": ";".join(classes),
            "evidence_source_count": str(len(sources)),
            "ann_consequence": best.get("ann_consequence", ""),
            "ann_tier": best.get("ann_tier", ""),
            "ann_topology_or_domain": best.get("ann_topology_or_domain", ""),
            "ann_axis_or_evidence": best.get("ann_axis_or_evidence", ""),
            "AFR_AF": best.get("AFR_AF", ""),
            "AMR_AF": best.get("AMR_AF", ""),
            "EAS_AF": best.get("EAS_AF", ""),
            "EUR_AF": best.get("EUR_AF", ""),
            "SAS_AF": best.get("SAS_AF", ""),
            "max_population": best.get("max_population", ""),
            "min_population": best.get("min_population", ""),
            "max_minus_min_AF": best.get("max_minus_min_AF", ""),
            "top_pairwise_comparison": best.get("top_pairwise_comparison", ""),
            "top_pairwise_abs_delta_af": best.get("top_pairwise_abs_delta_af", ""),
            "top_pairwise_fdr_q": best.get("top_pairwise_fdr_q", ""),
            "significant_pairwise_count_fdr_0_05_delta_0_10": best.get("significant_pairwise_count_fdr_0_05_delta_0_10", ""),
            "population_differentiation_class": best.get("population_differentiation_class", ""),
            "redocking_priority_class": action_class,
            "significant_pairwise_details": " | ".join(sig_details),
            "best_source_excerpt": best.get("source_row_excerpt", "")[:1000],
        })

    collapsed = sorted(
        collapsed,
        key=lambda r: (
            r["redocking_priority_class"][:1],
            as_int(r["significant_pairwise_count_fdr_0_05_delta_0_10"]),
            as_float(r["top_pairwise_abs_delta_af"]),
            as_int(r["best_contact_rank"]),
            as_float(r["best_contact_priority_score"]),
        ),
        reverse=False,
    )

    fields = [
        "variant_key", "gene", "ligand", "protein_change_detected",
        "best_contact_class", "best_contact_rank", "best_contact_priority_score",
        "best_source_path", "evidence_row_count", "evidence_contact_classes", "evidence_source_count",
        "ann_consequence", "ann_tier", "ann_topology_or_domain", "ann_axis_or_evidence",
        "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF",
        "max_population", "min_population", "max_minus_min_AF",
        "top_pairwise_comparison", "top_pairwise_abs_delta_af", "top_pairwise_fdr_q",
        "significant_pairwise_count_fdr_0_05_delta_0_10",
        "population_differentiation_class",
        "redocking_priority_class",
        "significant_pairwise_details",
        "best_source_excerpt",
    ]
    write_tsv(OUT_COLLAPSED, collapsed, fields)

    # Redocking panel: A/B/C candidates, prioritizing population support, but keep direct-contact low-population candidates too.
    redock = []
    for r in collapsed:
        cls = r["redocking_priority_class"]
        if cls.startswith("A_") or cls.startswith("B_") or cls.startswith("C_structural_context_population_supported"):
            redock.append(r)

    # Limit to a focused panel for next rescue module.
    redock = sorted(
        redock,
        key=lambda r: (
            1 if r["gene"] == "HRH4" else 0,
            1 if r["ligand"].lower() in {"toreforant", "adriforant"} else 0,
            as_int(r["significant_pairwise_count_fdr_0_05_delta_0_10"]),
            as_float(r["top_pairwise_abs_delta_af"]),
            as_int(r["best_contact_rank"]),
        ),
        reverse=True,
    )[:20]
    write_tsv(OUT_REDOCK_PANEL, redock, fields)

    # Summary by gene-ligand and redocking priority class.
    summary_map = defaultdict(lambda: defaultdict(int))
    for r in collapsed:
        key = (r["gene"], r["ligand"])
        summary_map[key]["collapsed_variant_ligand_pairs"] += 1
        if as_int(r["significant_pairwise_count_fdr_0_05_delta_0_10"]) > 0:
            summary_map[key]["population_significant_pairs"] += 1
        summary_map[key][r["redocking_priority_class"]] += 1
        summary_map[key]["max_abs_delta_af"] = max(summary_map[key]["max_abs_delta_af"], as_float(r["top_pairwise_abs_delta_af"]))

    summary_rows = []
    for (gene, ligand), d in sorted(summary_map.items()):
        summary_rows.append({
            "gene": gene,
            "ligand": ligand,
            "collapsed_variant_ligand_pairs": str(d["collapsed_variant_ligand_pairs"]),
            "population_significant_pairs": str(d["population_significant_pairs"]),
            "A_direct_contact_population_context": str(d["A_direct_contact_population_context"]),
            "A_direct_contact_low_population_signal": str(d["A_direct_contact_low_population_signal"]),
            "B_near_pocket_population_supported": str(d["B_near_pocket_population_supported"]),
            "B_near_pocket_low_population_signal": str(d["B_near_pocket_low_population_signal"]),
            "C_structural_context_population_supported": str(d["C_structural_context_population_supported"]),
            "C_structural_context_low_population_signal": str(d["C_structural_context_low_population_signal"]),
            "D_distal_or_unresolved_do_not_redock_first": str(d["D_distal_or_unresolved_do_not_redock_first"]),
            "D_population_signal_but_weak_structural_contact": str(d["D_population_signal_but_weak_structural_contact"]),
            "E_background_or_annotation_only": str(d["E_background_or_annotation_only"]),
            "max_abs_delta_af": f"{d['max_abs_delta_af']:.6g}",
        })

    summary_fields = [
        "gene", "ligand", "collapsed_variant_ligand_pairs", "population_significant_pairs",
        "A_direct_contact_population_context", "A_direct_contact_low_population_signal",
        "B_near_pocket_population_supported", "B_near_pocket_low_population_signal",
        "C_structural_context_population_supported", "C_structural_context_low_population_signal",
        "D_distal_or_unresolved_do_not_redock_first", "D_population_signal_but_weak_structural_contact",
        "E_background_or_annotation_only", "max_abs_delta_af",
    ]
    write_tsv(OUT_SUMMARY, summary_rows, summary_fields)

    checklist = [
        {"check": "input_matrix_found", "status": "PASS" if IN_MATRIX.exists() else "REVIEW", "observed": str(IN_MATRIX)},
        {"check": "raw_rows_read", "status": "PASS" if len(rows) > 0 else "REVIEW", "observed": str(len(rows))},
        {"check": "collapsed_variant_ligand_pairs", "status": "PASS" if len(collapsed) > 0 else "REVIEW", "observed": str(len(collapsed))},
        {"check": "redocking_panel_created", "status": "PASS" if len(redock) > 0 else "REVIEW", "observed": str(len(redock))},
        {"check": "direct_or_near_pocket_candidates_in_panel", "status": "PASS" if any(r["redocking_priority_class"].startswith(("A_", "B_")) for r in redock) else "REVIEW", "observed": str(sum(1 for r in redock if r["redocking_priority_class"].startswith(("A_", "B_"))))},
        {"check": "hrh4_toreforant_present_in_panel", "status": "PASS" if any(r["gene"] == "HRH4" and r["ligand"].lower() == "toreforant" for r in redock) else "REVIEW", "observed": str(sum(1 for r in redock if r["gene"] == "HRH4" and r["ligand"].lower() == "toreforant"))},
        {"check": "collapsed_summary_created", "status": "PASS" if OUT_SUMMARY.exists() else "REVIEW", "observed": str(OUT_SUMMARY)},
    ]
    write_tsv(OUT_CHECKLIST, checklist, ["check", "status", "observed"])

    report = []
    report.append("# HRH ligand-pocket collapsed population evidence report v6")
    report.append("")
    report.append(f"Generated: {datetime.datetime.now().isoformat(timespec='seconds')}")
    report.append("")
    report.append("## Headline results")
    report.append("")
    report.append(f"- Raw evidence rows: {len(rows)}")
    report.append(f"- Collapsed variant-ligand pairs: {len(collapsed)}")
    report.append(f"- Focused redocking candidate panel: {len(redock)}")
    report.append("")
    report.append("## Why this collapse was necessary")
    report.append("")
    report.append("The repaired 71C matrix intentionally captured evidence from many docking/contact/manuscript source files, so the same variant-ligand pair appeared repeatedly. This module collapses the table to one best evidence row per variant-ligand pair and separates direct/near-pocket candidates from distal, unresolved or annotation-inferred structural-context rows.")
    report.append("")
    report.append("## Top redocking candidates")
    report.append("")
    report.append("| Gene | Ligand | Variant | Protein change | Class | AF max-min | Top comparison | Significant pairwise count |")
    report.append("|---|---|---|---|---|---:|---|---:|")
    for r in redock[:20]:
        report.append(
            f"| {r['gene']} | {r['ligand']} | {r['variant_key']} | {r['protein_change_detected']} | "
            f"{r['redocking_priority_class']} | {r['top_pairwise_abs_delta_af']} | {r['top_pairwise_comparison']} | "
            f"{r['significant_pairwise_count_fdr_0_05_delta_0_10']} |"
        )
    report.append("")
    report.append("## Main outputs")
    report.append("")
    report.append(f"- Collapsed evidence table: `{OUT_COLLAPSED}`")
    report.append(f"- Redocking candidate panel: `{OUT_REDOCK_PANEL}`")
    report.append(f"- Summary: `{OUT_SUMMARY}`")
    report.append(f"- Checklist: `{OUT_CHECKLIST}`")
    OUT_REPORT.write_text("\n".join(report) + "\n", encoding="utf-8")

    print("\n".join(report))

if __name__ == "__main__":
    main()
