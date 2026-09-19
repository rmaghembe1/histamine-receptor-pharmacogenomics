#!/usr/bin/env python3

from pathlib import Path
import csv
import math
import datetime
from collections import defaultdict

PROJECT = Path("/mnt/d/histamine_receptor_pharmacogenomics")
OUTDIR = PROJECT / "11_results/monitoring_evaluation_v6_global_functional_diversity"
FIGDIR = PROJECT / "12_figures/monitoring_evaluation_v6_global_functional_diversity"
TABLEDIR = PROJECT / "13_tables/manuscript_ready"

VARIANT_ANNOT = OUTDIR / "HRH_receptor_fate_CMA_proteostasis_variant_annotation_v6.tsv"
REDOCK_RESULTS = OUTDIR / "HRH_mutant_vs_wt_redocking_results_pdbfixer_v6.tsv"
REDOCK_FAILURES = OUTDIR / "HRH_mutant_vs_wt_redocking_failures_pdbfixer_v6.tsv"
CLEAN_PANEL = OUTDIR / "HRH_ligand_pocket_population_redocking_candidate_panel_clean_v6.tsv"

OUT_LIGAND_AWARE = OUTDIR / "HRH_receptor_fate_CMA_proteostasis_ligand_aware_v6.tsv"
OUT_PRIORITY = OUTDIR / "HRH_receptor_fate_CMA_proteostasis_ligand_aware_prioritized_v6.tsv"
OUT_SUMMARY = OUTDIR / "HRH_receptor_fate_CMA_proteostasis_ligand_aware_summary_v6.tsv"
OUT_CHECKLIST = OUTDIR / "HRH_receptor_fate_CMA_proteostasis_ligand_aware_checklist_v6.tsv"
OUT_REPORT = OUTDIR / "HRH_receptor_fate_CMA_proteostasis_ligand_aware_report_v6.md"
OUT_TABLE = TABLEDIR / "Table_71E2_ligand_aware_receptor_fate_CMA_proteostasis_v6.tsv"

OUT_FIG_PNG = FIGDIR / "Figure_6B_ligand_aware_receptor_fate_CMA_proteostasis_priority_v6_900dpi.png"
OUT_FIG_SVG = FIGDIR / "Figure_6B_ligand_aware_receptor_fate_CMA_proteostasis_priority_v6.svg"

def read_tsv(path):
    if not path.exists():
        return []
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

def pf(x):
    try:
        if x is None or str(x).strip() == "":
            return None
        return float(str(x).strip())
    except Exception:
        return None

def variant_key(r):
    return (r.get("gene", ""), r.get("variant_key", ""), r.get("protein_change", ""))

def pair_key(r):
    return (
        str(r.get("rank", "")),
        r.get("gene", ""),
        r.get("ligand", ""),
        r.get("variant_key", ""),
        r.get("protein_change", ""),
    )

def pair_key_no_rank(r):
    return (
        r.get("gene", ""),
        r.get("ligand", ""),
        r.get("variant_key", ""),
        r.get("protein_change", ""),
    )

def score_row(r):
    score = 0.0

    motif = r.get("motif_consequence", "")
    if "N_glycosylation_lost" in motif or "N_glycosylation_created" in motif:
        score += 4.0
    if "CMA_like_motif_lost" in motif or "CMA_like_motif_created" in motif:
        score += 3.0

    effects = r.get("protein_fate_effect_classes", "")
    if "cysteine_disulfide_or_palmitoylation_proxy" in effects:
        score += 2.5
    if "lysine_ubiquitination" in effects:
        score += 2.0
    if "phosphorylation_site_proxy" in effects:
        score += 1.5
    if "proline_helix_kink_proxy" in effects:
        score += 1.5
    if "glycine_flexibility_proxy" in effects:
        score += 1.0
    if "local_charge_rewiring_proxy" in effects:
        score += 1.0

    delta = pf(r.get("delta_mutant_minus_wt_kcal_mol"))
    if delta is not None:
        ad = abs(delta)
        if ad >= 1.0:
            score += 4.0
        elif ad >= 0.5:
            score += 3.0
        elif ad >= 0.2:
            score += 2.0
        elif ad > 0:
            score += 0.5

    daf = pf(r.get("top_pairwise_abs_delta_af"))
    q = pf(r.get("top_pairwise_fdr_q"))
    if daf is not None:
        if daf >= 0.10:
            score += 3.0
        elif daf >= 0.05:
            score += 2.0
        elif daf >= 0.01:
            score += 1.0
    if q is not None:
        if q < 0.05:
            score += 1.5
        elif q < 0.10:
            score += 0.5

    if r.get("in_71c3_clean_redocking_panel") == "yes":
        score += 1.0
    if r.get("successful_71d3_redocking") == "yes":
        score += 1.5

    return score

def main():
    timestamp = datetime.datetime.now().isoformat(timespec="seconds")

    variant_annot = {variant_key(r): r for r in read_tsv(VARIANT_ANNOT)}
    redock = read_tsv(REDOCK_RESULTS)
    failures = read_tsv(REDOCK_FAILURES)
    clean_rows = read_tsv(CLEAN_PANEL)

    clean_by_pair = {}
    clean_by_pair_no_rank = {}
    clean_by_variant_ligand = defaultdict(list)
    for r in clean_rows:
        clean_by_pair[pair_key(r)] = r
        clean_by_pair_no_rank[pair_key_no_rank(r)] = r
        clean_by_variant_ligand[(r.get("gene", ""), r.get("variant_key", ""), r.get("protein_change", ""), r.get("ligand", ""))].append(r)

    rows = []

    def build_row(base, status):
        vk = variant_key(base)
        va = variant_annot.get(vk, {})
        cp = clean_by_pair.get(pair_key(base), {})
        if not cp:
            cp = clean_by_pair_no_rank.get(pair_key_no_rank(base), {})

        out = {}

        # Pair identity first.
        for col in ["rank", "gene", "ligand", "variant_key", "protein_change"]:
            out[col] = base.get(col, "")

        # Sequence-aware receptor-fate annotation from 71E after sequence recovery.
        for col in [
            "ref1", "alt1", "position", "sequence_found", "sequence_length",
            "sequence_source", "sequence_header", "reference_match", "position_zone",
            "wt_local_window", "mut_local_window", "protein_fate_effect_classes",
            "motif_consequence", "ngly_lost", "ngly_created", "cma_like_lost",
            "cma_like_created",
        ]:
            out[col] = va.get(col, "")

        # Population and panel metadata from clean panel where available.
        for col in [
            "redocking_priority_class", "top_pairwise_comparison",
            "top_pairwise_abs_delta_af", "top_pairwise_fdr_q",
        ]:
            out[col] = base.get(col, "") or cp.get(col, "") or va.get(col, "")

        out["in_71c3_clean_redocking_panel"] = "yes" if cp else base.get("in_71c3_clean_redocking_panel", "no")
        out["successful_71d3_redocking"] = "yes" if status == "success" else "no"
        out["failed_71d3_redocking"] = "yes" if status == "failure" else "no"
        out["71d3_failure_reason"] = base.get("reason", "") if status == "failure" else ""

        # Docking evidence remains ligand-pair specific. Never collapse across ligands.
        for col in [
            "wt_vina_affinity_kcal_mol", "mutant_vina_affinity_kcal_mol",
            "delta_mutant_minus_wt_kcal_mol", "predicted_direction",
            "receptor_pdb", "docked_pose_source", "pose_discovery_method",
            "ligand_input_pdbqt", "wt_receptor_pdbqt", "mutant_receptor_pdb",
            "mutant_receptor_pdbqt", "wt_vina_log", "mutant_vina_log", "run_dir",
        ]:
            out[col] = base.get(col, "")

        out["receptor_fate_ligand_aware_priority_score"] = f"{score_row(out):.2f}"
        rows.append(out)

    for r in redock:
        build_row(r, "success")

    for r in failures:
        build_row(r, "failure")

    rows.sort(
        key=lambda r: (
            -float(r.get("receptor_fate_ligand_aware_priority_score", "0") or 0),
            r.get("gene", ""),
            r.get("ligand", ""),
            r.get("protein_change", ""),
            r.get("rank", ""),
        )
    )

    fields = [
        "rank", "gene", "ligand", "variant_key", "protein_change",
        "ref1", "alt1", "position", "sequence_found", "sequence_length",
        "sequence_source", "sequence_header", "reference_match", "position_zone",
        "wt_local_window", "mut_local_window",
        "protein_fate_effect_classes", "motif_consequence",
        "ngly_lost", "ngly_created", "cma_like_lost", "cma_like_created",
        "redocking_priority_class", "top_pairwise_comparison",
        "top_pairwise_abs_delta_af", "top_pairwise_fdr_q",
        "in_71c3_clean_redocking_panel", "successful_71d3_redocking",
        "failed_71d3_redocking", "71d3_failure_reason",
        "wt_vina_affinity_kcal_mol", "mutant_vina_affinity_kcal_mol",
        "delta_mutant_minus_wt_kcal_mol", "predicted_direction",
        "receptor_fate_ligand_aware_priority_score",
        "receptor_pdb", "docked_pose_source", "pose_discovery_method",
        "ligand_input_pdbqt", "wt_receptor_pdbqt", "mutant_receptor_pdb",
        "mutant_receptor_pdbqt", "wt_vina_log", "mutant_vina_log", "run_dir",
    ]

    write_tsv(OUT_LIGAND_AWARE, rows, fields)

    prioritized = [r for r in rows if pf(r.get("receptor_fate_ligand_aware_priority_score")) and pf(r.get("receptor_fate_ligand_aware_priority_score")) > 0]
    write_tsv(OUT_PRIORITY, prioritized, fields)
    write_tsv(OUT_TABLE, prioritized, fields)

    # Summary.
    summary = []
    summary.append({"metric": "ligand_variant_rows_total", "value": str(len(rows))})
    summary.append({"metric": "successful_71d3_redocking_pairs", "value": str(sum(1 for r in rows if r.get("successful_71d3_redocking") == "yes"))})
    summary.append({"metric": "failed_71d3_redocking_pairs", "value": str(sum(1 for r in rows if r.get("failed_71d3_redocking") == "yes"))})
    summary.append({"metric": "sequence_found_rows", "value": str(sum(1 for r in rows if r.get("sequence_found") == "yes"))})
    summary.append({"metric": "reference_match_yes_rows", "value": str(sum(1 for r in rows if r.get("reference_match") == "yes"))})
    summary.append({"metric": "motif_change_rows", "value": str(sum(1 for r in rows if r.get("motif_consequence") != "no_local_ngly_or_CMA_like_change_detected"))})

    by_lig = defaultdict(list)
    for r in rows:
        if r.get("successful_71d3_redocking") == "yes":
            by_lig[r.get("ligand", "")].append(pf(r.get("delta_mutant_minus_wt_kcal_mol")))
    for lig, vals in sorted(by_lig.items()):
        vals = [v for v in vals if v is not None]
        if vals:
            summary.append({"metric": f"ligand:{lig}:successful_pairs", "value": str(len(vals))})
            summary.append({"metric": f"ligand:{lig}:mean_delta", "value": f"{sum(vals)/len(vals):.3f}"})
            summary.append({"metric": f"ligand:{lig}:max_abs_delta", "value": f"{max(abs(v) for v in vals):.3f}"})

    write_tsv(OUT_SUMMARY, summary, ["metric", "value"])

    # Figure.
    fig_status = "not_created"
    try:
        import matplotlib.pyplot as plt
        top = prioritized[:25]
        if top:
            labels = [f"{r['gene']} {r['ligand']} {r['protein_change']}" for r in top]
            vals = [float(r["receptor_fate_ligand_aware_priority_score"]) for r in top]
            height = max(5, 0.38 * len(top) + 1.5)
            plt.figure(figsize=(10, height))
            y = list(range(len(top)))
            plt.barh(y, vals)
            plt.yticks(y, labels, fontsize=7)
            plt.xlabel("Ligand-aware receptor-fate/proteostasis priority score")
            plt.gca().invert_yaxis()
            plt.tight_layout()
            plt.savefig(OUT_FIG_PNG, dpi=900)
            plt.savefig(OUT_FIG_SVG)
            plt.close()
            fig_status = "created"
        else:
            fig_status = "no_prioritized_rows"
    except Exception as e:
        fig_status = str(e)

    checklist = [
        {"check": "sequence_aware_variant_annotation_found", "status": "PASS" if VARIANT_ANNOT.exists() else "FAIL", "observed": str(VARIANT_ANNOT)},
        {"check": "71d3_redocking_results_found", "status": "PASS" if REDOCK_RESULTS.exists() else "FAIL", "observed": str(REDOCK_RESULTS)},
        {"check": "71d3_failure_table_found", "status": "PASS" if REDOCK_FAILURES.exists() else "REVIEW", "observed": str(REDOCK_FAILURES)},
        {"check": "clean_panel_found", "status": "PASS" if CLEAN_PANEL.exists() else "REVIEW", "observed": str(CLEAN_PANEL)},
        {"check": "ligand_variant_rows_total", "status": "PASS" if len(rows) >= 16 else "REVIEW", "observed": str(len(rows))},
        {"check": "successful_redocking_pairs", "status": "PASS" if sum(1 for r in rows if r.get("successful_71d3_redocking") == "yes") == 15 else "REVIEW", "observed": str(sum(1 for r in rows if r.get("successful_71d3_redocking") == "yes"))},
        {"check": "sequence_found_for_rows", "status": "PASS" if all(r.get("sequence_found") == "yes" for r in rows if r.get("gene") == "HRH4") else "REVIEW", "observed": str(sum(1 for r in rows if r.get("sequence_found") == "yes"))},
        {"check": "figure_created", "status": "PASS" if fig_status == "created" else "REVIEW", "observed": fig_status},
    ]
    write_tsv(OUT_CHECKLIST, checklist, ["check", "status", "observed"])

    report = []
    report.append("# Module 71E2 ligand-aware receptor-fate/CMA/proteostasis report v6")
    report.append("")
    report.append(f"Generated: {timestamp}")
    report.append("")
    report.append("## Headline")
    report.append("")
    report.append(f"- Ligand-variant rows: {len(rows)}")
    report.append(f"- Successful 71D3 redocking pairs represented: {sum(1 for r in rows if r.get('successful_71d3_redocking') == 'yes')}")
    report.append(f"- Failed 71D3 rows represented: {sum(1 for r in rows if r.get('failed_71d3_redocking') == 'yes')}")
    report.append(f"- Sequence-found rows: {sum(1 for r in rows if r.get('sequence_found') == 'yes')}")
    report.append(f"- Motif-change rows: {sum(1 for r in rows if r.get('motif_consequence') != 'no_local_ngly_or_CMA_like_change_detected')}")
    report.append(f"- Figure status: {fig_status}")
    report.append("")
    report.append("## Main outputs")
    report.append("")
    for label, path in [
        ("Ligand-aware annotation", OUT_LIGAND_AWARE),
        ("Prioritized ligand-aware table", OUT_PRIORITY),
        ("Summary", OUT_SUMMARY),
        ("Checklist", OUT_CHECKLIST),
        ("Manuscript table", OUT_TABLE),
        ("Figure PNG", OUT_FIG_PNG),
        ("Figure SVG", OUT_FIG_SVG),
    ]:
        report.append(f"- {label}: `{path}`")
    report.append("")
    report.append("## Top ligand-aware prioritized rows")
    report.append("")
    report.append("| Rank | Gene | Ligand | Protein change | Score | WT | Mutant | Delta | Motif consequence | Effect classes |")
    report.append("|---|---|---|---|---:|---:|---:|---:|---|---|")
    for r in prioritized[:25]:
        report.append(
            f"| {r.get('rank','')} | {r.get('gene','')} | {r.get('ligand','')} | {r.get('protein_change','')} | "
            f"{r.get('receptor_fate_ligand_aware_priority_score','')} | {r.get('wt_vina_affinity_kcal_mol','')} | "
            f"{r.get('mutant_vina_affinity_kcal_mol','')} | {r.get('delta_mutant_minus_wt_kcal_mol','')} | "
            f"{r.get('motif_consequence','')} | {r.get('protein_fate_effect_classes','')} |"
        )

    OUT_REPORT.write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))

if __name__ == "__main__":
    main()
