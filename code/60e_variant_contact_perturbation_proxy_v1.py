#!/usr/bin/env python3

from pathlib import Path
import csv
import re
from collections import Counter, defaultdict

PROJECT = Path("/mnt/d/histamine_receptor_pharmacogenomics")

PANEL = PROJECT / "10_docking/mutant_redocking_v1/hrh_mutant_redocking_candidate_panel_v1.tsv"
CONTACTS = PROJECT / "10_docking/contact_analysis_v1/hrh_representative_ligand_receptor_contact_residue_summary_v1.tsv"
VARIANT_CONTEXT = PROJECT / "10_docking/variant_docking_integration/final_interpretation/hrh_variant_docking_final_interpretive_table_v1.tsv"
PRECISION = PROJECT / "11_results/monitoring_evaluation_v4_precision/hrh_precision_validated_vs_exploratory_evidence_table_v1.tsv"

OUTDIR = PROJECT / "10_docking/mutant_redocking_v1/contact_perturbation_proxy"
FIGDIR = PROJECT / "12_figures/mutant_redocking_v1"

OUT_PROXY = OUTDIR / "hrh_variant_contact_perturbation_proxy_v1.tsv"
OUT_SUMMARY = OUTDIR / "hrh_variant_contact_perturbation_proxy_summary_v1.tsv"
OUT_FAILED_NOTE = OUTDIR / "hrh_mutant_redocking_technical_failure_note_v1.txt"

TABLE64C = PROJECT / "13_tables/manuscript_ready/Table_64C_variant_contact_perturbation_proxy_v1.tsv"
RESULTS_MD = PROJECT / "14_manuscript/results_sections/results_variant_contact_perturbation_proxy_v1.md"

FIG_PNG = FIGDIR / "Figure_13_variant_contact_perturbation_proxy_v1_900dpi.png"
FIG_SVG = FIGDIR / "Figure_13_variant_contact_perturbation_proxy_v1.svg"

AA = {
    "ALA": {"size":  89, "hydro":  1.8, "charge": 0,  "polar": 0, "aromatic": 0, "sulfur": 0},
    "ARG": {"size": 174, "hydro": -4.5, "charge": 1,  "polar": 1, "aromatic": 0, "sulfur": 0},
    "ASN": {"size": 132, "hydro": -3.5, "charge": 0,  "polar": 1, "aromatic": 0, "sulfur": 0},
    "ASP": {"size": 133, "hydro": -3.5, "charge": -1, "polar": 1, "aromatic": 0, "sulfur": 0},
    "CYS": {"size": 121, "hydro":  2.5, "charge": 0,  "polar": 1, "aromatic": 0, "sulfur": 1},
    "GLN": {"size": 146, "hydro": -3.5, "charge": 0,  "polar": 1, "aromatic": 0, "sulfur": 0},
    "GLU": {"size": 147, "hydro": -3.5, "charge": -1, "polar": 1, "aromatic": 0, "sulfur": 0},
    "GLY": {"size":  75, "hydro": -0.4, "charge": 0,  "polar": 0, "aromatic": 0, "sulfur": 0},
    "HIS": {"size": 155, "hydro": -3.2, "charge": 1,  "polar": 1, "aromatic": 1, "sulfur": 0},
    "ILE": {"size": 131, "hydro":  4.5, "charge": 0,  "polar": 0, "aromatic": 0, "sulfur": 0},
    "LEU": {"size": 131, "hydro":  3.8, "charge": 0,  "polar": 0, "aromatic": 0, "sulfur": 0},
    "LYS": {"size": 146, "hydro": -3.9, "charge": 1,  "polar": 1, "aromatic": 0, "sulfur": 0},
    "MET": {"size": 149, "hydro":  1.9, "charge": 0,  "polar": 0, "aromatic": 0, "sulfur": 1},
    "PHE": {"size": 165, "hydro":  2.8, "charge": 0,  "polar": 0, "aromatic": 1, "sulfur": 0},
    "PRO": {"size": 115, "hydro": -1.6, "charge": 0,  "polar": 0, "aromatic": 0, "sulfur": 0},
    "SER": {"size": 105, "hydro": -0.8, "charge": 0,  "polar": 1, "aromatic": 0, "sulfur": 0},
    "THR": {"size": 119, "hydro": -0.7, "charge": 0,  "polar": 1, "aromatic": 0, "sulfur": 0},
    "TRP": {"size": 204, "hydro": -0.9, "charge": 0,  "polar": 0, "aromatic": 1, "sulfur": 0},
    "TYR": {"size": 181, "hydro": -1.3, "charge": 0,  "polar": 1, "aromatic": 1, "sulfur": 0},
    "VAL": {"size": 117, "hydro":  4.2, "charge": 0,  "polar": 0, "aromatic": 0, "sulfur": 0},
}

GENE_COLORS = {
    "HRH1": "#1b9e77",
    "HRH2": "#d95f02",
    "HRH3": "#7570b3",
    "HRH4": "#e7298a",
}

def read_tsv(path):
    if not path.exists():
        return []
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

def row_text(row):
    return " ".join(norm(v) for v in row.values())

def load_contacts():
    rows = read_tsv(CONTACTS)
    by_gene_pos = defaultdict(list)
    for r in rows:
        gene = norm(r.get("gene"))
        pos = norm(r.get("receptor_resnum"))
        if gene and pos:
            by_gene_pos[(gene, pos)].append(r)
    return by_gene_pos

def load_variant_context():
    rows = read_tsv(VARIANT_CONTEXT)
    by_key_hgvsp = {}
    for r in rows:
        txt = row_text(r)
        key = ""
        m = re.search(r"\b(?:chr)?([0-9XYM]+):([0-9]+):([ACGTN]+):([ACGTN]+)\b", txt, flags=re.I)
        if m:
            chrom, pos, ref, alt = m.groups()
            key = f"{chrom.replace('chr','')}:{pos}:{ref.upper()}:{alt.upper()}"
        hgvsp = ""
        hm = re.search(r"p\.[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}", txt)
        if hm:
            hgvsp = hm.group(0)
        if key:
            by_key_hgvsp[key] = r
        if hgvsp:
            by_key_hgvsp[hgvsp] = r
    return by_key_hgvsp

def load_precision():
    rows = read_tsv(PRECISION)
    out = {}
    for r in rows:
        key = norm(r.get("variant_key"))
        if key:
            out[key] = r
    return out

def aa_change_metrics(ref, alt):
    if ref not in AA or alt not in AA:
        return {}, ["unsupported_amino_acid"]

    r = AA[ref]
    a = AA[alt]

    metrics = {
        "delta_size": a["size"] - r["size"],
        "abs_delta_size": abs(a["size"] - r["size"]),
        "delta_hydrophobicity": a["hydro"] - r["hydro"],
        "abs_delta_hydrophobicity": abs(a["hydro"] - r["hydro"]),
        "charge_change": int(a["charge"] != r["charge"]),
        "polarity_change": int(a["polar"] != r["polar"]),
        "aromatic_change": int(a["aromatic"] != r["aromatic"]),
        "aromatic_loss": int(r["aromatic"] == 1 and a["aromatic"] == 0),
        "aromatic_gain": int(r["aromatic"] == 0 and a["aromatic"] == 1),
        "sulfur_change": int(a["sulfur"] != r["sulfur"]),
        "cysteine_gain": int(alt == "CYS" and ref != "CYS"),
        "cysteine_loss": int(ref == "CYS" and alt != "CYS"),
    }

    labels = []
    if metrics["aromatic_loss"]:
        labels.append("aromatic_loss")
    if metrics["aromatic_gain"]:
        labels.append("aromatic_gain")
    if metrics["charge_change"]:
        labels.append("charge_change")
    if metrics["polarity_change"]:
        labels.append("polarity_change")
    if metrics["abs_delta_size"] >= 40:
        labels.append("large_side_chain_size_change")
    if metrics["abs_delta_hydrophobicity"] >= 3.0:
        labels.append("large_hydrophobicity_change")
    if metrics["cysteine_gain"]:
        labels.append("cysteine_gain")
    if metrics["cysteine_loss"]:
        labels.append("cysteine_loss")
    if metrics["sulfur_change"]:
        labels.append("sulfur_status_change")

    if not labels:
        labels.append("conservative_or_subtle_side_chain_change")

    return metrics, labels

def classify_score(score):
    if score >= 80:
        return "very_high_contact_perturbation_priority"
    if score >= 60:
        return "high_contact_perturbation_priority"
    if score >= 40:
        return "moderate_contact_perturbation_priority"
    return "contextual_or_low_contact_perturbation_priority"

def score_variant(row, contact_rows, context_row, precision_row, metrics, labels):
    score = 0
    reasons = []

    exact = len(contact_rows) > 0
    txt = (row_text(row) + " " + row_text(context_row)).lower()

    if exact:
        score += 45
        reasons.append("variant_residue_is_representative_ligand_contact_residue")
        max_contacts = max(int(float(c.get("contact_count_4A", 0) or 0)) for c in contact_rows)
        if max_contacts >= 10:
            score += 15
            reasons.append("high_contact_count_at_variant_residue")
        elif max_contacts >= 5:
            score += 8
            reasons.append("moderate_contact_count_at_variant_residue")

    if "ligand-contact-zone" in txt:
        score += 25
        reasons.append("module51_ligand_contact_zone")
    elif "pocket-proximal" in txt:
        score += 20
        reasons.append("module51_pocket_proximal")
    elif "near-pocket" in txt:
        score += 15
        reasons.append("module51_near_pocket")
    elif "structurally-near" in txt:
        score += 10
        reasons.append("module51_structurally_near")

    if metrics.get("aromatic_loss"):
        score += 20
        reasons.append("aromatic_side_chain_loss")
    if metrics.get("aromatic_gain"):
        score += 12
        reasons.append("aromatic_side_chain_gain")
    if metrics.get("charge_change"):
        score += 15
        reasons.append("charge_state_change")
    if metrics.get("polarity_change"):
        score += 8
        reasons.append("polarity_change")
    if metrics.get("abs_delta_size", 0) >= 40:
        score += 10
        reasons.append("large_side_chain_size_change")
    if metrics.get("abs_delta_hydrophobicity", 0) >= 3.0:
        score += 8
        reasons.append("large_hydrophobicity_change")
    if metrics.get("cysteine_gain"):
        score += 10
        reasons.append("cysteine_gain_possible_new_thiol_constraint")
    if metrics.get("cysteine_loss"):
        score += 10
        reasons.append("cysteine_loss_possible_disulfide_or_thiol_change")

    if precision_row:
        grade = norm(precision_row.get("interpretive_grade"))
        if grade.startswith(("A_", "B_", "C_", "D_")):
            score += 8
            reasons.append("precision_population_functional_support")

    if "tier_1_high_priority" in txt:
        score += 5
        reasons.append("tier1_candidate")

    return score, reasons

def make_interpretation(row, contact_rows, labels, score_class):
    gene = row["gene"]
    hgvsp = row["hgvsp"]
    ligand = row["representative_ligand"]

    if contact_rows:
        contact = contact_rows[0]
        residue_id = contact.get("residue_id", "")
        contact_count = contact.get("contact_count_4A", "")
        dominant = contact.get("dominant_contact_type", "")
        return (
            f"{gene} {hgvsp} directly overlaps representative {ligand} contact residue {residue_id} "
            f"with {contact_count} contacts dominated by {dominant}; side-chain change features include "
            f"{';'.join(labels)}. Classified as {score_class}."
        )

    return (
        f"{gene} {hgvsp} does not exactly overlap the representative ligand contact residues but remains a selected "
        f"structural-context candidate; side-chain change features include {';'.join(labels)}. "
        f"Classified as {score_class}."
    )

def make_figure(rows):
    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"WARNING: matplotlib unavailable; figure skipped: {type(e).__name__}: {e}")
        return

    if not rows:
        return

    rows = sorted(rows, key=lambda r: (-int(r["contact_perturbation_score"]), r["gene"], r["hgvsp"]))
    labels = [f"{r['gene']}\n{r['hgvsp']}" for r in rows]
    scores = [int(r["contact_perturbation_score"]) for r in rows]
    colors = [GENE_COLORS.get(r["gene"], "#666666") for r in rows]

    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(max(10, len(rows) * 0.72), 6.5))
    ax.bar(range(len(rows)), scores, color=colors)
    ax.axhline(40, color="#777777", linewidth=0.8, linestyle="--")
    ax.axhline(60, color="#777777", linewidth=0.8, linestyle=":")
    ax.axhline(80, color="#222222", linewidth=0.9, linestyle=":")
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Contact-perturbation priority score")
    ax.set_title("Selected HRH coding variants ranked by ligand-contact perturbation proxy")
    fig.tight_layout()
    fig.savefig(FIG_PNG, dpi=900)
    fig.savefig(FIG_SVG)
    plt.close(fig)

def write_results_md(rows, summary):
    top = sorted(rows, key=lambda r: -int(r["contact_perturbation_score"]))[:8]
    top_lines = [
        f"- {r['gene']} {r['hgvsp']} ({r['representative_ligand']}): score {r['contact_perturbation_score']}, "
        f"{r['contact_perturbation_class']}; {r['short_interpretation']}"
        for r in top
    ]

    summary_lines = [
        f"- {r['gene']}: {r['selected_variants']} selected variants; {r['exact_contact_variants']} exact contact variants; "
        f"{r['very_high_or_high_priority']} high/very-high proxy-priority variants."
        for r in summary
    ]

    md = f"""# Variant contact-perturbation proxy analysis

Mutant-versus-wild-type redocking was attempted but failed technically because the local PyMOL installation could not generate point-mutant receptor models reproducibly. Therefore, the failed redocking outputs should not be interpreted biologically. As a conservative replacement, selected coding variants were ranked using a ligand-contact perturbation proxy that integrates representative ligand-contact residues, final variant-pocket context and amino-acid physicochemical change.

## Gene-level summary

{chr(10).join(summary_lines)}

## Highest-priority contact-perturbation candidates

{chr(10).join(top_lines)}

## Interpretation

This analysis does not estimate mutant binding affinity. Instead, it identifies variants whose residue position and side-chain change plausibly perturb a ligand-contact environment or nearby pocket architecture. Exact representative contact-residue variants are the strongest candidates for future mutant modelling, molecular dynamics or experimental pharmacology. Structurally near variants are retained as indirect pocket-context candidates and should not be claimed to alter binding without further validation.

## Output files

- Proxy table: `10_docking/mutant_redocking_v1/contact_perturbation_proxy/hrh_variant_contact_perturbation_proxy_v1.tsv`
- Proxy summary: `10_docking/mutant_redocking_v1/contact_perturbation_proxy/hrh_variant_contact_perturbation_proxy_summary_v1.tsv`
- Technical failure note: `10_docking/mutant_redocking_v1/contact_perturbation_proxy/hrh_mutant_redocking_technical_failure_note_v1.txt`
- Figure 13 replacement: `12_figures/mutant_redocking_v1/Figure_13_variant_contact_perturbation_proxy_v1_900dpi.png`
- Manuscript table: `13_tables/manuscript_ready/Table_64C_variant_contact_perturbation_proxy_v1.tsv`
"""
    RESULTS_MD.write_text(md, encoding="utf-8")

def main():
    if not PANEL.exists():
        raise SystemExit(f"Missing selected mutant-redocking panel: {PANEL}")
    if not CONTACTS.exists():
        raise SystemExit(f"Missing contact residue table: {CONTACTS}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)

    panel = read_tsv(PANEL)
    contacts = load_contacts()
    context = load_variant_context()
    precision = load_precision()

    rows = []

    for row in panel:
        gene = norm(row.get("gene"))
        hgvsp = norm(row.get("hgvsp"))
        pos = norm(row.get("protein_position"))
        key = norm(row.get("variant_key"))
        ref = norm(row.get("reference_residue_3letter")).upper()
        alt = norm(row.get("alternate_residue_3letter")).upper()

        contact_rows = contacts.get((gene, pos), [])
        context_row = context.get(key, context.get(hgvsp, {}))
        precision_row = precision.get(key, {})

        metrics, labels = aa_change_metrics(ref, alt)
        score, reasons = score_variant(row, contact_rows, context_row, precision_row, metrics, labels)
        score_class = classify_score(score)

        if contact_rows:
            best_contact = sorted(
                contact_rows,
                key=lambda c: (-int(float(c.get("contact_count_4A", 0) or 0)), float(c.get("min_distance_angstrom", 99) or 99))
            )[0]
        else:
            best_contact = {}

        out = {
            "gene": gene,
            "variant_key": key,
            "hgvsp": hgvsp,
            "protein_position": pos,
            "reference_residue_3letter": ref,
            "alternate_residue_3letter": alt,
            "representative_pdb_id": norm(row.get("representative_pdb_id")),
            "representative_ligand": norm(row.get("representative_ligand")),
            "exact_representative_contact_residue": str(bool(contact_rows)),
            "contact_residue_id": best_contact.get("residue_id", ""),
            "contact_count_4A": best_contact.get("contact_count_4A", ""),
            "min_contact_distance_angstrom": best_contact.get("min_distance_angstrom", ""),
            "dominant_contact_type": best_contact.get("dominant_contact_type", ""),
            "delta_size": metrics.get("delta_size", ""),
            "abs_delta_size": metrics.get("abs_delta_size", ""),
            "delta_hydrophobicity": metrics.get("delta_hydrophobicity", ""),
            "abs_delta_hydrophobicity": metrics.get("abs_delta_hydrophobicity", ""),
            "charge_change": metrics.get("charge_change", ""),
            "polarity_change": metrics.get("polarity_change", ""),
            "aromatic_change": metrics.get("aromatic_change", ""),
            "aromatic_loss": metrics.get("aromatic_loss", ""),
            "aromatic_gain": metrics.get("aromatic_gain", ""),
            "cysteine_gain": metrics.get("cysteine_gain", ""),
            "cysteine_loss": metrics.get("cysteine_loss", ""),
            "side_chain_change_labels": ";".join(labels),
            "contact_perturbation_score": score,
            "contact_perturbation_class": score_class,
            "score_reasons": ";".join(reasons),
            "short_interpretation": "",
            "manuscript_claim": "contact_perturbation_candidate_not_mutant_binding_affinity",
        }

        out["short_interpretation"] = make_interpretation(out, contact_rows, labels, score_class)
        rows.append(out)

    rows = sorted(rows, key=lambda r: (-int(r["contact_perturbation_score"]), r["gene"], r["hgvsp"]))

    fields = [
        "gene", "variant_key", "hgvsp", "protein_position",
        "reference_residue_3letter", "alternate_residue_3letter",
        "representative_pdb_id", "representative_ligand",
        "exact_representative_contact_residue", "contact_residue_id",
        "contact_count_4A", "min_contact_distance_angstrom",
        "dominant_contact_type", "delta_size", "abs_delta_size",
        "delta_hydrophobicity", "abs_delta_hydrophobicity",
        "charge_change", "polarity_change", "aromatic_change",
        "aromatic_loss", "aromatic_gain", "cysteine_gain", "cysteine_loss",
        "side_chain_change_labels", "contact_perturbation_score",
        "contact_perturbation_class", "score_reasons",
        "short_interpretation", "manuscript_claim",
    ]

    write_tsv(OUT_PROXY, rows, fields)
    write_tsv(TABLE64C, rows, fields)

    summary = []
    for gene in ["HRH1", "HRH2", "HRH3", "HRH4"]:
        sub = [r for r in rows if r["gene"] == gene]
        high = [r for r in sub if r["contact_perturbation_class"] in {"very_high_contact_perturbation_priority", "high_contact_perturbation_priority"}]
        exact = [r for r in sub if r["exact_representative_contact_residue"] == "True"]
        summary.append({
            "gene": gene,
            "selected_variants": len(sub),
            "exact_contact_variants": len(exact),
            "very_high_or_high_priority": len(high),
            "top_variants": ";".join(f"{r['hgvsp']}({r['contact_perturbation_score']})" for r in sub[:8]),
        })

    write_tsv(OUT_SUMMARY, summary, ["gene", "selected_variants", "exact_contact_variants", "very_high_or_high_priority", "top_variants"])

    OUT_FAILED_NOTE.write_text(
        "Mutant-versus-wild-type redocking was attempted in Modules 60B and 60D but failed at PyMOL point-mutant model generation. "
        "The failure is technical and should not be interpreted as absence of variant effects. "
        "The contact-perturbation proxy table replaces the failed redocking figure/table for manuscript interpretation unless a reproducible modelling tool such as FoldX, Modeller or a working PyMOL API environment is later installed.\n",
        encoding="utf-8",
    )

    make_figure(rows)
    write_results_md(rows, summary)

    print(f"Wrote contact-perturbation proxy table: {OUT_PROXY}")
    print(f"Wrote contact-perturbation summary: {OUT_SUMMARY}")
    print(f"Wrote technical failure note: {OUT_FAILED_NOTE}")
    print(f"Wrote Table 64C: {TABLE64C}")
    print(f"Wrote Figure 13 replacement: {FIG_PNG}, {FIG_SVG}")
    print(f"Wrote Results subsection: {RESULTS_MD}")
    print()
    print(f"Selected variants scored: {len(rows)}")
    print(f"Exact representative contact variants: {sum(1 for r in rows if r['exact_representative_contact_residue'] == 'True')}")
    print("Top variants:")
    for r in rows[:10]:
        print(f"  {r['gene']} {r['hgvsp']} {r['representative_ligand']}: score={r['contact_perturbation_score']} class={r['contact_perturbation_class']} exact_contact={r['exact_representative_contact_residue']}")
    print("Gene-level summary:")
    for r in summary:
        print(f"  {r['gene']}: selected={r['selected_variants']}; exact_contact={r['exact_contact_variants']}; high_or_very_high={r['very_high_or_high_priority']}; top={r['top_variants']}")

if __name__ == "__main__":
    main()
