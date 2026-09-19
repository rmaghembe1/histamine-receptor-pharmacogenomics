#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import datetime
from collections import defaultdict

PROJECT = Path("/mnt/d/histamine_receptor_pharmacogenomics")
OUTDIR = PROJECT / "11_results/monitoring_evaluation_v6_global_functional_diversity"
OUTDIR.mkdir(parents=True, exist_ok=True)

IN_COLLAPSED = OUTDIR / "HRH_ligand_pocket_population_collapsed_variant_ligand_evidence_v6.tsv"
IN_REDOCK = OUTDIR / "HRH_ligand_pocket_population_redocking_candidate_panel_v6.tsv"

OUT_CLEAN_COLLAPSED = OUTDIR / "HRH_ligand_pocket_population_collapsed_variant_ligand_evidence_clean_v6.tsv"
OUT_CLEAN_REDOCK = OUTDIR / "HRH_ligand_pocket_population_redocking_candidate_panel_clean_v6.tsv"
OUT_MAPPING_AUDIT = OUTDIR / "HRH_variant_to_protein_mapping_audit_for_redocking_v6.tsv"
OUT_SUMMARY = OUTDIR / "HRH_redocking_candidate_panel_clean_summary_v6.tsv"
OUT_REPORT = OUTDIR / "HRH_redocking_candidate_panel_clean_report_v6.md"
OUT_CHECKLIST = OUTDIR / "HRH_redocking_candidate_panel_clean_checklist_v6.tsv"

RELIABLE_MAPPING_SOURCES = [
    PROJECT / "11_results/candidate_prioritization/hrh_protein_altering_candidates_canonical_safe_gpcr_topology_v2_ncbi.tsv",
    PROJECT / "11_results/candidate_prioritization/hrh_exact_candidate_prioritization_v2_ncbi.tsv",
    PROJECT / "10_docking/variant_docking_integration/final_interpretation/hrh_variant_docking_final_interpretive_table_v1.tsv",
    PROJECT / "13_tables/manuscript_ready/Table_41_HRH_variant_docking_final_interpretive_table_v1.tsv",
    PROJECT / "13_tables/manuscript_ready/Table_39_HRH_sequence_mapped_variant_docking_context_v2.tsv",
    PROJECT / "13_tables/manuscript_ready/Table_37_HRH_variant_docking_pocket_context_v1.tsv",
]

CONTACT_CLASS_RANK = {
    "direct_or_exact_contact": 5,
    "near_or_pocket_proximal": 4,
    "structural_context": 3,
    "ligand_or_contact_context": 2,
    "annotation_inferred_structural_context": 1,
}

def read_table(path):
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        delim = "\t" if sample.count("\t") >= sample.count(",") else ","
        return list(csv.DictReader(f, delimiter=delim))

def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

def norm_chrom(x):
    return str(x).replace("chr", "").replace("NC_000003.12", "3").replace("NC_000005.10", "5").replace("NC_000020.11", "20").replace("NC_000018.10", "18")

def cell_text(row):
    vals = []
    for k, v in row.items():
        if k is None:
            continue
        if v is not None:
            vals.append(str(v))
    return " | ".join(vals)

def extract_variant_keys(text):
    keys = set()
    if not text:
        return keys
    for m in re.finditer(r"\b(?:chr)?([0-9XYM]+)[:_](\d+)[:_]([ACGT]+)[:_]([ACGT]+)\b", str(text)):
        c, pos, ref, alt = m.groups()
        keys.add(f"{norm_chrom(c)}:{pos}:{ref}:{alt}")
    return keys

def extract_protein_changes(text):
    if not text:
        return set()
    return set(re.findall(r"\b(p\.[A-Za-z]{3}\d+[A-Za-z]{3})\b", str(text)))

def as_float(x, default=0.0):
    try:
        if x in {None, ""}:
            return default
        return float(str(x).replace(",", ""))
    except Exception:
        return default

def as_int(x, default=0):
    try:
        if x in {None, ""}:
            return default
        return int(float(str(x).replace(",", "")))
    except Exception:
        return default

def source_priority(path):
    s = str(path)
    if "hrh_protein_altering_candidates_canonical_safe" in s:
        return 100
    if "hrh_exact_candidate_prioritization" in s:
        return 90
    if "final_interpretive" in s:
        return 80
    if "sequence_mapped" in s:
        return 70
    return 50

def build_variant_to_protein_map():
    mapping = {}
    audit = []

    for source in RELIABLE_MAPPING_SOURCES:
        if not source.exists():
            audit.append({
                "source_path": str(source),
                "status": "MISSING",
                "rows_read": "",
                "mapped_rows": "",
                "note": "",
            })
            continue

        try:
            rows = read_table(source)
        except Exception as e:
            audit.append({
                "source_path": str(source),
                "status": "READ_FAIL",
                "rows_read": "",
                "mapped_rows": "",
                "note": str(e),
            })
            continue

        mapped = 0
        prio = source_priority(source)

        for row in rows:
            text = cell_text(row)
            keys = extract_variant_keys(text)
            proteins = extract_protein_changes(text)

            # Reliable source rows should usually contain one variant and one protein change.
            # Skip multi-protein rows to avoid contaminating the mapping.
            if len(keys) != 1 or len(proteins) != 1:
                continue

            vk = next(iter(keys))
            protein = next(iter(proteins))

            old = mapping.get(vk)
            if old is None or prio > old["priority"]:
                mapping[vk] = {
                    "variant_key": vk,
                    "canonical_protein_change": protein,
                    "mapping_source": str(source.relative_to(PROJECT)),
                    "priority": prio,
                }
            mapped += 1

        audit.append({
            "source_path": str(source.relative_to(PROJECT)),
            "status": "READ_OK",
            "rows_read": str(len(rows)),
            "mapped_rows": str(mapped),
            "note": "",
        })

    return mapping, audit

def clean_row(row, mapping):
    vk = row.get("variant_key", "")
    raw_proteins = sorted(extract_protein_changes(row.get("protein_change_detected", "")) | extract_protein_changes(row.get("best_source_excerpt", "")))
    raw_protein_field = row.get("protein_change_detected", "")

    map_entry = mapping.get(vk)
    if map_entry:
        canonical = map_entry["canonical_protein_change"]
        map_source = map_entry["mapping_source"]
        map_status = "EXACT_FROM_RELIABLE_SOURCE"
    else:
        only_raw = sorted(extract_protein_changes(raw_protein_field))
        if len(only_raw) == 1:
            canonical = only_raw[0]
            map_source = "existing_single_protein_field"
            map_status = "SINGLE_EXISTING_FIELD_USED"
        else:
            canonical = ""
            map_source = ""
            map_status = "NO_UNAMBIGUOUS_PROTEIN_MAPPING"

    new = dict(row)
    new["raw_protein_change_detected"] = raw_protein_field
    new["raw_protein_change_count_in_text"] = str(len(raw_proteins))
    new["canonical_protein_change_clean"] = canonical
    new["protein_mapping_status"] = map_status
    new["protein_mapping_source"] = map_source

    cls = row.get("redocking_priority_class", "")
    contact = row.get("best_contact_class", "")
    consequence = row.get("ann_consequence", "").lower()
    sig_count = as_int(row.get("significant_pairwise_count_fdr_0_05_delta_0_10"))
    abs_delta = as_float(row.get("top_pairwise_abs_delta_af"))

    is_missense = "missense" in consequence
    has_clean_protein = bool(canonical)

    if cls.startswith("D_") or cls.startswith("E_"):
        eligible = "NO"
        reason = "distal_unresolved_or_background"
    elif not is_missense:
        eligible = "NO"
        reason = "not_missense"
    elif not has_clean_protein:
        eligible = "NO"
        reason = "no_unambiguous_protein_mapping"
    elif contact in {"direct_or_exact_contact", "near_or_pocket_proximal", "structural_context"}:
        eligible = "YES"
        reason = "clean_protein_and_structural_contact_context"
    else:
        eligible = "NO"
        reason = "weak_contact_class"

    new["redocking_eligible_clean"] = eligible
    new["redocking_exclusion_reason_clean"] = reason

    rank = 0
    if eligible == "YES":
        rank += CONTACT_CLASS_RANK.get(contact, 0) * 100
        rank += min(sig_count, 5) * 50
        rank += min(abs_delta, 0.30) * 100
        if row.get("gene") == "HRH4":
            rank += 30
        if row.get("ligand", "").lower() in {"toreforant", "adriforant"}:
            rank += 20
        if cls.startswith("A_"):
            rank += 40
        if "population_supported" in cls:
            rank += 40
    new["redocking_clean_rank_score"] = f"{rank:.4g}"

    return new

def main():
    if not IN_COLLAPSED.exists():
        raise SystemExit(f"Missing input collapsed table: {IN_COLLAPSED}")
    if not IN_REDOCK.exists():
        raise SystemExit(f"Missing input redocking panel: {IN_REDOCK}")

    collapsed = read_table(IN_COLLAPSED)
    redock_raw = read_table(IN_REDOCK)

    mapping, mapping_audit = build_variant_to_protein_map()

    cleaned_collapsed = [clean_row(r, mapping) for r in collapsed]
    cleaned_redock_all = [clean_row(r, mapping) for r in redock_raw]

    # Keep only clean eligible rows for the next redocking rescue module.
    clean_panel = [r for r in cleaned_redock_all if r.get("redocking_eligible_clean") == "YES"]

    # Deduplicate same variant-ligand-protein.
    dedup = {}
    for r in clean_panel:
        key = (r.get("variant_key", ""), r.get("gene", ""), r.get("ligand", ""), r.get("canonical_protein_change_clean", ""))
        old = dedup.get(key)
        if old is None or as_float(r.get("redocking_clean_rank_score")) > as_float(old.get("redocking_clean_rank_score")):
            dedup[key] = r
    clean_panel = list(dedup.values())

    clean_panel = sorted(
        clean_panel,
        key=lambda r: (
            as_float(r.get("redocking_clean_rank_score")),
            as_int(r.get("significant_pairwise_count_fdr_0_05_delta_0_10")),
            as_float(r.get("top_pairwise_abs_delta_af")),
        ),
        reverse=True,
    )

    # Keep a focused but sufficiently broad panel for 71D.
    clean_panel = clean_panel[:16]

    base_fields = list(cleaned_collapsed[0].keys()) if cleaned_collapsed else []
    extra_fields = [
        "raw_protein_change_detected",
        "raw_protein_change_count_in_text",
        "canonical_protein_change_clean",
        "protein_mapping_status",
        "protein_mapping_source",
        "redocking_eligible_clean",
        "redocking_exclusion_reason_clean",
        "redocking_clean_rank_score",
    ]
    fields = []
    for f in base_fields + extra_fields:
        if f not in fields:
            fields.append(f)

    write_tsv(OUT_CLEAN_COLLAPSED, cleaned_collapsed, fields)
    write_tsv(OUT_CLEAN_REDOCK, clean_panel, fields)

    audit_fields = ["source_path", "status", "rows_read", "mapped_rows", "note"]
    write_tsv(OUT_MAPPING_AUDIT, mapping_audit, audit_fields)

    # Summary.
    summary_map = defaultdict(lambda: defaultdict(int))
    for r in cleaned_collapsed:
        key = (r.get("gene", ""), r.get("ligand", ""))
        summary_map[key]["collapsed_pairs"] += 1
        if r.get("redocking_eligible_clean") == "YES":
            summary_map[key]["redocking_eligible_clean"] += 1
        if r.get("protein_mapping_status") == "EXACT_FROM_RELIABLE_SOURCE":
            summary_map[key]["exact_protein_mapped"] += 1
        if as_int(r.get("significant_pairwise_count_fdr_0_05_delta_0_10")) > 0:
            summary_map[key]["population_significant_pairs"] += 1

    summary_rows = []
    for (gene, ligand), d in sorted(summary_map.items()):
        summary_rows.append({
            "gene": gene,
            "ligand": ligand,
            "collapsed_pairs": str(d["collapsed_pairs"]),
            "exact_protein_mapped": str(d["exact_protein_mapped"]),
            "redocking_eligible_clean": str(d["redocking_eligible_clean"]),
            "population_significant_pairs": str(d["population_significant_pairs"]),
        })
    write_tsv(
        OUT_SUMMARY,
        summary_rows,
        ["gene", "ligand", "collapsed_pairs", "exact_protein_mapped", "redocking_eligible_clean", "population_significant_pairs"],
    )

    # Targeted checks for known contamination.
    ala_rows = [r for r in clean_panel + cleaned_collapsed if r.get("variant_key") == "18:24476802:C:T"]
    tyr_rows = [r for r in clean_panel + cleaned_collapsed if r.get("variant_key") == "18:24477345:A:G"]

    ala_clean = any(r.get("canonical_protein_change_clean") == "p.Ala138Val" for r in ala_rows)
    tyr_clean = any(r.get("canonical_protein_change_clean") == "p.Tyr319Cys" for r in tyr_rows)

    checklist = [
        {"check": "collapsed_input_found", "status": "PASS" if IN_COLLAPSED.exists() else "REVIEW", "observed": str(IN_COLLAPSED)},
        {"check": "raw_redocking_input_found", "status": "PASS" if IN_REDOCK.exists() else "REVIEW", "observed": str(IN_REDOCK)},
        {"check": "variant_protein_mappings_built", "status": "PASS" if len(mapping) > 0 else "REVIEW", "observed": str(len(mapping))},
        {"check": "clean_collapsed_rows_written", "status": "PASS" if len(cleaned_collapsed) > 0 else "REVIEW", "observed": str(len(cleaned_collapsed))},
        {"check": "clean_redocking_panel_written", "status": "PASS" if len(clean_panel) > 0 else "REVIEW", "observed": str(len(clean_panel))},
        {"check": "ala138_mapping_clean", "status": "PASS" if ala_clean else "REVIEW", "observed": "18:24476802:C:T should be p.Ala138Val"},
        {"check": "tyr319_mapping_clean", "status": "PASS" if tyr_clean else "REVIEW", "observed": "18:24477345:A:G should be p.Tyr319Cys"},
        {"check": "no_semicolon_in_clean_panel_protein_change", "status": "PASS" if not any(";" in r.get("canonical_protein_change_clean", "") for r in clean_panel) else "REVIEW", "observed": str(sum(1 for r in clean_panel if ';' in r.get('canonical_protein_change_clean','')))},
        {"check": "hrh4_toreforant_in_clean_panel", "status": "PASS" if any(r.get("gene") == "HRH4" and r.get("ligand", "").lower() == "toreforant" for r in clean_panel) else "REVIEW", "observed": str(sum(1 for r in clean_panel if r.get('gene') == 'HRH4' and r.get('ligand','').lower() == 'toreforant'))},
    ]
    write_tsv(OUT_CHECKLIST, checklist, ["check", "status", "observed"])

    report = []
    report.append("# HRH cleaned redocking candidate panel report v6")
    report.append("")
    report.append(f"Generated: {datetime.datetime.now().isoformat(timespec='seconds')}")
    report.append("")
    report.append("## Headline results")
    report.append("")
    report.append(f"- Collapsed variant-ligand pairs cleaned: {len(cleaned_collapsed)}")
    report.append(f"- Reliable variant-to-protein mappings built: {len(mapping)}")
    report.append(f"- Clean redocking panel size: {len(clean_panel)}")
    report.append("")
    report.append("## Why this cleaning was necessary")
    report.append("")
    report.append("The 71C2 collapse correctly reduced repeated evidence rows, but some source excerpts contained multiple protein-change names in a single text field. This module rebuilds the protein-change assignment from reliable one-variant/one-protein mapping sources and prevents multi-protein strings from entering the redocking rescue panel.")
    report.append("")
    report.append("## Clean 71D redocking panel")
    report.append("")
    report.append("| Rank | Gene | Ligand | Variant | Clean protein change | Class | AF delta | Top comparison | Eligibility reason |")
    report.append("|---:|---|---|---|---|---|---:|---|---|")
    for i, r in enumerate(clean_panel, start=1):
        report.append(
            f"| {i} | {r.get('gene','')} | {r.get('ligand','')} | {r.get('variant_key','')} | "
            f"{r.get('canonical_protein_change_clean','')} | {r.get('redocking_priority_class','')} | "
            f"{r.get('top_pairwise_abs_delta_af','')} | {r.get('top_pairwise_comparison','')} | "
            f"{r.get('redocking_exclusion_reason_clean','')} |"
        )
    report.append("")
    report.append("## Main outputs")
    report.append("")
    report.append(f"- Clean collapsed table: `{OUT_CLEAN_COLLAPSED}`")
    report.append(f"- Clean redocking panel: `{OUT_CLEAN_REDOCK}`")
    report.append(f"- Mapping audit: `{OUT_MAPPING_AUDIT}`")
    report.append(f"- Summary: `{OUT_SUMMARY}`")
    report.append(f"- Checklist: `{OUT_CHECKLIST}`")
    OUT_REPORT.write_text("\n".join(report) + "\n", encoding="utf-8")

    print("\n".join(report))

if __name__ == "__main__":
    main()
