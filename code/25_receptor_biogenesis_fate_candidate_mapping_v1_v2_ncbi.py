#!/usr/bin/env python3

from pathlib import Path
import csv
import json
import re
from collections import Counter, defaultdict

PROJECT = Path("/mnt/d/histamine_receptor_pharmacogenomics")

PROTEIN_CANDIDATES = PROJECT / "11_results/candidate_prioritization/hrh_protein_altering_candidates_canonical_safe_gpcr_topology_v2_ncbi.tsv"
UNIPROT_JSON = PROJECT / "07_protein_topology/uniprot_hrh_receptor_records_v1.json"

OUTDIR = PROJECT / "11_results/receptor_fate"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_CANDIDATES = OUTDIR / "hrh_receptor_biogenesis_fate_candidates_v1_v2_ncbi.tsv"
OUT_SUMMARY = OUTDIR / "hrh_receptor_biogenesis_fate_summary_v1_v2_ncbi.tsv"
OUT_TOP = OUTDIR / "hrh_top_receptor_fate_candidates_v1_v2_ncbi.tsv"

TABLE8 = PROJECT / "13_tables/manuscript_ready/Table_8_HRH_receptor_biogenesis_fate_candidates_v1_v2_ncbi.tsv"
RESULTS_MD = PROJECT / "14_manuscript/results_sections/results_receptor_biogenesis_fate_v1_v2_ncbi.md"

AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Stop": "*"
}

GENE_TO_UNIPROT = {
    "HRH1": "P35367",
    "HRH2": "P25021",
    "HRH3": "Q9Y5N1",
    "HRH4": "Q9H3N8",
}

def read_tsv(path):
    with path.open("r", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

def parse_hgvsp(text):
    text = str(text or "")
    # Handles ENSP...:p.His206Arg, p.His206Arg, p.Val389Ile, p.Tyr319Cys
    m = re.search(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|Stop|\*)", text)
    if not m:
        return "", "", ""
    ref3, pos, alt3 = m.group(1), int(m.group(2)), m.group(3)
    ref = AA3.get(ref3, ref3)
    alt = AA3.get(alt3, alt3)
    return ref, pos, alt

def find_sequence_in_obj(obj):
    if isinstance(obj, dict):
        if "sequence" in obj:
            seq = obj["sequence"]
            if isinstance(seq, str) and len(seq) > 50:
                return seq
            if isinstance(seq, dict):
                for key in ["value", "sequence"]:
                    if isinstance(seq.get(key), str) and len(seq.get(key)) > 50:
                        return seq.get(key)
        for val in obj.values():
            found = find_sequence_in_obj(val)
            if found:
                return found
    elif isinstance(obj, list):
        for item in obj:
            found = find_sequence_in_obj(item)
            if found:
                return found
    return ""

def collect_records_by_gene(obj):
    records = {}

    if isinstance(obj, dict):
        for gene in GENE_TO_UNIPROT:
            if gene in obj:
                records[gene] = obj[gene]

        # Also search by accession.
        for gene, acc in GENE_TO_UNIPROT.items():
            for key, val in obj.items():
                if acc in str(key):
                    records[gene] = val

        # Some JSONs are accession keyed.
        for gene, acc in GENE_TO_UNIPROT.items():
            if acc in obj:
                records[gene] = obj[acc]

    return records

def load_sequences(path):
    if not path.exists():
        return {}

    data = json.loads(path.read_text(encoding="utf-8"))
    records = collect_records_by_gene(data)
    seqs = {}

    for gene, rec in records.items():
        seq = find_sequence_in_obj(rec)
        if seq:
            seqs[gene] = seq

    # Fallback: recursive search for accessions in records.
    if len(seqs) < 4:
        text = path.read_text(encoding="utf-8", errors="replace")
        for gene, acc in GENE_TO_UNIPROT.items():
            if gene in seqs:
                continue
            # crude but safe fallback around accession is hard; leave missing.
            pass

    return seqs

def region_value(row):
    for key in ["gpcr_region_broad", "candidate_gpcr_region_broad", "topology_region", "region"]:
        if row.get(key):
            return row.get(key)
    return ""

def motif_value(row):
    for key in ["gpcr_motif_neighborhood", "candidate_gpcr_motif_neighborhood", "nearest_motif", "motif"]:
        if row.get(key):
            return row.get(key)
    return ""

def gene_value(row):
    for key in ["uploaded_gene", "gene", "candidate_gene"]:
        if row.get(key):
            return row.get(key)
    return ""

def hgvsp_value(row):
    for key in ["selected_hgvsp", "hgvsp", "candidate_selected_hgvsp"]:
        if row.get(key):
            return row.get(key)
    return ""

def variant_key_value(row):
    for key in ["variant_key", "candidate_variant_key"]:
        if row.get(key):
            return row.get(key)
    return ""

def window(seq, pos, flank=8):
    if not seq or not pos:
        return ""
    start = max(1, pos - flank)
    end = min(len(seq), pos + flank)
    return f"{start}-{end}:{seq[start-1:end]}"

def has_nglyc(seq, pos):
    if not seq or not pos:
        return False, ""
    hits = []
    for i in range(max(0, pos - 4), min(len(seq) - 2, pos + 2)):
        tri = seq[i:i+3]
        if len(tri) == 3 and tri[0] == "N" and tri[1] != "P" and tri[2] in {"S", "T"}:
            if i + 1 <= pos <= i + 3:
                hits.append(f"N-glyc:{i+1}-{i+3}:{tri}")
    return bool(hits), ";".join(hits)

def mutate(seq, pos, alt):
    if not seq or not pos or pos < 1 or pos > len(seq) or len(alt) != 1:
        return seq
    return seq[:pos-1] + alt + seq[pos:]

def nglyc_change(seq, pos, alt):
    if not seq or not pos:
        return "", ""
    wt_has, wt_hits = has_nglyc(seq, pos)
    mut_seq = mutate(seq, pos, alt)
    mut_has, mut_hits = has_nglyc(mut_seq, pos)
    if wt_has and not mut_has:
        return "loss_predicted", wt_hits
    if not wt_has and mut_has:
        return "gain_predicted", mut_hits
    if wt_has and mut_has:
        return "within_existing_motif", wt_hits
    return "none", ""

def classify_fate(row, seqs):
    gene = gene_value(row)
    hgvsp = hgvsp_value(row)
    ref, pos, alt = parse_hgvsp(hgvsp)
    region = region_value(row)
    motif = motif_value(row)
    seq = seqs.get(gene, "")

    classes = []
    evidence = []
    score = 0

    region_l = region.lower()
    motif_l = motif.lower()

    wt_window = window(seq, pos, flank=10)

    # N-terminal / extracellular biogenesis
    if "extracellular" in region_l or (pos and pos <= 60):
        classes.append("extracellular_or_N_terminal_biogenesis_context")
        evidence.append("Variant lies in extracellular or N-terminal receptor-biogenesis context.")
        score += 1

    # N-glycosylation
    ng_change, ng_detail = nglyc_change(seq, pos, alt)
    if ng_change != "none":
        classes.append("N_glycosylation_motif_candidate")
        evidence.append(f"N-glycosylation motif context: {ng_change}; {ng_detail}")
        score += 4 if "loss" in ng_change or "gain" in ng_change else 2

    # Cysteine/disulfide
    if ref == "C" or alt == "C" or "disulfide" in motif_l:
        classes.append("cysteine_or_disulfide_support_candidate")
        evidence.append("Variant changes or introduces cysteine, or lies near disulfide-support annotation.")
        score += 3

    # TM insertion / folding
    if "transmembrane" in region_l or "tm" in motif_l:
        classes.append("transmembrane_folding_or_membrane_insertion_candidate")
        evidence.append("Variant lies in transmembrane or TM-neighborhood context.")
        score += 2

    # Conserved GPCR motif
    if any(x in motif_l for x in ["dry", "npxxy", "cwxp", "motif"]):
        classes.append("conserved_GPCR_motif_neighborhood_candidate")
        evidence.append(f"Variant lies near conserved GPCR motif annotation: {motif}.")
        score += 3

    # Phosphorylation / arrestin / desensitization
    if "intracellular" in region_l or "c-terminal" in region_l or "cytoplasmic" in region_l:
        if ref in {"S", "T", "Y"} or alt in {"S", "T", "Y"}:
            classes.append("phosphorylation_desensitization_candidate")
            evidence.append("Variant removes or introduces Ser/Thr/Tyr in intracellular context.")
            score += 4
        else:
            classes.append("intracellular_signaling_or_arrestin_context")
            evidence.append("Variant lies in intracellular region relevant to signaling/desensitization context.")
            score += 1

    # Ubiquitination/degradation
    if "intracellular" in region_l or "c-terminal" in region_l or "cytoplasmic" in region_l:
        if ref == "K" or alt == "K":
            classes.append("ubiquitination_or_receptor_turnover_candidate")
            evidence.append("Variant removes or introduces Lys in intracellular context.")
            score += 4

    # Palmitoylation/cytoplasmic cysteine
    if ("intracellular" in region_l or "c-terminal" in region_l or "cytoplasmic" in region_l) and (ref == "C" or alt == "C"):
        classes.append("palmitoylation_or_membrane_proximal_cysteine_candidate")
        evidence.append("Variant changes cysteine in intracellular/membrane-proximal context.")
        score += 4

    # Internalization/degron motif simple screen
    if seq and pos:
        local = seq[max(0, pos-8):min(len(seq), pos+8)]
        mut_local = mutate(seq, pos, alt)[max(0, pos-8):min(len(seq), pos+8)]
        if any(x in local for x in ["LL", "Y"]) and ("intracellular" in region_l or "c-terminal" in region_l):
            classes.append("internalization_motif_context")
            evidence.append("Variant lies near simple intracellular Y/di-leucine internalization motif context.")
            score += 1
        if local != mut_local and ("PP" in local or "PEST" in local):
            classes.append("degron_or_turnover_context")
            evidence.append("Variant lies near simple Pro/PEST-like turnover context.")
            score += 1

    if not classes:
        classes.append("no_specific_receptor_fate_signal_detected")
        evidence.append("No receptor-fate motif signal detected by this rule-based scan.")

    priority = "background"
    if score >= 8:
        priority = "FateTier_1_high_priority"
    elif score >= 5:
        priority = "FateTier_2_moderate_priority"
    elif score >= 2:
        priority = "FateTier_3_context_candidate"

    return {
        "gene": gene,
        "variant_key": variant_key_value(row),
        "selected_hgvsp": hgvsp,
        "aa_ref": ref,
        "protein_position": pos,
        "aa_alt": alt,
        "gpcr_region_broad": region,
        "gpcr_motif_neighborhood": motif,
        "wt_local_sequence_window": wt_window,
        "receptor_fate_classes": ";".join(sorted(set(classes))),
        "receptor_fate_evidence": " | ".join(evidence),
        "receptor_fate_score": score,
        "receptor_fate_priority": priority,
    }

def main():
    rows = read_tsv(PROTEIN_CANDIDATES)
    seqs = load_sequences(UNIPROT_JSON)

    annotated = []
    for row in rows:
        rec = classify_fate(row, seqs)

        # Preserve important upstream columns if present.
        for key in [
            "global_AF", "superpop_af_range", "candidate_score",
            "candidate_tier", "selected_hgvsc", "gpcr_region_fine",
            "nearest_tm_label", "nearest_tm_distance"
        ]:
            if key in row:
                rec[key] = row.get(key, "")

        annotated.append(rec)

    annotated_sorted = sorted(
        annotated,
        key=lambda r: (int(r["receptor_fate_score"]), str(r["gene"]), str(r["variant_key"])),
        reverse=True
    )

    fields = [
        "gene", "variant_key", "selected_hgvsp", "aa_ref", "protein_position", "aa_alt",
        "global_AF", "superpop_af_range", "selected_hgvsc",
        "gpcr_region_broad", "gpcr_region_fine", "nearest_tm_label",
        "gpcr_motif_neighborhood", "wt_local_sequence_window",
        "receptor_fate_classes", "receptor_fate_evidence",
        "receptor_fate_score", "receptor_fate_priority",
        "candidate_score", "candidate_tier"
    ]

    write_tsv(OUT_CANDIDATES, annotated_sorted, fields)
    write_tsv(TABLE8, annotated_sorted, fields)
    write_tsv(OUT_TOP, annotated_sorted[:40], fields)

    summary_rows = []

    def add_summary(group, subset):
        priorities = Counter(r["receptor_fate_priority"] for r in subset)
        class_counter = Counter()
        for r in subset:
            for c in r["receptor_fate_classes"].split(";"):
                class_counter[c] += 1
        summary_rows.append({
            "group": group,
            "candidates": len(subset),
            "FateTier_1_high_priority": priorities.get("FateTier_1_high_priority", 0),
            "FateTier_2_moderate_priority": priorities.get("FateTier_2_moderate_priority", 0),
            "FateTier_3_context_candidate": priorities.get("FateTier_3_context_candidate", 0),
            "background": priorities.get("background", 0),
            "class_counts": ";".join(f"{k}:{v}" for k, v in sorted(class_counter.items())),
        })

    add_summary("ALL", annotated)
    for gene in sorted(set(r["gene"] for r in annotated)):
        add_summary(f"gene={gene}", [r for r in annotated if r["gene"] == gene])
    for region in sorted(set(r["gpcr_region_broad"] for r in annotated)):
        add_summary(f"region={region}", [r for r in annotated if r["gpcr_region_broad"] == region])

    summary_fields = [
        "group", "candidates", "FateTier_1_high_priority",
        "FateTier_2_moderate_priority", "FateTier_3_context_candidate",
        "background", "class_counts"
    ]
    write_tsv(OUT_SUMMARY, summary_rows, summary_fields)

    priority_counts = Counter(r["receptor_fate_priority"] for r in annotated)
    gene_counts = Counter(r["gene"] for r in annotated)
    class_counts = Counter()
    for r in annotated:
        for c in r["receptor_fate_classes"].split(";"):
            class_counts[c] += 1

    top_lines = []
    for r in annotated_sorted[:12]:
        top_lines.append(
            f"- {r['gene']} {r['variant_key']} {r['selected_hgvsp']}: "
            f"{r['receptor_fate_priority']}; classes={r['receptor_fate_classes']}; "
            f"score={r['receptor_fate_score']}."
        )

    md = f"""# Receptor biogenesis and receptor-fate mapping of HRH protein candidates

To address receptor maturation, turnover and signaling-fate mechanisms, canonical-safe HRH protein-altering candidates were screened against GPCR topology context and sequence-based receptor-fate rules. The scan evaluated {len(annotated)} canonical-safe protein candidates for N-terminal or extracellular biogenesis context, N-glycosylation motif disruption or creation, cysteine/disulfide support, transmembrane folding context, conserved GPCR motif neighborhoods, intracellular phosphorylation/desensitization motifs, lysine-dependent ubiquitination or receptor-turnover context, palmitoylation-related cysteine context and simple internalization/degron motif context.

The receptor-fate priority distribution was: {', '.join(f'{k}: {v}' for k, v in sorted(priority_counts.items()))}. Gene-level candidate counts were: {', '.join(f'{k}: {v}' for k, v in sorted(gene_counts.items()))}. Rule-based fate classes were distributed as: {', '.join(f'{k}: {v}' for k, v in sorted(class_counts.items()))}.

These annotations extend the HRH pharmacogenomic model beyond coding consequence and static topology. Variants in extracellular or N-terminal regions may influence receptor maturation, folding or ligand-accessible extracellular architecture; variants in transmembrane helices may alter receptor folding or conformational stability; and intracellular Ser/Thr/Tyr, Lys or Cys changes may plausibly affect phosphorylation, beta-arrestin recruitment, ubiquitination, palmitoylation, internalization or receptor turnover. These predictions are mechanistic prioritization signals and should not be interpreted as confirmed biochemical effects without experimental validation.

## Top receptor-fate candidates

{chr(10).join(top_lines)}

## Output files

- Receptor-fate candidate table: `11_results/receptor_fate/hrh_receptor_biogenesis_fate_candidates_v1_v2_ncbi.tsv`
- Manuscript-ready Table 8: `13_tables/manuscript_ready/Table_8_HRH_receptor_biogenesis_fate_candidates_v1_v2_ncbi.tsv`
- Receptor-fate summary: `11_results/receptor_fate/hrh_receptor_biogenesis_fate_summary_v1_v2_ncbi.tsv`
"""
    RESULTS_MD.write_text(md, encoding="utf-8")

    print(f"Wrote receptor-fate candidates: {OUT_CANDIDATES}")
    print(f"Wrote receptor-fate summary: {OUT_SUMMARY}")
    print(f"Wrote top receptor-fate candidates: {OUT_TOP}")
    print(f"Wrote Table 8: {TABLE8}")
    print(f"Wrote Results subsection: {RESULTS_MD}")
    print()
    print(f"Canonical-safe protein candidates analyzed: {len(annotated)}")
    print("Sequences loaded:")
    for gene in sorted(GENE_TO_UNIPROT):
        print(f"  {gene}: {'yes' if gene in seqs else 'no'}")
    print("Priority counts:")
    for k, v in sorted(priority_counts.items()):
        print(f"  {k}: {v}")
    print("Top fate classes:")
    for k, v in class_counts.most_common(12):
        print(f"  {k}: {v}")

if __name__ == "__main__":
    main()
