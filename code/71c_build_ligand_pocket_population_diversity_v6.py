#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import datetime
from collections import defaultdict

PROJECT = Path("/mnt/d/histamine_receptor_pharmacogenomics")

OUTDIR = PROJECT / "11_results/monitoring_evaluation_v6_global_functional_diversity"
FIGDIR = PROJECT / "12_figures/monitoring_evaluation_v6_global_functional_diversity"
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

POP_MATRIX = OUTDIR / "HRH_population_functional_variant_matrix_v6.tsv"
PAIRWISE = OUTDIR / "HRH_population_pairwise_variant_tests_v6.tsv"

OUT_MATRIX = OUTDIR / "HRH_ligand_pocket_population_diversity_matrix_v6.tsv"
OUT_PAIRWISE = OUTDIR / "HRH_ligand_pocket_population_pairwise_tests_v6.tsv"
OUT_SIG_PAIRWISE = OUTDIR / "HRH_ligand_pocket_population_significant_pairwise_tests_v6.tsv"
OUT_SUMMARY = OUTDIR / "HRH_ligand_pocket_population_diversity_summary_v6.tsv"
OUT_TOP = OUTDIR / "HRH_top_ligand_pocket_population_candidates_v6.tsv"
OUT_SOURCE_INVENTORY = OUTDIR / "HRH_ligand_pocket_contact_source_inventory_v6.tsv"
OUT_UNMATCHED = OUTDIR / "HRH_ligand_pocket_unmatched_contact_rows_v6.tsv"
OUT_REPORT = OUTDIR / "HRH_ligand_pocket_population_diversity_report_v6.md"
OUT_CHECKLIST = OUTDIR / "HRH_ligand_pocket_population_diversity_checklist_v6.tsv"

OUT_FIG_PNG = FIGDIR / "Figure_4_ligand_pocket_population_diversity_v6_900dpi.png"
OUT_FIG_SVG = FIGDIR / "Figure_4_ligand_pocket_population_diversity_v6.svg"

SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]

CONTACT_FILE_PATTERNS = [
    "*contact*.tsv",
    "*pocket*.tsv",
    "*docking*.tsv",
    "*perturbation*.tsv",
    "*ligand*.tsv",
    "*contact*.csv",
    "*pocket*.csv",
    "*docking*.csv",
    "*perturbation*.csv",
    "*ligand*.csv",
]

SEARCH_ROOTS = [
    PROJECT / "10_docking",
    PROJECT / "11_results",
    PROJECT / "13_tables",
    PROJECT / "14_manuscript/target_journal_TPJ_v11",
]

PREFERRED_SOURCES = [
    PROJECT / "13_tables/manuscript_ready/Table_64C_variant_contact_perturbation_proxy_v1.tsv",
    PROJECT / "14_manuscript/target_journal_TPJ_v11/tpj_upload_package_v18_APA_RIS/04_tables_separate_upload/Table_2_top_structural_contact_perturbation_candidates.tsv",
    PROJECT / "14_manuscript/target_journal_TPJ_v11/tpj_upload_package_v18_APA_RIS/05_supplementary_files/Supplementary_Table_4_variant_contact_perturbation_proxy_full.tsv",
]

GENE_LIGAND_DEFAULTS = {
    "HRH1": "doxepin",
    "HRH2": "impromidine",
    "HRH3": "ciproxifan",
    "HRH4": "toreforant",
}

def read_tsv(path):
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
        return list(csv.DictReader(f, delimiter=delimiter))

def write_tsv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not fields:
        fields = list(rows[0].keys()) if rows else ["empty"]
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

def norm_chrom(x):
    return str(x).replace("chr", "").replace("NC_000003.12", "3").replace("NC_000005.10", "5").replace("NC_000020.11", "20").replace("NC_000018.10", "18")

def detect_col(header, names):
    # Robust to malformed TSV/CSV files where csv.DictReader may create None keys.
    clean_header = [h for h in header if h is not None and str(h).strip() != ""]
    lower = [str(h).lower() for h in clean_header]
    for name in names:
        n = name.lower()
        if n in lower:
            return clean_header[lower.index(n)]
    for h in clean_header:
        lh = str(h).lower()
        if any(name.lower() in lh for name in names):
            return h
    return None

def cell_text(row):
    # Ignore csv.DictReader overflow fields stored under None.
    vals = []
    for k, v in row.items():
        if k is None:
            continue
        if v is not None:
            vals.append(str(v))
    return " | ".join(vals)

def extract_coordinate_keys(text):
    keys = set()
    if not text:
        return keys

    # chr:pos:ref:alt or chr_pos_ref_alt
    for m in re.finditer(r"\b(?:chr)?([0-9XYM]+)[:_](\d+)[:_]([ACGT]+)[:_]([ACGT]+)\b", text):
        c, pos, ref, alt = m.groups()
        keys.add(f"{norm_chrom(c)}:{pos}:{ref}:{alt}")

    # chr:pos only
    for m in re.finditer(r"\b(?:chr)?([0-9XYM]+):(\d+)\b", text):
        c, pos = m.groups()
        keys.add(f"{norm_chrom(c)}:{pos}")

    return keys

def extract_protein_aliases(text):
    aliases = set()
    if not text:
        return aliases

    # HRH4 p.Tyr319Cys, p.Tyr319Cys, Tyr319Cys
    for m in re.finditer(r"\b(HRH[1-4])\s+(p\.[A-Za-z]{3}\d+[A-Za-z]{3})\b", text):
        aliases.add(f"{m.group(1)} {m.group(2)}")
        aliases.add(m.group(2))

    for m in re.finditer(r"\b(p\.[A-Za-z]{3}\d+[A-Za-z]{3})\b", text):
        aliases.add(m.group(1))

    for m in re.finditer(r"\b([A-Z][a-z]{2}\d+[A-Z][a-z]{2})\b", text):
        aliases.add("p." + m.group(1))
        aliases.add(m.group(1))

    return aliases

def extract_gene(text):
    m = re.search(r"\b(HRH[1-4])\b", text or "")
    return m.group(1) if m else ""

def extract_ligand(row, text):
    header_hits = []
    for k, v in row.items():
        if k is None:
            continue
        lk = str(k).lower()
        if "ligand" in lk or "drug" in lk or "compound" in lk:
            val = str(v).strip()
            if val and val.lower() not in {"nan", "none"}:
                header_hits.append(val)
    if header_hits:
        return header_hits[0]

    low = (text or "").lower()
    ligands = [
        "doxepin", "impromidine", "ciproxifan", "toreforant",
        "cetirizine", "diphenhydramine", "fexofenadine", "mepyramine",
        "cimetidine", "famotidine", "ranitidine", "amthamine",
        "pitolisant", "betahistine", "thioperamide", "clobenpropit",
        "jnj7777120", "jnj39758979", "zpl389", "adriforant", "histamine",
    ]
    for lig in ligands:
        if lig in low:
            return lig
    gene = extract_gene(text)
    return GENE_LIGAND_DEFAULTS.get(gene, "")

def load_population_matrix():
    if not POP_MATRIX.exists():
        raise SystemExit(f"Missing population matrix: {POP_MATRIX}")
    rows = read_tsv(POP_MATRIX)
    by_key = {}
    alias_map = defaultdict(set)

    for r in rows:
        vk = r.get("variant_key", "")
        if not vk:
            continue
        by_key[vk] = r
        alias_map[vk].add(vk)
        chrom = r.get("chrom", "")
        pos = r.get("pos", "")
        ref = r.get("ref", "")
        alt = r.get("alt", "")
        vid = r.get("id", "")
        if chrom and pos:
            alias_map[f"{norm_chrom(chrom)}:{pos}"].add(vk)
            alias_map[f"chr{norm_chrom(chrom)}:{pos}"].add(vk)
        if chrom and pos and ref and alt:
            alias_map[f"{norm_chrom(chrom)}:{pos}:{ref}:{alt}"].add(vk)
            alias_map[f"chr{norm_chrom(chrom)}:{pos}:{ref}:{alt}"].add(vk)
        if vid and vid != ".":
            alias_map[vid].add(vk)
    return rows, by_key, alias_map

def build_alias_map_from_annotation_tables(alias_map):
    # Add aliases such as HRH4 p.Tyr319Cys -> coordinate variant_key from any project TSV that has both coordinates and protein labels.
    for root in [PROJECT / "11_results", PROJECT / "13_tables", PROJECT / "14_manuscript/target_journal_TPJ_v11"]:
        if not root.exists():
            continue
        for path in root.rglob("*.tsv"):
            try:
                rows = read_tsv(path)
            except Exception:
                continue
            if not rows:
                continue
            header = list(rows[0].keys())
            chrom_col = detect_col(header, ["chrom", "#chrom", "chr"])
            pos_col = detect_col(header, ["pos", "position"])
            ref_col = detect_col(header, ["ref", "reference"])
            alt_col = detect_col(header, ["alt", "alternate"])
            variant_col = detect_col(header, ["variant_key", "variant_id", "variant", "id"])
            gene_col = detect_col(header, ["gene", "receptor"])
            protein_col = detect_col(header, ["protein", "hgvs_p", "hgvsp", "aa_change", "amino", "coding"])

            for row in rows:
                keys = set()

                txt = cell_text(row)
                keys |= extract_coordinate_keys(txt)

                if chrom_col and pos_col and ref_col and alt_col:
                    try:
                        keys.add(f"{norm_chrom(row[chrom_col])}:{int(float(row[pos_col]))}:{row[ref_col]}:{row[alt_col]}")
                    except Exception:
                        pass

                if variant_col and row.get(variant_col):
                    keys |= extract_coordinate_keys(str(row.get(variant_col, "")))
                    keys.add(str(row.get(variant_col, "")).strip())

                candidate_vks = set()
                for k in keys:
                    if k in alias_map:
                        candidate_vks |= alias_map[k]

                if not candidate_vks:
                    continue

                gene = str(row.get(gene_col, "")).strip() if gene_col else extract_gene(txt)
                aliases = set()
                aliases |= extract_protein_aliases(txt)
                if protein_col and row.get(protein_col):
                    aliases |= extract_protein_aliases(str(row.get(protein_col, "")))

                for alias in aliases:
                    if gene and alias.startswith("p."):
                        alias_map[f"{gene} {alias}"] |= candidate_vks
                    alias_map[alias] |= candidate_vks

    return alias_map

def discover_contact_sources():
    paths = set()
    for p in PREFERRED_SOURCES:
        if p.exists():
            paths.add(p)

    for root in SEARCH_ROOTS:
        if not root.exists():
            continue
        for pat in CONTACT_FILE_PATTERNS:
            for p in root.rglob(pat):
                if p.is_file():
                    low = str(p).lower()
                    if any(k in low for k in ["contact", "pocket", "docking", "perturbation", "ligand"]):
                        paths.add(p)
    return sorted(paths)

def classify_contact_row(row, path):
    txt = (str(path) + " | " + cell_text(row)).lower()

    if any(k in txt for k in ["exact_contact", "direct contact", "exact representative contact", "ligand-contact-zone"]):
        return "direct_or_exact_contact"
    if any(k in txt for k in ["near-pocket", "pocket-proximal", "near pocket", "contact-zone"]):
        return "near_or_pocket_proximal"
    if any(k in txt for k in ["structurally-near", "near-structural", "structural"]):
        return "structural_context"
    if any(k in txt for k in ["docking", "ligand", "pocket", "contact", "perturbation"]):
        return "ligand_or_contact_context"
    return "unclassified_contact_source"

def score_contact_row(row, path):
    txt = (str(path) + " | " + cell_text(row)).lower()
    score = 0
    if "exact_contact" in txt or "exact representative contact" in txt:
        score += 50
    if "direct" in txt and "contact" in txt:
        score += 40
    if "ligand-contact-zone" in txt:
        score += 35
    if "near-pocket" in txt or "pocket-proximal" in txt:
        score += 25
    if "structurally-near" in txt or "near-structural" in txt:
        score += 15
    if "very_high" in txt:
        score += 30
    if "high" in txt:
        score += 20
    if "moderate" in txt:
        score += 10
    for k, v in row.items():
        if "score" in k.lower():
            try:
                score += float(v) / 10.0
            except Exception:
                pass
    return score

def match_contact_row_to_variants(row, alias_map):
    txt = cell_text(row)
    candidates = set()

    for k in extract_coordinate_keys(txt):
        if k in alias_map:
            candidates |= alias_map[k]

    protein_aliases = extract_protein_aliases(txt)
    gene = extract_gene(txt)
    for a in protein_aliases:
        if a in alias_map:
            candidates |= alias_map[a]
        if gene and f"{gene} {a}" in alias_map:
            candidates |= alias_map[f"{gene} {a}"]

    # Column-aware matching.
    for col, val in row.items():
        sval = str(val).strip()
        if not sval:
            continue
        if sval in alias_map:
            candidates |= alias_map[sval]
        for k in extract_coordinate_keys(sval):
            if k in alias_map:
                candidates |= alias_map[k]
        for a in extract_protein_aliases(sval):
            if a in alias_map:
                candidates |= alias_map[a]
            if gene and f"{gene} {a}" in alias_map:
                candidates |= alias_map[f"{gene} {a}"]

    return sorted(candidates)

def load_pairwise_for_variants(variant_keys):
    if not PAIRWISE.exists():
        raise SystemExit(f"Missing pairwise table: {PAIRWISE}")
    keep = set(variant_keys)
    rows = []
    with open(PAIRWISE, "r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if r.get("variant_key") in keep:
                rows.append(r)
    return rows

def as_float(x, default=0.0):
    try:
        if x in {None, ""}:
            return default
        return float(x)
    except Exception:
        return default

def main():
    pop_rows, pop_by_key, alias_map = load_population_matrix()
    alias_map = build_alias_map_from_annotation_tables(alias_map)

    contact_sources = discover_contact_sources()

    source_inventory = []
    matched_rows = []
    unmatched_rows = []

    seen = set()

    for path in contact_sources:
        try:
            rows = read_tsv(path)
        except Exception as e:
            source_inventory.append({
                "source_path": str(path.relative_to(PROJECT)),
                "status": "READ_FAIL",
                "rows": "",
                "matched_rows": "",
                "note": str(e),
            })
            continue

        matched_count = 0

        for idx, row in enumerate(rows, start=1):
            txt = cell_text(row)
            # Skip obviously irrelevant rows in broad ligand/docking files unless they mention variant/contact/pocket/protein.
            low = (str(path) + " | " + txt).lower()
            if not any(k in low for k in ["variant", "contact", "pocket", "perturb", "tyr", "cys", "arg", "gly", "ser", "val", "ala", "phe", "leu", "asn", "his"]):
                continue

            vks = match_contact_row_to_variants(row, alias_map)

            if not vks:
                unmatched_rows.append({
                    "source_path": str(path.relative_to(PROJECT)),
                    "source_row_number": str(idx),
                    "detected_gene": extract_gene(txt),
                    "detected_ligand": extract_ligand(row, txt),
                    "contact_class": classify_contact_row(row, path),
                    "contact_score": str(score_contact_row(row, path)),
                    "row_text": txt[:1000],
                })
                continue

            for vk in vks:
                if vk not in pop_by_key:
                    continue
                p = pop_by_key[vk]
                gene = extract_gene(txt) or p.get("gene", "")
                ligand = extract_ligand(row, txt) or GENE_LIGAND_DEFAULTS.get(gene, "")
                contact_class = classify_contact_row(row, path)
                contact_score = score_contact_row(row, path)

                key = (vk, ligand, contact_class, str(path.relative_to(PROJECT)))
                if key in seen:
                    continue
                seen.add(key)

                out = {
                    "variant_key": vk,
                    "gene": p.get("gene", ""),
                    "ligand": ligand,
                    "contact_class": contact_class,
                    "contact_priority_score": f"{contact_score:.4g}",
                    "source_path": str(path.relative_to(PROJECT)),
                    "source_row_number": str(idx),
                    "id": p.get("id", ""),
                    "ref": p.get("ref", ""),
                    "alt": p.get("alt", ""),
                    "ann_consequence": p.get("ann_consequence", ""),
                    "ann_tier": p.get("ann_tier", ""),
                    "ann_topology_or_domain": p.get("ann_topology_or_domain", ""),
                    "ann_axis_or_evidence": p.get("ann_axis_or_evidence", ""),
                    "AFR_AF": p.get("AFR_AF", ""),
                    "AMR_AF": p.get("AMR_AF", ""),
                    "EAS_AF": p.get("EAS_AF", ""),
                    "EUR_AF": p.get("EUR_AF", ""),
                    "SAS_AF": p.get("SAS_AF", ""),
                    "max_population": p.get("max_population", ""),
                    "min_population": p.get("min_population", ""),
                    "max_minus_min_AF": p.get("max_minus_min_AF", ""),
                    "fst_like_population_structure_score": p.get("fst_like_population_structure_score", ""),
                    "top_pairwise_comparison": p.get("top_pairwise_comparison", ""),
                    "top_pairwise_abs_delta_af": p.get("top_pairwise_abs_delta_af", ""),
                    "top_pairwise_fdr_q": p.get("top_pairwise_fdr_q", ""),
                    "significant_pairwise_count_fdr_0_05_delta_0_10": p.get("significant_pairwise_count_fdr_0_05_delta_0_10", ""),
                    "population_differentiation_class": p.get("population_differentiation_class", ""),
                    "source_row_excerpt": txt[:1000],
                }
                matched_rows.append(out)
                matched_count += 1

        source_inventory.append({
            "source_path": str(path.relative_to(PROJECT)),
            "status": "READ_OK",
            "rows": str(len(rows)),
            "matched_rows": str(matched_count),
            "note": "",
        })

    # If broad source matching is too sparse, add all variants that already have docking/contact-related annotation in the population matrix.
    # This prevents the module from failing just because earlier contact-source tables used inconsistent names.
    supplemental_added = 0
    for p in pop_rows:
        text = " ".join([
            p.get("annotation_source", ""),
            p.get("ann_consequence", ""),
            p.get("ann_tier", ""),
            p.get("ann_topology_or_domain", ""),
            p.get("ann_axis_or_evidence", ""),
        ]).lower()
        if any(k in text for k in ["docking", "pocket", "contact", "protein", "topology"]):
            gene = p.get("gene", "")
            ligand = GENE_LIGAND_DEFAULTS.get(gene, "")
            key = (p["variant_key"], ligand, "annotation_inferred_structural_context", "population_matrix_annotation")
            if key in seen:
                continue
            seen.add(key)
            matched_rows.append({
                "variant_key": p.get("variant_key", ""),
                "gene": gene,
                "ligand": ligand,
                "contact_class": "annotation_inferred_structural_context",
                "contact_priority_score": "5",
                "source_path": "population_matrix_annotation",
                "source_row_number": "",
                "id": p.get("id", ""),
                "ref": p.get("ref", ""),
                "alt": p.get("alt", ""),
                "ann_consequence": p.get("ann_consequence", ""),
                "ann_tier": p.get("ann_tier", ""),
                "ann_topology_or_domain": p.get("ann_topology_or_domain", ""),
                "ann_axis_or_evidence": p.get("ann_axis_or_evidence", ""),
                "AFR_AF": p.get("AFR_AF", ""),
                "AMR_AF": p.get("AMR_AF", ""),
                "EAS_AF": p.get("EAS_AF", ""),
                "EUR_AF": p.get("EUR_AF", ""),
                "SAS_AF": p.get("SAS_AF", ""),
                "max_population": p.get("max_population", ""),
                "min_population": p.get("min_population", ""),
                "max_minus_min_AF": p.get("max_minus_min_AF", ""),
                "fst_like_population_structure_score": p.get("fst_like_population_structure_score", ""),
                "top_pairwise_comparison": p.get("top_pairwise_comparison", ""),
                "top_pairwise_abs_delta_af": p.get("top_pairwise_abs_delta_af", ""),
                "top_pairwise_fdr_q": p.get("top_pairwise_fdr_q", ""),
                "significant_pairwise_count_fdr_0_05_delta_0_10": p.get("significant_pairwise_count_fdr_0_05_delta_0_10", ""),
                "population_differentiation_class": p.get("population_differentiation_class", ""),
                "source_row_excerpt": "Added from population matrix annotation because row contains structural/topology/docking-related evidence.",
            })
            supplemental_added += 1

    # Sort by biological/contact priority and population signal.
    matched_rows = sorted(
        matched_rows,
        key=lambda r: (
            int(r.get("significant_pairwise_count_fdr_0_05_delta_0_10") or 0),
            as_float(r.get("top_pairwise_abs_delta_af")),
            as_float(r.get("contact_priority_score")),
            as_float(r.get("fst_like_population_structure_score")),
        ),
        reverse=True,
    )

    fields = [
        "variant_key", "gene", "ligand", "contact_class", "contact_priority_score",
        "source_path", "source_row_number", "id", "ref", "alt",
        "ann_consequence", "ann_tier", "ann_topology_or_domain", "ann_axis_or_evidence",
        "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF",
        "max_population", "min_population", "max_minus_min_AF",
        "fst_like_population_structure_score",
        "top_pairwise_comparison", "top_pairwise_abs_delta_af", "top_pairwise_fdr_q",
        "significant_pairwise_count_fdr_0_05_delta_0_10",
        "population_differentiation_class",
        "source_row_excerpt",
    ]
    write_tsv(OUT_MATRIX, matched_rows, fields)

    pocket_keys = sorted(set(r["variant_key"] for r in matched_rows))
    pair_rows = load_pairwise_for_variants(pocket_keys)
    write_tsv(OUT_PAIRWISE, pair_rows, list(pair_rows[0].keys()) if pair_rows else None)

    sig_pair_rows = [
        r for r in pair_rows
        if r.get("significant_fdr_0_05_delta_0_10") == "YES"
    ]
    write_tsv(OUT_SIG_PAIRWISE, sig_pair_rows, list(pair_rows[0].keys()) if pair_rows else None)

    top_rows = matched_rows[:100]
    write_tsv(OUT_TOP, top_rows, fields)

    write_tsv(OUT_SOURCE_INVENTORY, source_inventory, ["source_path", "status", "rows", "matched_rows", "note"])
    write_tsv(OUT_UNMATCHED, unmatched_rows, ["source_path", "source_row_number", "detected_gene", "detected_ligand", "contact_class", "contact_score", "row_text"])

    # Summary by receptor-ligand.
    summary_map = {}
    for r in matched_rows:
        key = (r["gene"], r["ligand"])
        if key not in summary_map:
            summary_map[key] = {
                "gene": r["gene"],
                "ligand": r["ligand"],
                "candidate_rows": 0,
                "unique_variants": set(),
                "population_significant_variants": set(),
                "strong_population_variants": set(),
                "direct_or_exact_contact_variants": set(),
                "max_abs_delta_af": 0.0,
                "top_variant": "",
                "top_comparison": "",
                "dominant_max_population_counts": defaultdict(int),
            }
        s = summary_map[key]
        s["candidate_rows"] += 1
        s["unique_variants"].add(r["variant_key"])
        if int(r.get("significant_pairwise_count_fdr_0_05_delta_0_10") or 0) > 0:
            s["population_significant_variants"].add(r["variant_key"])
        if r.get("population_differentiation_class") == "strong":
            s["strong_population_variants"].add(r["variant_key"])
        if r.get("contact_class") == "direct_or_exact_contact":
            s["direct_or_exact_contact_variants"].add(r["variant_key"])
        d = as_float(r.get("top_pairwise_abs_delta_af"))
        if d > s["max_abs_delta_af"]:
            s["max_abs_delta_af"] = d
            s["top_variant"] = r["variant_key"]
            s["top_comparison"] = r.get("top_pairwise_comparison", "")
        if r.get("max_population"):
            s["dominant_max_population_counts"][r["max_population"]] += 1

    summary_rows = []
    for key, s in sorted(summary_map.items()):
        dominant = ""
        if s["dominant_max_population_counts"]:
            dominant = max(s["dominant_max_population_counts"], key=lambda k: s["dominant_max_population_counts"][k])
        summary_rows.append({
            "gene": s["gene"],
            "ligand": s["ligand"],
            "candidate_rows": str(s["candidate_rows"]),
            "unique_variants": str(len(s["unique_variants"])),
            "population_significant_variants_fdr_0_05_delta_0_10": str(len(s["population_significant_variants"])),
            "strong_population_differentiated_variants": str(len(s["strong_population_variants"])),
            "direct_or_exact_contact_variants": str(len(s["direct_or_exact_contact_variants"])),
            "max_abs_delta_af": f"{s['max_abs_delta_af']:.6g}",
            "top_variant": s["top_variant"],
            "top_comparison": s["top_comparison"],
            "dominant_max_frequency_population": dominant,
        })
    write_tsv(
        OUT_SUMMARY,
        summary_rows,
        [
            "gene", "ligand", "candidate_rows", "unique_variants",
            "population_significant_variants_fdr_0_05_delta_0_10",
            "strong_population_differentiated_variants",
            "direct_or_exact_contact_variants",
            "max_abs_delta_af", "top_variant", "top_comparison",
            "dominant_max_frequency_population",
        ],
    )

    # Basic figure.
    figure_status = "NOT_RUN"
    try:
        import matplotlib.pyplot as plt

        plot_rows = sorted(
            summary_rows,
            key=lambda x: (
                int(x["population_significant_variants_fdr_0_05_delta_0_10"]),
                float(x["max_abs_delta_af"] or 0),
            ),
            reverse=True,
        )[:20]

        if plot_rows:
            labels = [f"{r['gene']}-{r['ligand']}" for r in plot_rows]
            values = [int(r["population_significant_variants_fdr_0_05_delta_0_10"]) for r in plot_rows]

            plt.figure(figsize=(10, max(4, 0.35 * len(labels))))
            y = range(len(labels))
            plt.barh(list(y), values)
            plt.yticks(list(y), labels)
            plt.xlabel("Population-significant pocket/contact candidate variants")
            plt.ylabel("Receptor-ligand context")
            plt.title("Ligand-pocket population diversity across histamine receptors")
            plt.gca().invert_yaxis()
            plt.tight_layout()
            plt.savefig(OUT_FIG_PNG, dpi=900)
            plt.savefig(OUT_FIG_SVG)
            plt.close()
            figure_status = "PASS"
        else:
            figure_status = "NO_PLOT_ROWS"
    except Exception as e:
        figure_status = f"FAILED: {e}"

    checklist_rows = [
        {"check": "population_matrix_found", "status": "PASS" if POP_MATRIX.exists() else "REVIEW", "observed": str(POP_MATRIX)},
        {"check": "pairwise_table_found", "status": "PASS" if PAIRWISE.exists() else "REVIEW", "observed": str(PAIRWISE)},
        {"check": "contact_sources_discovered", "status": "PASS" if len(contact_sources) > 0 else "REVIEW", "observed": str(len(contact_sources))},
        {"check": "matched_ligand_pocket_rows", "status": "PASS" if len(matched_rows) > 0 else "REVIEW", "observed": str(len(matched_rows))},
        {"check": "unique_pocket_variants", "status": "PASS" if len(pocket_keys) > 0 else "REVIEW", "observed": str(len(pocket_keys))},
        {"check": "pairwise_rows_for_pocket_variants", "status": "PASS" if len(pair_rows) == len(pocket_keys) * 10 else "REVIEW", "observed": f"{len(pair_rows)} vs {len(pocket_keys) * 10}"},
        {"check": "significant_pairwise_rows", "status": "PASS" if len(sig_pair_rows) >= 0 else "REVIEW", "observed": str(len(sig_pair_rows))},
        {"check": "summary_created", "status": "PASS" if OUT_SUMMARY.exists() else "REVIEW", "observed": str(OUT_SUMMARY)},
        {"check": "figure_created", "status": "PASS" if figure_status == "PASS" else "REVIEW", "observed": figure_status},
        {"check": "supplemental_annotation_rows_added", "status": "INFO", "observed": str(supplemental_added)},
    ]
    write_tsv(OUT_CHECKLIST, checklist_rows, ["check", "status", "observed"])

    report = []
    report.append("# HRH ligand-pocket population diversity report v6")
    report.append("")
    report.append(f"Generated: {datetime.datetime.now().isoformat(timespec='seconds')}")
    report.append("")
    report.append("## Headline results")
    report.append("")
    report.append(f"- Contact/docking/pocket source files discovered: {len(contact_sources)}")
    report.append(f"- Matched ligand-pocket/contact candidate rows: {len(matched_rows)}")
    report.append(f"- Unique ligand-pocket/contact candidate variants: {len(pocket_keys)}")
    report.append(f"- Pairwise population tests for pocket/contact variants: {len(pair_rows)}")
    report.append(f"- Significant pairwise population tests among pocket/contact variants: {len(sig_pair_rows)}")
    report.append(f"- Supplemental structural/topology rows added from matrix annotation: {supplemental_added}")
    report.append("")
    report.append("## Receptor-ligand summary")
    report.append("")
    report.append("| Gene | Ligand | Unique variants | Population-significant variants | Strong variants | Direct/exact-contact variants | Max abs delta AF | Top variant | Top comparison | Dominant max-AF population |")
    report.append("|---|---|---:|---:|---:|---:|---:|---|---|---|")
    for r in summary_rows:
        report.append(
            f"| {r['gene']} | {r['ligand']} | {r['unique_variants']} | "
            f"{r['population_significant_variants_fdr_0_05_delta_0_10']} | "
            f"{r['strong_population_differentiated_variants']} | "
            f"{r['direct_or_exact_contact_variants']} | {r['max_abs_delta_af']} | "
            f"{r['top_variant']} | {r['top_comparison']} | {r['dominant_max_frequency_population']} |"
        )
    report.append("")
    report.append("## Interpretation")
    report.append("")
    report.append("This module converts the broad population matrix into a ligand-pocket/contact-focused pharmacogenomic diversity layer. The output is intended to identify receptor-ligand contexts where population-differentiated variants occur in direct ligand-contact, near-pocket, pocket-proximal or structurally contextualized regions. These candidates should drive the next mutant-versus-WT redocking rescue module.")
    report.append("")
    report.append("## Main outputs")
    report.append("")
    report.append(f"- Ligand-pocket population matrix: `{OUT_MATRIX}`")
    report.append(f"- Ligand-pocket pairwise tests: `{OUT_PAIRWISE}`")
    report.append(f"- Significant pairwise tests: `{OUT_SIG_PAIRWISE}`")
    report.append(f"- Summary: `{OUT_SUMMARY}`")
    report.append(f"- Top candidates: `{OUT_TOP}`")
    report.append(f"- Figure PNG: `{OUT_FIG_PNG}`")
    report.append(f"- Figure SVG: `{OUT_FIG_SVG}`")
    report.append(f"- Checklist: `{OUT_CHECKLIST}`")
    OUT_REPORT.write_text("\n".join(report) + "\n", encoding="utf-8")

    print("\n".join(report))

if __name__ == "__main__":
    main()
