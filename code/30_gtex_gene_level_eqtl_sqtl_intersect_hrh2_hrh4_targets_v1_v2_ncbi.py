#!/usr/bin/env python3

from pathlib import Path
import csv
import json
import time
import urllib.parse
import urllib.request
import urllib.error
from collections import Counter, defaultdict

PROJECT = Path("/mnt/d/histamine_receptor_pharmacogenomics")

TARGETS = PROJECT / "06_regulatory_genomics/regulatory_deepening_targets/HRH2_HRH4_regulatory_deepening_targets_v1_v2_ncbi.tsv"

OUTDIR = PROJECT / "06_regulatory_genomics/gtex_eqtl_sqtl_lookup"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_GENE_IDS = OUTDIR / "HRH2_HRH4_GTEx_reference_gene_ids_v1_v2_ncbi.tsv"
OUT_GENE_LEVEL = OUTDIR / "HRH2_HRH4_GTEx_gene_level_eQTL_sQTL_records_v1_v2_ncbi.tsv"
OUT_TARGET_INTERSECT = OUTDIR / "HRH2_HRH4_GTEx_target_intersections_v1_v2_ncbi.tsv"
OUT_AUDIT = OUTDIR / "HRH2_HRH4_GTEx_gene_level_query_audit_v1_v2_ncbi.tsv"
OUT_SUMMARY = OUTDIR / "HRH2_HRH4_GTEx_gene_level_intersection_summary_v1_v2_ncbi.tsv"
OUT_RAW = OUTDIR / "HRH2_HRH4_GTEx_gene_level_raw_responses_v1_v2_ncbi.jsonl"

TABLE12 = PROJECT / "13_tables/manuscript_ready/Table_12_HRH2_HRH4_GTEx_gene_level_eQTL_sQTL_intersections_v1_v2_ncbi.tsv"
RESULTS_MD = PROJECT / "14_manuscript/results_sections/results_hrh2_hrh4_gtex_gene_level_eqtl_sqtl_intersection_v1_v2_ncbi.md"

GENES = ["HRH2", "HRH4"]
DATASETS = ["gtex_v10", "gtex_v8"]
MAX_PAGES = 20
ITEMS_PER_PAGE = 100000
TIMEOUT = 30
SLEEP = 0.2

API = "https://gtexportal.org/api/v2"

def read_tsv(path):
    with path.open("r", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

def request_json(path, params):
    query = urllib.parse.urlencode(params, doseq=True)
    url = f"{API}{path}?{query}"
    req = urllib.request.Request(url, headers={"User-Agent": "HRH-pharmacogenomics/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
            txt = resp.read().decode("utf-8", errors="replace")
            status = getattr(resp, "status", 200)
        try:
            data = json.loads(txt)
        except Exception:
            data = {"raw_text": txt[:1000]}
        return status, url, data, ""
    except urllib.error.HTTPError as e:
        txt = ""
        try:
            txt = e.read().decode("utf-8", errors="replace")
        except Exception:
            pass
        return e.code, url, {"error_text": txt[:1000]}, f"HTTPError:{e.code}"
    except Exception as e:
        return "", url, {}, f"{type(e).__name__}:{e}"

def extract_list(payload):
    if not isinstance(payload, dict):
        return []

    for key in ["data", "results", "result"]:
        val = payload.get(key)
        if isinstance(val, list):
            return val
        if isinstance(val, dict):
            for subkey in ["data", "results", "result"]:
                subv = val.get(subkey)
                if isinstance(subv, list):
                    return subv

    for val in payload.values():
        if isinstance(val, list) and (not val or isinstance(val[0], dict)):
            return val

    return []

def flatten(d):
    out = {}
    if not isinstance(d, dict):
        return out
    for k, v in d.items():
        if isinstance(v, (dict, list)):
            out[k] = json.dumps(v, sort_keys=True)
        else:
            out[k] = v
    return out

def find_field(d, candidates):
    for k in candidates:
        if k in d and d[k] not in ["", None]:
            return d[k]
    lower = {str(k).lower(): v for k, v in d.items()}
    for k in candidates:
        lk = str(k).lower()
        if lk in lower and lower[lk] not in ["", None]:
            return lower[lk]
    return ""

def get_versioned_gene_ids(rawf):
    rows = []

    for gene in GENES:
        for dataset in DATASETS:
            # The reference/gene endpoint uses latest-release defaults unless dataset-independent
            # gencodeVersion/genomeBuild are specified. Here we search by symbol under GRCh38.
            status, url, payload, error = request_json(
                "/reference/gene",
                {
                    "geneId": [gene],
                    "genomeBuild": "GRCh38/hg38",
                    "page": 0,
                    "itemsPerPage": 100,
                },
            )

            rawf.write(json.dumps({
                "query_type": "reference_gene",
                "gene": gene,
                "dataset": dataset,
                "status": status,
                "error": error,
                "url": url,
                "payload": payload,
            }) + "\n")

            records = extract_list(payload)
            for rec in records:
                flat = flatten(rec)
                gid = find_field(flat, ["gencodeId", "geneId", "id"])
                symbol = find_field(flat, ["geneSymbol", "symbol", "geneName"])
                chrom = find_field(flat, ["chromosome", "chrom"])
                start = find_field(flat, ["start", "startPos"])
                end = find_field(flat, ["end", "endPos"])

                if symbol and str(symbol).upper() != gene:
                    continue

                rows.append({
                    "gene": gene,
                    "dataset": dataset,
                    "query_symbol": gene,
                    "gtex_gencode_id": gid,
                    "gtex_gene_symbol": symbol or gene,
                    "chromosome": chrom,
                    "start": start,
                    "end": end,
                    "http_status": status,
                    "error": error,
                    "url": url,
                    "raw_json": json.dumps(flat, sort_keys=True),
                })

            time.sleep(SLEEP)

    # Fallback: if the API record parser does not capture IDs, use unversioned IDs as last resort.
    have = {r["gene"] for r in rows if r.get("gtex_gencode_id")}
    fallback = {
        "HRH2": "ENSG00000113749",
        "HRH4": "ENSG00000134489",
    }
    for gene in GENES:
        if gene not in have:
            rows.append({
                "gene": gene,
                "dataset": "fallback_unversioned",
                "query_symbol": gene,
                "gtex_gencode_id": fallback[gene],
                "gtex_gene_symbol": gene,
                "chromosome": "",
                "start": "",
                "end": "",
                "http_status": "",
                "error": "fallback_unversioned_id_used",
                "url": "",
                "raw_json": "",
            })

    # Deduplicate by gene and gencode id.
    seen = set()
    out = []
    for r in rows:
        key = (r["gene"], r["gtex_gencode_id"])
        if key in seen or not r["gtex_gencode_id"]:
            continue
        seen.add(key)
        out.append(r)

    return out

def fetch_associations(gene_id_rows, rawf):
    association_rows = []
    audit_rows = []

    endpoint_specs = [
        ("eQTL", "/association/singleTissueEqtl"),
        ("sQTL", "/association/singleTissueSqtl"),
        ("iQTL", "/association/singleTissueIEqtl"),
        ("isQTL", "/association/singleTissueISqtl"),
    ]

    for gene_row in gene_id_rows:
        gene = gene_row["gene"]
        gencode_id = gene_row["gtex_gencode_id"]

        for dataset in DATASETS:
            for kind, path in endpoint_specs:
                for page in range(MAX_PAGES):
                    params = {
                        "gencodeId": [gencode_id],
                        "datasetId": dataset,
                        "page": page,
                        "itemsPerPage": ITEMS_PER_PAGE,
                    }

                    status, url, payload, error = request_json(path, params)
                    records = extract_list(payload)

                    rawf.write(json.dumps({
                        "query_type": "association_gene_level",
                        "gene": gene,
                        "gencode_id": gencode_id,
                        "dataset": dataset,
                        "kind": kind,
                        "page": page,
                        "status": status,
                        "error": error,
                        "url": url,
                        "n_records": len(records),
                        "payload": payload,
                    }) + "\n")

                    audit_rows.append({
                        "gene": gene,
                        "gencode_id": gencode_id,
                        "dataset": dataset,
                        "kind": kind,
                        "page": page,
                        "http_status": status,
                        "error": error,
                        "n_records": len(records),
                        "url": url,
                    })

                    for rec in records:
                        flat = flatten(rec)
                        variant_id = find_field(flat, ["variantId", "variant_id", "snpId", "snp_id"])
                        tissue = find_field(flat, ["tissueSiteDetailId", "tissueSiteDetail", "tissueSiteDetailName", "tissue"])
                        pvalue = find_field(flat, ["pValue", "pvalue", "pval", "p_value"])
                        qvalue = find_field(flat, ["qValue", "qvalue", "q_value"])
                        nes = find_field(flat, ["nes", "slope", "effectSize", "beta"])
                        phenotype = find_field(flat, ["phenotypeId", "phenotype_id", "geneSymbol", "gencodeId"])

                        association_rows.append({
                            "gene": gene,
                            "query_gencode_id": gencode_id,
                            "dataset": dataset,
                            "association_kind": kind,
                            "variantId": variant_id,
                            "tissue": tissue,
                            "pValue": pvalue,
                            "qValue": qvalue,
                            "effect_or_nes": nes,
                            "phenotype_or_gene": phenotype,
                            "raw_json": json.dumps(flat, sort_keys=True),
                        })

                    # Stop pagination when no records or page was not successful.
                    if status != 200 or not records:
                        break

                    time.sleep(SLEEP)

    return association_rows, audit_rows

def main():
    targets = read_tsv(TARGETS)
    target_by_id = defaultdict(list)

    for t in targets:
        vid = t.get("gtex_variant_id_b38", "")
        if vid:
            target_by_id[vid].append(t)

    gene_id_rows = []
    association_rows = []
    audit_rows = []

    with OUT_RAW.open("w", encoding="utf-8") as rawf:
        gene_id_rows = get_versioned_gene_ids(rawf)
        association_rows, audit_rows = fetch_associations(gene_id_rows, rawf)

    gene_fields = [
        "gene", "dataset", "query_symbol", "gtex_gencode_id", "gtex_gene_symbol",
        "chromosome", "start", "end", "http_status", "error", "url", "raw_json"
    ]
    write_tsv(OUT_GENE_IDS, gene_id_rows, gene_fields)

    assoc_fields = [
        "gene", "query_gencode_id", "dataset", "association_kind",
        "variantId", "tissue", "pValue", "qValue", "effect_or_nes",
        "phenotype_or_gene", "raw_json"
    ]
    write_tsv(OUT_GENE_LEVEL, association_rows, assoc_fields)

    intersect_rows = []

    for a in association_rows:
        vid = str(a.get("variantId", ""))
        matched_targets = target_by_id.get(vid, [])

        # Also tolerate chrless or b37/b38-like minor differences.
        if not matched_targets and vid:
            vid_norm = vid.replace("_b37", "_b38")
            matched_targets = target_by_id.get(vid_norm, [])

        for t in matched_targets:
            row = dict(a)
            row.update({
                "target_variant_key": t.get("variant_key", ""),
                "target_gtex_variant_id_b38": t.get("gtex_variant_id_b38", ""),
                "regulatory_deepening_priority": t.get("regulatory_deepening_priority", ""),
                "regulatory_deepening_score": t.get("regulatory_deepening_score", ""),
                "evidence_layers": t.get("evidence_layers", ""),
                "labels": t.get("labels", ""),
                "candidate_tier": t.get("candidate_tier", ""),
                "selected_hgvsc": t.get("selected_hgvsc", ""),
                "ccre_class": t.get("ccre_class", ""),
                "population_structure_class": t.get("population_structure_class", ""),
                "max_pairwise_fst_hudson": t.get("max_pairwise_fst_hudson", ""),
            })
            intersect_rows.append(row)

    intersect_fields = [
        "gene", "query_gencode_id", "dataset", "association_kind",
        "variantId", "target_variant_key", "target_gtex_variant_id_b38",
        "tissue", "pValue", "qValue", "effect_or_nes", "phenotype_or_gene",
        "regulatory_deepening_priority", "regulatory_deepening_score",
        "evidence_layers", "labels", "candidate_tier", "selected_hgvsc",
        "ccre_class", "population_structure_class", "max_pairwise_fst_hudson",
        "raw_json"
    ]

    write_tsv(OUT_TARGET_INTERSECT, intersect_rows, intersect_fields)
    write_tsv(TABLE12, intersect_rows[:100], intersect_fields)

    audit_fields = [
        "gene", "gencode_id", "dataset", "kind", "page",
        "http_status", "error", "n_records", "url"
    ]
    write_tsv(OUT_AUDIT, audit_rows, audit_fields)

    summary_rows = []
    summary_rows.append({
        "group": "ALL",
        "target_variants": len(targets),
        "gene_ids_resolved": len(gene_id_rows),
        "gene_level_association_records": len(association_rows),
        "target_intersection_records": len(intersect_rows),
        "association_kind_counts": ";".join(f"{k}:{v}" for k, v in sorted(Counter(r["association_kind"] for r in association_rows).items())),
        "intersection_kind_counts": ";".join(f"{k}:{v}" for k, v in sorted(Counter(r["association_kind"] for r in intersect_rows).items())),
        "audit_status_counts": ";".join(f"{k}:{v}" for k, v in sorted(Counter(str(r["http_status"]) for r in audit_rows).items())),
    })

    for gene in GENES:
        summary_rows.append({
            "group": f"gene={gene}",
            "target_variants": sum(1 for t in targets if t.get("gene") == gene),
            "gene_ids_resolved": sum(1 for r in gene_id_rows if r.get("gene") == gene),
            "gene_level_association_records": sum(1 for r in association_rows if r.get("gene") == gene),
            "target_intersection_records": sum(1 for r in intersect_rows if r.get("gene") == gene),
            "association_kind_counts": ";".join(f"{k}:{v}" for k, v in sorted(Counter(r["association_kind"] for r in association_rows if r.get("gene") == gene).items())),
            "intersection_kind_counts": ";".join(f"{k}:{v}" for k, v in sorted(Counter(r["association_kind"] for r in intersect_rows if r.get("gene") == gene).items())),
            "audit_status_counts": ";".join(f"{k}:{v}" for k, v in sorted(Counter(str(r["http_status"]) for r in audit_rows if r.get("gene") == gene).items())),
        })

    summary_fields = [
        "group", "target_variants", "gene_ids_resolved",
        "gene_level_association_records", "target_intersection_records",
        "association_kind_counts", "intersection_kind_counts", "audit_status_counts"
    ]
    write_tsv(OUT_SUMMARY, summary_rows, summary_fields)

    if intersect_rows:
        top_lines = []
        for r in intersect_rows[:15]:
            top_lines.append(
                f"- {r['gene']} {r['target_variant_key']} matched {r['association_kind']} in {r.get('tissue') or 'unknown_tissue'} "
                f"(dataset={r['dataset']}, p={r.get('pValue', '')})."
            )
        result_text = f"""The gene-level GTEx retrieval produced {len(association_rows)} HRH2/HRH4 association records, of which {len(intersect_rows)} intersected the regulatory-deepening target list.

## Target intersections

{chr(10).join(top_lines)}
"""
    else:
        result_text = f"""The gene-level GTEx retrieval produced {len(association_rows)} HRH2/HRH4 association records, but none directly matched the current HRH2/HRH4 regulatory-deepening target IDs. This should be interpreted carefully: it means that the current prioritized variants were not recovered as precomputed significant GTEx single-tissue eQTL/sQTL/iQTL/isQTL hits through this API workflow, not that the loci lack regulatory biology. Gene-level association records, query audit files and raw API responses were retained for traceability.
"""

    md = f"""# Gene-level GTEx eQTL/sQTL retrieval and intersection with HRH2/HRH4 targets

A gene-first GTEx strategy was used after the initial variant-first probe returned no extracted target records. GTEx reference gene queries were used to resolve HRH2/HRH4 identifiers, then significant single-tissue eQTL, sQTL, interaction eQTL and interaction sQTL records were queried by gene and intersected locally with the HRH2/HRH4 regulatory-deepening target list.

{result_text}

## Output files

- Resolved GTEx gene identifiers: `06_regulatory_genomics/gtex_eqtl_sqtl_lookup/HRH2_HRH4_GTEx_reference_gene_ids_v1_v2_ncbi.tsv`
- Gene-level GTEx association records: `06_regulatory_genomics/gtex_eqtl_sqtl_lookup/HRH2_HRH4_GTEx_gene_level_eQTL_sQTL_records_v1_v2_ncbi.tsv`
- Target intersections: `06_regulatory_genomics/gtex_eqtl_sqtl_lookup/HRH2_HRH4_GTEx_target_intersections_v1_v2_ncbi.tsv`
- Query audit: `06_regulatory_genomics/gtex_eqtl_sqtl_lookup/HRH2_HRH4_GTEx_gene_level_query_audit_v1_v2_ncbi.tsv`
- Manuscript-ready Table 12: `13_tables/manuscript_ready/Table_12_HRH2_HRH4_GTEx_gene_level_eQTL_sQTL_intersections_v1_v2_ncbi.tsv`
"""
    RESULTS_MD.write_text(md, encoding="utf-8")

    print(f"Wrote GTEx gene IDs: {OUT_GENE_IDS}")
    print(f"Wrote gene-level associations: {OUT_GENE_LEVEL}")
    print(f"Wrote target intersections: {OUT_TARGET_INTERSECT}")
    print(f"Wrote query audit: {OUT_AUDIT}")
    print(f"Wrote summary: {OUT_SUMMARY}")
    print(f"Wrote Table 12: {TABLE12}")
    print(f"Wrote Results subsection: {RESULTS_MD}")
    print()
    print(f"Targets: {len(targets)}")
    print(f"Resolved gene ID rows: {len(gene_id_rows)}")
    print(f"Gene-level association records: {len(association_rows)}")
    print(f"Target intersection records: {len(intersect_rows)}")
    print("Association kind counts:")
    for k, v in sorted(Counter(r["association_kind"] for r in association_rows).items()):
        print(f"  {k}: {v}")
    print("Audit status counts:")
    for k, v in sorted(Counter(str(r["http_status"]) for r in audit_rows).items()):
        print(f"  {k}: {v}")

if __name__ == "__main__":
    main()
