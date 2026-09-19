#!/usr/bin/env python3

from pathlib import Path
import gzip
import csv
import re
import math
import itertools
import datetime
from collections import defaultdict

PROJECT = Path("/mnt/d/histamine_receptor_pharmacogenomics")
OUTDIR = PROJECT / "11_results/monitoring_evaluation_v6_global_functional_diversity"
OUTDIR.mkdir(parents=True, exist_ok=True)

VCF_CANDIDATES = [
    PROJECT / "04_variants/raw_vcf/1000genomes_20190312_GRCh38_v2_ncbi_nochr/hrh_exact_1000g_GRCh38_v2_ncbi_allchr.vcf.gz",
]

GENE_BED = PROJECT / "04_variants/gene_slices/hrh_genes_exact_GRCh38_v2_ncbi_nochr.bed"

ANNOTATION_CANDIDATES = [
    PROJECT / "11_results/candidate_prioritization/hrh_integrated_master_candidate_table_v2_ncbi.tsv",
    PROJECT / "11_results/candidate_prioritization/hrh_exact_candidate_prioritization_v2_ncbi.tsv",
    PROJECT / "11_results/candidate_prioritization/hrh_protein_altering_candidates_canonical_safe_gpcr_topology_v2_ncbi.tsv",
    PROJECT / "11_results/candidate_prioritization/hrh_regulatory_utr_splice_candidates_v2_ncbi.tsv",
    PROJECT / "11_results/monitoring_evaluation_v4_precision/hrh_precision_receptor_axis_summary_v4.tsv",
]

OUT_MATRIX = OUTDIR / "HRH_population_functional_variant_matrix_v6.tsv"
OUT_PAIRWISE = OUTDIR / "HRH_population_pairwise_variant_tests_v6.tsv"
OUT_TOP = OUTDIR / "HRH_top_population_differentiated_functional_variants_v6.tsv"
OUT_RECEPTOR_SUMMARY = OUTDIR / "HRH_receptor_population_diversity_summary_v6.tsv"
OUT_SUPERPOP_SUMMARY = OUTDIR / "HRH_superpopulation_variant_burden_summary_v6.tsv"
OUT_METHODS = OUTDIR / "HRH_population_functional_variant_matrix_methods_note_v6.md"
OUT_REPORT = OUTDIR / "HRH_population_functional_variant_matrix_report_v6.md"
OUT_CHECKLIST = OUTDIR / "HRH_population_functional_variant_matrix_checklist_v6.tsv"
OUT_PANEL_USED = OUTDIR / "HRH_1000G_sample_panel_used_v6.tsv"

SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
PAIRWISE = list(itertools.combinations(SUPERPOPS, 2))

GENE_INTERVALS_FALLBACK = [
    # chrom, start, end, gene, strand
    ("3", 11137238, 11263557, "HRH1", "+"),
    ("5", 175658071, 175710756, "HRH2", "+"),
    ("20", 62214960, 62220278, "HRH3", "-"),
    ("18", 24460637, 24479974, "HRH4", "+"),
]

def norm_chrom(x):
    return str(x).replace("chr", "").replace("NC_000003.12", "3").replace("NC_000005.10", "5").replace("NC_000020.11", "20").replace("NC_000018.10", "18")

def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")

def find_vcf():
    for p in VCF_CANDIDATES:
        if p.exists():
            return p
    vcfs = sorted(PROJECT.rglob("*hrh*exact*1000g*.vcf.gz"))
    if vcfs:
        return vcfs[0]
    vcfs = sorted(PROJECT.rglob("*.vcf.gz"))
    for p in vcfs:
        if "1000" in str(p).lower() and "hrh" in str(p).lower():
            return p
    raise SystemExit("No HRH 1000G VCF found.")

def find_sample_panel():
    candidates = []
    patterns = [
        "*integrated_call_samples*.panel*",
        "*1000*panel*",
        "*sample*panel*",
        "*samples*.tsv",
        "*samples*.csv",
        "*panel*.tsv",
        "*panel*.txt",
    ]
    for pat in patterns:
        candidates.extend(PROJECT.rglob(pat))
    scored = []
    for p in sorted(set(candidates)):
        if not p.is_file():
            continue
        name = str(p).lower()
        score = 0
        if "integrated_call_samples" in name:
            score += 20
        if "panel" in name:
            score += 10
        if "1000" in name or "1kg" in name:
            score += 5
        if p.suffix.lower() in {".panel", ".tsv", ".txt", ".csv"}:
            score += 2
        scored.append((score, p))
    scored.sort(reverse=True)
    for _, p in scored:
        panel = try_read_panel(p)
        if panel and len(panel) > 500:
            return p, panel
    raise SystemExit("Could not locate a usable 1000G sample panel with sample-to-superpopulation mapping.")

def split_line(line):
    if "\t" in line:
        return line.rstrip("\n").split("\t")
    if "," in line:
        return [x.strip() for x in line.rstrip("\n").split(",")]
    return line.strip().split()

def try_read_panel(path):
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except Exception:
        return None
    if not lines:
        return None

    first = split_line(lines[0])
    lower = [x.lower() for x in first]
    has_header = any(x in lower for x in ["sample", "sample_id", "super_pop", "superpopulation", "superpopulation code", "pop", "population"])

    rows = []
    if has_header:
        header = lower
        data_lines = lines[1:]
        def idx(options):
            for o in options:
                if o in header:
                    return header.index(o)
            return None
        sample_i = idx(["sample", "sample_id", "sample name", "id"])
        pop_i = idx(["pop", "population", "population code"])
        super_i = idx(["super_pop", "superpopulation", "superpopulation code", "super population"])
    else:
        data_lines = lines
        sample_i = 0
        pop_i = 1 if len(first) > 1 else None
        super_i = 2 if len(first) > 2 else None

    if sample_i is None or super_i is None:
        return None

    panel = {}
    for line in data_lines:
        if not line.strip() or line.startswith("#"):
            continue
        parts = split_line(line)
        if len(parts) <= max(sample_i, super_i):
            continue
        sample = parts[sample_i]
        pop = parts[pop_i] if pop_i is not None and len(parts) > pop_i else ""
        superpop = parts[super_i].upper()
        if superpop in SUPERPOPS:
            panel[sample] = {"pop": pop, "superpop": superpop}
    return panel

def load_gene_intervals():
    intervals = []
    if GENE_BED.exists():
        with open(GENE_BED, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                chrom = norm_chrom(parts[0])
                start = int(parts[1]) + 1  # BED 0-based to 1-based
                end = int(parts[2])
                gene = parts[3].split("|")[0]
                strand = parts[5] if len(parts) >= 6 else ""
                intervals.append((chrom, start, end, gene, strand))
    if not intervals:
        intervals = GENE_INTERVALS_FALLBACK
    return intervals

def gene_for_variant(chrom, pos, intervals):
    chrom = norm_chrom(chrom)
    for c, start, end, gene, strand in intervals:
        if norm_chrom(c) == chrom and start <= pos <= end:
            return gene
    return "UNKNOWN"

def variant_keys(chrom, pos, ref, alt, vid=""):
    c = norm_chrom(chrom)
    keys = {
        f"{c}:{pos}:{ref}:{alt}",
        f"chr{c}:{pos}:{ref}:{alt}",
        f"{c}_{pos}_{ref}_{alt}",
        f"chr{c}_{pos}_{ref}_{alt}",
        f"{c}:{pos}",
        f"chr{c}:{pos}",
    }
    if vid and vid != ".":
        keys.add(vid)
    return keys

def detect_col(header, options):
    lower = [h.lower() for h in header]
    for opt in options:
        if opt.lower() in lower:
            return header[lower.index(opt.lower())]
    for h in header:
        lh = h.lower()
        if any(opt.lower() in lh for opt in options):
            return h
    return None

def load_annotations():
    annotations = {}
    sources = []
    for path in ANNOTATION_CANDIDATES:
        if not path.exists():
            continue
        try:
            with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                if not reader.fieldnames:
                    continue
                header = reader.fieldnames
                chrom_col = detect_col(header, ["chrom", "chr", "#chrom"])
                pos_col = detect_col(header, ["pos", "position", "start"])
                ref_col = detect_col(header, ["ref", "reference"])
                alt_col = detect_col(header, ["alt", "alternate"])
                var_col = detect_col(header, ["variant_id", "variant", "variant_key", "id"])
                gene_col = detect_col(header, ["gene", "receptor"])
                consequence_col = detect_col(header, ["consequence", "variant_class", "annotation"])
                tier_col = detect_col(header, ["tier", "priority", "candidate_tier", "integrated_tier", "regulatory_tier"])
                topology_col = detect_col(header, ["topology", "domain", "gpcr_region"])
                axis_col = detect_col(header, ["axis", "evidence", "evidence_axis", "candidate_axis"])
                for row in reader:
                    keys = set()
                    if chrom_col and pos_col and ref_col and alt_col:
                        try:
                            keys |= variant_keys(row[chrom_col], int(float(row[pos_col])), row[ref_col], row[alt_col], row.get(var_col, "") if var_col else "")
                        except Exception:
                            pass
                    if var_col and row.get(var_col):
                        keys.add(row[var_col])
                    if not keys:
                        continue
                    ann = {
                        "annotation_source": str(path.relative_to(PROJECT)),
                        "ann_gene": row.get(gene_col, "") if gene_col else "",
                        "ann_consequence": row.get(consequence_col, "") if consequence_col else "",
                        "ann_tier": row.get(tier_col, "") if tier_col else "",
                        "ann_topology_or_domain": row.get(topology_col, "") if topology_col else "",
                        "ann_axis_or_evidence": row.get(axis_col, "") if axis_col else "",
                    }
                    for k in keys:
                        if k not in annotations:
                            annotations[k] = ann
                sources.append(str(path.relative_to(PROJECT)))
        except Exception:
            continue
    return annotations, sources

def parse_gt(gt_field):
    gt = gt_field.split(":", 1)[0]
    if gt in {".", "./.", ".|."}:
        return None
    alleles = re.split(r"[\/|]", gt)
    alt_count = 0
    called = 0
    for a in alleles:
        if a == ".":
            continue
        called += 1
        try:
            if int(a) > 0:
                alt_count += 1
        except Exception:
            pass
    if called == 0:
        return None
    return alt_count, called

def fisher_exact_two_sided(a, b, c, d):
    # 2x2 table:
    # group1 alt=a ref=b; group2 alt=c ref=d
    # exact two-sided probability using hypergeometric enumeration.
    n = a + b + c + d
    if n == 0:
        return 1.0
    row1 = a + b
    col1 = a + c
    min_x = max(0, row1 - (n - col1))
    max_x = min(row1, col1)

    def logchoose(n, k):
        if k < 0 or k > n:
            return float("-inf")
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    def hyperprob(x):
        return math.exp(logchoose(col1, x) + logchoose(n - col1, row1 - x) - logchoose(n, row1))

    p_obs = hyperprob(a)
    p = 0.0
    for x in range(min_x, max_x + 1):
        px = hyperprob(x)
        if px <= p_obs + 1e-15:
            p += px
    return min(1.0, p)

def bh_fdr(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [1.0] * n
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = n - rank + 1
        val = pvals[idx] * n / i
        prev = min(prev, val)
        q[idx] = min(prev, 1.0)
    return q

def fmt(x, nd=6):
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.{nd}g}"
    return str(x)

def main():
    vcf = find_vcf()
    panel_path, panel = find_sample_panel()
    intervals = load_gene_intervals()
    annotations, annotation_sources = load_annotations()

    sample_names = []
    sample_superpop = {}
    rows = []
    pair_rows = []

    with open_text(vcf) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_names = header[9:]
                for s in sample_names:
                    if s in panel:
                        sample_superpop[s] = panel[s]["superpop"]
                break

        if not sample_names:
            raise SystemExit("Could not read VCF #CHROM header.")

        pop_sample_counts = {sp: sum(1 for s in sample_names if sample_superpop.get(s) == sp) for sp in SUPERPOPS}

        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos_s, vid, ref, alt_s = parts[0], parts[1], parts[2], parts[3], parts[4]
            try:
                pos = int(pos_s)
            except Exception:
                continue
            # Only first ALT if multiallelic. The HRH VCF should mostly be normalized; flag if multi.
            alts = alt_s.split(",")
            alt = alts[0]
            multiallelic = "YES" if len(alts) > 1 else "NO"
            gene = gene_for_variant(chrom, pos, intervals)

            counts = {sp: {"alt": 0, "called": 0, "samples_called": 0} for sp in SUPERPOPS}
            for s, gtfield in zip(sample_names, parts[9:]):
                sp = sample_superpop.get(s)
                if sp not in counts:
                    continue
                parsed = parse_gt(gtfield)
                if parsed is None:
                    continue
                alt_count, called = parsed
                counts[sp]["alt"] += alt_count
                counts[sp]["called"] += called
                counts[sp]["samples_called"] += 1

            af = {}
            for sp in SUPERPOPS:
                af[sp] = counts[sp]["alt"] / counts[sp]["called"] if counts[sp]["called"] else None

            valid_af = [x for x in af.values() if x is not None]
            max_af = max(valid_af) if valid_af else None
            min_af = min(valid_af) if valid_af else None
            delta = (max_af - min_af) if valid_af else None
            dominant_pop = max((sp for sp in SUPERPOPS if af[sp] is not None), key=lambda sp: af[sp], default="")
            lowest_pop = min((sp for sp in SUPERPOPS if af[sp] is not None), key=lambda sp: af[sp], default="")

            # Approximate global population-structure score, not a formal Weir-Cockerham estimator.
            pbar = sum(valid_af) / len(valid_af) if valid_af else None
            fst_like = None
            if pbar is not None and 0 < pbar < 1 and valid_af:
                fst_like = sum((p - pbar) ** 2 for p in valid_af) / len(valid_af) / (pbar * (1 - pbar))

            keyset = variant_keys(chrom, pos, ref, alt, vid)
            ann = {}
            for k in keyset:
                if k in annotations:
                    ann = annotations[k]
                    break

            variant_key = f"{norm_chrom(chrom)}:{pos}:{ref}:{alt}"

            row = {
                "variant_key": variant_key,
                "chrom": norm_chrom(chrom),
                "pos": str(pos),
                "id": vid,
                "ref": ref,
                "alt": alt,
                "multiallelic_flag": multiallelic,
                "gene": gene,
                "annotation_source": ann.get("annotation_source", ""),
                "ann_gene": ann.get("ann_gene", ""),
                "ann_consequence": ann.get("ann_consequence", ""),
                "ann_tier": ann.get("ann_tier", ""),
                "ann_topology_or_domain": ann.get("ann_topology_or_domain", ""),
                "ann_axis_or_evidence": ann.get("ann_axis_or_evidence", ""),
                "AFR_AF": fmt(af["AFR"]),
                "AMR_AF": fmt(af["AMR"]),
                "EAS_AF": fmt(af["EAS"]),
                "EUR_AF": fmt(af["EUR"]),
                "SAS_AF": fmt(af["SAS"]),
                "AFR_alt_alleles": str(counts["AFR"]["alt"]),
                "AMR_alt_alleles": str(counts["AMR"]["alt"]),
                "EAS_alt_alleles": str(counts["EAS"]["alt"]),
                "EUR_alt_alleles": str(counts["EUR"]["alt"]),
                "SAS_alt_alleles": str(counts["SAS"]["alt"]),
                "AFR_called_alleles": str(counts["AFR"]["called"]),
                "AMR_called_alleles": str(counts["AMR"]["called"]),
                "EAS_called_alleles": str(counts["EAS"]["called"]),
                "EUR_called_alleles": str(counts["EUR"]["called"]),
                "SAS_called_alleles": str(counts["SAS"]["called"]),
                "max_population": dominant_pop,
                "min_population": lowest_pop,
                "max_minus_min_AF": fmt(delta),
                "fst_like_population_structure_score": fmt(fst_like),
            }
            rows.append(row)

            for sp1, sp2 in PAIRWISE:
                a = counts[sp1]["alt"]
                b = counts[sp1]["called"] - a
                c = counts[sp2]["alt"]
                d = counts[sp2]["called"] - c
                p1 = af[sp1]
                p2 = af[sp2]
                if p1 is None or p2 is None:
                    pval = 1.0
                    d_af = None
                else:
                    pval = fisher_exact_two_sided(a, b, c, d)
                    d_af = p1 - p2
                pair_rows.append({
                    "variant_key": variant_key,
                    "gene": gene,
                    "id": vid,
                    "ref": ref,
                    "alt": alt,
                    "comparison": f"{sp1}_vs_{sp2}",
                    "pop1": sp1,
                    "pop2": sp2,
                    "pop1_af": fmt(p1),
                    "pop2_af": fmt(p2),
                    "delta_af_pop1_minus_pop2": fmt(d_af),
                    "abs_delta_af": fmt(abs(d_af) if d_af is not None else None),
                    "pop1_alt": str(a),
                    "pop1_ref": str(b),
                    "pop2_alt": str(c),
                    "pop2_ref": str(d),
                    "fisher_p": pval,
                })

    # FDR correct all pairwise tests across all variants and comparisons.
    pvals = [float(r["fisher_p"]) for r in pair_rows]
    qvals = bh_fdr(pvals)
    for r, q in zip(pair_rows, qvals):
        r["fisher_p"] = fmt(float(r["fisher_p"]))
        r["fdr_q"] = fmt(q)
        abs_delta = float(r["abs_delta_af"]) if r["abs_delta_af"] else 0.0
        r["significant_fdr_0_05_delta_0_10"] = "YES" if q <= 0.05 and abs_delta >= 0.10 else "NO"
        r["significant_fdr_0_05_delta_0_20"] = "YES" if q <= 0.05 and abs_delta >= 0.20 else "NO"

    # Add best pairwise stats back to matrix.
    best_by_variant = {}
    sig_counts = defaultdict(int)
    for r in pair_rows:
        key = r["variant_key"]
        abs_delta = float(r["abs_delta_af"]) if r["abs_delta_af"] else -1
        q = float(r["fdr_q"]) if r["fdr_q"] else 1
        if key not in best_by_variant:
            best_by_variant[key] = r
        else:
            old = best_by_variant[key]
            old_abs = float(old["abs_delta_af"]) if old["abs_delta_af"] else -1
            if abs_delta > old_abs:
                best_by_variant[key] = r
        if r["significant_fdr_0_05_delta_0_10"] == "YES":
            sig_counts[key] += 1

    for row in rows:
        b = best_by_variant.get(row["variant_key"], {})
        row["top_pairwise_comparison"] = b.get("comparison", "")
        row["top_pairwise_abs_delta_af"] = b.get("abs_delta_af", "")
        row["top_pairwise_fdr_q"] = b.get("fdr_q", "")
        row["significant_pairwise_count_fdr_0_05_delta_0_10"] = str(sig_counts.get(row["variant_key"], 0))
        row["population_differentiation_class"] = (
            "strong" if sig_counts.get(row["variant_key"], 0) and float(row.get("top_pairwise_abs_delta_af") or 0) >= 0.30 else
            "moderate" if sig_counts.get(row["variant_key"], 0) and float(row.get("top_pairwise_abs_delta_af") or 0) >= 0.20 else
            "signal" if sig_counts.get(row["variant_key"], 0) else
            "background"
        )

    # Write outputs.
    matrix_fields = list(rows[0].keys()) if rows else []
    pair_fields = list(pair_rows[0].keys()) if pair_rows else []
    with open(OUT_MATRIX, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=matrix_fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    with open(OUT_PAIRWISE, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=pair_fields)
        w.writeheader()
        for r in pair_rows:
            w.writerow(r)

    top_rows = sorted(
        rows,
        key=lambda r: (
            int(r.get("significant_pairwise_count_fdr_0_05_delta_0_10") or 0),
            float(r.get("top_pairwise_abs_delta_af") or 0),
            float(r.get("fst_like_population_structure_score") or 0),
        ),
        reverse=True,
    )[:250]
    with open(OUT_TOP, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=matrix_fields)
        w.writeheader()
        for r in top_rows:
            w.writerow(r)

    # Receptor summary.
    receptor_rows = []
    for gene in ["HRH1", "HRH2", "HRH3", "HRH4", "UNKNOWN"]:
        subset = [r for r in rows if r["gene"] == gene]
        if not subset:
            continue
        sig = [r for r in subset if int(r["significant_pairwise_count_fdr_0_05_delta_0_10"]) > 0]
        strong = [r for r in subset if r["population_differentiation_class"] == "strong"]
        receptor_rows.append({
            "gene": gene,
            "variant_count": str(len(subset)),
            "population_significant_variants_fdr_0_05_delta_0_10": str(len(sig)),
            "strong_population_differentiated_variants": str(len(strong)),
            "max_abs_delta_af": fmt(max(float(r["top_pairwise_abs_delta_af"] or 0) for r in subset)),
            "top_variant": max(subset, key=lambda r: float(r["top_pairwise_abs_delta_af"] or 0))["variant_key"],
            "top_comparison": max(subset, key=lambda r: float(r["top_pairwise_abs_delta_af"] or 0))["top_pairwise_comparison"],
        })
    write_simple_tsv(OUT_RECEPTOR_SUMMARY, receptor_rows)

    # Superpopulation burden summary: number of variants where pop is max/min.
    super_rows = []
    for sp in SUPERPOPS:
        super_rows.append({
            "superpopulation": sp,
            "sample_count": str(pop_sample_counts.get(sp, 0)),
            "variants_where_population_has_max_af": str(sum(1 for r in rows if r["max_population"] == sp)),
            "variants_where_population_has_min_af": str(sum(1 for r in rows if r["min_population"] == sp)),
            "population_significant_variants_where_max_af": str(sum(1 for r in rows if r["max_population"] == sp and int(r["significant_pairwise_count_fdr_0_05_delta_0_10"]) > 0)),
        })
    write_simple_tsv(OUT_SUPERPOP_SUMMARY, super_rows)

    with open(OUT_PANEL_USED, "w", encoding="utf-8", newline="") as f:
        fields = ["sample", "pop", "superpop", "in_vcf"]
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        vcf_samples = set(sample_names)
        for s, meta in sorted(panel.items()):
            w.writerow({"sample": s, "pop": meta["pop"], "superpop": meta["superpop"], "in_vcf": "YES" if s in vcf_samples else "NO"})

    checklist_rows = [
        {"check": "vcf_found", "status": "PASS" if vcf.exists() else "REVIEW", "observed": str(vcf)},
        {"check": "sample_panel_found", "status": "PASS" if panel_path.exists() else "REVIEW", "observed": str(panel_path)},
        {"check": "vcf_samples_mapped_to_superpop", "status": "PASS" if len(sample_superpop) > 2000 else "REVIEW", "observed": str(len(sample_superpop))},
        {"check": "variant_rows_ge_5000", "status": "PASS" if len(rows) >= 5000 else "REVIEW", "observed": str(len(rows))},
        {"check": "pairwise_rows_expected", "status": "PASS" if len(pair_rows) == len(rows) * len(PAIRWISE) else "REVIEW", "observed": f"{len(pair_rows)} vs {len(rows)*len(PAIRWISE)}"},
        {"check": "receptor_summary_created", "status": "PASS" if OUT_RECEPTOR_SUMMARY.exists() else "REVIEW", "observed": str(OUT_RECEPTOR_SUMMARY)},
        {"check": "top_variant_table_created", "status": "PASS" if OUT_TOP.exists() else "REVIEW", "observed": str(OUT_TOP)},
    ]
    write_simple_tsv(OUT_CHECKLIST, checklist_rows)

    methods = f"""# HRH population-functional variant matrix methods note v6

Generated: {datetime.datetime.now().isoformat(timespec='seconds')}

## Input data

- VCF: `{vcf}`
- Sample panel: `{panel_path}`
- Gene intervals: `{GENE_BED if GENE_BED.exists() else 'fallback built-in GRCh38 intervals'}`
- Annotation sources detected:
{chr(10).join('- `' + s + '`' for s in annotation_sources) if annotation_sources else '- No candidate annotation tables detected'}

## Analysis

For each HRH1-HRH4 variant in the 1000 Genomes-derived VCF, the script computed alternate-allele frequency in AFR, AMR, EAS, EUR and SAS superpopulations. Pairwise superpopulation contrasts were evaluated using a two-sided Fisher exact test on alternate and reference allele counts. P values were corrected across all variant-comparison tests using Benjamini-Hochberg FDR.

A variant was flagged as population-significant when pairwise FDR q <= 0.05 and absolute allele-frequency difference >= 0.10. A stronger threshold of absolute difference >= 0.20 was also recorded. A simple FST-like population-structure score was calculated as between-superpopulation allele-frequency variance divided by pbar(1-pbar). This is used as a descriptive prioritization score, not as a substitute for a formal Weir-Cockerham estimator.

## Output tables

- Full matrix: `{OUT_MATRIX}`
- Pairwise tests: `{OUT_PAIRWISE}`
- Top population-differentiated variants: `{OUT_TOP}`
- Receptor summary: `{OUT_RECEPTOR_SUMMARY}`
- Superpopulation burden summary: `{OUT_SUPERPOP_SUMMARY}`
"""
    OUT_METHODS.write_text(methods, encoding="utf-8")

    report = []
    report.append("# HRH population-functional variant matrix report v6")
    report.append("")
    report.append(f"Generated: {datetime.datetime.now().isoformat(timespec='seconds')}")
    report.append("")
    report.append("## Headline results")
    report.append("")
    report.append(f"- Variants analyzed: {len(rows)}")
    report.append(f"- Pairwise population tests: {len(pair_rows)}")
    report.append(f"- VCF samples mapped to superpopulations: {len(sample_superpop)}")
    report.append(f"- Annotation sources detected: {len(annotation_sources)}")
    report.append("")
    report.append("## Receptor summary")
    report.append("")
    report.append("| Gene | Variants | Population-significant | Strong population-differentiated | Max abs delta AF | Top variant | Top comparison |")
    report.append("|---|---:|---:|---:|---:|---|---|")
    for r in receptor_rows:
        report.append(f"| {r['gene']} | {r['variant_count']} | {r['population_significant_variants_fdr_0_05_delta_0_10']} | {r['strong_population_differentiated_variants']} | {r['max_abs_delta_af']} | {r['top_variant']} | {r['top_comparison']} |")
    report.append("")
    report.append("## Superpopulation burden summary")
    report.append("")
    report.append("| Superpopulation | Sample count | Max-AF variants | Min-AF variants | Significant variants where max AF |")
    report.append("|---|---:|---:|---:|---:|")
    for r in super_rows:
        report.append(f"| {r['superpopulation']} | {r['sample_count']} | {r['variants_where_population_has_max_af']} | {r['variants_where_population_has_min_af']} | {r['population_significant_variants_where_max_af']} |")
    report.append("")
    report.append("## Interpretation")
    report.append("")
    report.append("This table provides the statistical backbone for the global pharmacogenomic diversity extension of the manuscript. It should be used to identify which receptor-specific functional candidates are population-differentiated and to prioritize ligand-pocket, regulatory, topology and receptor-fate subsets for the next modules.")
    report.append("")
    report.append(f"Full matrix: `{OUT_MATRIX}`")
    report.append(f"Pairwise tests: `{OUT_PAIRWISE}`")
    report.append(f"Top variants: `{OUT_TOP}`")
    report.append(f"Checklist: `{OUT_CHECKLIST}`")
    OUT_REPORT.write_text("\n".join(report) + "\n", encoding="utf-8")

    print("\n".join(report))

def write_simple_tsv(path, rows):
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)

if __name__ == "__main__":
    main()
