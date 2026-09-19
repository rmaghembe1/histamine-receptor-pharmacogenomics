#!/usr/bin/env python3

from pathlib import Path
import gzip
import csv
import math
from collections import defaultdict, Counter
from itertools import combinations

PROJECT = Path("/mnt/d/histamine_receptor_pharmacogenomics")

VCF = PROJECT / "04_variants/raw_vcf/1000genomes_20190312_GRCh38_v2_ncbi_nochr/hrh_exact_1000g_GRCh38_v2_ncbi_allchr.vcf.gz"
BED = PROJECT / "04_variants/gene_slices/hrh_genes_exact_GRCh38_v2_ncbi_nochr.bed"
PANEL = PROJECT / "01_data_raw/1000genomes_panel/integrated_call_samples_v3.20130502.ALL.panel"
INTEGRATED = PROJECT / "11_results/candidate_prioritization/hrh_integrated_master_candidate_table_v2_ncbi.tsv"

OUTDIR = PROJECT / "05_population_genetics/formal_population_structure_v1_v2_ncbi"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_VARIANTS = OUTDIR / "hrh_exact_variant_population_structure_v1_v2_ncbi.tsv"
OUT_GENE_SUMMARY = OUTDIR / "hrh_gene_population_structure_summary_v1_v2_ncbi.tsv"
OUT_PAIRWISE_SUMMARY = OUTDIR / "hrh_pairwise_superpopulation_fst_summary_v1_v2_ncbi.tsv"
OUT_ENRICHED = OUTDIR / "hrh_population_enriched_candidate_variants_v1_v2_ncbi.tsv"
OUT_TOP = OUTDIR / "hrh_top_population_structure_candidates_v1_v2_ncbi.tsv"
OUT_OVERLAY = OUTDIR / "hrh_population_genetics_candidate_overlay_v1_v2_ncbi.tsv"

TABLE6 = PROJECT / "13_tables/manuscript_ready/Table_6_HRH_formal_population_structure_summary_v1_v2_ncbi.tsv"
TABLE7 = PROJECT / "13_tables/manuscript_ready/Table_7_HRH_population_enriched_candidates_v1_v2_ncbi.tsv"
RESULTS_MD = PROJECT / "14_manuscript/results_sections/results_formal_population_genetics_fst_entropy_v1_v2_ncbi.md"

SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]

def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")

def load_bed(path):
    intervals = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, name = line.rstrip("\n").split("\t")[:4]
            gene = name.split("|")[0]
            intervals.append((chrom.replace("chr", ""), int(start), int(end), gene))
    return intervals

def assign_gene(chrom, pos1, intervals):
    c = chrom.replace("chr", "")
    pos0 = int(pos1) - 1
    for b_chrom, start, end, gene in intervals:
        if b_chrom == c and start <= pos0 < end:
            return gene
    return ""

def load_panel(path):
    sample_to_super = {}
    sample_to_pop = {}

    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split()
        lower = [h.lower() for h in header]

        # Known format: sample pop super_pop gender
        if "sample" in lower and ("super_pop" in lower or "superpop" in lower):
            sample_i = lower.index("sample")
            pop_i = lower.index("pop") if "pop" in lower else None
            sp_i = lower.index("super_pop") if "super_pop" in lower else lower.index("superpop")

            for line in f:
                if not line.strip():
                    continue
                parts = line.strip().split()
                if len(parts) <= max(sample_i, sp_i):
                    continue
                sample = parts[sample_i]
                sp = parts[sp_i]
                pop = parts[pop_i] if pop_i is not None and len(parts) > pop_i else ""
                if sp in SUPERPOPS:
                    sample_to_super[sample] = sp
                    sample_to_pop[sample] = pop
        else:
            # fallback: treat first line as data if no header was detected
            rows = [header] + [line.strip().split() for line in f if line.strip()]
            for parts in rows:
                if len(parts) < 3:
                    continue
                sample, pop, sp = parts[0], parts[1], parts[2]
                if sp in SUPERPOPS:
                    sample_to_super[sample] = sp
                    sample_to_pop[sample] = pop

    return sample_to_super, sample_to_pop

def read_tsv(path):
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def load_candidate_overlay(path):
    rows = read_tsv(path)
    store = defaultdict(lambda: defaultdict(set))

    for r in rows:
        vkey = r.get("variant_key") or r.get("flanking_variant_key") or r.get("best_exact_anchor_variant_key")
        if not vkey:
            continue
        for field in [
            "gene", "uploaded_gene", "source_layer", "functional_axis",
            "manuscript_status", "candidate_tier", "selected_hgvsc",
            "selected_hgvsp", "gpcr_region_broad", "gpcr_motif_neighborhood"
        ]:
            val = r.get(field, "")
            if val not in ["", None]:
                store[vkey][field].add(str(val))

    out = {}
    for vkey, d in store.items():
        out[vkey] = {k: ";".join(sorted(v)) for k, v in d.items()}
    return out

def parse_gt(gt):
    if gt in [".", "./.", ".|."]:
        return None
    alleles = gt.replace("|", "/").split("/")
    vals = []
    for a in alleles:
        if a == ".":
            continue
        try:
            vals.append(int(a))
        except ValueError:
            continue
    if not vals:
        return None
    return vals

def entropy_from_values(values):
    vals = [v for v in values if v > 0]
    if not vals:
        return 0.0
    total = sum(vals)
    probs = [v / total for v in vals if total > 0]
    if len(probs) <= 1:
        return 0.0
    h = -sum(p * math.log(p) for p in probs)
    return h / math.log(len(values))

def hudson_fst(p1, n1, p2, n2):
    if n1 <= 1 or n2 <= 1:
        return ""
    pi_between = p1 * (1 - p2) + p2 * (1 - p1)
    if pi_between <= 0:
        return 0.0
    pi_within1 = 2 * p1 * (1 - p1) * n1 / (n1 - 1)
    pi_within2 = 2 * p2 * (1 - p2) * n2 / (n2 - 1)
    fst = (pi_between - ((pi_within1 + pi_within2) / 2)) / pi_between
    if fst < 0:
        fst = 0.0
    return fst

def classify_population_signal(global_af, afs, max_fst, entropy, top_sp, second_af):
    max_af = max(afs.values()) if afs else 0.0
    min_af = min(afs.values()) if afs else 0.0
    af_range = max_af - min_af

    if max_af >= 0.01 and second_af <= 0.001:
        return "private_or_near_private_superpopulation_signal"
    if max_fst >= 0.25 or af_range >= 0.60:
        return "very_strong_population_structure"
    if max_fst >= 0.10 or af_range >= 0.35:
        return "strong_population_structure"
    if max_fst >= 0.05 or af_range >= 0.20:
        return "moderate_population_structure"
    if global_af >= 0.05 and entropy >= 0.85 and af_range < 0.10:
        return "broadly_shared_common_variant"
    return "low_or_background_population_structure"

def fmt(x):
    if x == "":
        return ""
    if isinstance(x, int):
        return str(x)
    return f"{float(x):.6g}"

def main():
    if not VCF.exists():
        raise FileNotFoundError(f"Missing VCF: {VCF}")
    if not BED.exists():
        raise FileNotFoundError(f"Missing BED: {BED}")
    if not PANEL.exists():
        raise FileNotFoundError(f"Missing panel: {PANEL}")

    intervals = load_bed(BED)
    sample_to_super, sample_to_pop = load_panel(PANEL)
    overlay = load_candidate_overlay(INTEGRATED)

    variant_rows = []
    pair_rows = []
    processed = 0

    with open_text(VCF) as f:
        sample_names = []
        sample_superpops = []
        sample_indices = []

        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_names = header[9:]
                for idx, s in enumerate(sample_names):
                    sp = sample_to_super.get(s)
                    if sp in SUPERPOPS:
                        sample_indices.append(idx)
                        sample_superpops.append(sp)
                print(f"VCF samples: {len(sample_names)}")
                print(f"Samples with superpopulation metadata: {len(sample_indices)}")
                continue

            if not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
            fmt_fields = parts[8].split(":")
            if "GT" not in fmt_fields:
                continue
            gt_i = fmt_fields.index("GT")
            sample_fields = parts[9:]

            gene = assign_gene(chrom, pos, intervals)
            variant_key = f"{chrom.replace('chr', '')}:{pos}:{ref}:{alt}"

            ac = {sp: 0 for sp in SUPERPOPS}
            an = {sp: 0 for sp in SUPERPOPS}

            for sample_idx, sp in zip(sample_indices, sample_superpops):
                if sample_idx >= len(sample_fields):
                    continue
                gt = sample_fields[sample_idx].split(":")[gt_i]
                vals = parse_gt(gt)
                if vals is None:
                    continue
                for a in vals:
                    if a == 0:
                        an[sp] += 1
                    elif a == 1:
                        ac[sp] += 1
                        an[sp] += 1
                    else:
                        # biallelic resource expected; ignore non-biallelic alleles if encountered
                        pass

            af = {sp: (ac[sp] / an[sp] if an[sp] else 0.0) for sp in SUPERPOPS}
            global_ac = sum(ac.values())
            global_an = sum(an.values())
            global_af = global_ac / global_an if global_an else 0.0

            af_values = list(af.values())
            ac_values = [ac[sp] for sp in SUPERPOPS]

            af_entropy = entropy_from_values(af_values)
            ac_entropy = entropy_from_values(ac_values)
            af_range = max(af_values) - min(af_values)

            sorted_af = sorted(af.items(), key=lambda x: x[1], reverse=True)
            top_sp, top_af = sorted_af[0]
            second_sp, second_af = sorted_af[1] if len(sorted_af) > 1 else ("", 0.0)

            pair_fsts = {}
            for sp1, sp2 in combinations(SUPERPOPS, 2):
                fst = hudson_fst(af[sp1], an[sp1], af[sp2], an[sp2])
                if fst == "":
                    continue
                pair = f"{sp1}_vs_{sp2}"
                pair_fsts[pair] = fst
                pair_rows.append({
                    "gene": gene,
                    "variant_key": variant_key,
                    "pair": pair,
                    "fst_hudson": fst,
                    "af_1": af[sp1],
                    "af_2": af[sp2],
                    "n_alleles_1": an[sp1],
                    "n_alleles_2": an[sp2],
                })

            if pair_fsts:
                max_pair, max_fst = max(pair_fsts.items(), key=lambda x: x[1])
                mean_fst = sum(pair_fsts.values()) / len(pair_fsts)
            else:
                max_pair, max_fst, mean_fst = "", 0.0, 0.0

            pop_class = classify_population_signal(global_af, af, max_fst, af_entropy, top_sp, second_af)

            row = {
                "gene": gene,
                "variant_key": variant_key,
                "chrom": chrom.replace("chr", ""),
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "global_AC": global_ac,
                "global_AN": global_an,
                "global_AF": global_af,
                "AFR_AC": ac["AFR"], "AFR_AN": an["AFR"], "AFR_AF": af["AFR"],
                "AMR_AC": ac["AMR"], "AMR_AN": an["AMR"], "AMR_AF": af["AMR"],
                "EAS_AC": ac["EAS"], "EAS_AN": an["EAS"], "EAS_AF": af["EAS"],
                "EUR_AC": ac["EUR"], "EUR_AN": an["EUR"], "EUR_AF": af["EUR"],
                "SAS_AC": ac["SAS"], "SAS_AN": an["SAS"], "SAS_AF": af["SAS"],
                "superpop_AF_range": af_range,
                "top_superpopulation": top_sp,
                "top_superpopulation_AF": top_af,
                "second_superpopulation": second_sp,
                "second_superpopulation_AF": second_af,
                "frequency_entropy_superpop": af_entropy,
                "alt_count_entropy_superpop": ac_entropy,
                "max_pairwise_fst_hudson": max_fst,
                "max_pairwise_fst_pair": max_pair,
                "mean_pairwise_fst_hudson": mean_fst,
                "population_structure_class": pop_class,
            }

            if variant_key in overlay:
                for k, v in overlay[variant_key].items():
                    row[f"candidate_{k}"] = v

            variant_rows.append(row)
            processed += 1
            if processed % 1000 == 0:
                print(f"Processed {processed} variants")

    base_fields = [
        "gene", "variant_key", "chrom", "pos", "ref", "alt",
        "global_AC", "global_AN", "global_AF",
        "AFR_AC", "AFR_AN", "AFR_AF",
        "AMR_AC", "AMR_AN", "AMR_AF",
        "EAS_AC", "EAS_AN", "EAS_AF",
        "EUR_AC", "EUR_AN", "EUR_AF",
        "SAS_AC", "SAS_AN", "SAS_AF",
        "superpop_AF_range", "top_superpopulation", "top_superpopulation_AF",
        "second_superpopulation", "second_superpopulation_AF",
        "frequency_entropy_superpop", "alt_count_entropy_superpop",
        "max_pairwise_fst_hudson", "max_pairwise_fst_pair", "mean_pairwise_fst_hudson",
        "population_structure_class",
    ]

    overlay_fields = sorted({k for r in variant_rows for k in r.keys() if k.startswith("candidate_")})
    variant_fields = base_fields + overlay_fields

    with OUT_VARIANTS.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=variant_fields)
        w.writeheader()
        for r in variant_rows:
            out = dict(r)
            for k in [
                "global_AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF",
                "superpop_AF_range", "top_superpopulation_AF", "second_superpopulation_AF",
                "frequency_entropy_superpop", "alt_count_entropy_superpop",
                "max_pairwise_fst_hudson", "mean_pairwise_fst_hudson"
            ]:
                out[k] = fmt(out.get(k, ""))
            w.writerow({k: out.get(k, "") for k in variant_fields})

    pair_fields = ["gene", "variant_key", "pair", "fst_hudson", "af_1", "af_2", "n_alleles_1", "n_alleles_2"]
    with OUT_PAIRWISE_SUMMARY.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=pair_fields)
        w.writeheader()
        for r in pair_rows:
            out = dict(r)
            for k in ["fst_hudson", "af_1", "af_2"]:
                out[k] = fmt(out[k])
            w.writerow(out)

    # Gene-level summary
    gene_rows = []
    for gene in sorted(set(r["gene"] for r in variant_rows)):
        rows = [r for r in variant_rows if r["gene"] == gene]
        if not rows:
            continue
        classes = Counter(r["population_structure_class"] for r in rows)
        gene_rows.append({
            "gene": gene,
            "variants": len(rows),
            "mean_global_AF": sum(r["global_AF"] for r in rows) / len(rows),
            "mean_superpop_AF_range": sum(r["superpop_AF_range"] for r in rows) / len(rows),
            "max_superpop_AF_range": max(r["superpop_AF_range"] for r in rows),
            "mean_pairwise_fst_hudson": sum(r["mean_pairwise_fst_hudson"] for r in rows) / len(rows),
            "max_pairwise_fst_hudson": max(r["max_pairwise_fst_hudson"] for r in rows),
            "mean_frequency_entropy_superpop": sum(r["frequency_entropy_superpop"] for r in rows) / len(rows),
            "private_or_near_private": classes.get("private_or_near_private_superpopulation_signal", 0),
            "very_strong_population_structure": classes.get("very_strong_population_structure", 0),
            "strong_population_structure": classes.get("strong_population_structure", 0),
            "moderate_population_structure": classes.get("moderate_population_structure", 0),
            "broadly_shared_common_variant": classes.get("broadly_shared_common_variant", 0),
            "low_or_background_population_structure": classes.get("low_or_background_population_structure", 0),
        })

    gene_fields = [
        "gene", "variants", "mean_global_AF", "mean_superpop_AF_range", "max_superpop_AF_range",
        "mean_pairwise_fst_hudson", "max_pairwise_fst_hudson", "mean_frequency_entropy_superpop",
        "private_or_near_private", "very_strong_population_structure", "strong_population_structure",
        "moderate_population_structure", "broadly_shared_common_variant", "low_or_background_population_structure"
    ]

    for path in [OUT_GENE_SUMMARY, TABLE6]:
        with path.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, delimiter="\t", fieldnames=gene_fields)
            w.writeheader()
            for r in gene_rows:
                out = dict(r)
                for k in [
                    "mean_global_AF", "mean_superpop_AF_range", "max_superpop_AF_range",
                    "mean_pairwise_fst_hudson", "max_pairwise_fst_hudson", "mean_frequency_entropy_superpop"
                ]:
                    out[k] = fmt(out[k])
                w.writerow(out)

    enriched_classes = {
        "private_or_near_private_superpopulation_signal",
        "very_strong_population_structure",
        "strong_population_structure",
        "moderate_population_structure",
    }

    enriched = [r for r in variant_rows if r["population_structure_class"] in enriched_classes]
    enriched_sorted = sorted(
        enriched,
        key=lambda r: (r["max_pairwise_fst_hudson"], r["superpop_AF_range"], r["global_AF"]),
        reverse=True
    )

    for path, rows, limit in [(OUT_ENRICHED, enriched_sorted, None), (TABLE7, enriched_sorted, 80), (OUT_TOP, enriched_sorted, 200)]:
        with path.open("w", encoding="utf-8", newline="") as f:
            fields = variant_fields
            w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
            w.writeheader()
            selected = rows[:limit] if limit else rows
            for r in selected:
                out = dict(r)
                for k in [
                    "global_AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF",
                    "superpop_AF_range", "top_superpopulation_AF", "second_superpopulation_AF",
                    "frequency_entropy_superpop", "alt_count_entropy_superpop",
                    "max_pairwise_fst_hudson", "mean_pairwise_fst_hudson"
                ]:
                    out[k] = fmt(out.get(k, ""))
                w.writerow({k: out.get(k, "") for k in fields})

    # Candidate overlay subset: variants that are both population-structured and already have candidate annotations.
    overlay_rows = [r for r in enriched_sorted if any(k.startswith("candidate_") for k in r)]
    with OUT_OVERLAY.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=variant_fields)
        w.writeheader()
        for r in overlay_rows:
            out = dict(r)
            for k in [
                "global_AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF",
                "superpop_AF_range", "top_superpopulation_AF", "second_superpopulation_AF",
                "frequency_entropy_superpop", "alt_count_entropy_superpop",
                "max_pairwise_fst_hudson", "mean_pairwise_fst_hudson"
            ]:
                out[k] = fmt(out.get(k, ""))
            w.writerow({k: out.get(k, "") for k in variant_fields})

    # Results subsection
    total = len(variant_rows)
    class_counts = Counter(r["population_structure_class"] for r in variant_rows)
    top_gene_fst = sorted(gene_rows, key=lambda r: r["max_pairwise_fst_hudson"], reverse=True)
    top_variants = enriched_sorted[:10]

    md_lines = []
    md_lines.append("# Formal population-genetic structure of HRH1-HRH4 variation")
    md_lines.append("")
    md_lines.append(
        "To move beyond descriptive superpopulation allele-frequency range, exact-window HRH variants were analyzed using "
        "superpopulation-resolved allele counts, normalized allele-frequency entropy and a Hudson-style pairwise FST estimator. "
        "This formal population-genetic layer provides a more explicit distinction between broadly shared variants, "
        "population-enriched contextual markers and variants that may warrant follow-up as population-relevant pharmacogenomic candidates."
    )
    md_lines.append("")
    md_lines.append(f"The analysis evaluated {total} exact-window HRH variants with 1000 Genomes superpopulation metadata.")
    md_lines.append("Population-structure classes were distributed as follows: " + ", ".join(f"{k}: {v}" for k, v in sorted(class_counts.items())) + ".")
    md_lines.append("")
    md_lines.append("At receptor level, maximum pairwise FST and allele-frequency entropy highlighted receptor-specific population architecture:")
    for r in top_gene_fst:
        md_lines.append(
            f"- {r['gene']}: variants={r['variants']}; max pairwise FST={r['max_pairwise_fst_hudson']:.4f}; "
            f"mean pairwise FST={r['mean_pairwise_fst_hudson']:.4f}; "
            f"max AF range={r['max_superpop_AF_range']:.4f}; mean frequency entropy={r['mean_frequency_entropy_superpop']:.4f}."
        )
    md_lines.append("")
    md_lines.append("The top population-structured variants were:")
    for r in top_variants:
        md_lines.append(
            f"- {r['gene']} {r['variant_key']}: class={r['population_structure_class']}; "
            f"top superpopulation={r['top_superpopulation']} AF={r['top_superpopulation_AF']:.4f}; "
            f"AF range={r['superpop_AF_range']:.4f}; "
            f"max pairwise FST={r['max_pairwise_fst_hudson']:.4f} ({r['max_pairwise_fst_pair']})."
        )
    md_lines.append("")
    md_lines.append(
        "These results should be interpreted as population-genetic prioritization rather than clinical association. "
        "High FST or population enrichment does not demonstrate altered antihistamine response by itself; instead, it identifies variants and haplotypes "
        "whose population distribution may be important when combined with consequence, GPCR topology, regulatory, receptor-fate or ligand-context evidence."
    )
    md_lines.append("")
    md_lines.append("## Output files")
    md_lines.append(f"- Variant-level population structure table: `{OUT_VARIANTS.relative_to(PROJECT)}`")
    md_lines.append(f"- Gene-level summary/Table 6: `{TABLE6.relative_to(PROJECT)}`")
    md_lines.append(f"- Population-enriched candidate/Table 7: `{TABLE7.relative_to(PROJECT)}`")
    md_lines.append(f"- Candidate overlay: `{OUT_OVERLAY.relative_to(PROJECT)}`")

    RESULTS_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(f"Wrote variant-level population structure: {OUT_VARIANTS}")
    print(f"Wrote gene summary: {OUT_GENE_SUMMARY}")
    print(f"Wrote pairwise FST table: {OUT_PAIRWISE_SUMMARY}")
    print(f"Wrote enriched candidates: {OUT_ENRICHED}")
    print(f"Wrote candidate overlay: {OUT_OVERLAY}")
    print(f"Wrote Table 6: {TABLE6}")
    print(f"Wrote Table 7: {TABLE7}")
    print(f"Wrote Results subsection: {RESULTS_MD}")
    print()
    print(f"Variants processed: {len(variant_rows)}")
    print("Population-structure class counts:")
    for k, v in sorted(class_counts.items()):
        print(f"  {k}: {v}")
    print()
    print("Gene-level max pairwise FST:")
    for r in top_gene_fst:
        print(f"  {r['gene']}: {r['max_pairwise_fst_hudson']:.4f}")

if __name__ == "__main__":
    main()
