#!/usr/bin/env python3

from pathlib import Path
from collections import Counter, defaultdict
import csv
import runpy

PROJECT = Path("/mnt/d/histamine_receptor_pharmacogenomics")

ORIGINAL_SCRIPT = PROJECT / (
    "16_scripts/python/"
    "28_hrh2_hrh4_mirna_rbp_tf_motif_disruption_v1_v2_ncbi.py"
)

INPUT_SCAN = PROJECT / (
    "06_regulatory_genomics/motif_disruption/"
    "HRH2_HRH4_regulatory_motif_disruption_scan_v1_v2_ncbi.tsv"
)

OUTDIR = PROJECT / (
    "14_manuscript/PLOS_ONE_revision_PONE-D-26-38698_20260918/"
    "07_corrected_analysis/R2-04_motif_indel_aware_v1/01_corrected_motif"
)

OUT_SCAN = OUTDIR / "HRH2_HRH4_regulatory_motif_disruption_scan_indel_aware_v1.tsv"
OUT_CHANGED = OUTDIR / "HRH2_HRH4_regulatory_motif_disruption_changed_targets_indel_aware_v1.tsv"
OUT_TRUE = OUTDIR / "HRH2_HRH4_true_motif_disruption_candidates_indel_aware_v1.tsv"
OUT_SUMMARY = OUTDIR / "R2-04_indel_aware_motif_correction_summary.tsv"

defs = runpy.run_path(str(ORIGINAL_SCRIPT), run_name="r2_04_original_defs")

MOTIFS = defs["MOTIFS"]
find_motifs = defs["find_motifs"]
parse_variant_key = defs["parse_variant_key"]
overlaps = defs["overlaps"]
score_row = defs["score_row"]
motif_names = defs["motif_names"]
motif_summary = defs["motif_summary"]
changed_overlapping_variant = defs["changed_overlapping_variant"]


def read_tsv(path):
    with path.open("r", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})


def normal_key(hit, start=None, end=None):
    return (
        hit["motif_group"],
        hit["motif_id"],
        hit["start"] if start is None else start,
        hit["end"] if end is None else end,
        hit["match"],
    )


def overlap_identity(hit):
    return (
        hit["motif_group"],
        hit["motif_id"],
        hit["match"],
    )


def multiset_excess(left, right):
    pool = defaultdict(list)
    for hit in right:
        pool[overlap_identity(hit)].append(hit)

    excess = []
    for hit in left:
        k = overlap_identity(hit)
        if pool[k]:
            pool[k].pop()
        else:
            excess.append(hit)
    return excess


def classify_direction_indel_aware(ref_hits, alt_hits, row):
    _, pos, ref_allele, alt_allele = parse_variant_key(row["variant_key"])

    win_start = int(row.get("ucsc_window_start", 0) or 0)
    local_start = pos - 1 - win_start

    ref_end = local_start + max(1, len(ref_allele))
    alt_end = local_start + max(1, len(alt_allele))

    ref_overlap = [
        h for h in ref_hits
        if overlaps(h["start"], h["end"], local_start, ref_end)
    ]

    alt_overlap = [
        h for h in alt_hits
        if overlaps(h["start"], h["end"], local_start, alt_end)
    ]

    ref_non = {}
    for h in ref_hits:
        if h in ref_overlap:
            continue
        ref_non[normal_key(h)] = h

    alt_non = {}
    shift = len(ref_allele) - len(alt_allele)

    for h in alt_hits:
        if h in alt_overlap:
            continue

        s = h["start"]
        e = h["end"]

        # Hits wholly downstream of the altered ALT allele are
        # projected back onto REF-window coordinates.
        if s >= alt_end:
            s += shift
            e += shift

        alt_non[normal_key(h, s, e)] = h

    lost_non = [
        ref_non[k]
        for k in sorted(set(ref_non) - set(alt_non))
    ]

    gained_non = [
        alt_non[k]
        for k in sorted(set(alt_non) - set(ref_non))
    ]

    # Within the altered allele span, compare motif identity/matched
    # sequence rather than raw local coordinates.
    lost_overlap = multiset_excess(ref_overlap, alt_overlap)
    gained_overlap = multiset_excess(alt_overlap, ref_overlap)

    return (
        lost_non + lost_overlap,
        gained_non + gained_overlap,
    )


rows = read_tsv(INPUT_SCAN)
fields = list(rows[0].keys())

corrected = []

for original in rows:
    row = dict(original)

    ref_seq = row.get("ref_sequence_window", "").upper()
    alt_seq = row.get("alt_sequence_window", "").upper()

    ref_hits = []
    alt_hits = []

    for motif in MOTIFS:
        ref_hits.extend(find_motifs(ref_seq, motif))
        alt_hits.extend(find_motifs(alt_seq, motif))

    lost, gained = classify_direction_indel_aware(
        ref_hits, alt_hits, row
    )

    lost_var = changed_overlapping_variant(
        lost, row, ref=True
    )
    gained_var = changed_overlapping_variant(
        gained, row, ref=False
    )

    score, priority = score_row(
        row, lost, gained, lost_var, gained_var
    )

    changed_groups = Counter(
        h["motif_group"] for h in lost + gained
    )

    row["motif_disruption_priority"] = priority
    row["motif_disruption_score"] = str(score)
    row["has_any_motif_change"] = str(bool(lost or gained))
    row["has_variant_overlapping_motif_change"] = str(
        bool(lost_var or gained_var)
    )
    row["changed_motif_group_counts"] = ";".join(
        f"{k}:{v}" for k, v in sorted(changed_groups.items())
    )
    row["lost_motifs"] = motif_names(lost)
    row["gained_motifs"] = motif_names(gained)
    row["lost_variant_overlapping_motifs"] = motif_names(lost_var)
    row["gained_variant_overlapping_motifs"] = motif_names(gained_var)
    row["ref_motif_summary"] = motif_summary(ref_hits)
    row["alt_motif_summary"] = motif_summary(alt_hits)

    corrected.append(row)

changed = [
    r for r in corrected
    if r["has_any_motif_change"].lower() == "true"
]

counts = Counter(r["gene"] for r in changed)

HRH2_TRUE = "5:175709501:TCTGCAGCTGCGTGC:T"
HRH4_ARTIFACT_1 = "18:24479080:GC:G"
HRH4_ARTIFACT_2 = "18:24456680:A:AT"

by_key = {r["variant_key"]: r for r in corrected}

assert len(changed) == 14, f"Expected 14 true motif rows, observed {len(changed)}"
assert counts["HRH2"] == 8, f"Expected HRH2=8, observed {counts['HRH2']}"
assert counts["HRH4"] == 6, f"Expected HRH4=6, observed {counts['HRH4']}"

assert by_key[HRH2_TRUE]["has_any_motif_change"] == "True"
assert by_key[HRH2_TRUE]["has_variant_overlapping_motif_change"] == "True"

assert by_key[HRH4_ARTIFACT_1]["has_any_motif_change"] == "False"
assert by_key[HRH4_ARTIFACT_2]["has_any_motif_change"] == "False"

write_tsv(OUT_SCAN, corrected, fields)
write_tsv(OUT_CHANGED, changed, fields)

true_fields = [
    "gene",
    "variant_key",
    "motif_disruption_score",
    "motif_disruption_priority",
    "has_variant_overlapping_motif_change",
    "changed_motif_group_counts",
    "lost_motifs",
    "gained_motifs",
    "lost_variant_overlapping_motifs",
    "gained_variant_overlapping_motifs",
]

write_tsv(OUT_TRUE, changed, true_fields)

summary = [
    {"metric": "input_rows", "value": len(corrected)},
    {"metric": "true_motif_total", "value": len(changed)},
    {"metric": "true_motif_HRH2", "value": counts["HRH2"]},
    {"metric": "true_motif_HRH4", "value": counts["HRH4"]},
    {
        "metric": "HRH2_deletion_variant_overlap",
        "value": by_key[HRH2_TRUE]["has_variant_overlapping_motif_change"],
    },
    {
        "metric": "HRH4_24479080_true_change",
        "value": by_key[HRH4_ARTIFACT_1]["has_any_motif_change"],
    },
    {
        "metric": "HRH4_24456680_true_change",
        "value": by_key[HRH4_ARTIFACT_2]["has_any_motif_change"],
    },
]

write_tsv(OUT_SUMMARY, summary, ["metric", "value"])

print("STATUS=PASS")
print(f"INPUT_ROWS={len(corrected)}")
print(f"TRUE_MOTIF_TOTAL={len(changed)}")
print(f"TRUE_MOTIF_HRH2={counts['HRH2']}")
print(f"TRUE_MOTIF_HRH4={counts['HRH4']}")
print(
    "HRH2_DELETION_VARIANT_OVERLAP="
    + by_key[HRH2_TRUE]["has_variant_overlapping_motif_change"]
)
print(
    "HRH4_24479080_TRUE_CHANGE="
    + by_key[HRH4_ARTIFACT_1]["has_any_motif_change"]
)
print(
    "HRH4_24456680_TRUE_CHANGE="
    + by_key[HRH4_ARTIFACT_2]["has_any_motif_change"]
)
print(f"OUT_SCAN={OUT_SCAN}")
print(f"OUT_CHANGED={OUT_CHANGED}")
print(f"OUT_TRUE={OUT_TRUE}")
print(f"OUT_SUMMARY={OUT_SUMMARY}")
