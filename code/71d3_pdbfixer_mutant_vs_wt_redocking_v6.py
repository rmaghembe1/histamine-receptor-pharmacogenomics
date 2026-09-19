#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import os
import shutil
import subprocess
import datetime
from collections import defaultdict

PROJECT = Path("/mnt/d/histamine_receptor_pharmacogenomics")

DISCOVERY = PROJECT / "11_results/monitoring_evaluation_v6_global_functional_diversity/HRH_redocking_input_discovery_repair_v6.tsv"
PANEL = PROJECT / "11_results/monitoring_evaluation_v6_global_functional_diversity/HRH_ligand_pocket_population_redocking_candidate_panel_clean_v6.tsv"

WORKDIR = PROJECT / "10_docking/mutant_redocking_v6_population_panel_pdbfixer"
OUTDIR = PROJECT / "11_results/monitoring_evaluation_v6_global_functional_diversity"
FIGDIR = PROJECT / "12_figures/monitoring_evaluation_v6_global_functional_diversity"

WORKDIR.mkdir(parents=True, exist_ok=True)
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

OUT_RESULTS = OUTDIR / "HRH_mutant_vs_wt_redocking_results_pdbfixer_v6.tsv"
OUT_FAILURES = OUTDIR / "HRH_mutant_vs_wt_redocking_failures_pdbfixer_v6.tsv"
OUT_SUMMARY = OUTDIR / "HRH_mutant_vs_wt_redocking_summary_pdbfixer_v6.tsv"
OUT_CHECKLIST = OUTDIR / "HRH_mutant_vs_wt_redocking_pdbfixer_checklist_v6.tsv"
OUT_REPORT = OUTDIR / "HRH_mutant_vs_wt_redocking_pdbfixer_report_v6.md"
OUT_METHODS = OUTDIR / "HRH_mutant_vs_wt_redocking_pdbfixer_methods_note_v6.md"

OUT_FIG_PNG = FIGDIR / "Figure_5_mutant_vs_wt_redocking_delta_affinity_pdbfixer_v6_900dpi.png"
OUT_FIG_SVG = FIGDIR / "Figure_5_mutant_vs_wt_redocking_delta_affinity_pdbfixer_v6.svg"

MAX_CANDIDATES = int(os.environ.get("MAX_71D3_CANDIDATES", "16"))

AA3 = {
    "Ala": "ALA", "Arg": "ARG", "Asn": "ASN", "Asp": "ASP", "Cys": "CYS",
    "Gln": "GLN", "Glu": "GLU", "Gly": "GLY", "His": "HIS", "Ile": "ILE",
    "Leu": "LEU", "Lys": "LYS", "Met": "MET", "Phe": "PHE", "Pro": "PRO",
    "Ser": "SER", "Thr": "THR", "Trp": "TRP", "Tyr": "TYR", "Val": "VAL",
}

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

def which(cmd):
    return shutil.which(cmd) or ""

def run_cmd(cmd, cwd=None, timeout=1800):
    try:
        p = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
        )
        return {
            "returncode": p.returncode,
            "stdout": p.stdout,
            "stderr": p.stderr,
            "cmd": " ".join(str(x) for x in cmd),
        }
    except Exception as e:
        return {
            "returncode": 999,
            "stdout": "",
            "stderr": str(e),
            "cmd": " ".join(str(x) for x in cmd),
        }

def safe_name(s):
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s))
    return s.strip("_")[:170]

def parse_protein_change(p):
    m = re.match(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})$", p or "")
    if not m:
        return None
    ref3, pos, alt3 = m.groups()
    return {
        "ref3": AA3.get(ref3, ref3.upper()),
        "pos": int(pos),
        "alt3": AA3.get(alt3, alt3.upper()),
    }

def to_abs_path(s):
    if not s:
        return None
    p = Path(s)
    if not p.is_absolute():
        p = PROJECT / p
    return p

def extract_first_model_pdbqt(in_pdbqt, out_pdbqt):
    """
    Extract the first docked ligand pose into a Vina-safe single-model PDBQT.

    Vina can reject ligand PDBQT files that contain MODEL/ENDMDL tags, even
    when only one model is retained. Therefore this function strips MODEL and
    ENDMDL records and writes only the first pose contents.
    """
    lines = in_pdbqt.read_text(encoding="utf-8", errors="replace").splitlines()

    has_model = any(line.startswith("MODEL") for line in lines)
    out = []

    if has_model:
        in_first_model = False
        saw_first_model = False

        for line in lines:
            if line.startswith("MODEL"):
                if saw_first_model:
                    break
                saw_first_model = True
                in_first_model = True
                continue

            if line.startswith("ENDMDL") and in_first_model:
                break

            if in_first_model:
                if not line.startswith(("MODEL", "ENDMDL")):
                    out.append(line)
    else:
        out = [
            line for line in lines
            if not line.startswith(("MODEL", "ENDMDL"))
        ]

    cleaned = []
    allowed_prefixes = (
        "REMARK", "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH",
        "ATOM", "HETATM", "TORSDOF"
    )
    for line in out:
        if line.startswith(allowed_prefixes):
            cleaned.append(line)

    if not cleaned:
        return False, "empty_ligand_pose_after_model_tag_stripping"

    out_pdbqt.write_text(chr(10).join(cleaned) + chr(10), encoding="utf-8")
    return True, ""

def ligand_box_from_pdbqt(pdbqt, padding=10.0, min_size=22.5):
    xs, ys, zs = [], [], []
    for line in pdbqt.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                parts = line.split()
                try:
                    x, y, z = float(parts[6]), float(parts[7]), float(parts[8])
                except Exception:
                    continue
            xs.append(x)
            ys.append(y)
            zs.append(z)

    if not xs:
        return None

    def center_size(vals):
        mn, mx = min(vals), max(vals)
        return (mn + mx) / 2.0, max(mx - mn + padding, min_size)

    cx, sx = center_size(xs)
    cy, sy = center_size(ys)
    cz, sz = center_size(zs)
    return {
        "center_x": cx, "center_y": cy, "center_z": cz,
        "size_x": sx, "size_y": sy, "size_z": sz,
    }

def convert_receptor_to_pdbqt(obabel, in_pdb, out_pdbqt):
    cmd = [obabel, "-ipdb", str(in_pdb), "-opdbqt", "-O", str(out_pdbqt), "-xr"]
    res = run_cmd(cmd, timeout=1200)
    if res["returncode"] != 0 or not out_pdbqt.exists() or out_pdbqt.stat().st_size == 0:
        return False, res
    return True, res

def parse_vina_affinity(stdout, out_pdbqt=None):
    for line in stdout.splitlines():
        m = re.match(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+", line)
        if m:
            return float(m.group(1))
    if out_pdbqt and out_pdbqt.exists():
        text = out_pdbqt.read_text(encoding="utf-8", errors="replace")
        m = re.search(r"REMARK VINA RESULT:\s*(-?\d+(?:\.\d+)?)", text)
        if m:
            return float(m.group(1))
    return None

def run_vina(vina, receptor_pdbqt, ligand_pdbqt, out_pdbqt, log_txt, box, seed):
    cmd = [
        vina,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", f"{box['center_x']:.3f}",
        "--center_y", f"{box['center_y']:.3f}",
        "--center_z", f"{box['center_z']:.3f}",
        "--size_x", f"{box['size_x']:.3f}",
        "--size_y", f"{box['size_y']:.3f}",
        "--size_z", f"{box['size_z']:.3f}",
        "--exhaustiveness", "8",
        "--num_modes", "9",
        "--seed", str(seed),
        "--out", str(out_pdbqt),
    ]
    res = run_cmd(cmd, timeout=2400)
    log_txt.write_text(
        "CMD: " + res["cmd"] + "\n\nSTDOUT:\n" + res["stdout"] + "\n\nSTDERR:\n" + res["stderr"],
        encoding="utf-8",
    )
    aff = parse_vina_affinity(res["stdout"], out_pdbqt)
    return res, aff

def load_pdbfixer():
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
        return True, "", PDBFixer, PDBFile
    except Exception as e:
        return False, str(e), None, None

def residue_id_number(resid):
    m = re.search(r"\d+", str(resid))
    return int(m.group(0)) if m else None

def find_target_residue(PDBFile, receptor_pdb, chain_hint, pos, ref3):
    pdb = PDBFile(str(receptor_pdb))
    candidates = []
    for chain in pdb.topology.chains():
        for res in chain.residues():
            num = residue_id_number(res.id)
            if num == pos:
                candidates.append({
                    "chain": chain.id,
                    "resid": res.id,
                    "name": res.name,
                    "ref_match": res.name.upper() == ref3.upper(),
                    "hint_match": str(chain.id) == str(chain_hint),
                })

    if not candidates:
        return None, "target_residue_number_not_found"

    ref_matches = [c for c in candidates if c["ref_match"]]
    if chain_hint:
        hinted_ref = [c for c in ref_matches if c["hint_match"]]
        if hinted_ref:
            return hinted_ref[0], ""
        hinted_any = [c for c in candidates if c["hint_match"]]
        if hinted_any and not hinted_any[0]["ref_match"]:
            return None, f"chain_hint_residue_ref_mismatch_found_{hinted_any[0]['name']}_expected_{ref3}"

    if len(ref_matches) == 1:
        return ref_matches[0], ""

    if len(ref_matches) > 1:
        return ref_matches[0], "multiple_ref_matching_residues_used_first"

    observed = ",".join(sorted(set(c["name"] for c in candidates)))
    return None, f"residue_ref_mismatch_observed_{observed}_expected_{ref3}"

def make_pdbfixer_mutant(PDBFixer, PDBFile, receptor_pdb, chain_hint, ref3, pos, alt3, out_pdb, run_dir):
    target, warning = find_target_residue(PDBFile, receptor_pdb, chain_hint, pos, ref3)
    if target is None:
        return False, warning, "", "", ""

    chain_id = target["chain"]
    resid = target["resid"]
    mutation = f"{ref3}-{resid}-{alt3}"

    try:
        fixer = PDBFixer(filename=str(receptor_pdb))

        # PDBFixer 1.12 expects missingResidues to exist before
        # findMissingAtoms/addMissingAtoms. Initialize it explicitly before
        # and after mutation to avoid the missingResidues attribute error.
        try:
            fixer.findMissingResidues()
        except Exception:
            fixer.missingResidues = {}
        if not hasattr(fixer, "missingResidues"):
            fixer.missingResidues = {}

        fixer.applyMutations([mutation], chain_id)

        if not hasattr(fixer, "missingResidues"):
            fixer.missingResidues = {}

        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        with open(out_pdb, "w", encoding="utf-8") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
    except Exception as e:
        return False, str(e), chain_id, resid, mutation

    if not out_pdb.exists() or out_pdb.stat().st_size == 0:
        return False, "mutant_pdb_not_created", chain_id, resid, mutation

    note = warning if warning else "ok"
    return True, note, chain_id, resid, mutation

def meta_key(row):
    return (
        str(row.get("rank", "")),
        row.get("gene", ""),
        row.get("ligand", ""),
        row.get("variant_key", ""),
        row.get("protein_change", ""),
    )

def build_panel_metadata():
    meta = {}
    for r in read_tsv(PANEL):
        k = meta_key(r)
        meta[k] = r
        k2 = ("", r.get("gene", ""), r.get("ligand", ""), r.get("variant_key", ""), r.get("protein_change", ""))
        meta.setdefault(k2, r)
    return meta

def pick_meta(row, meta):
    k = meta_key(row)
    if k in meta:
        return meta[k]
    k2 = ("", row.get("gene", ""), row.get("ligand", ""), row.get("variant_key", ""), row.get("protein_change", ""))
    return meta.get(k2, {})

def make_figure(rows):
    if not rows:
        return False, "no_successful_rows"
    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        return False, str(e)

    labels = [
        f"{r['rank']} {r['gene']} {r['ligand']} {r['protein_change']}"
        for r in rows
    ]
    vals = [float(r["delta_mutant_minus_wt_kcal_mol"]) for r in rows]

    height = max(4.5, 0.42 * len(rows) + 1.5)
    plt.figure(figsize=(9, height))
    y = list(range(len(vals)))
    plt.barh(y, vals)
    plt.yticks(y, labels, fontsize=7)
    plt.axvline(0, linewidth=1)
    plt.xlabel("Delta Vina affinity, mutant minus WT (kcal/mol)")
    plt.tight_layout()
    plt.savefig(OUT_FIG_PNG, dpi=900)
    plt.savefig(OUT_FIG_SVG)
    plt.close()
    return True, "created"

def main():
    timestamp = datetime.datetime.now().isoformat(timespec="seconds")
    discovery_rows = read_tsv(DISCOVERY)
    meta = build_panel_metadata()

    ok_pdbfixer, pdbfixer_err, PDBFixer, PDBFile = load_pdbfixer()
    obabel = which("obabel")
    vina = which("vina")

    results = []
    failures = []

    rows = discovery_rows[:MAX_CANDIDATES]

    for row in rows:
        mrow = pick_meta(row, meta)

        common = {
            "rank": row.get("rank", ""),
            "gene": row.get("gene", ""),
            "ligand": row.get("ligand", ""),
            "variant_key": row.get("variant_key", ""),
            "protein_change": row.get("protein_change", ""),
            "redocking_priority_class": mrow.get("redocking_priority_class", row.get("redocking_priority_class", "")),
            "top_pairwise_comparison": mrow.get("top_pairwise_comparison", row.get("top_pairwise_comparison", "")),
            "top_pairwise_abs_delta_af": mrow.get("top_pairwise_abs_delta_af", row.get("top_pairwise_abs_delta_af", "")),
            "top_pairwise_fdr_q": mrow.get("top_pairwise_fdr_q", row.get("top_pairwise_fdr_q", "")),
        }

        run_dir = WORKDIR / safe_name(
            f"{common['rank']}_{common['gene']}_{common['ligand']}_{common['variant_key']}_{common['protein_change']}"
        )
        run_dir.mkdir(parents=True, exist_ok=True)

        if not ok_pdbfixer:
            failures.append({**common, "stage": "tool_check", "reason": "pdbfixer_openmm_import_failed_" + pdbfixer_err, "run_dir": str(run_dir.relative_to(PROJECT))})
            continue
        if not obabel:
            failures.append({**common, "stage": "tool_check", "reason": "obabel_not_found", "run_dir": str(run_dir.relative_to(PROJECT))})
            continue
        if not vina:
            failures.append({**common, "stage": "tool_check", "reason": "vina_not_found", "run_dir": str(run_dir.relative_to(PROJECT))})
            continue

        pc = parse_protein_change(common["protein_change"])
        if not pc:
            failures.append({**common, "stage": "input_parse", "reason": "protein_change_parse_failed", "run_dir": str(run_dir.relative_to(PROJECT))})
            continue

        receptor_pdb = to_abs_path(row.get("receptor_pdb", ""))
        docked_pose = to_abs_path(row.get("docked_pose_pdbqt", ""))
        chain_hint = row.get("chain", "")

        if not receptor_pdb or not receptor_pdb.exists():
            failures.append({**common, "stage": "input_discovery", "reason": "receptor_pdb_not_found", "run_dir": str(run_dir.relative_to(PROJECT))})
            continue
        if not docked_pose or not docked_pose.exists():
            failures.append({**common, "stage": "input_discovery", "reason": "docked_pose_pdbqt_not_found", "run_dir": str(run_dir.relative_to(PROJECT))})
            continue

        ligand_pdbqt = run_dir / "ligand_first_pose_input.pdbqt"
        ok_lig, lig_reason = extract_first_model_pdbqt(docked_pose, ligand_pdbqt)
        if not ok_lig:
            failures.append({**common, "stage": "ligand_pose_extract", "reason": lig_reason, "run_dir": str(run_dir.relative_to(PROJECT))})
            continue

        box = ligand_box_from_pdbqt(ligand_pdbqt)
        if box is None:
            failures.append({**common, "stage": "vina_box", "reason": "ligand_coordinates_not_parsed", "run_dir": str(run_dir.relative_to(PROJECT))})
            continue

        wt_receptor_pdbqt = run_dir / "wt_receptor.pdbqt"
        mutant_receptor_pdb = run_dir / "mutant_receptor_pdbfixer.pdb"
        mutant_receptor_pdbqt = run_dir / "mutant_receptor_pdbfixer.pdbqt"

        wt_ok, wt_conv = convert_receptor_to_pdbqt(obabel, receptor_pdb, wt_receptor_pdbqt)
        (run_dir / "obabel_wt_to_pdbqt.log").write_text(
            "CMD: " + wt_conv["cmd"] + "\n\nSTDOUT:\n" + wt_conv["stdout"] + "\n\nSTDERR:\n" + wt_conv["stderr"],
            encoding="utf-8",
        )
        if not wt_ok:
            failures.append({**common, "stage": "wt_pdbqt_conversion", "reason": wt_conv["stderr"][:500] or "wt_conversion_failed", "run_dir": str(run_dir.relative_to(PROJECT))})
            continue

        mut_ok, mut_reason, chain_used, residue_id_used, mutation_string = make_pdbfixer_mutant(
            PDBFixer=PDBFixer,
            PDBFile=PDBFile,
            receptor_pdb=receptor_pdb,
            chain_hint=chain_hint,
            ref3=pc["ref3"],
            pos=pc["pos"],
            alt3=pc["alt3"],
            out_pdb=mutant_receptor_pdb,
            run_dir=run_dir,
        )
        (run_dir / "pdbfixer_mutation_status.txt").write_text(
            f"ok={mut_ok}\nreason={mut_reason}\nchain_used={chain_used}\nresidue_id_used={residue_id_used}\nmutation_string={mutation_string}\n",
            encoding="utf-8",
        )
        if not mut_ok:
            failures.append({**common, "stage": "pdbfixer_mutagenesis", "reason": mut_reason, "run_dir": str(run_dir.relative_to(PROJECT))})
            continue

        mut_conv_ok, mut_conv = convert_receptor_to_pdbqt(obabel, mutant_receptor_pdb, mutant_receptor_pdbqt)
        (run_dir / "obabel_mutant_to_pdbqt.log").write_text(
            "CMD: " + mut_conv["cmd"] + "\n\nSTDOUT:\n" + mut_conv["stdout"] + "\n\nSTDERR:\n" + mut_conv["stderr"],
            encoding="utf-8",
        )
        if not mut_conv_ok:
            failures.append({**common, "stage": "mutant_pdbqt_conversion", "reason": mut_conv["stderr"][:500] or "mutant_conversion_failed", "run_dir": str(run_dir.relative_to(PROJECT))})
            continue

        wt_out = run_dir / "wt_redocked_out.pdbqt"
        mut_out = run_dir / "mutant_redocked_out.pdbqt"
        wt_log = run_dir / "wt_redocking_vina.log"
        mut_log = run_dir / "mutant_redocking_vina.log"

        seed_base = 20260701 + int(common["rank"] or 0)
        wt_vina, wt_aff = run_vina(vina, wt_receptor_pdbqt, ligand_pdbqt, wt_out, wt_log, box, seed_base)
        mut_vina, mut_aff = run_vina(vina, mutant_receptor_pdbqt, ligand_pdbqt, mut_out, mut_log, box, seed_base)

        if wt_vina["returncode"] != 0 or wt_aff is None:
            failures.append({**common, "stage": "wt_vina", "reason": wt_vina["stderr"][:500] or "wt_affinity_not_parsed", "run_dir": str(run_dir.relative_to(PROJECT))})
            continue
        if mut_vina["returncode"] != 0 or mut_aff is None:
            failures.append({**common, "stage": "mutant_vina", "reason": mut_vina["stderr"][:500] or "mutant_affinity_not_parsed", "run_dir": str(run_dir.relative_to(PROJECT))})
            continue

        delta = mut_aff - wt_aff
        if delta > 0:
            direction = "mutant_weaker_predicted_binding"
        elif delta < 0:
            direction = "mutant_stronger_predicted_binding"
        else:
            direction = "no_predicted_affinity_difference"

        results.append({
            **common,
            "ref3": pc["ref3"],
            "alt3": pc["alt3"],
            "position": pc["pos"],
            "chain_hint": chain_hint,
            "chain_used": chain_used,
            "residue_id_used": residue_id_used,
            "pdbfixer_mutation_string": mutation_string,
            "pdbfixer_note": mut_reason,
            "receptor_pdb": str(receptor_pdb.relative_to(PROJECT)),
            "docked_pose_source": str(docked_pose.relative_to(PROJECT)),
            "pose_discovery_method": row.get("pose_discovery_method", ""),
            "wt_vina_affinity_kcal_mol": f"{wt_aff:.3f}",
            "mutant_vina_affinity_kcal_mol": f"{mut_aff:.3f}",
            "delta_mutant_minus_wt_kcal_mol": f"{delta:.3f}",
            "predicted_direction": direction,
            "ligand_input_pdbqt": str(ligand_pdbqt.relative_to(PROJECT)),
            "wt_receptor_pdbqt": str(wt_receptor_pdbqt.relative_to(PROJECT)),
            "mutant_receptor_pdb": str(mutant_receptor_pdb.relative_to(PROJECT)),
            "mutant_receptor_pdbqt": str(mutant_receptor_pdbqt.relative_to(PROJECT)),
            "wt_vina_log": str(wt_log.relative_to(PROJECT)),
            "mutant_vina_log": str(mut_log.relative_to(PROJECT)),
            "run_dir": str(run_dir.relative_to(PROJECT)),
        })

    result_fields = [
        "rank", "gene", "ligand", "variant_key", "protein_change", "ref3", "alt3", "position",
        "chain_hint", "chain_used", "residue_id_used", "pdbfixer_mutation_string", "pdbfixer_note",
        "redocking_priority_class", "top_pairwise_comparison", "top_pairwise_abs_delta_af", "top_pairwise_fdr_q",
        "receptor_pdb", "docked_pose_source", "pose_discovery_method",
        "wt_vina_affinity_kcal_mol", "mutant_vina_affinity_kcal_mol", "delta_mutant_minus_wt_kcal_mol",
        "predicted_direction", "ligand_input_pdbqt", "wt_receptor_pdbqt", "mutant_receptor_pdb",
        "mutant_receptor_pdbqt", "wt_vina_log", "mutant_vina_log", "run_dir",
    ]
    failure_fields = [
        "rank", "gene", "ligand", "variant_key", "protein_change",
        "redocking_priority_class", "top_pairwise_comparison", "top_pairwise_abs_delta_af", "top_pairwise_fdr_q",
        "stage", "reason", "run_dir",
    ]

    write_tsv(OUT_RESULTS, results, result_fields)
    write_tsv(OUT_FAILURES, failures, failure_fields)

    groups = defaultdict(list)
    for r in results:
        groups[(r["gene"], r["ligand"])].append(float(r["delta_mutant_minus_wt_kcal_mol"]))

    summary = []
    for (gene, ligand), vals in sorted(groups.items()):
        summary.append({
            "gene": gene,
            "ligand": ligand,
            "successful_redocking_pairs": len(vals),
            "mean_delta_mutant_minus_wt_kcal_mol": f"{sum(vals) / len(vals):.3f}",
            "max_abs_delta_kcal_mol": f"{max(abs(v) for v in vals):.3f}",
            "mutant_weaker_count": sum(1 for v in vals if v > 0),
            "mutant_stronger_count": sum(1 for v in vals if v < 0),
        })

    write_tsv(
        OUT_SUMMARY,
        summary,
        ["gene", "ligand", "successful_redocking_pairs", "mean_delta_mutant_minus_wt_kcal_mol", "max_abs_delta_kcal_mol", "mutant_weaker_count", "mutant_stronger_count"],
    )

    fig_ok, fig_reason = make_figure(results)

    checklist = [
        {"check": "discovery_table_found", "status": "PASS" if DISCOVERY.exists() else "FAIL", "observed": str(DISCOVERY)},
        {"check": "pdbfixer_openmm_import", "status": "PASS" if ok_pdbfixer else "FAIL", "observed": "OK" if ok_pdbfixer else pdbfixer_err},
        {"check": "obabel_found", "status": "PASS" if obabel else "FAIL", "observed": obabel or "not_found"},
        {"check": "vina_found", "status": "PASS" if vina else "FAIL", "observed": vina or "not_found"},
        {"check": "candidates_attempted", "status": "PASS" if rows else "FAIL", "observed": str(len(rows))},
        {"check": "successful_redocking_pairs", "status": "PASS" if results else "REVIEW", "observed": str(len(results))},
        {"check": "failure_rows_recorded", "status": "INFO", "observed": str(len(failures))},
        {"check": "figure_created_if_successful", "status": "PASS" if (fig_ok or not results) else "REVIEW", "observed": fig_reason},
    ]
    write_tsv(OUT_CHECKLIST, checklist, ["check", "status", "observed"])

    OUT_METHODS.write_text(
        "# Module 71D3 PDBFixer/OpenMM mutant-versus-WT redocking methods note\n\n"
        "PyMOL-dependent mutant modeling was bypassed. Module 71D3 used the repaired 71D2 input discovery table, "
        "which contains resolved receptor PDB files and docked ligand pose PDBQT files. For each candidate, the first "
        "docked ligand pose was extracted and used to define a Vina search box. Single-residue mutant receptor models "
        "were generated with PDBFixer/OpenMM using the residue identity and numbering encoded in the protein-change "
        "field. WT and mutant receptors were converted to PDBQT using the same Open Babel route, then the same ligand "
        "pose was redocked against WT and mutant receptors using AutoDock Vina. Delta affinity was calculated as mutant "
        "minus WT Vina affinity. Because Vina scores are negative, positive delta values indicate weaker predicted "
        "binding for the mutant, whereas negative values indicate stronger predicted binding. These outputs are "
        "computational docking hypotheses and are not measured binding affinities.\n",
        encoding="utf-8",
    )

    report = []
    report.append("# HRH mutant-versus-WT redocking report: PDBFixer/OpenMM rescue v6")
    report.append("")
    report.append(f"Generated: {timestamp}")
    report.append("")
    report.append("## Headline results")
    report.append("")
    report.append(f"- Candidates attempted: {len(rows)}")
    report.append(f"- Successful WT-versus-mutant redocking pairs: {len(results)}")
    report.append(f"- Failure rows recorded: {len(failures)}")
    report.append(f"- Figure status: {'CREATED' if fig_ok else 'NOT_CREATED'}")
    report.append("")
    report.append("## Tool status")
    report.append("")
    report.append(f"- PDBFixer/OpenMM import: {'OK' if ok_pdbfixer else 'FAILED: ' + pdbfixer_err}")
    report.append(f"- Open Babel: `{obabel or 'not found'}`")
    report.append(f"- AutoDock Vina: `{vina or 'not found'}`")
    report.append("")
    report.append("## Main outputs")
    report.append("")
    for label, path in [
        ("Results", OUT_RESULTS),
        ("Failures", OUT_FAILURES),
        ("Summary", OUT_SUMMARY),
        ("Checklist", OUT_CHECKLIST),
        ("Methods note", OUT_METHODS),
        ("Figure PNG", OUT_FIG_PNG),
        ("Figure SVG", OUT_FIG_SVG),
    ]:
        report.append(f"- {label}: `{path}`")

    if results:
        report.append("")
        report.append("## Successful pairs")
        report.append("")
        report.append("| Rank | Gene | Ligand | Variant | Protein change | WT | Mutant | Delta | Direction |")
        report.append("|---|---|---|---|---|---:|---:|---:|---|")
        for r in results:
            report.append(
                f"| {r['rank']} | {r['gene']} | {r['ligand']} | {r['variant_key']} | {r['protein_change']} | "
                f"{r['wt_vina_affinity_kcal_mol']} | {r['mutant_vina_affinity_kcal_mol']} | "
                f"{r['delta_mutant_minus_wt_kcal_mol']} | {r['predicted_direction']} |"
            )

    if failures:
        report.append("")
        report.append("## Failure summary")
        report.append("")
        counts = defaultdict(int)
        for f in failures:
            counts[f["stage"] + ":" + f["reason"]] += 1
        for k, n in sorted(counts.items()):
            report.append(f"- {k}: {n}")

    OUT_REPORT.write_text("\n".join(report) + "\n", encoding="utf-8")

    print("\n".join(report))

if __name__ == "__main__":
    main()
