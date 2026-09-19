#!/usr/bin/env python3
from pathlib import Path
import csv
import hashlib
import json
import math
import re
import shutil
import sys
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from PIL import Image

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 9,
    'axes.titlesize': 11,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.titlesize': 14,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

ROOT = Path('/mnt/d/histamine_receptor_pharmacogenomics')
REV = ROOT / '14_manuscript/PLOS_ONE_revision_PONE-D-26-38698_20260918'
CORR = REV / '07_corrected_analysis/R2-04_motif_indel_aware_v1/02_corrected_precision'
FORENSICS = REV / '06_methods_forensics'
OUT = REV / '08_final_figures_deterministic_v2'
OUT.mkdir(parents=True, exist_ok=True)

PREC_SUM = CORR / 'hrh_precision_evidence_layer_summary_indel_aware_v1.tsv'
AXIS_SUM = CORR / 'hrh_precision_receptor_axis_summary_indel_aware_v1.tsv'
CLAIMS = CORR / 'hrh_precision_validated_vs_exploratory_indel_aware_v1.tsv'
MATRIX = CORR / 'hrh_precision_functional_candidate_evidence_matrix_indel_aware_v1.tsv'
MOTIF_TRUE = CORR / 'hrh_true_motif_disruption_candidates_indel_aware_v1.tsv'
GTEX_TRUE = CORR / 'hrh_true_gtex_eqtl_intersection_candidates_indel_aware_v1.tsv'

BW = FORENSICS / 'R1_BW_GPCRdb_variant_mapping_v3_FINAL.tsv'

CONTACT = ROOT / '13_tables/manuscript_ready/Table_64C_variant_contact_perturbation_proxy_v1.tsv'
FIG12_SVG = ROOT / '12_figures/docking_contact_analysis_v1/Figure_12_representative_ligand_pocket_contact_analysis_v1.svg'
FIG12_PNG = ROOT / '12_figures/docking_contact_analysis_v1/Figure_12_representative_ligand_pocket_contact_analysis_v1_900dpi.png'

E2 = ROOT / ('11_results/monitoring_evaluation_v6_global_functional_diversity/'
             'LOCKED_71E2_LIGAND_AWARE_RECEPTOR_FATE_PASS/'
             'HRH_receptor_fate_CMA_proteostasis_ligand_aware_v6.tsv')

POP_TESTS = ROOT / '11_results/monitoring_evaluation_v6_global_functional_diversity/HRH_population_pairwise_variant_tests_v6.tsv'

EXPECTED_SHA256 = {
    POP_TESTS: 'de0862d3d0c22f7544435367e767490930d29d5c0852a7cebcfc511870626faf',
    E2: '590195d7cf2a7c5deda27e97edb9cd11dc2f72bb7295881da03b537a20a172f2',
    CONTACT: '9a997b858260887810adceac174966617ba49d5881f91479d78499fae98ff95e',
    MATRIX: '910070040124e4059f195ea1fbab25dcf7b989e42f2d4daf055adf758f7cef9d',
    PREC_SUM: '9c42ca993bdc1c682df0e57470c4781cb80b4248b4715886980201ab95147e31',
    AXIS_SUM: '4adef4af95aa3adde9f8af429c395661b8f2552780a33794f7652f70d0faea33',
    MOTIF_TRUE: 'cb81ffbe0359d81b77e4b9cef18e1ef4e91d8d567c21579bb2b7e1286f044633',
    GTEX_TRUE: 'f3bf4c0885bdc6cc0defd0cd449185aaa77f34997d45335190d515f32a7a71fc',
}

FIG12_EXPECTED_SHA = '08a2909ad41eafd03301d1ebf8307bff650ada510ef9081c8a94b2c275cd0cad'


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def require(path: Path):
    if not path.is_file():
        raise FileNotFoundError(f'MISSING REQUIRED FILE: {path}')


def check_inputs():
    for p in list(EXPECTED_SHA256) + [BW, FIG12_SVG, FIG12_PNG, CLAIMS]:
        require(p)
    failures = []
    for p, expected in EXPECTED_SHA256.items():
        actual = sha256(p)
        if actual != expected:
            failures.append((p, expected, actual))
    actual_fig12 = sha256(FIG12_PNG)
    if actual_fig12 != FIG12_EXPECTED_SHA:
        failures.append((FIG12_PNG, FIG12_EXPECTED_SHA, actual_fig12))
    if failures:
        msg = ['INPUT HASH MISMATCH: deterministic generation aborted.']
        for p, exp, act in failures:
            msg.append(f'{p}\n expected={exp}\n actual  ={act}')
        raise RuntimeError('\n'.join(msg))


def read_tsv(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def num(series):
    return pd.to_numeric(series, errors='coerce').fillna(0)


def truthy(x):
    return str(x).strip().lower() in {'1', 'true', 'yes', 'y', 'pass'}


def pretty_evidence(x):
    mapping = {
        'protein': 'Protein',
        'splice': 'Splice',
        'biogenesis_fate': 'Biogenesis/fate',
        'motif_disruption_true': 'True motif\ndisruption',
        'motif_context_only': 'Motif context\nonly',
        'ccre': 'cCRE',
        'gtex_true_hit': 'Direct GTEx/eQTL',
        'gtex_lookup_only': 'GTEx lookup\nonly',
        'ld_context': 'LD context',
        'docking_variant': 'Docking variant',
        'population': 'Population\ndifferentiation',
        'regulatory_priority': 'Regulatory priority',
    }
    return mapping.get(x, x.replace('_', ' '))


def save_all(fig, stem: str, dpi=900):
    svg = OUT / f'{stem}.svg'
    pdf = OUT / f'{stem}.pdf'
    png = OUT / f'{stem}_900dpi.png'
    tif = OUT / f'{stem}_900dpi.tiff'
    fig.savefig(svg, bbox_inches='tight')
    fig.savefig(pdf, bbox_inches='tight')
    fig.savefig(png, dpi=dpi, bbox_inches='tight')
    fig.savefig(tif, dpi=dpi, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)
    return [svg, pdf, png, tif]


def bw_lookup():
    df = read_tsv(BW)
    out = {}
    for _, r in df.iterrows():
        bw = r.get('ballesteros_weinstein_number', '').strip()
        seg = r.get('protein_segment', '').strip()
        suffix = bw if bw else seg
        out[(r.get('gene','').strip(), r.get('variant','').strip())] = suffix
    return out


def annotate_variant(gene, variant, lookup):
    suffix = lookup.get((gene, variant), '')
    if suffix:
        return f'{variant} [{suffix}]'
    return variant


def figure2():
    s = read_tsv(PREC_SUM)
    axis = read_tsv(AXIS_SUM)
    claims = read_tsv(CLAIMS)

    order = ['protein','splice','biogenesis_fate','motif_disruption_true','motif_context_only','ccre',
             'gtex_true_hit','gtex_lookup_only','ld_context','docking_variant','population','regulatory_priority']
    s = s.set_index('evidence_layer').reindex(order).reset_index()
    cand = num(s['candidate_variants']).to_numpy()
    pop = num(s['population_significant_variants_fdr_lt_0_05_delta_ge_0_10']).to_numpy()
    labels = [pretty_evidence(x) for x in order]

    fig, axs = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    fig.suptitle('Figure 2. Integrated population-functional evidence layers across HRH1-HRH4', fontweight='bold')

    ax = axs[0,0]
    x = np.arange(len(order))
    ax.bar(x, cand)
    ax.set_title('A. Evidence-layer distribution across the final precision-analysis set', fontweight='bold')
    ax.set_ylabel('Number of variants')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha='right')
    ax.spines[['top','right']].set_visible(False)

    ax = axs[0,1]
    ax.bar(x, pop)
    ax.set_title('B. Evidence layers among population-differentiated variants', fontweight='bold')
    ax.set_ylabel('Number of population-differentiated variants')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha='right')
    ax.spines[['top','right']].set_visible(False)

    # Panel C: receptor-specific composition. Keep evidence layers that are directly interpretable.
    fields = ['protein','splice','biogenesis_fate','motif_disruption_true','ccre','gtex_true_hit','gtex_lookup_only','ld_context','docking_variant']
    field_labels = [pretty_evidence(f).replace('\n',' ') for f in fields]
    genes = ['HRH1','HRH2','HRH3','HRH4']
    axis = axis.set_index('gene').reindex(genes)
    bottom = np.zeros(len(genes))
    ax = axs[1,0]
    for f, lab in zip(fields, field_labels):
        vals = num(axis[f]).to_numpy()
        ax.bar(genes, vals, bottom=bottom, label=lab)
        bottom += vals
    ax.set_title('C. Receptor-specific functional evidence composition', fontweight='bold')
    ax.set_ylabel('Evidence-row count')
    ax.legend(ncol=2, frameon=False, loc='upper left')
    ax.spines[['top','right']].set_visible(False)

    raw_counts = claims['manuscript_claim_strength'].value_counts()
    claim_order = [
        'I_single_layer_or_background_candidate',
        'H_multilayer_exploratory_candidate',
        'G_population_divergent_candidate',
        'F_true_motif_disruption_candidate',
        'E_expression_supported_candidate',
        'D_population_divergent_multilayer_candidate',
        'C_population_divergent_motif_cCRE_candidate',
        'B_population_divergent_expression_supported',
        'A_population_divergent_expression_motif_supported',
    ]
    claim_labels = {
        'I_single_layer_or_background_candidate': 'Single-layer or background',
        'H_multilayer_exploratory_candidate': 'Multilayer exploratory',
        'G_population_divergent_candidate': 'Population-differentiated',
        'F_true_motif_disruption_candidate': 'True motif-disruption',
        'E_expression_supported_candidate': 'Expression-supported',
        'D_population_divergent_multilayer_candidate': 'Population-differentiated multilayer',
        'C_population_divergent_motif_cCRE_candidate': 'Population-differentiated motif/cCRE',
        'B_population_divergent_expression_supported': 'Population-differentiated expression-supported',
        'A_population_divergent_expression_motif_supported': 'Population-differentiated expression+motif supported',
    }
    present = [c for c in claim_order if c in raw_counts.index]
    # Preserve any unforeseen class rather than silently discarding it.
    present += [c for c in raw_counts.index if c not in present]
    vals = [int(raw_counts.get(c, 0)) for c in present]
    labs = [claim_labels.get(c, re.sub(r'^[A-Z]_', '', c).replace('_',' ')) for c in present]
    ax = axs[1,1]
    yy = np.arange(len(present))
    ax.barh(yy, vals)
    ax.set_yticks(yy)
    ax.set_yticklabels(labs)
    ax.invert_yaxis()
    ax.set_xlabel('Number of variants')
    ax.set_title('D. Integrated evidence-class distribution', fontweight='bold')
    ax.spines[['top','right']].set_visible(False)

    return save_all(fig, 'Figure_2_integrated_population_functional_evidence_layers_v2')


def figure3():
    axis = read_tsv(AXIS_SUM).set_index('gene').reindex(['HRH1','HRH2','HRH3','HRH4'])
    fields = ['candidate_variants','population_significant_variants','motif_disruption_true','gtex_true_hit']
    labels = ['Variants retained in\nprecision-analysis set','Population-\ndifferentiated','True motif\ndisruption','Direct GTEx/eQTL\nsupport']
    vals = np.column_stack([num(axis[f]).to_numpy() for f in fields])

    fig, ax = plt.subplots(figsize=(10, 7.5), constrained_layout=True)
    fig.suptitle('Figure 3. Receptor-specific pharmacogenomic evidence layers', fontweight='bold')
    maxv = max(vals.max(), 1)
    for i, gene in enumerate(axis.index):
        for j, v in enumerate(vals[i]):
            size = 40 if v == 0 else 40 + 1700 * (float(v) / maxv)
            face = '0.85' if v == 0 else None
            ax.scatter(j, i, s=size, alpha=0.85, facecolor=face, edgecolor='0.25' if v == 0 else None)
            ax.text(j, i, str(int(v)), ha='center', va='center', fontsize=9 if v < 50 else 10,
                    color='black' if v < 50 else 'white', fontweight='bold')
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels)
    ax.set_yticks(range(4))
    ax.set_yticklabels(axis.index)
    ax.invert_yaxis()
    ax.set_xlabel('Evidence layer')
    ax.set_ylabel('Receptor')
    ax.grid(True, linestyle='--', linewidth=0.6, alpha=0.4)
    ax.set_axisbelow(True)
    return save_all(fig, 'Figure_3_receptor_specific_pharmacogenomic_evidence_layers_v2')


def make_panel_e_svg(lookup):
    df = read_tsv(CONTACT)
    df['score'] = num(df['contact_perturbation_score'])
    df = df.sort_values(['score','gene','hgvsp'], ascending=[False, True, True]).reset_index(drop=True)
    labels = []
    for _, r in df.iterrows():
        v = annotate_variant(r['gene'], r['hgvsp'], lookup)
        labels.append(f"{r['gene']}\n{v}")
    receptor_colors = {'HRH1':'#1b9e77','HRH2':'#d95f02','HRH3':'#7570b3','HRH4':'#e7298a'}
    colors = [receptor_colors.get(g, '#777777') for g in df['gene']]

    fig, ax = plt.subplots(figsize=(8, 6.2), constrained_layout=True)
    x = np.arange(len(df))
    ax.bar(x, df['score'], color=colors)
    for y in [40, 60, 80]:
        ax.axhline(y, color='0.5', linestyle='--' if y < 80 else ':', linewidth=0.8)
    ax.set_ylabel('Contact-perturbation priority score')
    ax.set_title('E. Selected HRH coding variants ranked by ligand-contact perturbation proxy', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=55, ha='right')
    ax.spines[['top','right']].set_visible(False)
    handles = [Patch(facecolor=receptor_colors[g], label=g) for g in ['HRH1','HRH2','HRH3','HRH4']]
    ax.legend(handles=handles, title='Receptor', ncol=4, frameon=False, loc='upper right')
    stem = OUT / 'Figure_4_panel_E_contact_perturbation_v2'
    fig.savefig(stem.with_suffix('.svg'), bbox_inches='tight')
    fig.savefig(stem.with_suffix('.pdf'), bbox_inches='tight')
    fig.savefig(Path(str(stem) + '_900dpi.png'), dpi=900, bbox_inches='tight')
    plt.close(fig)
    return stem.with_suffix('.svg'), Path(str(stem) + '_900dpi.png')


def parse_viewbox(root):
    vb = root.attrib.get('viewBox') or root.attrib.get('viewbox')
    if vb:
        x, y, w, h = map(float, vb.split())
        return x, y, w, h
    def n(v):
        return float(re.sub(r'[^0-9.+-]', '', v))
    return 0, 0, n(root.attrib['width']), n(root.attrib['height'])


def combine_svgs(left_path, right_path, out_path):
    ET.register_namespace('', 'http://www.w3.org/2000/svg')
    a = ET.parse(left_path).getroot()
    b = ET.parse(right_path).getroot()
    _, _, aw, ah = parse_viewbox(a)
    _, _, bw, bh = parse_viewbox(b)
    # The locked source A-D SVG already contains its scientific title. Do not add
    # another title in the vector composite; this prevents duplicate title text.
    title_h = 0.0
    target_h = max(ah, bh)
    a_scale = target_h / ah
    b_scale = target_h / bh
    a_w = aw * a_scale
    b_w = bw * b_scale
    gap = 25.0
    total_w = a_w + gap + b_w
    total_h = target_h
    ns = 'http://www.w3.org/2000/svg'
    root = ET.Element(f'{{{ns}}}svg', {
        'viewBox': f'0 0 {total_w:.3f} {total_h:.3f}',
        'width': f'{total_w:.3f}', 'height': f'{total_h:.3f}'
    })
    a.attrib.update({'x':'0', 'y':'0', 'width':str(a_w), 'height':str(target_h)})
    b.attrib.update({'x':str(a_w + gap), 'y':'0', 'width':str(b_w), 'height':str(target_h)})
    root.append(a)
    root.append(b)
    ET.ElementTree(root).write(out_path, encoding='utf-8', xml_declaration=True)


def figure4():
    lookup = bw_lookup()
    panel_e_svg, panel_e_png = make_panel_e_svg(lookup)

    combined_svg = OUT / 'Figure_4_structural_contact_composite_v2.svg'
    combine_svgs(FIG12_SVG, panel_e_svg, combined_svg)

    # Raster/PDF composite for journal upload and visual QC. Crop only the old internal title strip from source A-D.
    left = Image.open(FIG12_PNG).convert('RGB')
    crop_top = int(left.height * 0.055)
    left = left.crop((0, crop_top, left.width, left.height))
    right = Image.open(panel_e_png).convert('RGB')

    fig = plt.figure(figsize=(16, 8.5), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1])
    ax1 = fig.add_subplot(gs[0,0]); ax2 = fig.add_subplot(gs[0,1])
    ax1.imshow(left); ax1.axis('off')
    ax2.imshow(right); ax2.axis('off')
    fig.suptitle('Figure 4. Representative HRH ligand-pocket contacts and coding-variant structural context', fontweight='bold')
    # SVG already produced vector-wise; save other formats only here.
    pdf = OUT / 'Figure_4_structural_contact_composite_v2.pdf'
    png = OUT / 'Figure_4_structural_contact_composite_v2_900dpi.png'
    tif = OUT / 'Figure_4_structural_contact_composite_v2_900dpi.tiff'
    fig.savefig(pdf, bbox_inches='tight')
    fig.savefig(png, dpi=900, bbox_inches='tight')
    fig.savefig(tif, dpi=900, bbox_inches='tight', pil_kwargs={'compression':'tiff_lzw'})
    plt.close(fig)
    return [combined_svg, pdf, png, tif, panel_e_svg]


def figure5():
    lookup = bw_lookup()
    df = read_tsv(E2)
    df['rank_num'] = num(df['rank'])
    df['delta'] = pd.to_numeric(df['delta_mutant_minus_wt_kcal_mol'], errors='coerce')
    df['score'] = pd.to_numeric(df['receptor_fate_ligand_aware_priority_score'], errors='coerce')
    success = df['successful_71d3_redocking'].map(truthy)
    plot = df.loc[success & df['delta'].notna()].copy().sort_values(['score','rank_num'], ascending=[False, True])

    labels = []
    for _, r in plot.iterrows():
        v = annotate_variant(r['gene'], r['protein_change'], lookup)
        labels.append(f"{r['ligand']} | {v}")

    y = np.arange(len(plot))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8.5), sharey=True, constrained_layout=True)
    fig.suptitle('Figure 5. Integrated HRH4 ligand-variant prioritization', fontweight='bold')

    delta = plot['delta'].to_numpy(dtype=float)
    colors = np.where(delta < -0.05, '#2f9db5', np.where(delta > 0.05, '#f0642a', '#c8ced6'))
    ax1.barh(y, delta, color=colors, edgecolor='0.35', linewidth=0.4)
    ax1.axvline(0, color='black', linewidth=0.9)
    ax1.set_yticks(y); ax1.set_yticklabels(labels)
    ax1.invert_yaxis()
    ax1.set_xlabel('Delta Vina affinity, mutant minus WT (kcal/mol)')
    ax1.set_title('A. Paired WT-mutant redocking', fontweight='bold')
    for yi, v in zip(y, delta):
        ax1.text(v + (0.035 if v >= 0 else -0.035), yi, f'{v:+.2f}', va='center',
                 ha='left' if v >= 0 else 'right', fontsize=8)
    ax1.legend(handles=[Patch(facecolor='#2f9db5', label='Stronger predicted mutant binding'),
                        Patch(facecolor='#f0642a', label='Weaker predicted mutant binding'),
                        Patch(facecolor='#c8ced6', label='Negligible change')],
               frameon=False, loc='lower right')
    ax1.spines[['top','right']].set_visible(False)

    scores = plot['score'].to_numpy(dtype=float)
    ax2.barh(y, scores, color=plt.cm.Purples(np.linspace(0.35, 0.9, len(plot)))[::-1], edgecolor='0.45', linewidth=0.4)
    ax2.set_yticks(y); ax2.set_yticklabels(labels)
    ax2.set_xlabel('Ligand-aware receptor-biogenesis/proteostasis score')
    ax2.set_title('B. Ligand-aware receptor-biogenesis/proteostasis score', fontweight='bold')
    for yi, v in zip(y, scores):
        ax2.text(v + 0.08, yi, f'{v:.1f}', va='center', ha='left', fontsize=8)
    ax2.spines[['top','right']].set_visible(False)

    return save_all(fig, 'Figure_5_HRH4_ligand_variant_prioritization_v2')


def qc(outputs):
    matrix = read_tsv(MATRIX)
    axis = read_tsv(AXIS_SUM).set_index('gene')
    motif = read_tsv(MOTIF_TRUE)
    gtex = read_tsv(GTEX_TRUE)
    e2 = read_tsv(E2)
    contact = read_tsv(CONTACT)
    pop = read_tsv(POP_TESTS)

    checks = []
    def ck(name, observed, expected):
        status = 'PASS' if observed == expected else 'FAIL'
        checks.append((name, status, observed, expected))

    ck('precision_matrix_rows', len(matrix), 2011)
    ck('true_motif_rows', len(motif), 14)
    ck('direct_gtex_records', len(gtex), 22)
    ck('direct_gtex_unique_variants', gtex['variant_key'].nunique(), 11)
    ck('population_pairwise_rows', len(pop), 59410)
    ck('contact_proxy_rows', len(contact), 13)
    ck('ligand_aware_rows', len(e2), 16)
    ck('successful_71d3_rows', int(e2['successful_71d3_redocking'].map(truthy).sum()), 15)
    ck('failed_71d3_rows', int(e2['failed_71d3_redocking'].map(truthy).sum()), 1)
    for gene, expected in [('HRH1',678),('HRH2',400),('HRH3',409),('HRH4',524)]:
        ck(f'{gene}_precision_set', int(float(axis.loc[gene,'candidate_variants'])), expected)
    for gene, expected in [('HRH1',0),('HRH2',8),('HRH3',0),('HRH4',6)]:
        ck(f'{gene}_true_motif', int(float(axis.loc[gene,'motif_disruption_true'])), expected)
    for gene, expected in [('HRH1',0),('HRH2',0),('HRH3',0),('HRH4',11)]:
        ck(f'{gene}_direct_gtex', int(float(axis.loc[gene,'gtex_true_hit'])), expected)

    failed_variants = e2.loc[e2['failed_71d3_redocking'].map(truthy), 'protein_change'].tolist()
    ck('failed_variant_identity', ';'.join(failed_variants), 'p.Phe312Val')

    # Figure 5 is a prioritization figure: verify the highest ligand-aware score
    # is the first-ranked displayed hypothesis after sorting by score.
    e2_success = e2.loc[e2['successful_71d3_redocking'].map(truthy)].copy()
    e2_success['score_num'] = pd.to_numeric(e2_success['receptor_fate_ligand_aware_priority_score'], errors='coerce')
    e2_priority = e2_success.sort_values('score_num', ascending=False)
    top = e2_priority.iloc[0]
    ck('figure5_top_priority_identity', f"{top['ligand']}|{top['protein_change']}", 'immepip|p.Tyr319Cys')
    ck('figure5_top_priority_score', float(top['score_num']), 8.5)

    # Publication-facing SVGs must not contain historical figure-number references or internal monitoring terminology.
    prohibited = ['Figure 12', 'Figure 13', 'Precision monitoring', 'claim-strength grading']
    for p in outputs:
        if p.suffix.lower() == '.svg' and p.is_file():
            txt = p.read_text(encoding='utf-8', errors='ignore')
            for token in prohibited:
                ck(f'no_{token.replace(" ","_")}_{p.name}', token in txt, False)

    report = OUT / 'FINAL_FIGURE_QC_v2.txt'
    with report.open('w', encoding='utf-8') as f:
        f.write('PLOS ONE deterministic final-figure QC\n')
        f.write('======================================\n\n')
        for name, status, obs, exp in checks:
            f.write(f'{status}\t{name}\tobserved={obs}\texpected={exp}\n')
        f.write('\nOUTPUT HASHES\n')
        for p in sorted(set(outputs)):
            if p.is_file():
                f.write(f'{sha256(p)}  {p}\n')
    if any(s == 'FAIL' for _, s, _, _ in checks):
        print(report.read_text())
        raise RuntimeError(f'FINAL FIGURE QC FAILED: {report}')
    return report


def main():
    check_inputs()
    outputs = []
    outputs += figure2()
    outputs += figure3()
    outputs += figure4()
    outputs += figure5()
    report = qc(outputs)
    print('STATUS=PASS')
    print(f'OUTPUT_DIR={OUT}')
    print(f'QC_REPORT={report}')
    for p in sorted(set(outputs)):
        if p.is_file():
            print(f'{sha256(p)}  {p}')

if __name__ == '__main__':
    main()
