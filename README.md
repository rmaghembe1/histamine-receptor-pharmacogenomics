# Histamine receptor population pharmacogenomics

Public reproducibility release for:

**Population pharmacogenomics of human histamine receptors reveals receptor-specific architectures and ligand-dependent H4 receptor docking shifts**

## Scope

This repository contains code, derived result tables, final manuscript figure sources, and provenance/QC records for a harmonized analysis of HRH1, HRH2, HRH3, and HRH4.

The study integrates public population variation, regulatory evidence, expression/eQTL evidence, GPCR topology, ligand-pocket context, and focused mutant-versus-wild-type redocking.

## Frozen analytical results

- Final precision-analysis set: 2,011 variants
  - HRH1: 678
  - HRH2: 400
  - HRH3: 409
  - HRH4: 524
- Population pairwise tests: 59,410
- Explicit motif-disruption variants after indel-aware correction: 14
  - HRH2: 8
  - HRH4: 6
- Direct GTEx/eQTL-supported variants: 11, all HRH4
- Ligand-aware receptor-fate/proteostasis rows: 16
- Successful paired mutant-versus-WT redocking rows: 15
- Failed redocking candidate: p.Phe312Val

## Repository structure

- `code/` - analysis and final figure-generation scripts
- `data/precision/` - final precision-analysis matrix and summaries
- `data/regulatory/` - corrected motif and GTEx/eQTL evidence
- `data/population/` - population-functional matrix and pairwise tests
- `data/structural/` - redocking, contact-perturbation and receptor-fate outputs
- `figures/` - final manuscript figure sources
- `provenance/` - frozen QC and SHA-256 records

## Data provenance

The analysis uses public resources and derived summary data. No private participant-level dataset is distributed here.

## Interpretation boundary

The computational and population evidence prioritizes testable hypotheses. It does not establish clinical actionability, causal drug-response effects, or experimentally validated receptor phenotypes.

## Figure 1 disclosure

The illustrative workflow schematic in Figure 1 was generated from an author-defined scientific brief using OpenAI image-generation assistance and was reviewed against the documented study workflow. Quantitative Figures 2-5 were generated deterministically from frozen analysis tables.

## Version

v1.0.0 - 2026-09-19

## Citation

Please cite the associated manuscript and the archived repository DOI once assigned.
