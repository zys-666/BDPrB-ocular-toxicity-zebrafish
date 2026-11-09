# Analysis code — BDPrB-induced ocular injury in zebrafish

Reproducible analysis code accompanying the manuscript **"Gut microbiota-associated
retinoid deficiency and a decompensated antioxidant response drive BDPrB-induced
ocular injury in zebrafish"**.

This repository contains **custom code only**. Third-party tools (HISAT2, StringTie,
DESeq2, WGCNA, clusterProfiler, STRING/cytoHubba, Cytoscape, fastshap, Trimmomatic,
cutadapt, DADA2, SILVA, LEfSe, PICRUSt2, AutoDock Vina, GROMACS, gmx_MMPBSA) are
cited in the manuscript and are not redistributed here.

> **Read this first.** No data are stored in this repository (see *Data layout*).
> Sections marked **missing** below are analyses reported in the manuscript for
> which no script is included in this release.

---

## Coverage against the manuscript

| Manuscript item | Folder | Status |
|---|---|---|
| Figure 2A–J (DEG, heatmap, volcano, GO/KEGG circos, WGCNA) | `01_transcriptome/` | provided |
| Figure 2K–U (130 combinations, stability, atf3 model, SHAP) | `02_machine_learning/` | provided |
| Figure 2 panel assembly | `05_figures/31_figure2_panels.R` | provided |
| Figure 4A–T (docking, MD, MM-PBSA) | `04_docking_md/` | provided |
| Figure 5 (16S diversity, PCoA, PERMANOVA, LEfSe, PICRUSt2) | `03_microbiota/` | **missing** |
| Figure 5P–Q (genus–gene correlations) | `03_microbiota/` | provided |
| Figure 1B–R (behaviour, morphometry, histology) | — | **missing** |
| Figure 3B–R (ELISA, RT-qPCR) | — | **missing** |
| Figure 6B–J (recovery experiment) | — | **missing** |
| Figure S1 (DEG counts) | `01_transcriptome/` | provided |
| Figure S2 (blue module ∩ PPI hubs) | `02_machine_learning/40_build_107.R` | provided |
| Figure S3 (platform coverage) | `02_machine_learning/40_build_107.R` | provided |
| Figure S4 (atf3 in the BDPrB cohort) | `02_machine_learning/12_atf3_final.R` | provided |
| Figure S5 (effect sizes across cohorts) | `02_machine_learning/22_external_run.R` | provided |
| Figure S6 (bootstrap optimism) | `02_machine_learning/19_validation_final.R` | provided |
| Figure S7 (WGCNA soft threshold) | `01_transcriptome/05_wgcna.R` | provided |
| Figure S8 (pLDDT, N-terminal truncation) | `04_docking_md/` | partially |
| Figure S9 (MD supplementary) | `04_docking_md/` | provided |
| Figure S10 (16S supplementary) | `03_microbiota/` | **missing** |
| Table S1 (primers) | — | no code required |
| Table S2 (genus–gene correlations) | `03_microbiota/01_...R` | provided |

**Note on `04_docking_md/`.** Receptor and ligand structures used for docking
and for the MD starting complexes are included in `04_docking_md/vina/`
(`receptor_zf_ALDH1A2_26-518.pdb`, `ligand_BDPrB.pdb` / `.sdf`,
`KyA33_reference_pose.pdb`, the docked poses and `complex_BDPrB.pdb`). The GROMACS
`.mdp` inputs and the trajectory-analysis scripts are in `04_docking_md/mdp/` and
`04_docking_md/scripts/`; GROMACS selection files are in
`04_docking_md/gromacs_selections/`. The AutoDock Vina configuration files and
`.pdbqt` conversions are **not** included — regenerate them with
`04_docking_md/scripts/dock_prepare_receptor.py` and the parameters given in
Methods 2.9. Trajectories (`.xtc`, `.trr`, `.tpr`, `.edr`) are excluded by
`.gitignore` and are available from the corresponding author.

**Note on the missing wet-lab statistics.** Figures 1, 3 and 6 were analysed in
SPSS v27 by nested one-way ANOVA with tank as a random factor followed by Tukey
HSD (Methods 2.11); those analyses are not scripted here.

---

## Folder contents

```
01_transcriptome/        six R scripts: volcano, DEG filtering + heatmap, GO circos,
                         KEGG circos, WGCNA, module GO/KEGG
02_machine_learning/     preparation, the 107-gene candidate build, 130-combination
                         selection, repeated resampling, atf3 single-gene model,
                         Firth helpers, cohort and external validation, SHAP
03_microbiota/           genus abundance vs retinal transcript correlations
04_docking_md/           mdp/ (GROMACS inputs) · vina/ (receptor, ligand, poses)
                         gromacs_selections/ (selection files) · scripts/
05_figures/              Figure 2 panel assembly, supplementary figures
06_wetlab_stats/         README only - the statistics for Figures 1, 3 and 6 are not scripted
00_docs/                 note explaining which data are not distributed with this repository
archive/                 earlier iterations of the docking/MD scripts, kept for
                         provenance, not needed to reproduce any figure
```

---

## Requirements

R ≥ 4.6.0 with:

```
caret 7.0.1        glmnet 5.0          randomForest 4.7.1.2   e1071 1.7.17
xgboost 3.2.1.1    gbm 2.3.1           nnet 7.3.20            rpart 4.1.27
pROC 1.19.0.1      fastshap 0.1.1      readxl 1.4.5           WGCNA 1.74
clusterProfiler 4.20.0   org.Dr.eg.db 3.22.0   openxlsx   ggplot2 4.0.3
patchwork 1.3.2    svglite 2.2.2       geomtextpath  cowplot  corrplot
```

`logistf`, `dcurves` and `pmsampsize` are **not** required: Firth penalised logistic
regression, decision-curve analysis and the Riley sample-size criterion are
implemented from first principles in `02_machine_learning/15_firth_fun.R` and the
calling scripts, so the pipeline runs without outbound network access.

Python ≥ 3.10 with `MDAnalysis`, `numpy`, `pandas` is needed only for
`04_docking_md/scripts/`.

---

## Data layout

**Data are not included.** Every script reads from `DATA_ROOT` (set in `config.R`)
and writes to `./out`, `./tables` and `./figures`; `.gitignore` excludes `data/`,
`out/`, `tables/` and `figures/`. To run anything you must first restore the
project tree the scripts were written against:

```
data/
├─ 3.转录组/
│  ├─ 1.火山图/data.txt                       # volcano input (gene × comparison)
│  ├─ 1.火山图/*_all_geneSymbol.xls           # six DESeq2 result tables
│  ├─ 2.热图/*_all_geneSymbol.xls
│  ├─ 3.GO/多组GO前5通路环形图/Pathway_GO.xlsx    # GOtools export
│  ├─ 4.KEGG/多组KEGG前5通路环形图/Pathway_KEGG.xlsx
│  └─ 5.WGCNA/输入/{表达谱.csv, 样本信息.csv}      # 9 samples × top-25 % variable genes
├─ 9.菌群测序/菌群数量和核心基因表达量的相关性/     # genus / qPCR matrices
└─ 2.机器学习/ …
```

Public datasets:

| Accession | Content |
|---|---|
| **GSE331109** | RNA-seq of this study (whole eye, 12 samples) |
| **PRJNA1465755** | 16S rRNA amplicon sequences |
| GSE234646 | Müller glia ablation time series (screening platform) |
| GSE310406 | PFNA versus vehicle, zebrafish eye |
| GSE284728 | optic nerve crush versus uninjured retinal ganglion cells |

The BDPrB group means and the WGCNA network object that are not part of GSE331109
are available from the corresponding author on reasonable request.

---

## Running

```bash
cd <repo>

# 1. transcriptome
Rscript 01_transcriptome/02_deg_filter_heatmap.R
Rscript 01_transcriptome/01_volcano_six_comparisons.R
Rscript 01_transcriptome/03_go_circos.R
Rscript 01_transcriptome/04_kegg_circos.R
Rscript 01_transcriptome/05_wgcna.R
Rscript 01_transcriptome/06_module_go_kegg.R

# 2. machine learning (≈35 min total)
Rscript 02_machine_learning/01_prepare.R
Rscript 02_machine_learning/40_build_107.R
Rscript 02_machine_learning/10_ml130.R                              # ~15 min
ML130_K=5 ML130_REP=4 Rscript 02_machine_learning/11_ml130_rep.R    # ~14 min
Rscript 02_machine_learning/12_atf3_final.R
Rscript 02_machine_learning/19_validation_final.R
Rscript 02_machine_learning/22_external_run.R
Rscript 02_machine_learning/23_GSE284728.R
Rscript 02_machine_learning/50_shap.R

# 3. figures
Rscript 05_figures/31_figure2_panels.R
Rscript 05_figures/80_supplementary_figures.R
```

`03_microbiota/` and `04_docking_md/` are documented in their own README files.

---

## Design notes

1. **Feature selection is inside the training folds.** No test-fold information enters
   selection or scaling, which is why the honest cross-validated AUC of the multi-gene
   models is 0.56–0.80 rather than the ~1.0 obtained by resubstitution.
2. **Leave-one-run-out folds** keep all samples of a sequencing run on the same side
   of the split.
3. **107 → 63 is platform coverage, not gene filtering**: only 63 of the 107 candidate
   genes have probes on the Affymetrix GPL1319 array used by the screening dataset.
4. **Two resampling schemes are reported** (leave-one-run-out: *atf3* retained in
   388/390 runs; 4 × 5-fold: 2346/2600). The identity of the runner-up gene is not
   stable between them; only *atf3* is.
5. **Firth penalised logistic regression** is used for all small-sample contrasts
   because complete separation makes ordinary maximum-likelihood estimates unreliable.
6. **SHAP is computed on the best multi-gene combination**, where *atf3* ranks third by
   mean |SHAP|; its primacy rests on selection frequency, single-gene AUC and external
   replication, not on the SHAP ranking.
7. **WGCNA used an unsigned network with soft-threshold power = 1**, returned by the
   automatic routine; the module partition is reported as exploratory (see manuscript
   Methods 2.8 and Figure S7).

---

## License

MIT — see `LICENSE`.
