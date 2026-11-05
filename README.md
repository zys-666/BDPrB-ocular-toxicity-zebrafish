# Analysis code — BDPrB-induced ocular injury in zebrafish

Reproducible analysis code for the manuscript **"Gut microbiota-associated retinoid deficiency and a decompensated antioxidant response drive BDPrB-induced ocular injury in zebrafish"**.

The repository contains the custom code only. Third-party tools (DESeq2, HISAT2, StringTie, WGCNA, clusterProfiler, AutoDock Vina, GROMACS, gmx_MMPBSA, DADA2) are cited in the manuscript and are not redistributed here.

---

## Contents

| Folder | Analysis | Manuscript section |
|---|---|---|
| `01_transcriptome/` | Differential expression, heatmap, volcano plot, GO/KEGG circos plots, WGCNA, module enrichment | 2.8, Figure 2A–J, S1 |
| `02_machine_learning/` | 130 method–learner combinations, selection stability, atf3 single-gene model, cohort validation, external validation, SHAP | 2.8, Figure 2K–U, S2–S7 |
| `03_microbiota/` | Genus abundance vs retinal gene expression correlations | 2.10, Table S2 |
| `04_docking_md/` | GROMACS `.mdp` inputs, Vina/MMPBSA inputs, trajectory-analysis scripts | 2.9, Figure 4, S8–S9 |
| `05_figures/` | Figure 2 panel assembly and supplementary figures S1–S7 | Figure 2, S1–S7 |

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

`logistf`, `dcurves` and `pmsampsize` are **not** required: Firth penalised logistic regression, decision-curve analysis and the Riley sample-size criterion are implemented from first principles in `02_machine_learning/15_firth_fun.R` and the calling scripts, so the pipeline runs on a machine without outbound network access.

Python ≥ 3.10 with `MDAnalysis`, `numpy`, `pandas` is needed only for `04_docking_md/scripts/`.

---

## Data layout

All scripts read from `DATA_ROOT` (set by `config.R`) and write to `./out`, `./tables`, `./figures`. The expected layout mirrors the original project tree:

```
data/
├─ 3.转录组/
│  ├─ 1.火山图/data.txt                       # volcano input (one row per gene × comparison)
│  ├─ 1.火山图/*_all_geneSymbol.xls           # six DESeq2 result tables
│  ├─ 2.热图/*_all_geneSymbol.xls
│  ├─ 3.GO/多组GO前5通路环形图/Pathway_GO.xlsx    # GOtools export (Pathway, Abbre sheets)
│  ├─ 4.KEGG/多组KEGG前5通路环形图/Pathway_KEGG.xlsx
│  └─ 5.WGCNA/输入/{表达谱.csv, 样本信息.csv}      # 9 samples × top-25 % variable genes
├─ 9.菌群测序/菌群数量和核心基因表达量的相关性/     # genus / qPCR matrices
└─ 2.机器学习/ …                              # see 02_machine_learning/README.md
```

Public datasets:

- **GSE331109** — RNA-seq of this study (deposited; see Data availability)
- **PRJNA1465755** — 16S rRNA amplicon sequences
- **GSE234646** — Müller glia ablation time series (screening platform)
- **GSE310406** — PFNA versus vehicle, zebrafish eye
- **GSE284728** — optic nerve crush versus uninjured retinal ganglion cells

The BDPrB group means and the WGCNA network object that are not part of GSE331109 are available from the corresponding author on reasonable request.

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
Rscript 02_machine_learning/10_ml130.R                    # 130 combinations, ~15 min
ML130_K=5 ML130_REP=4 Rscript 02_machine_learning/11_ml130_rep.R   # ~14 min
Rscript 02_machine_learning/12_atf3_final.R
Rscript 02_machine_learning/19_validation_final.R
Rscript 02_machine_learning/22_external_run.R
Rscript 02_machine_learning/23_GSE284728.R
Rscript 02_machine_learning/50_shap.R

# 3. figures
Rscript 05_figures/31_figure2_panels.R
Rscript 05_figures/80_supplementary_figures.R
```

`04_docking_md/` and `03_microbiota/` are documented in their own README files.

---

## Design notes

1. **Feature selection is inside the training folds.** No test-fold information enters selection or scaling, which is why the honest cross-validated AUC of the multi-gene models is 0.56–0.80 rather than the ~1.0 obtained by resubstitution.
2. **Leave-one-run-out folds** keep all samples of a sequencing run on the same side of the split.
3. **107 → 63 is platform coverage, not gene filtering**: only 63 of the 107 candidate genes have probes on the Affymetrix GPL1319 array used by the screening dataset.
4. **Two resampling schemes are reported** (leave-one-run-out: atf3 retained in 388/390 runs; 4 × 5-fold: 2346/2600). The identity of the runner-up gene is not stable between them; only atf3 is.
5. **Firth penalised logistic regression** is used for all small-sample contrasts because complete separation makes ordinary maximum-likelihood estimates unreliable.
6. **SHAP is computed on the best multi-gene combination**, where atf3 ranks third by mean |SHAP|; its primacy rests on selection frequency, single-gene AUC and external replication, not on the SHAP ranking.

---

## License

MIT — see `LICENSE`.
