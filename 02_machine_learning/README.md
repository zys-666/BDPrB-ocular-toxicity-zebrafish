# 02_machine_learning

Feature selection and model evaluation for Figure 2K–U and Figures S2–S7.

## Pipeline

| Step | Script | Output |
|---|---|---|
| Build the analysis object (132 samples, 66 vs 66; 63 measurable genes) | `01_prepare.R` | `out/prep.rds` |
| Rebuild the 107-gene candidate panel (WGCNA ∩ PPI) and check platform coverage | `40_build_107.R` | `tables/univariate_auc_107.csv`, `tables/discovery_107_auc.csv` |
| **130 combinations × leave-one-run-out CV** | `10_ml130.R` | `out/ml130_preds.rds`, `tables/ml130_combo_performance.csv`, `tables/ml130_gene_frequency.csv` |
| **130 combinations × 4 repeats × 5 folds** (stability) | `11_ml130_rep.R` | `out/ml130_preds_rep.rds`, `tables/ml130rep_gene_frequency.csv` |
| atf3 single-gene model (AUC, permutation, calibration, DCA, OR) | `12_atf3_final.R` | `tables/atf3_single_gene_model.csv`, `out/atf3_model.rds` |
| BDPrB cohort validation (4 groups, Firth, exact permutation) | `19_validation_final.R` | `tables/validation_final_*.csv` |
| External validation GSE310406 | `22_external_run.R` | `tables/external_GSE310406_atf3.csv`, `tables/external_validation_stats.csv` |
| External validation GSE284728 | `23_GSE284728.R` | `tables/external_GSE284728_stats.csv` |
| SHAP on the best combination | `50_shap.R` | `tables/shap_importance_best.csv` |
| Firth penalised logistic regression (self-implemented) | `15_firth_fun.R` | — |
| 130-combination method library | `mime_functions.R` | — |
| Plotting theme | `style.R` | — |

Paths come from `config.R` at the repository root; run every script from the root:

```bash
Rscript 02_machine_learning/10_ml130.R
```

## Inputs expected under `DATA_ROOT`

```
2.机器学习/训练/输入/normalize.txt        # 9,324 genes × 138 samples, log2-normalised
2.机器学习/训练/输入/interGenes.txt       # genes represented on GPL1319 (63 of the 107)
2.机器学习/训练/blue.csv                  # WGCNA blue-module genes (6,453)
2.机器学习/训练/HC和LC的ppi的交集基因.xls   # 152 PPI hub genes exported from CytoHubba
1.归一化去批次/训练集/sample_info.csv      # sample → group / run / batch
130机器学习组合/refer.txt                  # the 130 combination definitions
2.machine_learning/validation/expression_norm_log2.csv   # BDPrB cohort (12 samples, corrected matrix)
ext/GSE310406_PFNA-fpkm.xlsx
ext/GSE284728_BulkRNAseq_raw_counts.tsv.gz
```

`GSE234646` is downloaded from GEO and pre-processed into `normalize.txt`; the probe→gene step used
`1.归一化去批次/探针矩阵转基因矩阵.R` (not redistributed, because it depends on the Affymetrix CDF file).

## Design notes

1. **Feature selection sits inside the training folds.** No test-fold information enters selection or scaling; this is why the honest cross-validated AUC is 0.56–0.80 rather than the ~1.0 obtained by resubstitution.
2. **Leave-one-run-out folds** keep all samples of a sequencing run on the same side of the split.
3. **107 → 63 is platform coverage**, not gene filtering.
4. **Two resampling schemes are reported**: leave-one-run-out (atf3 selected 388/390 = 99.5 %) and 4 × 5-fold (2346/2600 = 90.2 %). The identity of the runner-up gene is *not* stable between them; only atf3 is.
5. **Firth penalised logistic regression** is used for every small-sample contrast; `logistf` is not needed.
6. **SHAP** is computed on the best multi-gene combination, where atf3 ranks third by mean |SHAP|; its primacy rests on selection frequency, single-gene AUC and external replication, not on the SHAP ranking.
