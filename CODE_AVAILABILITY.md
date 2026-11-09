# Code availability — draft statements

The manuscript currently has a **Data availability** section (PRJNA1465755 for 16S, GSE331109 for
RNA-seq) but **no Code availability**. Add one of the following.

---

## Option A — repository on GitHub (recommended; replace the URL after creating the repo)

**English**

> Code availability
> All custom code used for the transcriptome, machine-learning, microbiota-correlation and
> molecular-dynamics analyses is available at https://github.com/zys-666/BDPrB-ocular-toxicity-zebrafish (archived at
> https://doi.org/10.5281/zenodo.<ID-PENDING>). Third-party software is cited in the Methods:
> HISAT2, StringTie, DESeq2, WGCNA, clusterProfiler, STRING/cytoHubba, Cytoscape, fastshap,
> Trimmomatic, cutadapt, DADA2, LEfSe, PICRUSt2, AutoDock Vina, GROMACS and gmx_MMPBSA.

**中文**

> 代码可用性
> 本研究转录组、机器学习、菌群相关性及分子动力学分析所用的全部自编代码见
> https://github.com/zys-666/BDPrB-ocular-toxicity-zebrafish（存档于 https://doi.org/10.5281/zenodo.<编号>）。
> 第三方软件已在方法部分引用：HISAT2、StringTie、DESeq2、WGCNA、clusterProfiler、
> STRING/cytoHubba、Cytoscape、fastshap、Trimmomatic、cutadapt、DADA2、LEfSe、
> PICRUSt2、AutoDock Vina、GROMACS 与 gmx_MMPBSA。

---

## Option B — no repository, on request

**English**

> Code availability
> The custom analysis scripts (differential-expression, WGCNA, the 130-combination feature
> selection, cohort validation, microbiota correlations and trajectory analysis) are available
> from the corresponding author on reasonable request.

---

## Before submission (status as of 2026-11-09)

1. **Repository created and populated** — https://github.com/zys-666/BDPrB-ocular-toxicity-zebrafish
2. **Mint a DOI** through Zenodo (link the GitHub account, then release a tag such as `v1.0`) and
   put the DOI in the statement — journals increasingly require an archived snapshot, not just a URL.
3. **Do not commit** `data/`, `out/`, `tables/`, `figures/`; `.gitignore` already excludes them, and
   the raw data belong to GSE331109 / PRJNA1465755.
4. **Outstanding before submission**: the 16S pipeline (Figure 5A-Q, Figure S10) and the
   wet-lab statistics for Figures 1, 3 and 6 are not yet in the repository. Either add them
   or use wording that does not imply full coverage of every figure.
5. Check that the author list in `LICENSE` and the repository `README.md` matches the manuscript.
