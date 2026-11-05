# 03_microbiota

Code for the gut-microbiota section (Methods 2.10).

## Provided here

| Script | Purpose | Manuscript |
|---|---|---|
| `01_genus_expression_correlation_heatmap.R` | Spearman correlation between genus relative abundance and retinal gene expression across the twelve pools; BH correction; heatmap | Table S2, Figure S10 |
| `02_genus_expression_butterfly.R` | Butterfly plot of the same correlations, used for the visual summary | Figure S10 |

Both scripts read the genus and qPCR matrices from `DATA_ROOT/9.菌群测序/菌群数量和核心基因表达量的相关性/`
and write to `WORK/out/`.

## Upstream 16S processing (not in this repository)

The amplicon pipeline was run with standard tools, and the exact parameters are given in Methods 2.10:

```
Trimmomatic v0.33  →  cutadapt 1.9.1  →  DADA2  →  SILVA 138 classification
```

Diversity and composition statistics (Methods 2.10) used:

- alpha diversity: Kruskal–Wallis with 9,999 label permutations, Jonckheere–Terpstra for monotonic
  dose responses, Dunn's test with Benjamini–Hochberg correction;
- beta diversity: Bray–Curtis, binary Jaccard, unweighted and weighted UniFrac; PCoA, PERMANOVA
  (9,999 permutations), ANOSIM, PERMDISP with Lingoes correction;
- biomarkers: LEfSe (LDA > 4, P < 0.05);
- predicted function: PICRUSt2, collapsed to KEGG level-3 pathways.

These steps are third-party software driven by command-line parameters rather than custom code.
If the review asks for the exact command lines, they can be added here verbatim; they are not
reconstructed in this repository because they were executed on the sequencing provider's
infrastructure.
