# RNA-seq Differential Expression Analysis

## Overview

This project presents an RNA-seq differential expression analysis performed on the human gene expression dataset GSE196006.

The analysis was conducted in R to identify differentially expressed genes and characterize the biological pathways and processes associated with the observed expression changes.

## Dataset

- **Dataset:** GSE196006
- **Organism:** Homo sapiens
- **Data type:** RNA-seq raw gene counts
- **Source:** NCBI Gene Expression Omnibus (GEO)

## Analysis Workflow

The analysis included the following steps:
  
  1. Preparation and quality assessment of the count matrix
2. Sample metadata organization
3. RNA-seq normalization using DESeq2
4. Variance-stabilizing transformation (VST)
5. Principal Component Analysis (PCA)
6. Differential expression analysis
7. Volcano plot visualization
8. Identification of upregulated and downregulated genes
9. Gene Ontology (GO) enrichment analysis
10. KEGG pathway enrichment analysis

## Differential Expression Results

A total of:
  
  - **1,837 genes** were identified as upregulated
- **2,715 genes** were identified as downregulated
- **17,101 genes** were not significantly differentially expressed

Differential expression was assessed using adjusted p-value and log2 fold-change thresholds.

## Functional Enrichment Analysis

### GO Biological Process

GO enrichment analysis revealed that the upregulated genes were strongly associated with biological processes related to:
  
  - Chromosome segregation
- Nuclear division
- Mitotic nuclear division
- Organelle fission
- DNA replication
- Sister chromatid segregation
- Cell cycle regulation

### KEGG Pathway Analysis

KEGG enrichment analysis identified several significantly enriched pathways.

Among the upregulated genes, prominent pathways included:
  
  - Cell cycle
- DNA replication
- Cytokine-cytokine receptor interaction
- IL-17 signaling
- Wnt signaling
- p53 signaling
- Homologous recombination

Among the downregulated genes, enriched pathways included:
  
  - Neuroactive ligand-receptor interaction
- Calcium signaling
- cAMP signaling
- Bile secretion
- Retinol metabolism
- Drug metabolism

## Visualizations

The project includes:
  
  - PCA plot
- Volcano plot
- GO enrichment plots
- KEGG enrichment plots

## Tools and Packages

The analysis was performed using R and Bioconductor packages, including:
  
  - DESeq2
- clusterProfiler
- org.Hs.eg.db
- enrichplot
- ggplot2
- pheatmap

## Project Structure

```text
RNA-seq-Differential-Expression/
  │
├── scripts/
  │   └── RNAseq_analysis.R
│
├── results/
  │   ├── DEG_results.csv
│   ├── GO_Upregulated.csv
│   ├── GO_Downregulated.csv
│   ├── KEGG_Upregulated.csv
│   └── KEGG_Downregulated.csv
│
└── figures/
  ├── PCA_Plot.png
├── Volcano_Plot.png
├── GO_Upregulated.png
├── GO_Downregulated.png
├── KEGG_Upregulated.png
└── KEGG_Downregulated.png
```