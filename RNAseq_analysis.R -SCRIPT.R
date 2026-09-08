# Read the RNA-seq dataset
data <- read.csv("GSE196006_raw_counts.csv")

# Check the first 6 rows
head(data)

# Check the dimensions of the dataset
dim(data)

# Check the structure of the dataset
str(data)
getwd()
list.files()
file.access("GSE196006_raw_counts.csv", 4)
data <- read.csv("GSE196006_raw_counts.csv")
data <- read.csv("C:/Users/YourName/Downloads/GSE196006_raw_counts.csv")
file.info("GSE196006_raw_counts.csv")
data <- read.csv(file.choose())
head(data)
dim(data)
colnames(data)
gene_ids <- data$X

counts <- data[, -1]

rownames(counts) <- gene_ids
str(counts[, 1:5])
head(counts[, 1:5])
colnames(counts)
gene_ids <- data$X

counts <- data[, -1]

rownames(counts) <- gene_ids
head(counts[, 1:5])
samples <- colnames(counts)
samples
sample_info <- data.frame(
  sample = samples
)

sample_info$individual <- sub("_.*", "", samples)
sample_info$cell <- sub(".*_", "", samples)
sample_info$individual <- sub("_.*", "", samples)

sample_info$cell <- sub(".*_([A-Z])[0-9]_.*", "\\1", samples)

sample_info$time <- sub(".*_([A-Z])([0-9])_.*", "\\2", samples)

sample_info
# Check the dimensions of the count matrix
dim(counts)

# See the first few column names
head(colnames(counts))

# Compare count matrix samples with sample_info
all(colnames(counts) == sample_info$sample)
sample_info$time <- factor(
  sample_info$time,
  levels = c("0", "7")
)
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ time
)
library(DESeq2)
keep <- rowSums(counts(dds) >= 10) >= 3

dds <- dds[keep, ]
dim(dds)
sample_info$time <- factor(
  sample_info$time,
  levels = c("0", "7")
)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dim(dds)
dds <- DESeq(dds)
res <- results(dds)
summary(res)
dds <- DESeq(dds)
head(res)
res_ordered <- res[order(res$padj), ]

head(res_ordered, 20)
order(res$padj)
deg <- subset(
  res,
  padj < 0.05 & abs(log2FoldChange) >= 1
)
dim(deg)
res_ordered <- res[order(res$padj), ]
head(res_ordered, 20)

deg <- subset(
  res,
  padj < 0.05 & abs(log2FoldChange) >= 1
)

dim(deg)
head(deg)
up_genes <- deg[deg$log2FoldChange > 0, ]

down_genes <- deg[deg$log2FoldChange < 0, ]
nrow(up_genes)
nrow(down_genes)
up_genes <- up_genes[order(up_genes$padj), ]
down_genes <- down_genes[order(down_genes$padj), ]
head(up_genes, 10)
head(down_genes, 10)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ individual + time
)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds <- DESeq(dds)

res <- results(dds, contrast = c("time", "7", "0"))

summary(res)
vsd <- vst(dds, blind = FALSE)

plotPCA(vsd, intgroup = "time")
pcaData <- plotPCA(vsd, intgroup = c("time", "individual"), returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

library(ggplot2)

ggplot(pcaData, aes(PC1, PC2, color = time)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_minimal()
res_df <- as.data.frame(res)

res_df$significant <- ifelse(
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1,
  "DEG",
  "Not DEG"
)
+library(ggplot2)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Day 7 vs Day 0",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )
library(ggplot2)

res_df$group <- ifelse(
  res_df$padj < 0.05 & res_df$log2FoldChange >= 1,
  "Upregulated",
  ifelse(
    res_df$padj < 0.05 & res_df$log2FoldChange <= -1,
    "Downregulated",
    "Not significant"
  )
)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = group)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(
    values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not significant" = "grey"
    )
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Day 7 vs Day 0",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Gene status"
  )
res_ordered <- res_df[order(res_df$padj), ]
top20 <- head(res_ordered, 20)

top20
library(biomaRt)
install.packages("biomaRt")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_names <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = rownames(res_df),
  mart = mart
)
BiocManager::install("biomaRt")
library(biomaRt)
mart <- useMart(
  "ensembl",
  dataset = "hsapiens_gene_ensembl"
)
gene_names <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = rownames(res_df),
  mart = mart
)
res_df$ensembl_gene_id <- rownames(res_df)

res_annotated <- merge(
  res_df,
  gene_names,
  by = "ensembl_gene_id",
  all.x = TRUE
)

res_annotated <- res_annotated[order(res_annotated$padj), ]

head(res_annotated, 20)
library(biomaRt)
mart <- useMart(
  "ensembl",
  dataset = "hsapiens_gene_ensembl"
)
gene_names <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = rownames(res_df),
  mart = mart
)
library(biomaRt)

mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)
gene_names <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = rownames(res_df),
  mart = mart
)
head(gene_names)
head(counts)
# Check sample names
colnames(counts)
sample_info$sample
# Check if they are in the same order
all(colnames(counts) == sample_info$sample)
library(DESeq2)
sample_info$individual <- factor(sample_info$individual)
sample_info$cell <- factor(sample_info$cell)
sample_info$time <- factor(sample_info$time)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ individual + cell + time
)
table(sample_info$individual, sample_info$cell)
table(sample_info$individual, sample_info$time)
table(sample_info$cell, sample_info$time)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ individual + time
)
dds
design(dds)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
nrow(dds)
summary(keep)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name = "time_7_vs_0")
summary(res)
head(res)
res_sig <- res[
  !is.na(res$padj) &
    res$padj < 0.05 &
    abs(res$log2FoldChange) >= 1,
]
nrow(res_sig)
head(res_sig[order(res_sig$padj), ])
res_df <- as.data.frame(res)
res_df$ensembl_gene_id <- rownames(res_df)

res_df <- merge(
  res_df,
  gene_names,
  by = "ensembl_gene_id",
  all.x = TRUE
)
head(res_df)
sum(!is.na(res_df$external_gene_name))
head(res_df)
sum(!is.na(res_df$external_gene_name))
res_df$significant <- "Not significant"

res_df$significant[
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    res_df$log2FoldChange >= 1
] <- "Upregulated"

res_df$significant[
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    res_df$log2FoldChange <= -1
] <- "Downregulated"
table(res_df$significant)
table(res_df$significant)
library(ggplot2)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  theme_minimal() +
  labs(
    title = "Differential Gene Expression: Time 7 vs Time 0",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Expression"
  )
vsd <- vst(dds, blind = FALSE)
library(ggplot2)

pca_data <- plotPCA(vsd, intgroup = c("time"), returnData = TRUE)

percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(x = PC1, y = PC2, color = time)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of RNA-seq Samples")
top_genes <- head(
  rownames(res[order(res$padj), ]),
  50
)
vsd_mat <- assay(vsd)[top_genes, ]
gene_symbols <- gene_names$external_gene_name[
  match(top_genes, gene_names$ensembl_gene_id)
]

gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- top_genes[
  is.na(gene_symbols) | gene_symbols == ""
]

rownames(vsd_mat) <- make.unique(gene_symbols)
annotation_col <- data.frame(
  time = sample_info$time
)

rownames(annotation_col) <- sample_info$sample
install.packages("pheatmap")
library(pheatmap)

pheatmap(
  vsd_mat,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "Top 50 Differentially Expressed Genes"
)
pcaData <- plotPCA(vsd, intgroup = c("cell", "time"), returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = cell, shape = factor(time))) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot") +
  theme_minimal()
pcaData <- plotPCA(
  vsd,
  intgroup = c("cell", "time"),
  returnData = TRUE
)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(
  PC1, PC2,
  color = cell,
  shape = factor(time)
)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot") +
  theme_minimal()
library(ggplot2)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(
    values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not significant" = "grey"
    )
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10 Adjusted p-value"
  ) +
  theme_minimal()
up_genes <- res_df %>%
  filter(significant == "Upregulated") %>%
  pull(gene)

down_genes <- res_df %>%
  filter(significant == "Downregulated") %>%
  pull(gene)

length(up_genes)
length(down_genes)
library(dplyr)
up_genes <- res_df %>%
  filter(significant == "Upregulated") %>%
  pull(gene)
colnames(res_df)
colnames(res_df)
[1] "baseMean" "log2FoldChange" "lfcSE" "stat" "pvalue" "padj" "significant"
head(rownames(res_df))
head(res_df)up_genes <- res_df %>%
  filter(
    significant == "Upregulated",
    grepl("^ENSG", ensembl_gene_id)
  ) %>%
  pull(ensembl_gene_id)

down_genes <- res_df %>%
  filter(
    significant == "Downregulated",
    grepl("^ENSG", ensembl_gene_id)
  ) %>%
  pull(ensembl_gene_id)
length(up_genes)
length(down_genes)
table(res_df$significant)
head(
  res_df[res_df$significant == "Upregulated",
         c("ensembl_gene_id", "external_gene_name", "log2FoldChange", "padj")],
  20
)
up_genes <- res_df %>%
  filter(significant == "Upregulated") %>%
  pull(ensembl_gene_id)

down_genes <- res_df %>%
  filter(significant == "Downregulated") %>%
  pull(ensembl_gene_id)
length(up_genes)
length(down_genes)
head(up_genes)
if (!requireNamespace("BiocManager", quietly = TRUE))
  

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
BiocManager::install("enrichplot")
library(enrichplot)
library(clusterProfiler)
up_entrez <- bitr(
  up_genes,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_entrez <- bitr(
  down_genes,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
head(up_entrez)
head(down_entrez)
ego_up <- enrichGO(
  gene = up_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
ego_down <- enrichGO(
  gene = down_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
head(as.data.frame(ego_up))
head(as.data.frame(ego_down))
dotplot(
  ego_up,
  showCategory = 15,
  title = "GO Biological Process - Upregulated Genes"
)
dotplot(
  ego_down,
  showCategory = 15,
  title = "GO Biological Process - Downregulated Genes"
)
head(as.data.frame(ego_up), 20)
write.csv(
  as.data.frame(ego_up),
  "GO_Upregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_down),
  "GO_Downregulated.csv",
  row.names = FALSE
)
ekegg_up <- enrichKEGG(
  gene = up_entrez$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
ekegg_down <- enrichKEGG(
  gene = down_entrez$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
head(as.data.frame(ekegg_up), 15)
head(as.data.frame(ekegg_down), 15)
dotplot(
  ekegg_up,
  showCategory = 15,
  title = "KEGG Pathways - Upregulated Genes"
)
dotplot(
  ekegg_down,
  showCategory = 15,
  title = "KEGG Pathways - Downregulated Genes"
)
write.csv(
  as.data.frame(ekegg_up),
  "KEGG_Upregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ekegg_down),
  "KEGG_Downregulated.csv",
  row.names = FALSE
)
list.files()
write.csv(
  res_df,
  "DEG_results.csv",
  row.names = FALSE
)
list.files()
ggsave(
  "Volcano_Plot.png",
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  "PCA_Plot.png",
  width = 8,
  height = 6,
  dpi = 300
)
p_go_up <- dotplot(
  ego_up,
  showCategory = 15,
  title = "GO Biological Process - Upregulated Genes"
)

p_go_up
ggsave(
  "GO_Upregulated.png",
  p_go_up,
  width = 8,
  height = 6,
  dpi = 300
)
list.files()
p_go_down <- dotplot(
  ego_down,
  showCategory = 15,
  title = "GO Biological Process - Downregulated Genes"
)

p_go_down

ggsave(
  "GO_Downregulated.png",
  p_go_down,
  width = 8,
  height = 6,
  dpi = 300
)
list.files()
p_kegg_up <- dotplot(
  ekegg_up,
  showCategory = 15,
  title = "KEGG Pathways - Upregulated Genes"
)

p_kegg_up

ggsave(
  "KEGG_Upregulated.png",
  p_kegg_up,
  width = 8,
  height = 6,
  dpi = 300
)
p_kegg_down <- dotplot(
  ekegg_down,
  showCategory = 15,
  title = "KEGG Pathways - Downregulated Genes"
)

p_kegg_down

ggsave(
  "KEGG_Downregulated.png",
  p_kegg_down,
  width = 8,
  height = 6,
  dpi = 300
)
list.files()
dir.create("RNA-seq-Differential-Expression")
dir.exists("RNA-seq-Differential-Expression")
dir.create("RNA-seq-Differential-Expression/scripts")
dir.create("RNA-seq-Differential-Expression/results")
dir.create("RNA-seq-Differential-Expression/figures")
file.copy(
  "DEG_results.csv",
  "RNA-seq-Differential-Expression/results/DEG_results.csv"
)
file.copy(
  c("GO_Upregulated.csv", "GO_Downregulated.csv"),
  "RNA-seq-Differential-Expression/results"
)
file.copy(
  c("KEGG_Upregulated.csv", "KEGG_Downregulated.csv"),
  "RNA-seq-Differential-Expression/results"
)
file.copy(
  "PCA_Plot.png",
  "RNA-seq-Differential-Expression/figures/PCA_Plot.png"
)
file.copy(
  "Volcano_Plot.png",
  "RNA-seq-Differential-Expression/figures/Volcano_Plot.png"
)
file.copy(
  "GO_Upregulated.png",
  "RNA-seq-Differential-Expression/figures/GO_Upregulated.png"
)
file.copy(
  "GO_Downregulated.png",
  "RNA-seq-Differential-Expression/figures/GO_Downregulated.png"
)
file.copy(
  "KEGG_Upregulated.png",
  "RNA-seq-Differential-Expression/figures/KEGG_Upregulated.png"
)

file.copy(
  "KEGG_Downregulated.png",
  "RNA-seq-Differential-Expression/figures/KEGG_Downregulated.png"
)
file.exists(
  "RNA-seq-Differential-Expression/figures/GO_Downregulated.png"
)
file.copy(
  "KEGG_Upregulated.png",
  "RNA-seq-Differential-Expression/figures/KEGG_Upregulated.png"
)

file.copy(
  "KEGG_Downregulated.png",
  "RNA-seq-Differential-Expression/figures/KEGG_Downregulated.png"
)
file.exists(
  "RNA-seq-Differential-Expression/figures/KEGG_Upregulated.png"
)

file.exists(
  "RNA-seq-Differential-Expression/figures/KEGG_Downregulated.png"
)
list.files(
  "RNA-seq-Differential-Expression",
  recursive = TRUE
)
file.exists(
  "RNA-seq-Differential-Expression/scripts/RNAseq_analysis.R"
)
list.files(
  "RNA-seq-Differential-Expression",
  recursive = TRUE
)
dir.create("RNA-seq-Differential-Expression/scripts")
dir.exists("RNA-seq-Differential-Expression/scripts")
list.files("RNA-seq-Differential-Expression/scripts")
list.files("RNA-seq-Differential-Expression/scripts")
list.files("RNA-seq-Differential-Expression/scripts")
